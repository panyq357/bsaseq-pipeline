config <- list(
    mc.cores = snakemake@threads,
    vcf = snakemake@input$vcf_anno,
    outdir = snakemake@output$outdir,
    vep_fields = snakemake@params$vep_fields,
    vline_list=NULL,  #' vline_list example: list(c(chr="chr08", pos=1234567))
    window_length=2E6,
    window_step=1E4,
    fisher_th=0.05
)

main <- function() {

    if (!dir.exists(config$outdir)) dir.create(config$outdir, recursive=T)

    log <- file(file.path(config$outdir, "log.txt"), "wt")
    sink(log, type="message")

    vcf <- readr::read_tsv(config$vcf, comment="##")

    # Remove multi allelic variants.
    if (length(grep(",", vcf$ALT)) > 0) {
        vcf <- vcf[-grep(",", vcf$ALT),]
    }

    var_table <- vcf[c("#CHROM", "POS", "REF", "ALT", "QUAL")]

    ad_table <- get_ad_table(vcf)

    csq_string <- get_csq_string(vcf)
    var_table <- cbind(var_table, CSQ=csq_string, ad_table)
    var_table$HIGH_OR_MODERATE <- FALSE
    var_table$HIGH_OR_MODERATE[grep("(HIGH|MODERATE)", var_table$CSQ)] <- TRUE

    sample_num <- length(grep("\\.REF", names(ad_table)))
    if (sample_num != 2) {
        readr::write_csv(var_table, file.path(config$outdir, "var_table.csv"))
        return()
    }

    var_table$FISHER <- get_ad_fisher(ad_table, mc.cores=config$mc.cores)

    window_info_df <- sliding_stats(var_table, window_length=config$window_length, window_step=config$window_step, fisher_th=config$fisher_th)

    pdf(file.path(config$outdir, "plot.pdf"), width=24, height=3)
    print(plot_sliding_window(
        window_info_df,
        title=sprintf("window_length=%d, window_step=%d, fisher_th=%s", config$window_length, config$window_step, config$fisher_th),
        vline_list=config$vline_list)
    )
    dev.off()

    readr::write_csv(window_info_df, file.path(config$outdir, "window_info_df.csv"))
    readr::write_csv(var_table, file.path(config$outdir, "var_table.csv"))
}

get_ad_table <- function(vcf) {
    bulks <- names(vcf)[10:ncol(vcf)]

    ad_pos <- sapply(strsplit(vcf$FORMAT, ":"), function(x) which(x == "AD"))

    ad <- list()
    for (bulk in bulks) {
        ad[[bulk]] <- vcf[[bulk]] |>
            strsplit(":") |>
            mapply(function(x, y) x[y], x=_, y=ad_pos) |>
            strsplit(",") |>
            lapply(as.numeric)
    }
    
    for (name in names(ad))
        ad[[name]] <- do.call(rbind, ad[[name]])

    ad_table <- do.call(cbind, ad) |>
        as.data.frame()

    names(ad_table) <- paste(rep(names(ad), each=2), c("REF", "ALT"), sep=".")

    for (name in names(ad)) {
        ref <- ad_table[[sprintf("%s.REF", name)]]
        alt <- ad_table[[sprintf("%s.ALT", name)]]
        ad_table[[sprintf("%s.AF", name)]] <- alt / (ref+alt)
    }
    return(ad_table)
}

get_csq_string <- function(vcf) {
    info <- strsplit(vcf$INFO, "[;=]")
    csq_pos <- sapply(info, function(x) which(x == "CSQ"))
    csq_string <- mapply(function(x, y) x[y], info, csq_pos + 1, SIMPLIFY=TRUE)
    return(csq_string)
}

get_ad_fisher <- function(ad_table, mc.cores) {

    ad_col_num <- grep("(\\.REF|\\.ALT)", names(ad_table))

    if (length(ad_col_num) != 4)
        stop("Error: There should be two and only two bulks!")

    ad_matrix_list <- apply(ad_table[ad_col_num], 1, function(x) matrix(x, nrow=2, byrow=T), simplify=FALSE)
    pvalue_vect <- parallel::mclapply(ad_matrix_list, function(x) fisher.test(x)$p.value, mc.cores=mc.cores) |>
        unlist()

    return(pvalue_vect)
}

remove_noise <- function(var_table) {
    ref_col_num <- grep(".REF", names(var_table))

    # Remove SNP with 0 REF in both bulk.
    var_table <- subset(var_table, !(var_table[[ref_col_num[[1]]]] == 0 & var_table[[ref_col_num[[2]]]] == 0))
    return(var_table)
}

sliding_stats <- function(var_table, window_length, window_step, fisher_th) {

    var_table <- remove_noise(var_table)

    chr_range_list <- split(var_table, var_table[["#CHROM"]]) |>
        lapply(function(x) range(x$POS))    

    chr_sliding <- function(chr) {

        chr_table <- subset(var_table, `#CHROM` == chr)
        
        window_info_list <- list()

        window_start <- seq(
            from = (chr_range_list[[chr]][[1]] %/% window_step) * window_step,
            to = (((chr_range_list[[chr]][[2]] - window_length) %/% window_step) + 1) * window_step,
            by = window_step
        )
        window_end <- window_start + window_length

        for (i in 1:length(window_start)) {
            window_fisher <- chr_table$FISHER[chr_table$POS > window_start[i] & chr_table$POS < window_end[i]]
            window_info_list[[length(window_info_list) + 1]] <- list(
                chr = chr,
                window_start = window_start[i],
                window_end = window_end[i],
                total_snp = length(window_fisher),
                s_snp = sum(window_fisher < fisher_th)
            )
        }

        return(window_info_list)
    }

    window_info_list <- parallel::mclapply(names(chr_range_list), chr_sliding, mc.cores=config$mc.cores) |>
        unlist(recursive=F)

    window_info_df <- data.table::rbindlist(window_info_list)
    window_info_df$window_midpoint <- with(window_info_df, (window_start + window_end)/2)
    window_info_df$ratio <- with(window_info_df, s_snp / total_snp)

    return(window_info_df)
}


#' vline_list example: list(c(chr="chr08", pos=1234567))
plot_sliding_window <- function(window_info_df, title, vline_list=NULL) {
    plot <- ggplot2::ggplot(window_info_df) +
        ggplot2::geom_line(ggplot2::aes(x=window_midpoint, y=ratio)) +
        ggplot2::facet_grid(~chr, scales="free_x") +
        ggplot2::scale_x_continuous(breaks=seq(0, 1E8, by=1E7), labels=sprintf("%dM", seq(0, 100, by=10))) +
        ggplot2::labs(x="", y="sSNP/totalSNP", title=title)

    if (length(vline_list) > 0) {
        for (vline in vline_list) {
            plot <- plot + 
                ggplot2::geom_vline(data=subset(window_info_df, chr %in% vline[["chr"]]), ggplot2::aes(xintercept=as.numeric(vline[["pos"]])), color="red")
        }
    }
    return(plot)
}

main()

