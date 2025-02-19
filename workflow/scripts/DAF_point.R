library(ggplot2)

main <- function() {
    tab <- readr::read_tsv(snakemake@input[[1]])

    samples <- names(tab)[grep(".AD", names(tab))] |> sub(".AD", "", x = _)

    target <- snakemake@wildcards$target

    # Remove multiallelic variants.
    tab <- tab[-grep(",", tab$ALT), ]

    target_af <- ad_to_af(tab[[sprintf("%s.AD", target)]])

    rest_af_mat <- tab[sprintf("%s.AD", samples[samples != target])] |>
        apply(MARGIN = 2, ad_to_af)

    rest_af_mat[is.na(rest_af_mat)] <- 0

    rest_mean_af <- rowSums(rest_af_mat) / ncol(rest_af_mat)

    df <- cbind(tab[c("CHROM", "POS", "REF", "ALT")], TargetAF = target_af, RestMeanAF = rest_mean_af, DeltaAF = target_af - rest_mean_af)

    df$POS <- as.integer(df$POS)
    df |>
        panlabr::shrink_digit(3) |>
        readr::write_excel_csv(gzfile(snakemake@output$csv_gz), quote = "needed")

    png(snakemake@output$png, width = 3000, height = 1500)
    print(
        ggplot(df) +
            geom_point(mapping = aes(x = POS, y = DeltaAF), alpha = 0.2) +
            facet_wrap(~CHROM, nrow = 3, scales = "free_x") +
            scale_y_continuous(limits = c(-1, +1))
    )
    dev.off()
}

ad_to_af <- function(ad) {
    ad_mat <- stringr::str_split(ad, ",", simplify = TRUE)
    mode(ad_mat) <- "integer"
    af <- ad_mat[, 2] / rowSums(ad_mat)
    return(af)
}

if (!interactive()) {
    main()
}
