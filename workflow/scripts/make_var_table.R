library(panlabr)

vcf <- read_vcf(snakemake@input$vcf)
anno <- read_vcf(snakemake@input$anno)

ad_str_df <- get_allele_depth_string_df(vcf)
ad_mat_list <- get_allele_depth_matrix_list(ad_str_df)
af_df <- get_allele_freqency_df(ad_mat_list)

vcf <- cbind(vcf[c(1,2,4,5,6)], ad_str_df, af_df)
vcf$CSQ <- sub("CSQ=", "", anno$INFO)
vcf$HIGH_OR_MODERATE <- FALSE
vcf$HIGH_OR_MODERATE[grep("(HIGH|MODERATE)", vcf$CSQ)] <- TRUE

readr::write_tsv(vcf, snakemake@output[[1]])
