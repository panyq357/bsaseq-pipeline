library(ggplot2)

if (!dir.exists(snakemake@output$outdir)) dir.create(snakemake@output$outdir, recursive = TRUE)

vcf <- panlabr::read_vcf(snakemake@input$vcf)

vcf[["#CHROM"]] <- factor(vcf[["#CHROM"]], levels=snakemake@params$chromosomes)

ad_string_mat <- panlabr::vcf_to_field_matrix(vcf, "AD")

ad_3d_array <- panlabr::ad_string_mat_to_ad_3d_array(ad_string_mat)

dp_mat <- panlabr::ad_3d_array_to_dp_mat(ad_3d_array)

snp_index_mat <- panlabr::ad_3d_array_to_snp_index_mat(ad_3d_array)


for (sample in colnames(snp_index_mat)) {

  excel <- readxl::read_excel(snakemake@input$xlsx, sheet = sample)
  excel[["#CHROM"]] <- factor(excel[["#CHROM"]], levels=snakemake@params$chromosomes)
  excel <- subset(excel, ANNOTATION != "")

  delta_snp_index <- snp_index_mat[, sample] - rowMeans(snp_index_mat[, !grepl(sample, colnames(snp_index_mat))])
  dp <- dp_mat[, sample]
  df <- cbind(
    vcf[, c("#CHROM", "POS")],
    `Delta SNP-Index` = delta_snp_index,
    `Sequencing Depth` = dp
  )

  excel$`Delta SNP-Index` <- df$`Delta SNP-Index`[
    match(sprintf("%s:%s", excel$`#CHROM`, excel$POS), sprintf("%s:%s", df$`#CHROM`, df$POS))
  ]

  plt_obj <- ggplot(df) +
    geom_point(aes(x=POS, y=`Delta SNP-Index`, alpha=`Sequencing Depth`)) +
    geom_point(data=excel, mapping=aes(x=POS, y=`Delta SNP-Index`), shape=1, color="red", size=3, stroke=1) +
    facet_wrap(~`#CHROM`, nrow=snakemake@params$nrow, scales="free_x") +
    scale_y_continuous(limits=c(0.5, 1)) +
    scale_x_continuous(breaks=c(0, 1e7, 2e7, 3e7, 4e7), labels = c("0", "1", "2", "3", "4")) +
    theme_bw() +
    labs(x = "Pos (x10MB)")

  ggsave(file.path(snakemake@output$outdir, sprintf("%s.png", sample)), plt_obj, width=snakemake@params$width, height=snakemake@params$height)
}
