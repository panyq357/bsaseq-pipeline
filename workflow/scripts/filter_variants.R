vcf <- panlabr::read_vcf(snakemake@input$vcf)

anno <- panlabr::read_vcf(snakemake@input$anno)

ad_string_mat <- panlabr::vcf_to_field_matrix(vcf, "AD")

ad_3d_array <- panlabr::ad_string_mat_to_ad_3d_array(ad_string_mat)

dp_mat <- panlabr::ad_3d_array_to_dp_mat(ad_3d_array)

snp_index_mat <- panlabr::ad_3d_array_to_snp_index_mat(ad_3d_array)

out_list <- list()

for (sample in colnames(snp_index_mat)) {
  keep_row <- list()
  for (i in seq_len(nrow(snp_index_mat))) {
    cat(sprintf("\rSample: %s, Row: %d", sample, i))
    sample_dp <- dp_mat[i, sample]
    sample_snp_index <- snp_index_mat[i, sample]
    all_snp_index <- snp_index_mat[i, ]
    if (eval(parse(text = snakemake@params$filter_expr))) {
      keep_row[[length(keep_row)+1]] <- i
    }
  }
  keep_row <- unlist(keep_row)

  out_list[[sample]] <- cbind(
    anno[keep_row, c(1, 2, 4, 5, 6, 8)],
    data.frame(ad_string_mat[keep_row, sample]) |> setNames(sprintf("Target.%s", sample)),
    ad_string_mat[keep_row, !grepl(sample, colnames(ad_string_mat))]
  )
}

writexl::write_xlsx(out_list, snakemake@output$xlsx)
