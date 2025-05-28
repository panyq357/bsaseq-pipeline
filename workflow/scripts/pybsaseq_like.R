library(panlabr)
library(GenomicRanges)

var_table <- readr::read_tsv(snakemake@input[[1]])

ad_mat_list <- panlabr::get_allele_depth_matrix_list(var_table[grep("\\.AD", names(var_table))])

target_mat <- ad_mat_list[[snakemake@wildcards$target]]

rest_mat_list <- ad_mat_list[names(ad_mat_list)[names(ad_mat_list) != snakemake@wildcards$target]]

rest_mat <- Reduce(`+`, rest_mat_list)

fisher_p <- numeric(nrow(target_mat))

for (i in 1:nrow(target_mat)) {
  cat(sprintf("\r%d", i))
  m <- rbind(target_mat[i,], rest_mat[i,])
  m <- m[1:2,!apply(m, 2, anyNA)]
  fisher_p[i] <- fisher.test(m)$p.value
}

var_table$FisherP <- fisher_p

s_snp_over_total_snp <- sliding_window(var_table, window = snakemake@params$window_length, step = snakemake@params$window_step, fun = function(df) {
  mean(df$FisherP < snakemake@params$fisher_th)
})

svg(snakemake@output$svg, width=38, height=4)
  plot_sliding_window(s_snp_over_total_snp)
dev.off()

readr::write_csv(s_snp_over_total_snp, snakemake@output$csv)
