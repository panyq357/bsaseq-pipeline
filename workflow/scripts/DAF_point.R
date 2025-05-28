library(panlabr)
library(ggplot2)

tab <- readr::read_tsv(snakemake@input[[1]])

samples <- names(tab)[grep(".AD", names(tab))] |> sub(".AD", "", x = _)

target <- snakemake@wildcards$target

target_af <- tab[[sprintf("%s.AF", target)]]

rest_af_mat <- tab[sprintf("%s.AF", samples[samples != target])]

rest_mean_af <- rowSums(rest_af_mat) / ncol(rest_af_mat)

df <- data.frame(
  chrom = tab[["#CHROM"]],
  pos = tab[["POS"]],
  delta_af = target_af - rest_mean_af
)
png(snakemake@output[[1]], width = 3000, height = 1500)
print(
  ggplot(df) +
    geom_point(mapping = aes(x = pos, y = delta_af), alpha = 0.2) +
    facet_wrap(~chrom, nrow = 3, scales = "free_x") +
    scale_y_continuous(limits = c(-1, +1))
)
dev.off()

