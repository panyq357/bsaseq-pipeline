id_converter <- snakemake@params$id_converter

if (is.null(id_converter)) {
  cat("Not using id_converter\n")
} else if (!file.exists(id_converter)) {
  cat("Not using id_converter\n")
  id_converter <- NULL
} else {
  cat("Using id_converter:", id_converter)
  id_converter <- readr::read_tsv(id_converter, col_names = c("From", "To")) |> with(data = _, setNames(To, From))
}

vcf <- readr::read_tsv(snakemake@input$vcf, comment = "##")

anno <- vroom::vroom(snakemake@input$anno)

splitted_vep_fields <- strsplit(snakemake@params$vep_fields, ",")[[1]]
impact_pos <- which(splitted_vep_fields == "IMPACT")
gene_pos <- which(splitted_vep_fields == "Gene")

fetch_anno_str_from_anno <- function(query_id) {
  anno[match(query_id, anno[[snakemake@params$anno_id_col]]), snakemake@params$anno_info_col] |> unlist() |> paste(collapse = "!")
}

## First, split INFO by commas.
## out_split <- list(
##   row_of_one_variant = c(annotation_of_one_gene, annotation_of_another_gene),
##   row_of_another_variant = c(annotation_of_one_gene, annotation_of_another_gene),
## )
out_split <- strsplit(vcf$INFO, ",")

for (i in seq_along(out_split)) {

  cat(sprintf("\r%d/%d", i, length(out_split)))

  ## Second, split annotation by vertical bar.
  ## inner_split <- list(
  ##   annotation_of_one_gene = c(IMPACT, Consequence, Gene, ...)  ( order depends on vep_fields in config)
  ##   annotation_of_another_gene = c(IMPACT, Consequence, Gene, ...)  
  ## )
  inner_split <- strsplit(out_split[[i]], "\\|")

  for (j in seq_along(inner_split)) {

    original_id <- inner_split[[j]][[gene_pos]]

    if (original_id == "") { next }

    is_aa_changing <- grepl("MODERATE|HIGH", inner_split[[j]][[impact_pos]])

    if (is.null(id_converter)) {
      if (is_aa_changing) {
        anno_str <- fetch_anno_str_from_anno(original_id)
        inner_split[[j]][[gene_pos]] <- sprintf("%s!%s", original_id, anno_str)
      } else {
        NULL
      }
    } else {
      query_id <- id_converter[original_id]
      if (is_aa_changing) {
        anno_str <- fetch_anno_str_from_anno(query_id)
        inner_split[[j]][[gene_pos]] <- sprintf("%s!%s!%s", original_id, query_id, anno_str)
      } else {
        inner_split[[j]][[gene_pos]] <- sprintf("%s!%s", original_id, query_id)
      }
    }
  }
  out_split[[i]] <- sapply(inner_split, paste, collapse = "|")
}
cat("\n")

vcf$INFO <- sapply(out_split, paste, collapse = ",")

readr::write_tsv(vcf, snakemake@output[[1]], col_names=FALSE)  # Add header after using tabix.
