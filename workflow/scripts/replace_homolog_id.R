id_converter <- readr::read_tsv(snakemake@input$id_converter, col_names = c("From", "To")) |> with(data = _, setNames(To, From))

vcf <- readr::read_tsv(snakemake@input$vcf, comment = "##")

anno <- vroom::vroom(snakemake@input$anno)

out_split <- strsplit(vcf$INFO, ",")

for (i in seq_along(out_split)) {
  cat(sprintf("\r%d/%d", i, length(out_split)))
  inner_split <- strsplit(out_split[[i]], "\\|")
  for (j in seq_along(inner_split)) {
    original_id <- inner_split[[j]][[3]]
    if (original_id == "") {
      next
    }
    replace_id <- id_converter[original_id]
    if (grepl("MODERATE|HIGH", inner_split[[j]][[1]])) {
      anno_str <- anno[match(replace_id, anno[[snakemake@params$anno_id_col]]), snakemake@params$anno_info_col] |>
        unlist() |>
        paste(collapse = ", ")
      inner_split[[j]][[3]] <- sprintf("%s, %s, %s", original_id, replace_id, anno_str)
    } else {
      inner_split[[j]][[3]] <- sprintf("%s, %s", original_id, replace_id)
    }
  }
  out_split[[i]] <- sapply(inner_split, paste, collapse = "|")
}

vcf$INFO <- sapply(out_split, paste, collapse = "$")

readr::write_tsv(vcf, snakemake@output[[1]])
