anno <- readr::read_tsv(snakemake@input$anno, col_names=c("GeneID", "Symbol", "Description"))

id_converter <- readr::read_tsv(snakemake@input$id_converter, col_names=c("From", "To"))
id_converter <- setNames(id_converter$To, id_converter$From)

id_to_anno <- setNames(sprintf("%s = %s", anno$Symbol, anno$Description), anno$GeneID)

anno_res_list <- list()

for (sheet in readxl::excel_sheets(snakemake@input$xlsx)) {
  df <- readxl::read_excel(snakemake@input$xlsx, sheet=sheet)

  id_list <- stringr::str_extract_all(df$INFO, "(MODERATE|HIGH)\\|[^|]+\\|([^|]+)\\|") |>
    lapply(stringr::str_extract, "\\|(OsZH11[^|]+)\\|", group=1)

  anno_list <- list()

  for (i in 1:length(id_list)) {
    if (length(id_list[[i]]) == 0) {
      anno_list[[i]] <- ""
    } else {
      anno_list[[i]] <- sprintf("%s = %s = %s", id_list[[i]], id_converter[id_list[[i]]], id_to_anno[id_converter[id_list[[i]]]])
    }
  }

  anno_res_list[[sheet]] <- tibble::add_column(df, ANNOTATION=anno_list |> sapply(paste, collapse="|"), .before="INFO")
}

writexl::write_xlsx(anno_res_list, snakemake@output$xlsx)
