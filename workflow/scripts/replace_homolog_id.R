
id_converter <- readr::read_tsv(snakemake@input$id_converter, col_names=c("From", "To")) |> with(data=_, setNames(To, From))

var_table <- readr::read_csv(snakemake@input$var_table)

replaced_CSQ <- strsplit(var_table$CSQ, ",") |>
    lapply(function(x) {
    y <- strsplit(x, "\\|")
    for (i in 1:length(y)) {
        y[[i]][[3]] <- sprintf("%s/%s", y[[i]][[3]], id_converter[y[[i]][[3]]])
    }
    sapply(y, paste, collapse="|")
    }) |>
    lapply(paste, collapse=",") |>
    unlist()

var_table$CSQ <- replaced_CSQ

readr::write_csv(var_table, snakemake@output[[1]])
