rule filter_variants:
    input:
        vcf = "results/var_calling/{grp}.filter.vcf.gz",
        anno = "results/var_calling/{grp}.filter.vep_anno.vcf.gz"
    output:
        xlsx = "results/mutant_only_tables/{grp}.filter.xlsx"
    params:
        filter_expr = config["filter_variants"]["R_code"]
    script:
        "../scripts/filter_variants.R"


rule add_gene_function_annotation:
    input:
        xlsx = "results/mutant_only_tables/{grp}.filter.xlsx",
        id_converter = config["id_converter"],
        anno = config["gene_function_annotation"]
    output:
        xlsx = "results/mutant_only_tables/{grp}.filter.anno.xlsx",
    script:
        "../scripts/add_gene_function_annotation.R"


rule plot_snp_index_dot_plot:
    input:
        vcf = "results/var_calling/{grp}.filter.vcf.gz",
        xlsx = "results/mutant_only_tables/{grp}.filter.anno.xlsx"
    output:
        outdir = directory("results/mutant_only_tables/{grp}.filter.snp_index_plot"),
    params:
        chromosomes = config["chromosomes"],
        nrow = 2,
        width = 12,
        height = 4
    script:
        "../scripts/plot_snp_index_dot_plot.R"


rule insert_plot_to_excel:
    input:
        xlsx = "results/mutant_only_tables/{grp}.filter.anno.xlsx",
        plot_dir = "results/mutant_only_tables/{grp}.filter.snp_index_plot"
    output:
        xlsx = "results/mutant_only_tables/{grp}.filter.anno.add_plot.xlsx",
    script:
        "../scripts/insert_plot_to_excel.py"
