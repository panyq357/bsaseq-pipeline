rule VariantsToTable:
    input:
        vcf = "{prefix}.vcf.gz",
        tbi = "{prefix}.vcf.gz.tbi",
        genome = config["genome"]
    output:
        temp(pipe("{prefix}.tsv"))
    shell:
        '''
        gatk VariantsToTable \
            -R {input.genome} \
            -V {input.vcf} \
            -F CHROM -F POS -F REF -F ALT -F QUAL -F CSQ \
            -GF AD -GF DP -GF GQ -GF PL \
            -O /dev/stdout >> {output}
        '''


rule gzip_table:
    input:
        "{prefix}.tsv"
    output:
        "{prefix}.tsv.gz"
    shell:
        "gzip -c {input} > {output}"


rule delta_allele_frequence_point:
    input:
        "{prefix}.tsv.gz"
    output:
        csv_gz = "{prefix}.DAF-point.target-{target}.csv.gz",
        png = "{prefix}.DAF-point.target-{target}.png"
    script:
        "../scripts/DAF_point.R"


# rule replace_homolog_id:
#     input:
#         id_converter = lambda w: config["id_converter"][w.to],
#         table = "{prefix}.tsv.gz"
#     output:
#         "{prefix}.replace_homolog_id_to_{to}.tsv.gz"
#     script:
#         "../scripts/replace_homolog_id.R"
# 
# 
# rule fisher_table:
#     input:
#         "results/table/{grp}.tsv.gz"
#     output:
#         outdir = directory("results/fisher_table/{grp}"),
#         var_table = "results/fisher_table/{grp}/var_table.csv"
#     threads:
#         10
#     params:
#         vep_fields = config["vep_fields"]
#     script:
#         "../scripts/fisher_table.R"
# 
# 
