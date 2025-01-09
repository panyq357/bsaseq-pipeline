rule VariantsToTable:
    input:
        vcf = "results/var_calling/{grp}.filter.anno.vcf.gz",
        tbi = "results/var_calling/{grp}.filter.anno.vcf.gz.tbi",
        genome = config["genome"]
    output:
        temp(pipe("results/table/{grp}.tsv"))
    shell:
        '''
        gatk VariantsToTable \
            -R {input.genome} \
            -V {input.vcf} \
            -F CHROM -F POS -F REF -F ALT \
            -GF AD -GF DP -GF GQ -GF PL -F CSQ \
            -O {output}
        '''


rule gzip_table:
    input:
        "results/table/{grp}.tsv"
    output:
        "results/table/{grp}.tsv.gz"
    shell:
        "cat {input} | gzip > {output}"


rule fisher_table:
    input:
        vcf = "results/var_calling/{grp}.filter.anno.vcf.gz",
    output:
        outdir = directory("results/fisher_table/{group}"),
        var_table = "results/fisher_table/{group}/var_table.csv"
    threads:
        10
    params:
        vep_fields = config["vep_fields"]
    script:
        "../scripts/fisher_table.R"


rule replace_homolog_id:
    input:
        id_converter = lambda w: config["id_converter"][w.to],
        var_table = "results/fisher_table/{group}/var_table.csv"
    output:
        "results/fisher_table/{group}/var_table.replace_homolog_id_to_{to}.csv"
    script:
        "../scripts/replace_homolog_id.R"

