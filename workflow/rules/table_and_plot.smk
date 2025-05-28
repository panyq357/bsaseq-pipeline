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
            -F CHROM -F POS -F REF -F ALT -F QUAL \
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


rule make_var_table:
    input:
        anno = "{prefix}.anno.replace_homolog_id_to_rap.vcf.gz",
        vcf= "{prefix}.vcf.gz"
    output:
        "{prefix}.var_table.tsv"
    script:
        "../scripts/make_var_table.R"


rule delta_allele_frequence_point:
    input:
        "{prefix}.var_table.tsv"
    output:
        "{prefix}.DAF-point.target-{target}.png"
    script:
        "../scripts/DAF_point.R"


rule pybsaseq_like:
    input:
        "{prefix}.var_table.tsv"
    output:
        svg = "{prefix}.pybsaseq_like.target-{target}.svg",
        csv = "{prefix}.pybsaseq_like.target-{target}.csv"
    params:
        window_length=2000000,
        window_step=10000,
        fisher_th=0.05
    script:
        "../scripts/pybsaseq_like.R"
