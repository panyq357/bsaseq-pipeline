
rule sort_gtf:
    input:
        gtf = config["gtf"]
    output:
        sorted_gtf = "resources/sorted.gtf.gz"
    shell:
        '''
        (zcat < {input.gtf} || cat < {input.gtf}) \
        | grep -v '^#' \
        | sort -k1,1 -k4,4n -k5,5n -t$'\\t' \
        | bgzip -c > {output.sorted_gtf}
        tabix -p gff {output.sorted_gtf}
        '''


rule vep_anno:
    input:
        vcf = "{prefix}.vcf.gz",
        genome = config["genome"],
        gtf = "resources/sorted.gtf.gz"
    output:
        vcf = "{prefix}.anno.vcf.gz",
        summary = "{prefix}.anno.html"
    log:
        "{prefix}.anno.log"
    threads:
        4
    params:
        fields = config["vep_fields"]
    shell:
        '''
        bcftools view \
            -G {input.vcf} \
        | bcftools annotate \
            -x INFO \
        | vep \
            --fork {threads} \
            --gtf {input.gtf} \
            --fasta {input.genome} \
            --format vcf \
            --output_file STDOUT \
            --vcf \
            --force_overwrite \
            --stats_file {output.summary} \
            --per_gene \
            --fields "{params.fields}" \
            2> {log} \
        | bgzip -c > {output.vcf} 2>> {log}
        '''


rule replace_homolog_id:
    input:
        id_converter = lambda w: config["id_converter"][w.to],
        vcf = "results/var_calling/{grp}.anno.vcf.gz",
        anno_file = config["anno"]["file"],
    output:
        vcf = temp("results/var_calling/{grp}.anno.replace_homolog_id_to_{to}.vcf.not_bgzip.gz"),
    params:
        anno_id_col = config["anno"]["id_col"],
        anno_info_col = config["anno"]["info_col"]
    script:
        "../scripts/replace_homolog_id.R"


rule bgzip_and_add_header:
    input:
        anno = "results/var_calling/{grp}.anno.vcf.gz",
        anno_tbi = "results/var_calling/{grp}.anno.vcf.gz.tbi",
        replace = "results/var_calling/{grp}.anno.replace_homolog_id_to_{to}.vcf.not_bgzip.gz",
    output:
        vcf = "results/var_calling/{grp}.anno.replace_homolog_id_to_{to}.vcf.gz"
    shell:
        "(tabix -H {input.anno}; gzip -dc {input.replace}) | bgzip > {output.vcf}"

