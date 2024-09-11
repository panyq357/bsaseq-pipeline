rule filter_vcf:
    input:
        vcf_norm = "results/concat_and_norm/{group}.norm.vcf.gz",
    output:
        vcf_filter = temp("resources/filter_vcf/{group}.norm.filter.vcf")
    log:
        "logs/filter_vcf/{group}.log"
    params:
        filter_expression = config["bcftools_filter_expression"]
    shell:
        '''
        bcftools view \
            -e {params.filter_expression:q} \
            {input.vcf_norm} \
            > {output.vcf_filter} \
            2> {log}
        '''


rule sort_gtf:
    input:
        gtf = config["gtf"]
    output:
        sorted_gtf = "resources/sorted.gtf.gz"
    log:
        "logs/sort_gtf.log"
    shell:
        '''
        (zcat < {input.gtf} || cat < {input.gtf}) \
        | grep -v '^#' \
        | sort -k1,1 -k4,4n -k5,5n -t$'\\t' \
        | bgzip -c > {output.sorted_gtf} 2> {log}
        tabix -p gff {output.sorted_gtf} 2>> {log}
        '''


rule vep_anno:
    input:
        vcf_filter = "resources/filter_vcf/{group}.norm.filter.vcf",
        ref_genome = config["bwa_index_prefix"],
        sorted_gtf = "resources/sorted.gtf.gz"
    output:
        vcf_anno = "results/vep_anno/{group}.norm.filter.anno.vcf.gz",
        vep_summary = "results/vep_anno/{group}.html"
    log:
        "logs/vep_anno/{group}.log"
    threads:
        4
    params:
        fields = config["vep_fields"]
    shell:
        '''
        vep --fork {threads} \
            --gtf {input.sorted_gtf} --fasta {input.ref_genome} \
            --input_file {input.vcf_filter} --format vcf \
            --output_file STDOUT --force_overwrite --vcf \
            --stats_file {output.vep_summary} \
            --per_gene --fields "{params.fields}" \
            2> {log} \
        | bgzip -c > {output.vcf_anno} 2>> {log}
        '''


rule fisher_table:
    input:
        vcf_anno = "results/vep_anno/{group}.norm.filter.anno.vcf.gz"
    output:
        outdir = directory("results/fisher_table/{group}"),
        var_table = "results/fisher_table/{group}/var_table.csv"
    threads:
        10
    params:
        vep_fields = config["vep_fields"]
    script:
        "../scripts/fisher_table.R"

