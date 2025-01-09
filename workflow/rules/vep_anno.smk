
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
        vcf = "results/var_calling/{grp}.filter.vcf.gz",
        genome = config["genome"],
        gtf = "resources/sorted.gtf.gz"
    output:
        vcf = "results/var_calling/{grp}.filter.anno.vcf.gz",
        summary = "results/var_calling/{grp}.filter.anno.html"
    log:
        "logs/vep_anno/{grp}.log"
    threads:
        4
    params:
        fields = config["vep_fields"]
    shell:
        '''
        vep --fork {threads} \
            --gtf {input.gtf} --fasta {input.genome} \
            --input_file {input.vcf} --format vcf \
            --output_file STDOUT --force_overwrite --vcf \
            --stats_file {output.summary} \
            --per_gene --fields "{params.fields}" \
            2> {log} \
        | bgzip -c > {output.vcf} 2>> {log}
        '''

