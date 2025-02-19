
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
        vep --fork {threads} \
            --gtf {input.gtf} --fasta {input.genome} \
            --input_file {input.vcf} --format vcf \
            --output_file STDOUT --force_overwrite --vcf \
            --stats_file {output.summary} \
            --per_gene --fields "{params.fields}" \
            2> {log} \
        | bgzip -c > {output.vcf} 2>> {log}
        '''

