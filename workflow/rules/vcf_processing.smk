
rule tabix:
    input:
        "{prefix}.gz"
    output:
        "{prefix}.gz.tbi"
    shell:
        "tabix {input}"


rule filter_vcf:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.filter.vcf.gz"
    params:
        filter_expression = config["bcftools_filter_expression"]
    shell:
        '''
        bcftools view -e {params.filter_expression:q} {input} | bgzip > {output}
        '''

