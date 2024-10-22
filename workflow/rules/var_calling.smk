from pathlib import Path


def get_bam_list(w):
    bam_list = []
    for x in config["groups"][w.group]:
        if Path(x).exists():
            bam_list.append(x)
        else:
            bam_list.append(f"results/bwa_mapping/{x}.bam")
    return bam_list


rule call_by_chr:
    input:
        bam_list = get_bam_list,
        ref_genome = config["bwa_index_prefix"]
    output:
        chr_vcf = temp("resources/call_by_chr/{group}.{chr}.vcf.gz")
    log:
        "logs/call_by_chr/{group}.{chr}.log"
    priority:
        10
    shell:
        '''
        bcftools mpileup \
            -O u -a AD,DP \
            -r {wildcards.chr} \
            -f {input.ref_genome} \
            {input.bam_list} \
            2> {log} \
        | bcftools call \
            -v -m -O z -a GQ \
            -o {output.chr_vcf} \
            2>> {log}
        '''

rule concat_and_norm:
    input:
        chr_vcf_list = expand("resources/call_by_chr/{{group}}.{chrom}.vcf.gz", chrom=config["chromosomes"]),
        ref_genome = config["bwa_index_prefix"]
    output:
        norm_vcf = "results/concat_and_norm/{group}.norm.vcf.gz",
        norm_vcf_tbi = "results/concat_and_norm/{group}.norm.vcf.gz.tbi"
    log:
        "logs/concat_and_norm/{group}.log"
    priority:
        100
    shell:
        '''
        bcftools concat \
            {input.chr_vcf_list} \
            2> {log} \
        | bcftools norm \
            -m-both \
            -f {input.ref_genome} \
            2>> {log} \
        | bgzip \
            -c \
            > {output.norm_vcf} \
            2>> {log}
        tabix {output.norm_vcf} 2>> {log}
        '''

