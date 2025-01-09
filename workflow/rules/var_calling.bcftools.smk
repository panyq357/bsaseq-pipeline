from pathlib import Path


def get_bam_list(w):
    bam_list = []
    for x in config["groups"][w.group]:
        if Path(x).exists():
            bam_list.append(x)
        else:
            bam_list.append(f"results/read_mapping/{x}.bam")
    return bam_list


rule bcftools_call_by_chr:
    input:
        bam_list = get_bam_list,
        bai_list = lambda w: [x + ".bai" for x in get_bam_list(w)],
        genome = config["genome"]
    output:
        temp("resources/call_by_chr/{group}.{chr}.vcf.gz")
    log:
        "logs/bcftools/call_by_chr/{group}.{chr}.log"
    priority:
        10
    shell:
        '''
        bcftools mpileup -O u -a AD,DP -r {wildcards.chr} -f {input.genome} {input.bam_list} 2> {log} \
        | bcftools call -v -m -O z -a GQ -o {output} 2>> {log}
        '''


rule bcftools_concat:
    input:
        chr_vcf_list = expand("resources/call_by_chr/{{grp}}.{chrom}.vcf.gz", chrom=config["chromosomes"]),
    output:
        "results/var_calling/{grp}.vcf.gz",
    wildcard_constraints:
        grp = "[^.]+"
    log:
        "logs/bcftools/concat/{grp}.log"
    priority:
        100
    shell:
        '''
        bcftools concat {input.chr_vcf_list} 2> {log} \
        | bgzip -c > {output} 2>> {log}
        '''


