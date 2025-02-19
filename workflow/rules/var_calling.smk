from pathlib import Path


def get_cram_list(w):
    cram_list = []
    for x in config["groups"][w.group]:
        if Path(x).exists():
            cram_list.append(x)
        else:
            cram_list.append(f"results/read_mapping/{x}.cram")
    return cram_list


def get_call_io(w):
    io = len(get_cram_list(w)) * 10
    if io < 100:
        return io
    else:
        return 100


rule bcftools_call_by_chr:
    input:
        cram_list = get_cram_list,
        crai_list = lambda w: [x + ".crai" for x in get_cram_list(w)],
        genome = config["genome"]
    output:
        temp("resources/call_by_chr/{group}.{chr}.vcf.gz")
    log:
        temp("resources/call_by_chr/{group}.{chr}.log")
    priority:
        10
    resources:
        io = get_call_io
    shell:
        '''
        bcftools mpileup -O u -a AD,DP -r {wildcards.chr} -f {input.genome} {input.cram_list} 2> {log} \
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
        "results/var_calling/{grp}.log",
    priority:
        100
    shell:
        '''
        bcftools concat {input.chr_vcf_list} 2> {log} \
        | bgzip -c > {output} 2>> {log}
        '''


