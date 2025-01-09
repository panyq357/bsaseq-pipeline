
rule HaplotypeCaller:
    input:
        bam = "results/read_mapping/{sample_id}.bam",
        bai = "results/read_mapping/{sample_id}.bam.bai",
        genome = config["genome"]
    output:
        temp("resources/GATK/HaplotypeCaller/{sample_id}.g.vcf.gz")
    resources:
        mem_gb = config["GATK"]["HaplotypeCaller"]["mem_gb"]
    threads:
        1
    log:
        "logs/GATK/HaplotypeCaller/{sample_id}.log"
    shell:
        '''
        gatk --java-options "-Xmx{resources.mem_gb}G" HaplotypeCaller \
            --reference {input.genome} \
            --input {input.bam} \
            --emit-ref-confidence GVCF \
            --native-pair-hmm-threads {threads} \
            --output {output} 2> {log}
        '''


def get_gvcf_list(w):
    bam_list = []
    for x in config["groups"][w.group]:
        if Path(x).exists():
            bam_list.append(x)
        else:
            bam_list.append(f"results/read_mapping/{x}.bam")
    return bam_list


