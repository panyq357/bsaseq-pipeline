
rule samtools_index:
    input:
        "results/read_mapping/{sample_id}.bam"
    output:
        "results/read_mapping/{sample_id}.bam.bai"
    log:
        "logs/samtools/index/{sample_id}.log"
    shell:
        "samtools index {input} 2> {log}"


rule samtools_flagstats:
    input:
        bam = "results/read_mapping/{sample_id}.bam",
        bai = "results/read_mapping/{sample_id}.bam.bai"
    output:
        "results/read_mapping/{sample_id}.flagstats.txt"
    log:
        "logs/samtools/flagstats/{sample_id}.log"
    priority:
        80
    shell:
        '''
        samtools flagstats {input.bam} > {output} 2> {log}
        '''

