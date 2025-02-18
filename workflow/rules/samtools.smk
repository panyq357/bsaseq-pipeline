
rule samtools_index:
    input:
        "results/read_mapping/{sample_id}.cram"
    output:
        "results/read_mapping/{sample_id}.cram.crai"
    log:
        "logs/samtools/index/{sample_id}.log"
    resources:
        io = 50
    shell:
        "samtools index {input} 2> {log}"


rule samtools_flagstat:
    input:
        cram = "results/read_mapping/{sample_id}.cram",
        crai = "results/read_mapping/{sample_id}.cram.crai"
    output:
        "results/read_mapping/{sample_id}.flagstat.txt"
    log:
        "logs/samtools/flagstat/{sample_id}.log"
    priority:
        80
    resources:
        io = 50
    shell:
        '''
        samtools flagstat {input.cram} > {output} 2> {log}
        '''

