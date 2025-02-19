
rule samtools_index:
    input:
        "{prefix}.cram"
    output:
        "{prefix}.cram.crai"
    resources:
        io = 50
    shell:
        "samtools index {input}"


rule samtools_flagstat:
    input:
        cram = "{prefix}.cram",
        crai = "{prefix}.cram.crai"
    output:
        "{prefix}.flagstat.txt"
    priority:
        80
    resources:
        io = 50
    shell:
        "samtools flagstat {input.cram} > {output}"

