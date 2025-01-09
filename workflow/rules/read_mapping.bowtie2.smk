
rule bowtie2_mapping:
    input:
        r1_fastp = rules.fastp_pe.output.r1_fastp,
        r2_fastp = rules.fastp_pe.output.r2_fastp,
        bowtie2_index = [config["genome"] + suffix for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]]
    output:
        bam = "results/read_mapping/{sample_id}.bam",
        bai = "results/read_mapping/{sample_id}.bam.bai"
    params:
        index_prefix = config["genome"],
        sort_threads = config["samtools"]["sort_threads"],
        sort_mem_per_thread = config["samtools"]["sort_mem_per_thread"]
