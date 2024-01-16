This is a pipeline for BSA-seq analysis.

## How to run this pipeline

Step 1. Filling `config/reads2bwa-config.yaml` with bwa index and raw FASTQ path.

Step 2. Filling `config/bams2vcf-config.yaml` with VCF groups.

Step 3. Filling `config/vcf2fisher-config.yaml` with ref genome and GTF for VEP annotation.

Finally, run this.

```bash
snakemake --cores 20
```

## Results explained

```txt
results/
├── bwa_mapping                   # read mapping results
│   ├── ES.bam
│   ├── ES.bam.bai
│   ├── ES.flagstats.txt
│   ├── ET.bam
│   ├── ET.bam.bai
│   └── ET.flagstats.txt
├── concat_and_norm               # normalized vcf.
│   ├── NIPxLPBG.norm.vcf.gz
│   └── NIPxLPBG.norm.vcf.gz.tbi
├── fastp                         # raw FASTQ QC reports.
│   ├── ES.html
│   ├── ES.json
│   ├── ET.html
│   └── ET.json
├── fisher_table                  # Fisher exact test results and plots.
│   └── NIPxLPBG
│       ├── log.txt
│       ├── plot.pdf
│       ├── var_table.csv
│       └── window_info_df.csv
└── vep_anno                      # VEP annotation reports.
    ├── NIPxLPBG.html
    └── NIPxLPBG.norm.filter.anno.vcf.gz
```

In this pipeline, `workflow/scripts/fisher_table.R` is a lightweight
implementation of [PyBSASeq](https://doi.org/10.1186/s12859-020-3435-8),
and it produces similar results.

