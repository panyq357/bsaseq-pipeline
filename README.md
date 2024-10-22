A pipeline for BSA-seq analysis.

## How to run this pipeline

First, run these.

```bash
cp workflow/Snakefile_template workflow/Snakefile
cp -r config_template config
```

Then filling config files in `config` directory.

Finally, run this.

```bash
snakemake --cores 20 --resources io=100
```

