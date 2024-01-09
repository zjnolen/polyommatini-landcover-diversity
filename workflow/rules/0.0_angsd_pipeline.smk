module angsd:
    snakefile:
        "https://raw.githubusercontent.com/zjnolen/angsd-snakemake-pipeline/13c5c8a0127180b1e36c99746650c2bbe7b25fc8/workflow/Snakefile"
    config:
        config


use rule * from angsd exclude all, samtools_subsample


use rule all from angsd as angsd_all
