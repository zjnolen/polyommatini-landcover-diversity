module angsd:
    snakefile:
        "https://raw.githubusercontent.com/zjnolen/angsd-snakemake-pipeline/3a0bbfbb7ec5bd36c0e0601d8013c957accf0722/workflow/Snakefile"
    config:
        config


use rule * from angsd exclude all, samtools_subsample


use rule all from angsd as angsd_all
