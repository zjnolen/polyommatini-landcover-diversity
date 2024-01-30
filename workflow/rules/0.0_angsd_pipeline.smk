# This section loads the main ANGSD pipeline used for most analyses

module angsd:
    snakefile:
        github(
            "zjnolen/angsd-snakemake-pipeline", path="workflow/Snakefile", tag="v0.2.0"
        )
    config:
        config


use rule * from angsd exclude all, samtools_subsample


use rule all from angsd as angsd_all
