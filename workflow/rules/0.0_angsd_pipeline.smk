module angsd:
    snakefile:
        "https://raw.githubusercontent.com/zjnolen/angsd-snakemake-pipeline/migrate_to_prune_graph/workflow/Snakefile"
    config:
        config


use rule * from angsd exclude all

use rule all from angsd as angsd_all
