rule around_site_cover:
    input:
        raster=config["landuse-raster"],
        sites=config["field-sites"],
    output:
        tsv="results/landscape/around_site_cover.tsv",
    log:
        "logs/landscape/around_site_cover.log",
    conda:
        "../envs/landscapemets.yaml"
    resources:
        runtime=360,
    retries: 0
    script:
        "../scripts/around-site-cover.R"

localrules: diversity_landscape_modeling
rule diversity_landscape_modeling:
    input:
        cover="results/landscape/around_site_cover.tsv",
        hz="results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all_{sites}-filts_heterozygosity.tsv",
    output:
        modelout="results/datasets/{dataset}/analyses/landscape_models/{dataset}.{ref}_all_{sites}-filts_heterozygosity-landuse.tsv",
    log:
        "logs/landscape/{dataset}.{ref}_all_{sites}-filts_heterozygosity-landuse.log",
    conda:
        "../envs/landscapemets.yaml"
    retries: 0
    script:
        "../scripts/diversity-landscape-modeling.R"
