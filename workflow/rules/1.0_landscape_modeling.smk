# These rules are used to perform the models examining genetic diversity and
# differentiation responses to land cover

rule around_site_cover:
    """
    Calculate land cover proportions in radii around all sample sites.
    """
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


localrules:
    diversity_landscape_modeling,
    differentiation_landscape_modeling,


rule diversity_landscape_modeling:
    """
    Fit models with genetic diversity as a response to land cover.
    """
    input:
        cover="results/landscape/around_site_cover.tsv",
        hz="results/datasets/{dataset}/analyses/heterozygosity/{dataset}.{ref}_all_{sites}-filts_heterozygosity.tsv",
    output:
        modelout="results/datasets/{dataset}/analyses/landscape_models/{dataset}.{ref}_all_{sites}-filts_heterozygosity-landuse.tsv",
    log:
        "logs/landscape/{dataset}.{ref}_all_{sites}-filts_heterozygosity-landuse.log",
    conda:
        "../envs/landscapemodels.yaml"
    retries: 0
    script:
        "../scripts/diversity-landscape-modeling.R"


rule between_site_cover:
    """
    Calculate land cover proportions between all pairs of sample sites.
    """
    input:
        raster=config["landuse-raster"],
        sites=config["field-sites"],
    output:
        tsv="results/landscape/{dataset}_between_site_cover.tsv",
    log:
        "logs/landscape/{dataset}_between_site_cover.log",
    conda:
        "../envs/landscapemets.yaml"
    resources:
        runtime=360,
    threads: 2
    retries: 0
    params:
        pops=angsd.pop_list,
    script:
        "../scripts/between-site-cover.R"


rule differentiation_landscape_modeling:
    """
    Fit models with genetic differentiation as a response to land cover.
    """
    input:
        cover="results/landscape/{dataset}_between_site_cover.tsv",
        fst="results/datasets/{dataset}/analyses/fst/{dataset}.{ref}_poppairs_{sites}-filts.fst.global.tsv",
    output:
        modelout="results/datasets/{dataset}/analyses/landscape_models/{dataset}.{ref}_all_{sites}-filts_fst-landuse.tsv",
    log:
        "logs/landscape/{dataset}.{ref}_all_{sites}-filts_fst-landuse.log",
    conda:
        "../envs/landscapemodels.yaml"
    retries: 0
    script:
        "../scripts/differentiation-landscape-modeling.R"
