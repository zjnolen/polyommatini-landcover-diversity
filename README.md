# Polyommatini Landscape Genomics

This repository contains the scripts utilized for processing the data in our
manuscript focusing on the effects of land use on genetic diversity and
differentiation in three blue-wing butterfly species.

All of the analyses were performed using a Snakemake based workflow, and can
be reproduced with this repository, along with a copy of the raw data and
external resources (land use raster, reference genomes, repeat libraries)
described in the manuscript.

The majority of the workflow was carried out using
[`zjnolen/angsd-snakemake-pipeline`](https://www.github.com/zjnolen/angsd-snakemake-pipeline)
v0.2.0. Additional analyses were added as an extension of this workflow, and
can be viewed in the [`workflow`](workflow) folder. The
sample lists and configuration files used to configure the combined workflow
can be found in the [`config`](config) folder. The table with sampling site
names and coordinates can be found in [`resources/gis`](resources/gis).

Figures for the manuscript were generated outside the workflow, with scripts
available in [`scripts/figures`](scripts/figures).

If you have any questions about the analyses we performed in this manuscript,
please feel free to reach out to the corresponding author or through the issues
section of this repository. If you would like to use the main genotype
likelihood based Snakemake workflow on your own WGS data, please check out
the [main repository for that](https://www.github.com/zjnolen/angsd-snakemake-pipeline).

## Reproducing our analyses

To reproduce our analyses for one of the species, you will first need to
[install Snakemake to your system](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html),
then clone this repository to your machine:

```bash
git clone
```

Download the raw data from NCBI, and change the paths in [`config/units.tsv`](config/units.tsv)
to correspond to the locations of the fastq files you download.

Download the additional external resources needed and place them in the correct
folder:

- [Swedish Land Cover Database 2018 v1.1](https://geodatakatalogen.naturvardsverket.se/geonetwork/srv/swe/catalog.search#/metadata/8853721d-a466-4c01-afcc-9eae57b17b39)
  -> `resources/gis/nmd2018bas_ogeneraliserad_v1_1.tif`
- Species reference genome:
  - [*P. icarus* GCA_937595015.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_937595015.1/)
    -> `resources/ref/GCA_937595015.1_ilPolIcar1.1_genomic.fna`
  - [*P. argus* GCA_905404155.3](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_905404155.3/)
    -> `resources/ref/GCA_905404155.3_ilPleArgu1.3_genomic.fna`
  - [*C. semiargus* GCA_905187585.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_905187585.1/)
    -> `resources/ref/GCA_905187585.1_ilCyaSemi1.1_genomic.fna`
- Species repeat libraries from Darwin Tree of Life:
  - [*P. icarus* GCA_937595015.1](https://ftp.ebi.ac.uk/pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/polyommatus_icarus/GCA_937595015.1.repeatmodeler.fa)
    -> `resources/repeatlibs/GCA_937595015.1.repeatmodeler.fa`
  - [*P. argus* GCA_905404155.1](https://ftp.ebi.ac.uk/pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/plebejus_argus/GCA_905404155.1.repeatmodeler.fa)
    -> `resources/repeatlibs/GCA_905404155.1.repeatmodeler.fa`
  - [*C. semiargus* GCA_905187585.1](https://ftp.ebi.ac.uk/pub/databases/ensembl/repeats/unfiltered_repeatmodeler/species/cyaniris_semiargus/GCA_905187585.1.repeatmodeler.fa)
    -> `resources/repeatlibs/GCA_905187585.1.repeatmodeler.fa`

Run the workflow from the cloned repository with Snakemake:

```bash
snakemake --configfile config/config_<species>.yaml
```

You will most likely need to use some additional options to adapt Snakemake to
your cluster's configuration. See [Snakemake's documentation](https://snakemake.readthedocs.io)
for more information on this.
