# Plotting scripts

This folder contains scripts for plotting the figures in the manuscript. Most
of these are run in R using the conda environment described in [`environment.yaml`](environment.yaml),
with the exception of [`supp_admix.R`](supp_admix.R), which uses simply base R
v3.6.

The scripts relate to the manuscript plots as follows:

Main:

- Figure 1: [`landscape_diversity.R`](landscape_diversity.R)
- Figure 2: [`pop_structure.R`](pop_structure.R) (for the PCA, the plots were
  taken directly from the angsd-snakemake-pipeline outputs and edited manually)
- Figure 3: [`landscape_differentiation.R`](landscape_differentiation.R)
- Figure 4: [`inbreeding_plots.R`](inbreeding_plots.R)
- Figure 5: [`downsampling_plots.R`](downsampling_plots.R)

Supplementary:

- Figure S1: Created with QGIS
- Figures S2-5: [`supp_diversity.R`](supp_diversity.R)
- Figure S6-7: [`supp_admix.R`](supp_admix.R)
- Figure S8: Created with QGIS
- Figure S9-10: [`inbreeding_plots.R`](inbreeding_plots.R)
- Figure S11: [`landscape_diversity.R`](landscape_diversity.R)
