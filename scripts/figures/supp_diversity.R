# Script to plot violin plots of various genetic diversity statistics produced
# by ANGSD

library(ggplot2)
library(ggridges)
library(Hmisc)
library(dplyr)

setwd("~/working/bioinfo/projects/landuse-manuscript")

species <- c("picarus","pargus","csemiargus")

# S2-5. Pi, Theta, Heterozygosity Tajima's

thetas <- c()

for (s in species) {
  df <- read.table(
    paste0(
      "results/datasets/landuse-",s,"/analyses/thetas/landuse-",s,"_thetas_combined.pestPG"
      ),
      header = FALSE
    )
  df$species <- s
  thetas <- rbind(thetas, df)
}

colnames(thetas) <- c("window", "chr", "win_center", "tW", "tP", "tF", "tH",
                      "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites",
                      "pop", "species")

thetas <- thetas[!thetas$nSites < 10000, ]

thetas$pi <- thetas$tP / thetas$nSites
thetas$watterson <- thetas$tW / thetas$nSites

popmeans <- thetas %>%
  group_by(species, pop) %>%
  summarize(pi = mean_cl_boot(pi, B = 10000),
            theta = mean_cl_boot(watterson, B = 10000),
            tajima = mean_cl_boot(Tajima, B = 10000))

hzcs <- read.table("results/datasets/landuse-csemiargus/analyses/heterozygosity/landuse-csemiargus.ilCyaSemi1.1_all_allsites-filts_heterozygosity.tsv", header = TRUE)
hzcs$species <- "csemiargus"
hzpa <- read.table("results/datasets/landuse-pargus/analyses/heterozygosity/landuse-pargus.ilPleArgu1.3_all_allsites-filts_heterozygosity.tsv", header = TRUE)
hzpa$species <- "pargus"
hzpi <- read.table("results/datasets/landuse-picarus/analyses/heterozygosity/landuse-picarus.ilPolIcar1.1_all_allsites-filts_heterozygosity.tsv", header = TRUE)
hzpi$species <- "picarus"

hz <- rbind(hzpi,hzpa,hzcs)


thetas$species <- factor(thetas$species, levels = c("picarus","pargus","csemiargus"), labels = c("P. icarus", "P. argus", "C. semiargus"))

hz$species <- factor(hz$species, levels = c("picarus","pargus","csemiargus"), labels = c("P. icarus", "P. argus", "C. semiargus"))
                                                  

## S2 Pi

ggplot(thetas, aes(x = pop, y = pi)) +
  geom_violin(fill = "grey70") +
  facet_grid(cols = vars(species), scales = "free_x") +
  ylab("Per-site Pairwise Nucleotide Diversity") + 
  xlab("Population") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text=element_text(face="italic",size=15),
        axis.title.y = element_text(vjust=3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  stat_summary(fun = "mean", geom = "point")
ggsave("results/figures/figureS2.svg", width = 8, height = 4)

## S3 Theta

ggplot(thetas, aes(x = pop, y = watterson)) +
  geom_violin(fill = "grey70") +
  facet_grid(cols = vars(species), scales = "free_x") + 
  ylab("Per-site Watterson's Theta") + 
  xlab("Population") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text=element_text(face="italic",size=15),
        axis.title.y = element_text(vjust=3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  stat_summary(fun = "mean", geom = "point")
ggsave("results/figures/figureS3.svg", width = 8, height = 4)

## S4 Heterozygosity

ggplot(hz, aes(x = pop, y = heterozygosity)) +
  geom_boxplot(fill = "grey70") +
  facet_grid(cols = vars(species), scales = "free_x") + 
  ylab("Individual heterozygosity per 1000bp") + 
  xlab("Population") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text=element_text(face="italic",size=15),
        axis.title.y = element_text(vjust=3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("results/figures/figureS4.svg", width = 8, height = 4)

## S5 Tajima

ggplot(thetas, aes(x = pop, y = Tajima)) +
  geom_violin(fill = "grey70") +
  facet_grid(cols = vars(species), scales = "free_x") + 
  ylab("Tajima's D") + 
  xlab("Population") +
  geom_hline(yintercept = 0, linetype='dashed') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text=element_text(face="italic",size=15),
        axis.title.y = element_text(vjust=3),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  stat_summary(fun = "mean", geom = "point")
ggsave("results/figures/figureS5.svg", width = 8, height = 4)
