library(ggplot2)
library(dplyr)

setwd("~/working/bioinfo/projects/landuse-manuscript")

# Compiles downsampled diversity statistics and calculates a normalized root
# mean squared error from the full sample size estimates.

# Individual heterozygosity

picarus <- read.table(
  "results/datasets/landuse-picarus/analyses/heterozygosity/landuse-picarus.ilPolIcar1.1_all_allsites-filts_heterozygosity.tsv",
  header = TRUE,
  sep = '\t')
pargus <- read.table(
  "results/datasets/landuse-pargus/analyses/heterozygosity/landuse-pargus.ilPleArgu1.3_all_allsites-filts_heterozygosity.tsv",
  header = TRUE,
  sep = '\t')
csemiargus <- read.table(
  "results/datasets/landuse-csemiargus/analyses/heterozygosity/landuse-csemiargus.ilCyaSemi1.1_all_allsites-filts_heterozygosity.tsv",
  header = TRUE,
  sep = '\t')

picarus$species <- "picarus"
pargus$species <- "pargus"
csemiargus$species <- "csemiargus"
hz <- rbind(picarus,pargus,csemiargus)

sample_sizes <- hz %>% count(pop, species)
colnames(sample_sizes) <- c("pop","species","full_n")

df <- c()
means <- c()

for (s in c("picarus", "pargus", "csemiargus")) {
  species_means <- hz[hz$species == s, ] %>%
    group_by(pop) %>%
    summarise(mean = mean(heterozygosity))
  species_means$species <- s
  colnames(species_means) <- c("pop", "dataset_mean", "species")
  means <- rbind(means, species_means)
  for (i in c(1:100)) {
    for (n in c(1,2,3,4,5,6,7,8,9,10)) {
      downsampled <- hz[hz$species == s, ] %>% 
        group_by(pop) %>% 
        slice_sample(n=n) %>% 
        summarise(mean = mean(heterozygosity))
      downsampled$species <- s
      downsampled$rep <- i
      downsampled$n <- n
      df <- rbind(df,downsampled)
    }
  }
}

df <- merge(df, means, by = c('species', 'pop'))

df$n <- as.factor(df$n)

df <- merge(df, sample_sizes, by = c("species","pop"))

# Fst

picarus <- read.table(
  "results/datasets/landuse-picarus/analyses/fst/downsampled/landuse-picarus.ilPolIcar1.1_poppairs.downsampled_allsites-filts.fst.global.tsv",
  header = TRUE,
  sep = '\t')
pargus <- read.table(
  "results/datasets/landuse-pargus/analyses/fst/downsampled/landuse-pargus.ilPleArgu1.3_poppairs.downsampled_allsites-filts.fst.global.tsv",
  header = TRUE,
  sep = '\t')
csemiargus <- read.table(
  "results/datasets/landuse-csemiargus/analyses/fst/downsampled/landuse-csemiargus.ilCyaSemi1.1_poppairs.downsampled_allsites-filts.fst.global.tsv",
  header = TRUE,
  sep = '\t')

picarus_full <- read.table(
  "results/datasets/landuse-picarus/analyses/fst/landuse-picarus.ilPolIcar1.1_poppairs_allsites-filts.fst.global.tsv",
  header = TRUE,
  sep = '\t')
pargus_full <- read.table(
  "results/datasets/landuse-pargus/analyses/fst/landuse-pargus.ilPleArgu1.3_poppairs_allsites-filts.fst.global.tsv",
  header = TRUE,
  sep = '\t')
csemiargus_full <- read.table(
  "results/datasets/landuse-csemiargus/analyses/fst/landuse-csemiargus.ilCyaSemi1.1_poppairs_allsites-filts.fst.global.tsv",
  header = TRUE,
  sep = '\t')

picarus$species <- "picarus"
pargus$species <- "pargus"
csemiargus$species <- "csemiargus"

picarus_full$species <- "picarus"
pargus_full$species <- "pargus"
csemiargus_full$species <- "csemiargus"

fst <- rbind(picarus, pargus, csemiargus)

fst_full <- rbind(picarus_full, pargus_full, csemiargus_full)

fst[fst$weight.fst < 0, "weight.fst"] <- 0

fst_full[fst_full$weight.fst < 0, "weight.fst"] <- 0

fst <- merge(fst, fst_full, by = c('species','pop1','pop2'))

fst$downsample.n <- as.factor(fst$downsample.n)

fst$poppair <- paste0(fst$pop1,"-",fst$pop2)

# theta, pi, tajima's

picarus <- read.table(
  "results/datasets/landuse-picarus/analyses/thetas/downsampled/landuse-picarus.ilPolIcar1.1_all.downsampled_allsites-filts.thetaMean.50000_10000.tsv",
  header = TRUE,
  sep = '\t')
pargus <- read.table(
  "results/datasets/landuse-pargus/analyses/thetas/downsampled/landuse-pargus.ilPleArgu1.3_all.downsampled_allsites-filts.thetaMean.50000_10000.tsv",
  header = TRUE,
  sep = '\t')
csemiargus <- read.table(
  "results/datasets/landuse-csemiargus/analyses/thetas/downsampled/landuse-csemiargus.ilCyaSemi1.1_all.downsampled_allsites-filts.thetaMean.50000_10000.tsv",
  header = TRUE,
  sep = '\t')

picarus_full <- read.table(
  "results/datasets/landuse-picarus/analyses/thetas/landuse-picarus_thetas_combined.pestPG",
  header = FALSE,
  sep = '\t')
pargus_full <- read.table(
  "results/datasets/landuse-pargus/analyses/thetas/landuse-pargus_thetas_combined.pestPG",
  header = FALSE,
  sep = '\t')
csemiargus_full <- read.table(
  "results/datasets/landuse-csemiargus/analyses/thetas/landuse-csemiargus_thetas_combined.pestPG",
  header = FALSE,
  sep = '\t')

picarus$species <- "picarus"
pargus$species <- "pargus"
csemiargus$species <- "csemiargus"

picarus_full$species <- "picarus"
pargus_full$species <- "pargus"
csemiargus_full$species <- "csemiargus"

thetas <- rbind(picarus, pargus, csemiargus)
thetas_full <- rbind(picarus_full, pargus_full, csemiargus_full)
colnames(thetas_full) <- c("coords","chr","WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites", "pop", "species")
thetas_full$pi.mean <- thetas_full$tP / thetas_full$nSites
thetas_full$watterson.mean <- thetas_full$tW / thetas_full$nSites

theta_means <- c()

for (s in c("picarus", "pargus", "csemiargus")) {
  species_means <- thetas_full[thetas_full$species == s & thetas_full$nSites > 1000, ] %>%
    group_by(pop) %>%
    summarise_at(c("pi.mean","watterson.mean","Tajima"), mean)
  species_means$species <- s
  theta_means <- rbind(theta_means, species_means)
}
colnames(theta_means) <- c("pop","pi.mean","watterson.mean","tajima.mean","species")

thetas <- merge(thetas, theta_means, by = c('species', "pop"))

thetas$downsample.n <- as.factor(thetas$downsample.n)

thetas <- merge(thetas, sample_sizes, by = c('species', 'pop'))

# Calculate upper and lower bounds of the middle 95% of each downsampled
# estimate distribution

min_pop_size <- 8

cis <- c()

for (s in c("picarus","pargus","csemiargus")) {
  for (n in c(2:(min_pop_size-1))) {
    subset <- df[df$n == n & df$species == s & df$full_n > (min_pop_size - 1), ]
    hz_ci <- quantile((subset$mean - subset$dataset_mean), probs = c(0.025, 0.975))
    subset <- fst[fst$downsample.n == n & fst$species == s & fst$pop1.full.samplesize > (min_pop_size - 1) & fst$pop2.full.samplesize > (min_pop_size - 1), ]
    fst_ci <- quantile((subset$weight.fst.x - subset$weight.fst.y), probs = c(0.025, 0.975))
    subset <- thetas[thetas$downsample.n == n & thetas$species == s & thetas$full_n > (min_pop_size - 1), ]
    pi_ci <- quantile((subset$pi - subset$pi.mean), probs = c(0.025, 0.975))
    wat_ci <- quantile((subset$watterson - subset$watterson.mean), probs = c(0.025, 0.975))
    taj_ci <- quantile((subset$tajima - subset$tajima.mean), probs = c(0.025, 0.975))
    rows <- rbind(
      c(s, n, "hz", hz_ci),
      c(s, n, "fst", fst_ci),
      c(s, n, "pi", pi_ci),
      c(s, n, "watterson", wat_ci),
      c(s, n, "tajima", taj_ci)
    )
    cis <- rbind(cis, data.frame(rows))
  }
}

# Calculate normalized root mean squared error for downsampled estimates

nrmse <- c()

for (s in c("picarus","pargus","csemiargus")) {
  for (n in c(2:(min_pop_size-1))) {
    subset <- df[df$n == n & df$species == s & df$full_n > (min_pop_size - 1), ]
    hz_nrmse <- sqrt(mean((subset$mean - subset$dataset_mean)^2))/sd(subset$mean)
    subset <- fst[fst$downsample.n == n & fst$species == s & fst$pop1.full.samplesize > (min_pop_size - 1) & fst$pop2.full.samplesize > (min_pop_size - 1), ]
    fst_nrmse <- sqrt(mean((subset$weight.fst.x - subset$weight.fst.y)^2))/sd(subset$weight.fst.x)
    subset <- thetas[thetas$downsample.n == n & thetas$species == s & thetas$full_n > (min_pop_size - 1), ]
    pi_nrmse <- sqrt(mean((subset$pi - subset$pi.mean)^2))/sd(subset$pi)
    wat_nrmse <- sqrt(mean((subset$watterson - subset$watterson.mean)^2))/sd(subset$watterson)
    taj_nrmse <- sqrt(mean((subset$tajima - subset$tajima.mean)^2))/sd(subset$tajima)
    rows <- rbind(
      c(s, n, "hz", hz_nrmse),
      c(s, n, "fst", fst_nrmse),
      c(s, n, "pi", pi_nrmse),
      c(s, n, "watterson", wat_nrmse),
      c(s, n, "tajima", taj_nrmse)
    )
    nrmse <- rbind(nrmse, data.frame(rows))
  }
}

colnames(nrmse) <- c("species","downsample.n","metric", "nrmse")

nrmse$downsample.n <- as.numeric(nrmse$downsample.n)
nrmse$nrmse <- as.numeric(nrmse$nrmse)
nrmse$species <- factor(nrmse$species, levels = c("picarus", "pargus", "csemiargus"),
                        labels = c("P. icarus", "P. argus", "C. semiargus"))

ggplot(nrmse, aes(x=downsample.n, y=nrmse, shape = metric, linetype = metric, color = metric, group = metric)) +
  facet_wrap(vars(species)) +
  geom_point() +
  geom_line() +
  xlab("Downsampled Sample Size") +
  ylab("NRMSE") +
  scale_color_grey(start = 0, end = 0.5) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic"),
        axis.title.y = element_text(face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
