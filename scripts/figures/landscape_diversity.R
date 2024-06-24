library(ggplot2)
library(reshape2)
library(dplyr)

setwd("~/working/bioinfo/projects/landuse-manuscript")

# Read in results from heterozygosity ~ land use models for all species

picarus <- read.table("results/datasets/landuse-picarus/analyses/landscape_models/landuse-picarus.ilPolIcar1.1_all_allsites-filts_heterozygosity-landuse.tsv", header = TRUE)
pargus <- read.table("results/datasets/landuse-pargus/analyses/landscape_models/landuse-pargus.ilPleArgu1.3_all_allsites-filts_heterozygosity-landuse.tsv", header = TRUE)
csemiargus <- read.table("results/datasets/landuse-csemiargus/analyses/landscape_models/landuse-csemiargus.ilCyaSemi1.1_all_allsites-filts_heterozygosity-landuse.tsv", header = TRUE)

picarus$species <- "picarus"
pargus$species <- "pargus"
csemiargus$species <- "csemiargus"

# Create a results object with results from all species. Set all columns to
# appropriate types

results <- rbind(picarus, pargus, csemiargus)

results$radius <- factor(
  results$radius,
  levels = as.character(c(500, 1000, 1500, 2000, 3000, 4000, 5000, 6000,
                          8000, 10000, 15000, 20000))
)
results$landuse <- factor(
  results$landuse,
  levels = rev(c(
    "grassland", "arable", "forest", "water"
  )),
  labels = rev(c(
    "Grassland", "Arable", "Forest", "Water"
  ))
)
results$species <- factor(
  results$species,
  levels = c("picarus", "pargus", "csemiargus"),
  labels = c("P. icarus", "P. argus", "C. semiargus")
)

# Set helper columns for direction of relationship and delta AIC
results$neg <- as.numeric(with(results, ifelse(t.val < 0, -1, 1)))
results <- results %>%
  group_by(radius, species) %>%
  mutate(deltaAIC = aicc - min(aicc))

# Plot results from landscape models

ggplot(
  data = results,
  aes(x = radius, y = landuse, col = neg, size = r2m,
    label = round(deltaAIC, 2))) +
  geom_point(alpha = ifelse(results$deltaAIC == 0, 0.9, 0.5)) +
  geom_point(stroke = ifelse(results$deltaAIC < 2, 1.25, 0), shape = 1, col = ifelse(results$deltaAIC == 0, "black", "grey40"), alpha = 1) +
  facet_grid(rows = vars(species), switch = "y") +
  scale_color_gradientn(colors = rev(c("#15607A", "#A63716"))) +
  scale_y_discrete(position = "right") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    strip.text = element_text(face = "italic"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(
    color = "none",
    alpha = "none"
  ) +
  labs(
    size = "Marginal R2"
  ) +
  xlab("Radius [m]")

ggsave("results/figures/landscape_diversity.svg", width = 4.5, height = 5)

# Assess correlations between land use variables at multiple scales


landuse <- read.table("results/landscape/around_site_cover_tuva.tsv", header = TRUE)

corr <- c()

for (s in c("picarus", "pargus", "csemiargus")) {
  
  samples <- read.table(paste0("config/samples_",s,".tsv"), header = TRUE)
  pops <- unique(samples$population)
  
  species_landuse <- landuse[landuse$site %in% pops, ]
  
  correlations <- c()
  
  for (r in unique(species_landuse$radius)) {
    radius <- species_landuse[species_landuse$radius == r, ]
    grasstuva <- cor(radius$grassland, radius$sng, method = "pearson")
    grassarable <- cor(radius$grassland, radius$arable, method = "pearson")
    grassforest <- cor(radius$grassland, radius$forest, method = "pearson")
    arableforest <- cor(radius$arable, radius$forest, method = "pearson")
    watergrass <- cor(radius$water, radius$grass, method = "pearson")
    waterforest <- cor(radius$water, radius$forest, method = "pearson")
    waterarable <- cor(radius$water, radius$arable, method = "pearson")
    row <- data.frame(r, grasstuva, grassarable, grassforest, arableforest, watergrass, waterforest, waterarable)
    correlations <- rbind(row, correlations)
  }
  
  correlations <- melt(correlations, id = "r")
  correlations$variable <- factor(correlations$variable, levels = c("grasstuva", "grassarable","grassforest","arableforest","watergrass","waterforest","waterarable"), labels = c("Grass:TUVA", "Grass:Arable", "Grass:Forest", "Arable:Forest", "Water:Grass","Water:Forest","Water:Arable"))
  
  correlations$species <- s
  corr <- rbind(corr, correlations)
}


corr$species <- factor(
  corr$species,
  levels = c("picarus", "pargus", "csemiargus"),
  labels = c("P. icarus", "P. argus", "C. semiargus")
)

# Plot correlations

ggplot(corr, aes(x = r, y = value, shape = variable, linetype = variable)) +
  facet_wrap( ~ species) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0) +
  ylab("Pearson correlation, r") +
  xlab("Radius (m)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "italic", size = 13),
    legend.title = element_blank()
  )
  
  
ggsave("results/figures/landscape_correlation.svg", width = 12, height = 5)
  