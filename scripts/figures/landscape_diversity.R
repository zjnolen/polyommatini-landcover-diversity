library(ggplot2)
library(dplyr)

setwd("~/working/bioinfo/projects/landuse-manuscript")

picarus <- read.table("results/datasets/landuse-picarus/analyses/landscape_models/landuse-picarus.ilPolIcar1.1_all_allsites-filts_heterozygosity-landuse.tsv", header = TRUE)
pargus <- read.table("results/datasets/landuse-pargus/analyses/landscape_models/landuse-pargus.ilPleArgu1.3_all_allsites-filts_heterozygosity-landuse.tsv", header = TRUE)
csemiargus <- read.table("results/datasets/landuse-csemiargus/analyses/landscape_models/landuse-csemiargus.ilCyaSemi1.1_all_allsites-filts_heterozygosity-landuse.tsv", header = TRUE)

picarus$species <- "picarus"
pargus$species <- "pargus"
csemiargus$species <- "csemiargus"

results <- rbind(picarus, pargus, csemiargus)

results$radius <- factor(
  results$radius,
  levels = as.character(c(100, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000,
                          8000, 10000, 15000, 20000))
)
results$landuse <- factor(
  results$landuse,
  levels = rev(c(
    "grassland", "arable", "forest", "grassland_forest",
    "arable_grassland", "arable_forest", "water"
  )),
  labels = rev(c(
    "Grassland", "Arable", "Forest", "Grassland & Forest",
    "Arable & Grassland", "Arable & Forest", "Water"
  ))
)
results$species <- factor(
  results$species,
  levels = c("picarus", "pargus", "csemiargus"),
  labels = c("P. icarus", "P. argus", "C. semiargus")
)
results$sig <- as.factor(with(results, ifelse(p.val < 0.05, 1, 0)))
results$neg <- as.numeric(with(results, ifelse(t.val < 0, -1, 1)))
results$r2sign <- as.numeric(results$r2m * results$neg)
results <- results %>%
  group_by(radius, species) %>%
  mutate(deltaAIC = aic - min(aic))

ggplot(
  data = results,
  aes(x = radius, y = landuse, col = neg, size = r2m, alpha = sig,
    label = round(deltaAIC, 2))) +
  geom_point() +
  geom_point(stroke = ifelse(results$deltaAIC < 2, 1, 0), shape = 1, col = ifelse(results$deltaAIC == 0, "black", "grey40"), alpha = 1) +
  # geom_text(color = "black", alpha = 1, size = 2, vjust = -2, fontface = ifelse(results$deltaAIC < 2, 2, 1)) +
  facet_grid(rows = vars(species), switch = "y") +
  scale_color_gradientn(colors = rev(c("#15607A", "#A63716"))) +
  scale_alpha_manual(values = c(0.5, 0.9)) +
  scale_y_discrete(position = "right") +
  # geom_text(color = "black", alpha = 1, size = 2, vjust = -2, fontface = ifelse(results$deltaAIC < 2, 2, 1)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(angle = 0),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    strip.text = element_text(face = "italic")
  ) +
  guides(
    color = "none",
    alpha = "none"
  ) +
  labs(
    size = "Marginal R2"
  ) +
  xlab("Radius [m]")

ggsave("results/figures/landscape_diversity.svg", width = 5, height = 5)
