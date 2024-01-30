# This script plots the results for the models of genetic differentiation ~ distance + land use

library(ggplot2)
library(stringr)
library(dplyr)

setwd("~/working/bioinfo/projects/landuse-manuscript")

# Read in results from models for all three species and compile them into a
# results object

picarus <- read.table("results/datasets/landuse-picarus/analyses/landscape_models/landuse-picarus.ilPolIcar1.1_all_allsites-filts_fst-landuse.tsv", header = TRUE)
pargus <- read.table("results/datasets/landuse-pargus/analyses/landscape_models/landuse-pargus.ilPleArgu1.3_all_allsites-filts_fst-landuse.tsv", header = TRUE)
csemiargus <- read.table("results/datasets/landuse-csemiargus/analyses/landscape_models/landuse-csemiargus.ilCyaSemi1.1_all_allsites-filts_fst-landuse.tsv", header = TRUE)

picarus$species <- "picarus"
pargus$species <- "pargus"
csemiargus$species <- "csemiargus"

results <- rbind(picarus, pargus, csemiargus)

results$landuse <- factor(
  results$landuse,
  levels = c(
    "distance", "grassland", "arable", "forest", "water"
  ),
  labels = c(
    "Distance", "+ Grassland", "+ Arable", "+ Forest", "+ Water"
  )
)

results$species <- factor(
  results$species,
  levels = c("picarus", "pargus", "csemiargus"),
  labels = c("P. icarus", "P. argus", "C. semiargus")
)

# Generate helper columns that have direction of relationship and delta AIC
results$neg <- as.numeric(with(results, ifelse(t.val < 0, -1, 1)))
results <- results %>%
  group_by(dataset) %>%
  mutate(deltaAIC = aicc - min(aicc))

# Plot land use relationships with genetic differentaiton per species

ggplot(
  data = results,
  aes(
    y = 1,
    x = landuse,
    col = neg,
    size = r2m,
    label = round(deltaAIC, 2)
  )
) +
  facet_grid(rows = vars(species), switch = "both") +
  geom_point(alpha = ifelse(results$deltaAIC == 0, 0.9, 0.5)) +
  geom_point(
    stroke = ifelse(results$deltaAIC < 2, 1, 0),
    shape = 1,
    col = ifelse(results$deltaAIC == 0, "black", "grey40"),
    alpha = 1
  ) +
  geom_text(
    color = "black",
    alpha = 1,
    size = 2.5,
    vjust = -2,
    fontface = ifelse(results$deltaAIC < 2, 2, 1)
  ) +
  scale_color_gradientn(colors = (c("#15607A", "#A63716"))) +
  theme_bw() +
  guides(
    color = "none",
    alpha = "none"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.8),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text = element_text(face = "italic"),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
labs(
  size = "Marginal R2"
)

# This portion plots the isolation by distance plots of the three species

library(lmerMultiMember)
library(ggeffects)
library(MuMIn)
library(ggplot2)

# Read in FST estimates and compile them into an object fstall

picarus <- read.table("results/datasets/landuse-picarus/analyses/fst/landuse-picarus.ilPolIcar1.1_poppairs_allsites-filts.fst.global.tsv", header = TRUE)
pargus <- read.table("results/datasets/landuse-pargus/analyses/fst/landuse-pargus.ilPleArgu1.3_poppairs_allsites-filts.fst.global.tsv", header = TRUE)
csemiargus <- read.table("results/datasets/landuse-csemiargus/analyses/fst/landuse-csemiargus.ilCyaSemi1.1_poppairs_allsites-filts.fst.global.tsv", header = TRUE)


picarus$species <- "picarus"
pargus$species <- "pargus"
csemiargus$species <- "csemiargus"

fstall <- rbind(picarus, pargus, csemiargus)


# Set FST values < 0 to 0

fstall$weight.fst <- ifelse(fstall$weight.fst < 0, 0, fstall$weight.fst)

# Read in table with distances betweeen sampling sites for all species and
# compile into landscapes object
picaruslandscapes <- data.frame(read.table("results/landscape/landuse-picarus_between_site_cover.tsv", sep = "\t", header = TRUE))
parguslandscapes <- data.frame(read.table("results/landscape/landuse-pargus_between_site_cover.tsv", sep = "\t", header = TRUE))
csemiarguslandscapes <- data.frame(read.table("results/landscape/landuse-csemiargus_between_site_cover.tsv", sep = "\t", header = TRUE))

landscapes <- rbind(picaruslandscapes, parguslandscapes, csemiarguslandscapes)

# Generate population combination variable to merge distances and fst tables
landscapes$combname <- ifelse(landscapes$pop1 < landscapes$pop2, paste0(landscapes$pop1," - ",landscapes$pop2),paste0(landscapes$pop2," - ",landscapes$pop1))

fstall$combname <- ifelse(fstall$pop1 < fstall$pop2, paste0(fstall$pop1," - ",fstall$pop2),paste0(fstall$pop2," - ",fstall$pop1))


# Merge fst and distance tables by population pair
fstlands <- merge(fstall,landscapes, by="combname")


# Generate a predictive model of isolation by distance to plot that matches
# the one tested in the manuscript (i.e. with a population pair membership
# random effect)
pred <- c()
stats <- c()
fst <- fstlands

species <- c("picarus", "pargus", "csemiargus")

for (s in species) {
  Wt <- weights_from_columns(fst[fst$species == s, c("pop1.x","pop2.x")])
  model <- lmer(weight.fst ~ scale(distance) + (1|pop), memberships = list(pop = Wt), data = fst[fst$species == s,])
  pr <- ggpredict(model,terms=c("distance"))
  pr$species <- s
  slope <- (max(pr$predicted)-min(pr$predicted)) / (max(pr$x)/1000-min(pr$x)/1000)
  r2m <- r.squaredGLMM(model)[1]
  stats <- rbind(stats, c(s, slope, r2m))
  pred <- rbind(pred,pr)
}

# Calculate stats to write onto plots - slope and R2m
stats <- data.frame(stats)
colnames(stats) <- c("species", "slope", "r2m")
pred$x <- pred$x / 1000
fst$distance <- fst$distance/1000

# Set variable types for nice plotting
pred$species <- factor(pred$species, levels = c("picarus","pargus","csemiargus"), labels = c("P. icarus", "P. argus", "C. semiargus"))
fst$species <- factor(fst$species, levels = c("picarus","pargus","csemiargus"), labels = c("P. icarus", "P. argus", "C. semiargus"))
stats$species <- factor(stats$species, levels = c("picarus","pargus","csemiargus"), labels = c("P. icarus", "P. argus", "C. semiargus"))
stats$slope <- as.numeric(stats$slope)
stats$r2m <- as.numeric(stats$r2m)
stats$r2m <- round(stats$r2m, digits = 3)

# Plot isolation by distance for the three species together
IsoByDist <- ggplot(pred,aes(x,predicted))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "grey60",
              linetype=0)+
  geom_point(data=fst, aes(y=weight.fst,x=distance))+
  facet_grid(cols = vars(species)) +
  geom_line()+
  xlab("Distance (km)") +
  ylab(expression("Weighted pairwise F"["ST"])) +
  scale_x_continuous(limits = c(0,300), labels = ~ format(.x, scientific = FALSE)) +
  geom_label(data = stats, aes(x = -Inf, y = Inf, label = paste("R2m",r2m,"\nSlope", format(slope, scientific = TRUE, digits = 3)), hjust = 0, vjust = 1)) +
  theme_bw() +
  ylim(c(0,0.15)) +
  theme(
    strip.text = element_text(face = "italic")
    
  )

IsoByDist

