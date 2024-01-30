sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(lme4)
library(dplyr)
library(car)
library(MuMIn)

# read in land cover and heterzozygosity data
land_cover <- data.frame(
  read.table(snakemake@input[["cover"]], sep = "\t", header = TRUE)
)
hz <- data.frame(
  read.table(snakemake@input[["hz"]], sep = "\t", header = TRUE)
)
outfile <- snakemake@output[["modelout"]]
dataset <- snakemake@wildcards[["dataset"]]

# define buffer radii around sample sites
buffer_radius <- c(
  500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 15000, 20000
)

# clean up column names to match land_cover table
hz$dataset <- snakemake@wildcards[["dataset"]]
colnames(hz) <- c("sample", "site", "heterozygosity", "dataset")

# combine land use and heterozygosity dataframes
ind_hz_land <- merge(hz, land_cover, by = "site")

results <- data.frame()

# set landscapes of interest
landscapes <- c(
  "arable", "grassland", "forest", "water"
)

# generate models of heterozygosity ~ landuse for each radii
for (r in buffer_radius) {
  for (m in landscapes) {
    if (any(c(ind_hz_land[ind_hz_land$radius == r, m] > 0))) {
      mod <- as.formula(sprintf("heterozygosity ~ %s + (1|site)", m))
      model <- lmer(
        mod,
        data = ind_hz_land[ind_hz_land$radius == r, ]
      )
      summ <- summary(model)
      aicc <- AICc(model)
      tval <- summ$coefficients[6]
      pval <- Anova(model)[[3]][1]
      r2m <- r.squaredGLMM(model)[1]
      r2c <- r.squaredGLMM(model)[2]
      row <- c(dataset, r, m, tval, pval, r2m, r2c, aicc)
      results <- rbind(results, row)
    }
  }
}

colnames(results) <- c("dataset", "radius", "landuse", "t-val", "p-val", "r2m",
                       "r2c", "aicc")

# write results to output table
write.table(
  results,
  file = outfile,
  quote = FALSE,
  append = FALSE,
  sep = "\t",
  row.names = FALSE
)