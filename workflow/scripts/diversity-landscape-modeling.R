sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(nlme)
library(MuMIn)
library(ggplot2)


land_cover <- data.frame(
  read.table(snakemake@input[["cover"]], sep = "\t", header = TRUE)
)
hz <- data.frame(
  read.table(snakemake@input[["hz"]], sep = "\t", header = TRUE)
)
outfile <- snakemake@output[["modelout"]]
buffer_radius <- c(
  100, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 15000, 20000
)


hz$dataset <- snakemake@wildcards[["dataset"]]
colnames(hz) <- c("sample", "site", "heterozygosity", "dataset")
land_cover$grassland_forest <- land_cover$grassland + land_cover$forest
land_cover$arable_grassland <- land_cover$arable + land_cover$grassland
land_cover$arable_forest <- land_cover$arable + land_cover$forest

ind_hz_land <- merge(hz, land_cover, by = "site")

results <- data.frame()

landscapes <- c(
  "arable", "grassland", "forest", "water", "arable_grassland",
  "arable_forest", "grassland_forest"
)

for (r in buffer_radius) {
  for (m in landscapes) {
    if (any(c(ind_hz_land[ind_hz_land$radius == r, m] > 0))) {
      mod <- as.formula(sprintf("heterozygosity ~ %s", m))
      model <- lme(
        mod,
        random = ~ 1 | site,
        data = ind_hz_land[ind_hz_land$radius == r, ]
      )
      summ <- summary(model)
      aic <- summ$AIC
      r2 <- r.squaredGLMM(model)[1]
      row <- c(
        snakemake@wildcards[["dataset"]],
        r,
        m,
        summ$tTable[8],
        summ$tTable[10],
        r2,
        aic
      )
      results <- rbind(results, row)
    }
  }
}

colnames(results) <- c(
  "dataset", "radius", "landuse", "t-val", "p-val", "R2m", "aic"
)

write.table(
  results,
  file = outfile,
  quote = FALSE,
  append = FALSE,
  sep = "\t",
  row.names = FALSE
)