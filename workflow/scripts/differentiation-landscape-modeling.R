library(lmerMultiMember)
library(car)
library(MuMIn)

# read in fst and land cover between populations
fst <- data.frame(
  read.table(snakemake@input[["fst"]], sep = "\t", header = TRUE)
)
land_cover <- data.frame(
  read.table(snakemake@input[["cover"]], sep = "\t", header = TRUE)
)
outfile <- snakemake@output[["modelout"]]
dataset <- snakemake@wildcards[["dataset"]]

# generate population combo names
land_cover$combname <- ifelse(
  land_cover$pop1 < land_cover$pop2,
  paste0(land_cover$pop1, "-", land_cover$pop2),
  paste0(land_cover$pop2, "-", land_cover$pop1)
)
fst$combname <- paste0(fst$pop1,"-",fst$pop2)


# merge fst and land cover matrices
fstland <- merge(fst, land_cover, by = "combname")

# define landscapes to test
landscapes <- c(
  "arable", "grassland", "forest", "water"
)

# create membership matrix for multi membership random effect
wt <- weights_from_columns(fstland[, c("pop1.x", "pop2.x")])

# generate output dataframe
results <- data.frame()

# model fst ~ distance + landscape (as well as just distance alone)
for (m in c("distance", landscapes)) {
  if (any(c(fstland[, m] > 0))) {
    if (m == "distance") {
      mod <- as.formula("weight.fst ~ scale(distance) + (1|pop)")
    } else {
      mod <- as.formula(
        sprintf("weight.fst ~ scale(distance) + scale(%s) + (1|pop)", m)
      )
    }
    model <- lmer(
      mod,
      memberships = list(pop = wt),
      data = fstland
    )
    summ <- summary(model)
    aicc <- AICc(model)
    r2m <- r.squaredGLMM(model)[1]
    r2c <- r.squaredGLMM(model)[2]
    if (m == "distance") {
      pval <- Anova(model)[[3]][1]
      row <- c(dataset, m, summ$coefficients[6], pval, r2m, r2c, aicc)
    } else {
      pval <- Anova(model)[[3]][2]
      row <- c(dataset, m, summ$coefficients[9], pval, r2m, r2c, aicc)
    }
    results <- rbind(results, row)
  }
}

colnames(results) <- c("dataset", "landuse", "t-val", "p-val", "r2m", "r2c", "aicc")

write.table(
  results,
  file = outfile,
  quote = FALSE,
  append = FALSE,
  sep = "\t",
  row.names = FALSE
)
