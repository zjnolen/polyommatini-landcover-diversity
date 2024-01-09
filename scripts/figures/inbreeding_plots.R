library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpmisc)
library(lme4)
library(car)
library(MuMIn)

# Script for plotting various plots surrounding inbreeding coefficients
# estimated for ROH. Very much specific to this specific manuscript, but should
# be relatively easily converted to run with other bcftools roh outputs. First 
# blocks of code should be run carefully, as they select what species/input
# to plot.

# Basic inputs

# Run this block first, commenting out the species you want to plot.

setwd("~/working/bioinfo/projects/landuse-manuscript")
options(scipen=999)

s <- "csemiargus"
r <- "ilCyaSemi1.1"
cml <- 1150
autol <- 416794061
popord <- c("Silvakra", "SESkane", "Fastan", "Agunnaryd", "Gotafors", "Oland")
poplab <- c("Silvåkra", "SE Skåne", "Fästan", "Agunnaryd", "Götafors", "Öland")
popcol <- c("#244f6e", "#658399", "#a7b8c5", "#c43c55", "#8a2a3b", "#ffcd00")

# s <- "pargus"
# r <- "ilPleArgu1.3"
# cml <- 1100
# autol <- 363250311
# popord <- c("Falsterbo", "Drakamollan", "Gotafors", "Aspo", "Branthalla",
#             "Gosslunda", "Jordtorpasen")
# poplab <- c("Falsterbo", "Drakamöllan", "Götafors", "Aspö", "Branthalla",
#             "Gosslunda", "Jordtorpåsen")

# s <- "picarus"
# r <- "ilPolIcar1.1"
# cml <- 1100
# autol <- 472204672
# popord <- c("Falsterbo", "Backakra", "MaxIV", "Tvedora", "Havang", "Ljungby",
#             "Blaberget", "Sturko", "Gotafors", "SandbyBorg", "Ismantorp")

inds <- read.table(
  paste0(
    "results/datasets/landuse-",
    s,
    "/poplists/landuse-",
    s,
    "_all.indiv.list"
  ),
  header = TRUE
)

poplist <- unique(inds$population)

# Set up autosomal size from VCF (total length of autosomes between first and
# last SNP on each chromosome) and read in ROHs

## For calls - Run this block if you want to plot call based results

autosize <- c()

for (i in poplist) {
  pos <- read.table(paste0(
    "results/datasets/landuse-",
    s,
    "/vcfs/landuse-",
    s,
    ".",
    r,
    "_",
    i,
    "_allsites-filts.filter.pos"), header = FALSE)
  maxs <- aggregate(V2 ~ V1, data = pos, max)
  mins <- aggregate(V2 ~ V1, data = pos, min)
  autolen <- sum(maxs$V2 - (mins$V2+1))
  autosize <- data.frame(rbind(autosize, c(i, autolen)))
}

colnames(autosize) <- c("population", "autosize")

rohs <- read.table(
  paste0(
    "results/datasets/landuse-",
    s,
    "/analyses/roh/bcftools/landuse-",
    s,
    ".",
    r,
    "_all_allsites-filts.regs.calls.roh"
  ),
  header = FALSE
)

rohs <- rohs[rohs$V8 >= 85,]
rohs <- rohs[rohs$V6 >= 100000,]

## For likes - Run this block if you want to plot likelihood based results

autosize <- c()

for (i in poplist) {
  pos <- read.table(paste0(
    "results/datasets/landuse-",
    s,
    "/bcfs/landuse-",
    s,
    ".",
    r,
    "_",
    i,
    "_allsites-filts.GLonly.pos"), header = FALSE)
  maxs <- aggregate(V2 ~ V1, data = pos, max)
  mins <- aggregate(V2 ~ V1, data = pos, min)
  autolen <- sum(maxs$V2 - (mins$V2+1))
  autosize <- data.frame(rbind(autosize, c(i, autolen)))
}

colnames(autosize) <- c("population", "autosize")

rohs <- read.table(
  paste0(
    "results/datasets/landuse-",
    s,
    "/analyses/roh/bcftools/landuse-",
    s,
    ".",
    r,
    "_all_allsites-filts.regs.GLonly.roh"
  ),
  header = FALSE
)



rohs <- rohs[rohs$V8 >= 85,]
rohs <- rohs[rohs$V6 >= 100000,]


# Remaining blocks are for plotting, run for each plot you want

## Bar plot of individual Froh

norun <- c()
norun$type <- rep("RG", nrow(inds))
norun$sample <- inds$sample
norun$chr <- rep(0, nrow(inds))
norun$start <- rep(0, nrow(inds))
norun$end <- rep(0, nrow(inds))
norun$length <- rep(0, nrow(inds))
norun$inform_sites <- rep(0, nrow(inds))
norun$phred <- rep(0, nrow(inds))
norun <- data.frame(norun)

## This can be adjusted to bin into whatever generation categories you want
generations <- c(1,10,20,40,80,100,100000)

bins <- (50/generations)/(cml/autol)

frohs <- c()
ranges <- c()

for (i in seq(length(bins))) {
  if (i <= length(bins)-1) {
    runs <- rohs[rohs$V6 < bins[i],]
    runs <- runs[runs$V6 >= bins[i+1],]
    colnames(runs) <- c("type", "sample", "chr", "start", "end", "length", "inform_sites", "phred")
    runs <- rbind(runs, norun)
    runs <- runs %>% group_by(sample) %>%
      summarize(length = sum(length))
    runs <- merge(runs, inds, by = "sample")
    runs <- merge(runs, autosize, by = "population")
    runs$froh <- runs$length / as.numeric(runs$autosize)
    runs$range <- paste0(generations[i],"-",generations[i+1])
    ranges <- c(ranges, paste0(generations[i],"-",generations[i+1]))
    frohs <- rbind(frohs, runs)
  }
}

range_labels <- c(ranges[1:5],">100")
frohs$range <- as.factor(frohs$range)
frohs$range <- factor(frohs$range, levels = ranges, labels = range_labels)

indfroh <- frohs %>% group_by(sample) %>%
  summarize(froh = sum(froh))

frohs$sample <- factor(frohs$sample, levels = indfroh$sample[order(indfroh$froh)])

frohs$population <- factor(frohs$population, levels = popord, labels = poplab)

ggplot(data = frohs, aes(x = sample, y = froh, fill = range)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(cols = vars(population), scales = "free_x") +
  theme_bw() +
  labs(y = expression(F[ROH]), fill = "Generations in past") +
  scale_fill_manual(values = brewer.pal(n=9, "YlGn")[4:9]) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )




## NROH ~ SumROH

rohsumlen <- aggregate(V6 ~ V2, data = rohs, FUN = function(x) c(sum = sum(x), count = length(x), med = median(x)))

rohsumlen <- do.call(data.frame, rohsumlen)

colnames(rohsumlen) <- c("sample", "rohsum", "rohcount", "rohmed")

rohsumlen <- merge(rohsumlen, inds, by = "sample")

popmeans <- rohsumlen %>%
  group_by(population) %>%
  summarize(rohcount = mean_cl_boot(rohcount),
            rohsum = mean_cl_boot(rohsum),
            rohmed = mean_cl_boot(rohmed))

rohsumlen$population <- factor(rohsumlen$population, level = popord)

ggplot(rohsumlen, aes(x = rohsum/1000, y = rohcount, color = population)) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linewidth = 0.5) +
  geom_point(alpha = 0.6) +
  geom_errorbarh(data = popmeans, aes(x = rohsum$y/1000, y = rohcount$y, xmin=rohsum$ymin/1000, xmax=rohsum$ymax/1000), linewidth = 0.6, height = 10) + 
  geom_errorbar(data = popmeans, aes(x = rohsum$y/1000, y = rohcount$y, ymin=rohcount$ymin, ymax=rohcount$ymax), linewidth = 0.6, width = 3000) +
  geom_point(data = popmeans, aes(x = rohsum$y/1000, y = rohcount$y), size = 3) +
  scale_color_manual(values = popcol) +
  ylab("Total number of runs of homozygosity") +
  xlab("Total length of runs of homozygosity (kbp)") +
  theme_classic() +
  guides(color = "none")




## ROH adjusted Heterozygosity

norun <- c()
norun$V1 <- rep("RG", nrow(inds))
norun$V2 <- inds$sample
norun$V3 <- rep(0, nrow(inds))
norun$V4 <- rep(0, nrow(inds))
norun$V5 <- rep(0, nrow(inds))
norun$V6 <- rep(0, nrow(inds))
norun$V7 <- rep(0, nrow(inds))
norun$V8 <- rep(0, nrow(inds))
norun <- data.frame(norun)

hz <- read.table(
  paste0(
    "results/datasets/landuse-",
    s,
    "/analyses/heterozygosity/landuse-",
    s,
    ".",
    r,
    "_all_allsites-filts_heterozygosity.tsv"
    ),
  header = TRUE)

f <- aggregate(V6 ~ V2, data = rbind(rohs, norun), sum)
colnames(f) <- c("sample", "cumroh")
f <- merge(f, inds, by = "sample")
f <- merge(f, autosize, by = "population")
df <- merge(f, hz, by = c("sample"))
df$froh <- df$cumroh/as.numeric(df$autosize)
df2 <- df
df2$heterozygosity <- df2$heterozygosity / ( 1 - df2$froh)
  
df$type <- "all"
df2$type <- "noroh"
  
df <- rbind(df, df2)

df$type <- factor(df$type, levels = c("all","noroh"), labels = c("Heterozygosity including ROHs", "Heterozygosity excluding ROHs"))

df$population <- factor(df$population, levels = popord, labels = poplab)

ggplot(df, aes(x = type, y = heterozygosity, fill = population, linetype = type, alpha = type)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(cols = vars(population)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  scale_alpha_manual(values = c(1, 0.4)) +
  scale_fill_manual(values = popcol) +
  guides(fill = "none") +
  ylab("Heterozygous sites per 1000bp") +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(linetype = NULL, alpha = NULL)


# Chromosomal roh cover - likes vs. calls

callroh <- read.table(
  paste0(
    "results/datasets/landuse-",
    s,
    "/analyses/roh/bcftools/landuse-",
    s,
    ".",
    r,
    "_all_allsites-filts.regs.calls.roh"
  ),
  header = FALSE
)

likeroh <- read.table(
  paste0(
    "results/datasets/landuse-",
    s,
    "/analyses/roh/bcftools/landuse-",
    s,
    ".",
    r,
    "_all_allsites-filts.regs.GLonly.roh"
  ),
  header = FALSE
)

likeroh$V9 <- "likes"
callroh$V9 <- "calls"

allroh <- rbind(likeroh,callroh)

allroh <- allroh[allroh$V8 >= 30,]
allroh <- allroh[allroh$V2 %in% c("CYSE0078","MZLU107500","MZLU107463", "CYSE0198", "CYSE0199"),]

ggplot(allroh[allroh$V3 == 2,], aes(y = V2, xmin = V4, xmax = V5, color = V9)) +
  geom_linerange(size = 5, position = position_dodge(0.8)) +
  theme_classic()

allroh <- allroh[allroh$V8 >= 85,]

ggplot(allroh[allroh$V3 == 2,], aes(y = V2, xmin = V4, xmax = V5, color = V9)) +
  geom_linerange(size = 5, position = position_dodge(0.8)) +
  theme_classic()
