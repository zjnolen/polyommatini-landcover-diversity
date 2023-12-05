library(pophelper)

setwd("~/working/bioinfo/projects/landuse-manuscript")

# Plot admix

# Set species to run for

#species <- "picarus"
#species <- "pargus"
species <- "csemiargus"

if (species == "picarus") {
  qopt <- readQ("results/datasets/landuse-picarus/analyses/ngsadmix/landuse-picarus.ilPolIcar1.1_all_excl_pca-admix_allsites-filts_K2.qopt", filetype="basic")
  pop <- read.table("results/datasets/landuse-picarus/poplists/landuse-picarus_all_excl_pca-admix.indiv.list", header = TRUE)
  popord <- c("Falsterbo", "Backakra", "MaxIV", "Tvedora", "Havang", "Ljungby",
              "Blaberget", "Sturko", "Gotafors", "SandbyBorg", "Ismantorp")
  cols <- c("#254f6e","#fecd03")
  sorting <- "Cluster2"
} else if (species == "pargus") {
  qopt <- readQ("results/datasets/landuse-pargus/analyses/ngsadmix/landuse-pargus.ilPleArgu1.3_all_allsites-filts_K3.qopt", filetype="basic")
  pop <- read.table("results/datasets/landuse-pargus/poplists/landuse-pargus_all.indiv.list", header = TRUE)
  popord <- c("Falsterbo", "Drakamollan", "Gotafors", "Aspo", "Branthalla", 
              "Gosslunda", "Jordtorpasen")
  cols <- c("#98AE6E","#fecd03","#254f6e")
  sorting <- "all"
} else if (species == "csemiargus") {
  qopt <- readQ("results/datasets/landuse-csemiargus/analyses/ngsadmix/landuse-csemiargus.ilCyaSemi1.1_all_excl_pca-admix_allsites-filts_K6.qopt", filetype="basic")
  pop <- read.table("results/datasets/landuse-csemiargus/poplists/landuse-csemiargus_all_excl_pca-admix.indiv.list", header = TRUE)
  popord <- c("Silvakra", "SESkane", "Fastan", "Agunnaryd", "Gotafors", 
              "Oland")
  qopt[[1]] <- qopt[[1]][,rev(c(2,5,3,1,6,4))]
  qopt <- as.qlist(qopt)
  colnames(qopt[[1]]) <- c("Cluster1","Cluster2","Cluster3","Cluster4","Cluster5","Cluster6")
  cols <- rev(c("#244F6E","#658399","#A7B8C5","#C43C55","#8A2A3B","#fecd03"))
  sorting <- "all"
}


row.names(qopt[[1]]) <- pop$sample

pops <- pop[2]

ord <- order(pops$population)

pops$population <- pops[ord, "population"]

qopt[[1]] <- qopt[[1]][ord,]

plotQ(qopt, grplab=pops, subsetgrp = popord, sortind = sorting, 
      clustercol=cols, barbordercolour="white", width = 10, height = 5, dpi = 600, imgtype = "png", exportpath=getwd())