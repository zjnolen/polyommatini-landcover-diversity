# Requires R 3.6, not the environment in the plotting folder

# Plots all values of K admixture plots, as well as evalAdmix correlation of
# residual plots for supplementary materials

source("https://raw.githubusercontent.com/GenisGE/evalAdmix/2a51aebaca70c9d3b5fb359d6c5c40145c58fce5/visFuns.R")

setwd("~/working/bioinfo/projects/landuse-manuscript")

dev.off()

par(mfrow=c(7,1), mar = c(1,2,1,1))

for (k in c(2:8)) {

  pop <- read.table("results/datasets/landuse-csemiargus/poplists/landuse-csemiargus_all_excl_pca-admix.indiv.list", header = TRUE)

  q <- read.table(paste0("results/datasets/landuse-csemiargus/analyses/ngsadmix/landuse-csemiargus.ilCyaSemi1.1_all_excl_pca-admix_allsites-filts_K",k,".qopt"),stringsAsFactors=T)
  
  popord <- c("Silvakra", "SESkane", "Fastan", "Agunnaryd", "Gotafors", 
              "Oland")
  
  ord <- orderInds(pop = as.vector(pop[,2]), q = q, popord = popord)

  plotAdmix(
    q,
    pop = as.vector(pop[,2]),
    ord=ord,
    main=""
  )
}
# 750x900
dev.off()

for (k in c(1:8)) {
  
  pop <- read.table("results/datasets/landuse-csemiargus/poplists/landuse-csemiargus_all_excl_pca-admix.indiv.list", header = TRUE)
  
  q <- read.table(paste0("results/datasets/landuse-csemiargus/analyses/ngsadmix/landuse-csemiargus.ilCyaSemi1.1_all_excl_pca-admix_allsites-filts_K",k,".qopt"),stringsAsFactors=T)
  
  popord <- c("Silvakra", "SESkane", "Fastan", "Agunnaryd", "Gotafors", 
              "Oland")
  
  ord <- orderInds(pop = as.vector(pop[,2]), q = q, popord = popord)
  
  r<-as.matrix(read.table(paste0("results/datasets/landuse-csemiargus/analyses/ngsadmix/landuse-csemiargus.ilCyaSemi1.1_all_excl_pca-admix_allsites-filts_K",k,".corres")))
  
  png(paste0("results/figures/csemiargus_evaladmix_k",k,".png"), height = 600, width = 800)
  
  plotCorRes(
    cor_mat = r,
    pop = as.vector(pop[,2]),
    ord=ord,
    title=paste0("Residuals K=",k),
    rotatelabpop=45,
    max_z=0.1,
    min_z=-0.1)
  
  dev.off()
  
}

dev.off()

par(mfrow=c(7,1), mar = c(1,2,1,1))

for (k in c(2:8)) {
  
  pop <- read.table("results/datasets/landuse-pargus/poplists/landuse-pargus_all.indiv.list", header = TRUE)
  
  q <- read.table(paste0("results/datasets/landuse-pargus/analyses/ngsadmix/landuse-pargus.ilPleArgu1.3_all_allsites-filts_K",k,".qopt"),stringsAsFactors=T)
  
  popord <- c("Falsterbo", "Drakamollan", "Gotafors", "Aspo", "Branthalla", 
              "Gosslunda", "Jordtorpasen")
  
  ord <- orderInds(pop = as.vector(pop[,2]), q = q, popord = popord)
  
  plotAdmix(
    q,
    pop = as.vector(pop[,2]),
    ord=ord,
    main=""
  )
}

dev.off()

for (k in c(1:8)) {
  
  pop <- read.table("results/datasets/landuse-pargus/poplists/landuse-pargus_all.indiv.list", header = TRUE)
  
  q <- read.table(paste0("results/datasets/landuse-pargus/analyses/ngsadmix/landuse-pargus.ilPleArgu1.3_all_allsites-filts_K",k,".qopt"),stringsAsFactors=T)
  
  popord <- c("Falsterbo", "Drakamollan", "Gotafors", "Aspo", "Branthalla", 
              "Gosslunda", "Jordtorpasen")
  
  ord <- orderInds(pop = as.vector(pop[,2]), q = q, popord = popord)
  
  r<-as.matrix(read.table(paste0("results/datasets/landuse-pargus/analyses/ngsadmix/landuse-pargus.ilPleArgu1.3_all_allsites-filts_K",k,".corres")))
  
  png(paste0("results/figures/pargus_evaladmix_k",k,".png"), height = 600, width = 800)
  
  plotCorRes(
    cor_mat = r,
    pop = as.vector(pop[,2]),
    ord=ord,
    title=paste0("Residuals K=",k),
    rotatelabpop=45,
    max_z=0.1,
    min_z=-0.1)
  
  dev.off()
  
}

dev.off()

par(mfrow=c(7,1), mar = c(1,2,1,1))

for (k in c(2:7)) {
  
  pop <- read.table("results/datasets/landuse-picarus/poplists/landuse-picarus_all_excl_pca-admix.indiv.list", header = TRUE)
  
  q <- read.table(paste0("results/datasets/landuse-picarus/analyses/ngsadmix/landuse-picarus.ilPolIcar1.1_all_excl_pca-admix_allsites-filts_K",k,".qopt"),stringsAsFactors=T)
  
  popord <- c("Falsterbo", "Backakra", "MaxIV", "Tvedora", "Havang", "Ljungby",
              "Blaberget", "Sturko", "Gotafors", "SandbyBorg", "Ismantorp")
  
  ord <- orderInds(pop = as.vector(pop[,2]), q = q, popord = popord)
  
  plotAdmix(
    q,
    pop = as.vector(pop[,2]),
    ord=ord,
    main=""
  )
  
}

dev.off()

for (k in c(1:8)) {
  
  pop <- read.table("results/datasets/landuse-picarus/poplists/landuse-picarus_all_excl_pca-admix.indiv.list", header = TRUE)
  
  q <- read.table(paste0("results/datasets/landuse-picarus/analyses/ngsadmix/landuse-picarus.ilPolIcar1.1_all_excl_pca-admix_allsites-filts_K",k,".qopt"),stringsAsFactors=T)
  
  popord <- c("Falsterbo", "Backakra", "MaxIV", "Tvedora", "Havang", "Ljungby",
              "Blaberget", "Sturko", "Gotafors", "SandbyBorg", "Ismantorp")
  
  ord <- orderInds(pop = as.vector(pop[,2]), q = q, popord = popord)
  
  r<-as.matrix(read.table(paste0("results/datasets/landuse-picarus/analyses/ngsadmix/landuse-picarus.ilPolIcar1.1_all_excl_pca-admix_allsites-filts_K",k,".corres")))

  png(paste0("results/figures/picarus_evaladmix_k",k,".png"), height = 600, width = 800)
  
  plotCorRes(
    cor_mat = r,
    pop = as.vector(pop[,2]),
    ord=ord,
    title=paste0("Residuals K=",k),
    rotatelabpop=45,
    adjlab=0.1,
    max_z=0.1,
    min_z=-0.1)
  
  dev.off()
  
}
 