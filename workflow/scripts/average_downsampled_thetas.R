sink(file(snakemake@log[[1]], open = "wt"), type = "message")

average_pestpg <- function(pestpg, popname, subsize, rep, minsites) {
  theta <- as.data.frame(read.table(pestpg, header = TRUE, comment.char = ""))
  theta <- theta[!theta$nSites < minsites, ]
  theta$watterson <- theta$tW / theta$nSites
  theta$pi <- theta$tP / theta$nSites
  return(
    cbind(
      popname,
      subsize,
      rep,
      mean(theta$pi),
      mean(theta$watterson),
      mean(theta$Tajima)
    )
  )
}

avg <- average_pestpg(
  snakemake@input[[1]],
  snakemake@wildcards[["population"]],
  snakemake@wildcards[["samplesize"]],
  snakemake@wildcards[["rep"]],
  snakemake@params[["minsites"]]
)

write.table(
  avg,
  file = snakemake@output[[1]],
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)
