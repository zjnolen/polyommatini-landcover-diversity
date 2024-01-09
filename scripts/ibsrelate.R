setwd("~/working/bioinfo/projects/landuse-manuscript")

source("scripts/read_IBS.R")

# P. icarus

df <- c()

pops <- c("Backakra","Blaberget","Falsterbo","Gotafors","Havang", "Ismantorp",
          "Ljungby", "MaxIV","SandbyBorg","Sturko","Tvedora")

for (n in c(1:17)) {
  for (p in pops) {
  kins <- do_derived_stats(
    read_ibspair_model0(
      paste0(
        "results/datasets/landuse-picarus/analyses/kinship/ibsrelate/landuse-picarus.ilPolIcar1.1_",
        p,
        "_chunk",
        n,
        "_allsites-filts.ibspair"
      )
    )
  )
  kins$chunk <- n
  kins$pop <- p
  df <- rbind(df, kins)
}
}

df_all <- df

df <- aggregate(cbind(A,B,C,D,E,F,G,H,I) ~ pair + pop, data = df_all, FUN = sum)
df['HETHET'] = df['E']
df['IBS0'] = df['C'] + df['G']
df['IBS1'] = df['B'] + df['D'] + df['F'] + df['H']
df['IBS2'] = df['A'] + df['E'] + df['I']
df['R0'] = df['IBS0'] / df['HETHET']
df['R1'] = df['HETHET'] / (df['IBS0'] +  df['IBS1'])
df['Kin'] = (df['HETHET'] - 2*(df['IBS0'])) / (df['IBS1'] + 2*df['HETHET'])
df$misschunk <- "none"

for (n in c(1:16)) {
  df2 <- aggregate(cbind(A,B,C,D,E,F,G,H,I) ~ pair + pop, data = df_all[!df_all$chunk == n,], FUN = sum)
  df2['HETHET'] = df2['E']
  df2['IBS0'] = df2['C'] + df2['G']
  df2['IBS1'] = df2['B'] + df2['D'] + df2['F'] + df2['H']
  df2['IBS2'] = df2['A'] + df2['E'] + df2['I']
  df2['R0'] = df2['IBS0'] / df2['HETHET']
  df2['R1'] = df2['HETHET'] / (df2['IBS0'] +  df2['IBS1'])
  df2['Kin'] = (df2['HETHET'] - 2*(df2['IBS0'])) / (df2['IBS1'] + 2*df2['HETHET'])
  df2$misschunk <- n
  df <- rbind(df, df2)
}

# P. argus

df <- c()

pops <- c("all")

for (n in c(1:17)) {
  for (p in pops) {
    kins <- do_derived_stats(
      read_ibspair_model0(
        paste0(
          "results/datasets/landuse-pargus/analyses/kinship/ibsrelate/landuse-pargus.ilPleArgu1.3_",
          p,
          "_chunk",
          n,
          "_allsites-filts.ibspair"
        )
      )
    )
    kins$chunk <- n
    kins$pop <- p
    df <- rbind(df, kins)
  }
}


df_all <- df

df <- aggregate(cbind(A,B,C,D,E,F,G,H,I) ~ pair + pop, data = df_all, FUN = sum)
df['HETHET'] = df['E']
df['IBS0'] = df['C'] + df['G']
df['IBS1'] = df['B'] + df['D'] + df['F'] + df['H']
df['IBS2'] = df['A'] + df['E'] + df['I']
df['R0'] = df['IBS0'] / df['HETHET']
df['R1'] = df['HETHET'] / (df['IBS0'] +  df['IBS1'])
df['Kin'] = (df['HETHET'] - 2*(df['IBS0'])) / (df['IBS1'] + 2*df['HETHET'])
df$misschunk <- "none"

for (n in c(1:17)) {
  df2 <- aggregate(cbind(A,B,C,D,E,F,G,H,I) ~ pair + pop, data = df_all[!df_all$chunk == n,], FUN = sum)
  df2['HETHET'] = df2['E']
  df2['IBS0'] = df2['C'] + df2['G']
  df2['IBS1'] = df2['B'] + df2['D'] + df2['F'] + df2['H']
  df2['IBS2'] = df2['A'] + df2['E'] + df2['I']
  df2['R0'] = df2['IBS0'] / df2['HETHET']
  df2['R1'] = df2['HETHET'] / (df2['IBS0'] +  df2['IBS1'])
  df2['Kin'] = (df2['HETHET'] - 2*(df2['IBS0'])) / (df2['IBS1'] + 2*df2['HETHET'])
  df2$misschunk <- n
  df <- rbind(df, df2)
}

# C. semiargus

s <- "csemiargus"

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

df <- c()

for (n in c(1:21)) {
    kins <- do_derived_stats(
      read_ibspair_model0(
        paste0(
          "results/datasets/landuse-csemiargus/analyses/kinship/ibsrelate/landuse-csemiargus.ilCyaSemi1.1_",
          "all",
          "_chunk",
          n,
          "_allsites-filts.ibspair"
        )
      )
    )
    kins$chunk <- n
    kins$pop <- p
    df <- rbind(df, kins)
}

df_all <- df

df <- aggregate(cbind(A,B,C,D,E,F,G,H,I) ~ pair + pop, data = df_all, FUN = sum)
df['HETHET'] = df['E']
df['IBS0'] = df['C'] + df['G']
df['IBS1'] = df['B'] + df['D'] + df['F'] + df['H']
df['IBS2'] = df['A'] + df['E'] + df['I']
df['R0'] = df['IBS0'] / df['HETHET']
df['R1'] = df['HETHET'] / (df['IBS0'] +  df['IBS1'])
df['Kin'] = (df['HETHET'] - 2*(df['IBS0'])) / (df['IBS1'] + 2*df['HETHET'])
df$misschunk <- "none"

for (n in c(1:16)) {
  df2 <- aggregate(cbind(A,B,C,D,E,F,G,H,I) ~ pair + pop, data = df_all[!df_all$chunk == n,], FUN = sum)
  df2['HETHET'] = df2['E']
  df2['IBS0'] = df2['C'] + df2['G']
  df2['IBS1'] = df2['B'] + df2['D'] + df2['F'] + df2['H']
  df2['IBS2'] = df2['A'] + df2['E'] + df2['I']
  df2['R0'] = df2['IBS0'] / df2['HETHET']
  df2['R1'] = df2['HETHET'] / (df2['IBS0'] +  df2['IBS1'])
  df2['Kin'] = (df2['HETHET'] - 2*(df2['IBS0'])) / (df2['IBS1'] + 2*df2['HETHET'])
  df2$misschunk <- n
  df <- rbind(df, df2)
}

library(ggplot2)

ggplot(df[df$pair == "3_5" & df$pop == "Gotafors",], aes(x = R1, y = R0)) +
  geom_point() +
  ylim(c(0,0.5)) +
  xlim(c(0, 0.8))
