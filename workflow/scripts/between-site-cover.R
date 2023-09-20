sink(file(snakemake@log[[1]], open = "wt"), type = "message")

library(landscapemetrics)
library(raster)
library(sf)

# read in land use raster and sample site coordinates
land_raster <- raster(snakemake@input[["raster"]])
field_sites <- read.table(snakemake@input[["sites"]], header = TRUE, sep = "\t")
outfile <- snakemake@output[["tsv"]]

# specify raster codes corresponding to each land use type
grass_codes <- c(42)
arable_codes <- c(3)
forest_codes <- c(
  111, 112, 113, 114, 115, 116, 117, 118, 121, 122, 123, 124, 125, 126, 127, 128
)
water_codes <- c(61, 62)

# generate all combinations of field sites
combos <- combn(field_sites$Site.name, 2)

sites_sf <- st_as_sf(
  field_sites, coords = c("SWEREF99.E","SWEREF99.N"), crs = "EPSG:3006"
)

# estimate proportion of each land use type in 10 km buffer between all sites
df <- c()
for (i in c(seq_len(ncol(combos)))) {
  pop1 <- combos[1, i]
  pop2 <- combos[2, i]
  line <- st_cast(
    st_union(
      sites_sf[sites_sf$Site.name == pop1, ]$geometry,
      sites_sf[sites_sf$Site.name == pop2, ]$geometry
    ),
    "LINESTRING"
  )
  distance <- st_length(line)
  bufferline <- st_buffer(line, dist = 5000, endCapStyle = "FLAT")
  landuse <- lsm_c_pland(crop(land_raster, st_as_sf(bufferline), mask = TRUE))
  arable <- sum(landuse[landuse$class %in% arable_codes, "value"])
  grassland <- sum(landuse[landuse$class %in% grass_codes, "value"])
  forest <- sum(landuse[landuse$class %in% forest_codes, "value"])
  water <- sum(landuse[landuse$class %in% water_codes, "value"])
  other <- 100 - sum(arable, grassland, forest, water)
  row <- data.frame(
    pop1, pop2, distance, arable, grassland, forest, water, other
  )
  df <- rbind(df, row)
}

# write results to table
write.table(
  df,
  file = outfile,
  quote = FALSE,
  append = FALSE,
  sep = "\t",
  row.names = FALSE
)
