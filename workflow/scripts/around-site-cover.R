sink(file(snakemake@log[[1]], open="wt"), type = "message")

library(landscapemetrics)
library(raster)
library(sf)
library(fasterize)

# read in land use raster and field site coordinates
land_raster <- raster(snakemake@input[["raster"]])
field_sites <- read.table(snakemake@input[["sites"]], header = TRUE, sep = "\t")
outfile <- snakemake@output[["tsv"]]

land_raster <- raster("~/working/bioinfo/projects/landuse-manuscript/resources/gis/nmd2018bas_ogeneraliserad_v1_1.tif")

field_sites <- read.table("~/working/bioinfo/projects/landuse-manuscript/resources/gis/polyommatini_land-use_sampling.tsv", header = TRUE, sep = "\t")

tuva <- read_sf('~/working/gis/tuva-nmd-corr/Arslager-naturtyper-2023/Naturtyper 2023.shp')


# tuva raster
buffered <- st_buffer(st_as_sf(
  field_sites,
  coords = c("SWEREF99.E", "SWEREF99.N"),
  crs = "EPSG:3006"
), dist = 25000)

temp_raster <- raster(ext = extent(buffered[,1]), res = res(land_raster))
tuvarast <- fasterize(tuva, temp_raster)
values(tuvarast)[is.na(values(tuvarast))] <- 0

# convert field site coordinates to spatial object
sites_sp <- as_Spatial(
  st_as_sf(
    field_sites,
    coords = c("SWEREF99.E", "SWEREF99.N"),
    crs = "EPSG:3006"
  )
)

# set codes for each land use type from raster
grass_codes <- c(42)
arable_codes <- c(3)
forest_codes <- c(
  111, 112, 113, 114, 115, 116, 117, 118, 121, 122, 123, 124, 125, 126, 127, 128
)
water_codes <- c(61, 62)

# define radii to estimate land use within
buffer_radius <- c(
  100, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 15000, 20000
)


df <- data.frame()

# calculate proportion of each land use type within the defined radii
for (r in buffer_radius) {
  classmets_df <- sample_lsm(
    landscape = land_raster,
    y = sites_sp,
    plot_id = sites_sp$Site.name,
    shape = "circle",
    size = r,
    what = "lsm_c_pland",
    return_raster = F
  )
  
  tuvamets_df <- sample_lsm(
    landscape = tuvarast,
    y = sites_sp,
    plot_id = sites_sp$Site.name,
    shape = "circle",
    size = r,
    what = "lsm_c_pland",
    return_raster = F
  )
  
  for (i in unique(classmets_df$plot_id)) {
    site <- i
    radius <- r
    arable <- sum(classmets_df[classmets_df$plot_id == i & classmets_df$class %in% arable_codes, ]$value)
    grassland <- sum(classmets_df[classmets_df$plot_id == i & classmets_df$class %in% grass_codes, ]$value)
    sng <- sum(tuvamets_df[tuvamets_df$plot_id == i & tuvamets_df$class == 1, ]$value)
    forest <- sum(classmets_df[classmets_df$plot_id == i & classmets_df$class %in% forest_codes, ]$value)
    water <- sum(classmets_df[classmets_df$plot_id == i & classmets_df$class %in% water_codes, ]$value)
    other <- 100 - sum(arable, grassland, forest, water)
    row <- data.frame(site, radius, arable, grassland, sng, forest, water, other)
    df <- rbind(df, row)
  }
}

colnames(df) <- c(
  "site", "radius", "arable", "grassland", "sng", "forest", "water", "other"
)

write.table(df, file = "results/landscape/around_site_cover_tuva.tsv", quote = FALSE, append = FALSE, sep = "\t", row.names = FALSE)
# write table to file
write.table(
  df,
  file = outfile,
  quote = FALSE,
  append = FALSE,
  sep = "\t",
  row.names = FALSE
)
