library(lidR)
library(future)
library(terra)
library(sf)

setwd("C:/0_Msc_Loek/M7_Fernerkundung/lidar")

plan(multisession, workers = availableCores()-2)

ctg <-  readLAScatalog("clouddc03becdccbf85e0.las")

opt_chunk_size(ctg) <- 40
opt_stop_early(ctg) <- TRUE
opt_chunk_buffer(ctg) <- 4

plot(classified_ctg, chunk=TRUE)


## quality checks

las_check(ctg)


## classify ground
opt_output_files(ctg) <- paste0(tempdir(), "{ID}_ground")
classified_ctg <- classify_ground(ctg, algorithm = csf())

#classified_ctg <- readLAScatalog("ground_points")

plot(classified_ctg, color = "Classification", size=3)

## calculate the dtm
opt_output_files(classified_ctg) <- paste0(tempdir(), "{ID}_dtm_aktuell")
opt_stop_early(classified_ctg) <- FALSE
opt_chunk_buffer(classified_ctg) <- 4

dtm_tin <- rasterize_terrain(classified_ctg, res = 1, algorithm = tin())

writeRaster(dtm_tin, "dtm.tif")

dtm_prod <- terra::terrain(dtm_tin, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- terra::shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)
plot(dtm_hillshade, col = gray(0:50/50), legend = FALSE)


## height normalization (not parralel!)
opt_output_files(classified_ctg) <- paste0(tempdir(), "{ID}_nlas1")
plan(sequential)

ctg_norm <- normalize_height(classified_ctg, dtm_tin)


## calculate digital height model 
plan(multisession, workers = availableCores()-2)

opt_output_files(classified_ctg) <- paste0(tempdir(), "{ID}_dhm_new")
dhm <- rasterize_canopy(classified_ctg, res = 1, algorithm = p2r())

writeRaster(dhm, "dhm.tif")

# post processing of dhm
fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
w <- matrix(1, 3, 3)
filled <- terra::focal(dhm, w, fun = fill.na)
smoothed <- terra::focal(dhm, w, fun = mean, na.rm = TRUE)

dhms <- c(dhm, filled, smoothed)
names(dhms) <- c("Base", "Filled", "Smoothed")
plot(dhms, col = col)

## calculate canopy height model

chm <- rasterize_canopy(ctg_norm, res = 1, algorithm = p2r())
chm <- dhm - dtm_tin
writeRaster(chm, "chm.tif")


## segment trees

opt_output_files(ctg_norm) <- ""
ttops <- locate_trees(ctg_norm, lmf(ws=10), uniqueness = "bitmerge")
ttops_chm <- locate_trees(chm, lmf(ws=10))

st_write(ttops, "ttops.shp")

# using lidar data
opt_output_files(ctg_norm) <- paste0(tempdir(), "/{*}_segmented")
algo <- dalponte2016(ctg_norm, ttops)
ctg_segmented <- segment_trees(ctg_norm, algo)

crowns <- crown_metrics(ctg_norm, func = .stdtreemetrics, geom = "convex")
plot(crowns["convhull_area"], main = "Crown area (convex hull)")

# using the chm only (Without lidar data)
algo <- dalponte2016(chm, ttops)
crowns <- algo()

plot(crowns, col = pastel.colors(200))


## complex metrices
ctg <- retrieve_pulses(ctg)
ctg <- filter_firstctgt(ctg)
ctg <- filter_poi(ctg, NumberOfReturns > 1)
ctg@data[, d := Z[1]-Z[2], by = pulseID]
ctg <- filter_first(ctg)
ctg <- filter_poi(ctg, !is.na(d))
plot(ctg, bg = "white", color = "d", colorPalette = viridis::magma(50), size = 4)
# distance to first return
d_first_ctgt <- pixel_metrics(ctg, ~mean(d), 5)
plot(d_first_ctgt, col = height.colors(50))
# average

## plots
plot_dtm3d(dtm_tin, bg = "white")
plot(nlas, size = 4, bg = "white")
plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

