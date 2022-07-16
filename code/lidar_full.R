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

plot(ctg, chunk=TRUE)


## quality checks

las_check(ctg)


## classify ground
opt_output_files(ctg) <- paste0(tempdir(), "{ID}_ground")
classified_ctg <- classify_ground(ctg, algorithm = csf())

#classified_ctg <- readLAScatalog("ground_points")

#plot(classified_ctg, color = "Classification", size=3)

## calculate the dtm
opt_output_files(classified_ctg) <- paste0(tempdir(), "{ID}_dtm_aktuell")
opt_stop_early(classified_ctg) <- FALSE
opt_chunk_buffer(classified_ctg) <- 4

dtm_tin <- rasterize_terrain(classified_ctg, res = 1, algorithm = tin())

writeRaster(dtm_tin, "dtm.tif")

#dtm_tin <- rast("dtm.tif")

dtm_prod <- terra::terrain(dtm_tin, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- terra::shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)
plot(dtm_hillshade, col = gray(0:50/50), legend = FALSE)


## height normalization (not in parralel!)
dir.create(paste0(getwd(), "/nlas"))
opt_output_files(classified_ctg) <- paste0(getwd(), "/nlas/{ID}_nlas")
plan(sequential)

ctg_norm <- normalize_height(classified_ctg, dtm_tin)


## calculate digital height model 
plan(multisession, workers = availableCores()-2)

opt_output_files(classified_ctg) <- paste0(tempdir(), "{ID}_dhm_new")
dhm <- rasterize_canopy(classified_ctg, res = 1, algorithm = p2r())

writeRaster(dhm, "dhm.tif")
#dhm <- rast("dhm.tif")


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
#chm <- dhm - dtm_tin
writeRaster(chm, "chm.tif")
#chm <- rast("chm.tif")

## segment trees

#ctg_norm@output_options[["drivers"]][["SpatRaster"]][["param"]][["overwrite"]] <- TRUE

opt_filter(ctg_norm) <- "-keep_random_fraction 0.01" # accelelerate processing
opt_output_files(ctg_norm) <- ""
opt_chunk_buffer(ctg_norm) <- 4
#opt_restart(ctg_norm) <- 1
ttops <- locate_trees(ctg_norm, lmf(ws=40), uniqueness = "bitmerge")
#ttops_chm <- locate_trees(chm, lmf(ws=10))

ttops <- do.call(rbind,ttops)

st_write(ttops, "ttops.gpkg")

plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

# using lidar data
plan(sequential)

opt_output_files(ctg_norm) <- paste0(tempdir(), "/{*}_segmented")
algo <- dalponte2016(chm, ttops)
ctg_segmented <- segment_trees(ctg_norm, algo)

plan(multisession, workers=availableCores()-2)

crowns <- crown_metrics(ctg_segmented, func = .stdtreemetrics, geom = "convex")
plot(crowns["convhull_area"], main = "Crown area (convex hull)")

 # using the chm only (Without lidar data)
algo <- dalponte2016(chm, ttops)
crowns <- algo()

plot(crowns, col = pastel.colors(200))


## complex metrices
las <- readLAS("clouddc03becdccbf85e0.las", filter = "-keep_random_fraction 0.1")
las <- retrieve_pulses(las)
las <- filter_firstlast(las)
las <- filter_poi(las, NumberOfReturns > 1)
las@data <- las@data %>% 
  
  mutate(d = Z[1]-Z[2])

las@data[, d := Z[1]-Z[2], by = pulseID]
las <- filter_first(las)
las <- filter_poi(las, !is.na(d))
plot(las, bg = "white", color = "d", pal = viridis::magma(50), size = 4)
# distance to first return
d_first_last <- pixel_metrics(las, ~mean(d), 5)
plot(d_first_last, col = height.colors(50))
# average

## plots
plot_dtm3d(dtm_tif, bg = "white")
plot(nlas, size = 4, bg = "white")
plot(dhm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)



### shadow calculation (W h m^-2); June 21

shfiles <- list.files("shadows", pattern = "*.tif")

shadows <- rast(paste0("shadows/", shfiles))

shadow <- app(shadows, "sum")
plot(shadow, col = pals::inferno(50), 
     plg = list(title = expression(paste('Solar Irradiation', " ", "[Wh", " ",  m^-2,']'))))
