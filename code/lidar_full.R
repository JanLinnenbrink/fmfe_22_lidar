library(lidR)
library(future)
library(terra)
library(sf)

setwd("C:/0_Msc_Loek/M7_Fernerkundung/lidar")

plan(multisession, workers = availableCores()-4)

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

plot(classified_ctg, color = "Classification", size=3)

## calculate the dtm
opt_output_files(classified_ctg) <- paste0(tempdir(), "{ID}_dtm_aktuell_neu")
opt_stop_early(classified_ctg) <- FALSE
opt_chunk_buffer(classified_ctg) <- 4

dtm_tin <- rasterize_terrain(classified_ctg, res = 0.1, algorithm = tin())

writeRaster(dtm_tin, "dtm_0-1.tif")

dtm_tin <- rast("dtm.tif")

plot_dtm3d(dtm_tin, bg = "white") 

dtm_prod <- terra::terrain(dtm_tin, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- terra::shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)
plot(dtm_hillshade, col = gray(0:50/50), legend = FALSE)
plot(dtm_tin, col = height.colors(50),plg = list(title = "height (m above NN)"))


## height normalization (not in parralel!)
dir.create(paste0(getwd(), "/nlas"))
opt_output_files(classified_ctg) <- paste0(getwd(), "/nlas/{ID}_nlas")
plan(sequential)

ctg_norm <- normalize_height(classified_ctg, dtm_tin)

#ctg_norm <- readLAScatalog(paste0(getwd(), "/nlas/"))

## calculate digital height model 
plan(multisession, workers = availableCores()-2)

opt_output_files(classified_ctg) <- paste0(tempdir(), "{ID}_dhm_new_highres")
dhm <- rasterize_canopy(classified_ctg, res = 0.1, algorithm = p2r())

writeRaster(dhm, "dhm_0-1.tif")
dhm <- rast("dhm_smoothed.tif")


# post processing of dhm
fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
filled <- terra::focal(dhm, w, fun = fill.na)

w <- matrix(1, 3, 3)
dhm_smoothed <- terra::focal(filled, w, fun = mean, na.rm = TRUE)

writeRaster(dhm_smoothed, "dhm_smoothed.tif", overwrite=TRUE)
dhm_smoothed <- rast("dhm_smoothed.tif")
plot(dhm_smoothed,col = height.colors(50),plg = list(title = "height (m above NN)"))

## calculate canopy height model

chm <- rasterize_canopy(ctg_norm, res = 0.1, algorithm = p2r())
filled <- terra::focal(chm, w, fun = fill.na)
chm_smoothed <- focal(filled, w, fun = mean, na.rm = TRUE)

#chm <- dhm - dtm_tin
writeRaster(chm_smoothed, "chm_smoothed_0-1.tif", overwrite=TRUE)
chm_smoothed <- rast("chm_smoothed.tif")

plot(chm_smoothed,col = height.colors(50),plg = list(title = "vegetation height (m)"))


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
ttops <- st_read("ttops.gpkg")



plot(chm_smoothed, col = height.colors(50),plg = list(title = "vegetation height (m)"))
plot(sf::st_geometry(ttops[ttops$Z > 5,]), add = TRUE, pch = 3)

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

x <- plot(las, bg = "white", size = 4)
add_treetops3d(x, ttops)


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
plot(dhm_smoothed, col = height.colors(50),plg = list(title = "Height (m above NN)"))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)



### shadow calculation (W h m^-2); June 21

shfiles <- list.files("shadows", pattern = "*.tif")

shadows <- rast(paste0("shadows/", shfiles))

shadow <- app(shadows, "sum")

w <- matrix(1, 11,11)
fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
shadow_filled <- terra::focal(shadow, w, fun = fill.na)

plot(shadow_filled, col = pals::inferno(50), 
     plg = list(title = expression(paste('Solar Irradiation', " ", "[Wh", " ",  m^-2,']'))))

writeRaster(shadow_filled, "shadow.tif")

thermal <- rast("thermal.tif")

ext(thermal) <- c(ext(thermal)[1]/10, ext(thermal)[2]/10, ext(thermal)[3], ext(thermal)[4] )

res(thermal) <- res(shadow)


plot(thermal)
plot(shadow)
