library(lidR)
library(future)
library(terra)
library(sf)
library(ggplot2)

setwd("C:/0_Msc_Loek/M7_Fernerkundung/fmfe_22_lidar")

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

dtm_tin <- rast("raster/allPoints_dtm.tif")

plot_dtm3d(dtm_tin, bg = "white") 

dtm_prod <- terra::terrain(dtm_tin, v = c("slope", "aspect"), unit = "radians")
dtm_hillshade <- terra::shade(slope = dtm_prod$slope, aspect = dtm_prod$aspect)
plot(dtm_hillshade, col = gray(0:50/50), legend = FALSE)

p <- recordPlot()
plot.new() ## clean up device
plot(dtm_tin, col = height.colors(50),plg = list(title = "height (m above NN)"))
plot(river_vec, add=TRUE, border="black")

d <- plot(dtm_tin, col = height.colors(50),plg = list(title = "height (m above NN)"))
d <- plot(river_vec, add=TRUE, border="black")

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
chm_smoothed <- rast("./raster/chm_smoothed.tif")

plot(chm_smoothed,col = height.colors(50),plg = list(title = "vegetation height (m)"))
plot(river_vec, add=TRUE)


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
plot(river_vec$geometry, add=TRUE)



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


### combine shadow and thermal data

dtm_tin <- rast("raster/allPoints_dtm.tif")

dtm_tin_low <- dtm_tin
dtm_tin_low[dtm_tin_low>100.5] <- NA
plot(dtm_tin_low)

river_vec <- terra::as.polygons(dtm_tin_low, dissolve = TRUE) |> 
  st_as_sf() |> 
  st_buffer(10) |> 
  st_buffer(-10) |> 
  nngeo::st_remove_holes()

river_vec <- st_union(river_vec)
plot(river_vec)

st_write(river_vec, "river.shp")

shadow <- rast(paste0(getwd(), "/raster/shadow.tif"))
thermal <- rast(paste0(substr(getwd(), start = 1, stop = nchar(getwd())-13),
                       "Thermal/Feldmethoden_Thermal_22062022_ETRS_modifiziert_p1_test2.tif")) |> 
  project(crs(shadow))
thermal_bbox <- st_bbox(st_as_sf(as.polygons(thermal, dissolve=TRUE)))

thermal <- crop(thermal,vect(st_as_sfc(thermal_bbox)))
shadow <- crop(shadow, thermal) |> 
  resample(thermal) |> 
  mask(thermal)

th_sh <- c(thermal, shadow)
names(th_sh) <- c("thermal", "shadow")

th_sh$thermal[th_sh$thermal > 20] <- NA

th_sh_r <- mask(th_sh, st_transform(river_vec, st_crs(th_sh$thermal)))

raster::writeRaster(raster::stack(th_sh), "./raster/thermal_shadow.tif")

th_sh <- rast("./raster/thermal_shadow.tif")
names(th_sh) <- c("thermal", "shadow")

thermal <- th_sh$thermal
shadow <- th_sh$shadow


writeRaster(thermal, "./raster/thermal_res.tif")
writeRaster(shadow, "./raster/shadow_res.tif")

thermal <- rast("./raster/thermal_res.tif")
shadow <- rast("./raster/shadow_res.tif")

plot(th_sh_r, col = pals::inferno(50))

plot(shadow)
plot(st_transform(river_vec, crs(shadow)), add=TRUE)
plot(thermal)
plot(st_transform(river_vec, crs(shadow)), add=TRUE)

# sample random points
source(paste0(substr(getwd(), start = 1, stop = nchar(getwd())-13),
              "sampleFromArea.R"))

rp <- spatSample(thermal_river, 10000, method="regular", replace=FALSE, na.rm=FALSE,
                 as.raster=FALSE, as.df=TRUE, as.points=TRUE, values=TRUE, 
                 cells=FALSE, xy=FALSE, ext=NULL, warn=TRUE, weights=NULL) |> 
  st_as_sf()

rp$thermal <- rp$Feldmethoden_Thermal_22062022_ETRS_modifiziert_p1_test2
rp$Feldmethoden_Thermal_22062022_ETRS_modifiziert_p1_test2 <- NULL
rp <- rp[!is.na(rp$thermal),]

plot(rp)

sh_th <- extract(shadow, rp, bind = TRUE)
sh_th$irradiance <- sh_th$sum
sh_th$sum <- NULL
sh_th <- sh_th[,!is.na(sh_th$irradiance)]
sh_df <- st_drop_geometry(st_as_sf(sh_th))
sh_df <- sh_df[!is.na(sh_df$irradiance) & sh_df$irradiance>0, ]


ggplot(sh_df, aes(x=irradiance,y=thermal)) +
  geom_bin2d(bins = 40) +
  geom_smooth(method = "lm") +
  xlab( expression(paste('Solar Irradiation', " ", "[W", "h", " ",  m^-2,']'))) +
  ylab("Temperature [°C]")

ggplot(sh_df, aes(x=irradiance,y=thermal)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab( expression(paste('Solar Irradiation', " ", "[W", "h", " ",  m^-2,']'))) +
  ylab("Temperature [°C]")

lm_pl <- ggplot(sh_df, aes(x=irradiance,y=thermal)) +
  geom_density_2d_filled(show.legend = FALSE, contour_var = "density") +
  geom_smooth(method = "lm", se = FALSE, colour = "white") +
  xlab( expression(paste('Solar Irradiation', " ", "[W", "h", " ",  m^-2,']'))) +
  ylab("Temperature [°C]") +
  theme(aspect.ratio = 1)

ggsave(paste0(getwd(), "/plots/temp_shadow_only_river.pdf"), lm_pl,
       width = 12, height = 12, units = "cm")





# plots -------------------------------------------------------------------

dtm_tin <- rast("raster/allPoints_dtm.tif")

chm_smoothed <- rast("./raster/chm_smoothed_0-1.tif")
ttops <- st_read("ttops.gpkg")

shadow_filled <- rast("./raster/shadow.tif")
river_vec <- st_read("river.shp") |> st_union()


adj = 0
font = 2

png("comb.png", width = 50, height = 20, units="cm", res=1000)
par(mfrow=c(1,3),oma=c(0,0,5,0), ann=TRUE, cex=1.2)

plot(dtm_tin, col = height.colors(50))
plot(river_vec, add=TRUE, border="black")
mtext("Height above NN [m]", side=3, line=1, cex=1.4, adj=adj, font = font)

plot(chm_smoothed, col = height.colors(50))
plot(river_vec, add=TRUE, border="black")
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)
mtext("Canopy height [m]",  side=3, line=1, cex=1.4, adj=adj, font = font)

plot(shadow_filled, col = pals::inferno(50))
plot(river_vec, add=TRUE, border="black")
mtext(expression(paste(bold('Solar Irradiation'), " ", bold("[Wh"), " ",  bold(m^-2),bold(']'))),  side=3, line=1,
      cex=1.4, adj=adj)

dev.off()


mat <- matrix(c(1, 2, 3), 
              nrow = 1, ncol = 3,
              byrow = TRUE)

l <- layout(mat = mat)

plot(dtm_tin, col = height.colors(50),plg = list(title = "height (m above NN)"))
plot(river_vec, add=TRUE, border="black")


plot(shadow_filled, col = pals::inferno(50), 
     plg = list(title = expression(paste('Solar Irradiation', " ", "[Wh", " ",  m^-2,']'))))

layout.show(l)
plot.new()



