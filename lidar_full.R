library(lidR)
library(future)
library(terra)

setwd("C:/0_Msc_Loek/M7_Fernerkundung/lidar")

plan(multisession, workers = 4L)

ctg <-  readLAScatalog("clouddc03becdccbf85e0.las")

opt_chunk_size(ctg) <- 120
opt_stop_early(ctg) <- FALSE

plot(ctg, chunk=TRUE)


## quality checks

las_check(ctg)


## classify ground
ctg <- classify_ground(ctg, algorithm = pmf(ws = 5, th = 3))


## calculate the dtm
opt_output_files(ctg) <- paste0(tempdir(), "{ID}_dtm")
dtm_tin <- rasterize_terrain(ctg, res = 1, algorithm = tin())


## height normalization
opt_output_files(ctg) <- paste0(tempdir(), "{ID}_nlas")
nlas <- ctg - dtm


## calculate canopy height model
opt_output_files(ctg) <- paste0(tempdir(), "{ID}_chm")
chm <- rasterize_canopy(ctg, res = 1, algorithm = p2r())


## calculate tree tops
opt_output_files(ctg) <- paste0(tempdir(), "{ID}_trees")
ttops <- locate_trees(las, lmf(ws = 5))


## c


## plots
plot_dtm3d(dtm_tin, bg = "white")
plot(nlas, size = 4, bg = "white")
plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

