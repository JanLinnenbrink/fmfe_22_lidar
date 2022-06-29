library(lidR)
library(future)
library(terra)

setwd("C:/0_Msc_Loek/M7_Fernerkundung/lidar")

plan(multisession, workers = 4L)

ctg <-  readLAScatalog("clouddc03becdccbf85e0.las")

opt_chunk_size(ctg) <- 120
opt_stop_early(ctg) <- FALSE

plot(ctg, chunk=TRUE)


## calculate tree tops
opt_output_files(ctg) <- paste0(tempdir(), "{ID}_trees")

ttops <- locate_trees(las, lmf(ws = 5))

#plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

