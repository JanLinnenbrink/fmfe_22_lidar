setwd("~/Desktop/Uni_Muenster/Master/Fernerkundung/LIDAR Aa/fmfe_22_lidar")

install.packages("lidR")
install.packages("RCSF")
install.packages("rgl")
install.packages("MBA")
library(lidR)
library(RCSF)
library(rgl)
library(MBA)
library(terra)

#unselected ist zu groß, kann eingeladen werden, aber nicht bearbeitet
#las <- readLAS("Punktwolke.las")

#Validating LIDAR Data
#las_check(las)


#Verkleinerten Datensatz einladen
#las_short <- readLAS("Punktwolke.las", select = "xyziatnr")
las_short <- readLAS("Punktwolke.las", filter = "-keep_random_fraction 0.1")
las_check(las_short)
plot(las_short, color="RGB")

#Querschnitt 2D rendering
p1 <- c(51.944803, 7.573000)
p2 <- c(51.944600, 7.573300)
las_tr <- clip_transect(las_short, p1, p2, width = 4, xz = TRUE)
#Plots kann ich durch Probleme mit rgl Package nicht darstellen


#Ground Classification
las_short <- classify_ground(las_short, algorithm = csf())

#Digital Terrain Model
dtm_tin <- rasterize_terrain(las_short, algorithm = tin())
writeRaster(dtm_tin, "Terrain_Model.tif")
plot(dtm_tin)
plot(las_short, size = 3, bg = "white")

#
#nlas <- las_short - dtm_tin
nlas <- normalize_height(las_short, tin(), dtm = dtm_tin)
plot(nlas, size = 4, bg = "wihte")



