# Even after thresholding the model output, some rasters were not binary. 
# To solve the issue, we had to repeat the process.

options(java.parameters = "-Xmx6g")

library(sp)
library(raster)
library(dismo)
library(dplyr)
library(sf)
library(rgeos)
library(data.table)
library(stringr)
library(rgdal)
library(parallel)
library(foreach)

# In parallel
n_threads <- 5

SpeciesList <- read.csv("S1.SpeciesList.csv", header = T)

sp <- unique(SpeciesList$species)

tifs <- list.files(path = "MosaicRasters", pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)

tifs <- tifs[stringr::str_detect(tifs, "S1")]
tifs <- tifs[stringr::str_detect(tifs, "bin")]

cl <- makeCluster(n_threads, "PSOCK") # create workers
clusterEvalQ(cl, { # load packages into workers
  library(sp)
  library(dplyr)
  library(sf)
  library(raster)
  library(dismo)
  library(rgeos)
  library(data.table)
  library(stringr)
  library(rgdal)
  library(rgdal)
  })
clusterExport(cl, c("tifs")) # send to cluster

# Main processing
result <- try(parLapply(cl, sp, function(i) {
  print(i)
  speciesname <- gsub(" ", "_", i)
  raster_files <- tifs[stringr::str_detect(tifs, as.character(i))]
  #raster_files <- r[(stringr::str_detect(r, "bin"))]
  
  if(length(raster_files) > 0) {
    world_layer <- raster(raster_files)
    snap <- raster(resolution = c(0.04166667, 0.04166667), 
                   xmn = -180, xmx = 180, ymn = -90, ymx = 90, 
                   crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    world_layer <- projectRaster(world_layer, snap, method = "bilinear")
    
    world_layer <- (world_layer > 0) * 1
    
    dir.create(path = speciesname)
    output_dir <- paste0(speciesname, "/")
    writeRaster(world_layer, paste0(output_dir, "S1bin", speciesname, ".tif"), NAflag=-9999)
    
    # world <- getMap(resolution = "low")
    # png(world_layer, file = paste0(output_dir, "S1bin", speciesname, ".png"), width = 11, height = 8.5,res=110, units="in")
    # plot <- plot(world_layer)
    # plot <- plot(world, add = TRUE)
    # print(plot)
    # dev.off()
    
  }
}), silent = FALSE)

# Stop cluster
cl <- stopCluster(cl)

