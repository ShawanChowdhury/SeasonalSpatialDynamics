options(java.parameters = "-Xmx6g")

library(sp)
library(dplyr)
library(sf)
library(raster)
library(dismo)
library(rgeos)
library(data.table)
library(stringr)
library(rgdal)
library(rworldmap)
library(tmap)

# Mosaic multiple rasters and saving as a single raster
# Repeat the same procedure for other seasonal maps
SpeciesList <- read.csv("S1.SpeciesList.csv", header = T)

sp <- unique(SpeciesList$species11)

tifs <- list.files(path = "Rasters_SeasonalEcoregional", pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)

tifs <- tifs[stringr::str_detect(tifs, "S1")]
tifs <- tifs[stringr::str_detect(tifs, "binary")]


for (i in sp) {
  print(i)
  speciesname <- gsub(" ", "_", i)
  r <- tifs[stringr::str_detect(tifs, as.character(i))]
  raster_files <- r[stringr::str_detect(r, "binary")]
  
  snap <- raster(resolution = c(0.04166667, 0.04166667), 
                 xmn = -180, xmx = 180, ymn = -90, ymx = 90, 
                 crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
  
  mosaicList <- function(rasList){
    
    #Internal function to make a list of raster objects from list of files.
    ListRasters <- function(list_names) {
      raster_list <- list() # initialise the list of rasters
      
      for (j in 1:(length(list_names))){ 
        grd_name <- list_names[j] # list_names contains all the names of the images in .grd format
        raster_file <- raster::raster(grd_name)
        
        raster_file <- projectRaster(raster_file, snap, method = "bilinear")
        
      }
      raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
      
      
    }
    
    #convert every raster path to a raster object and create list of the results
    raster.list <- sapply(rasList, FUN = ListRasters)
    
    # edit settings of the raster list for use in do.call and mosaic
    names(raster.list) <- NULL
    #####This function deals with overlapping areas
    raster.list$fun <- sum
    # raster.list$tolerance <- 0.1
    
    #run do call to implement mosaic over the list of raster objects.
    mos <- do.call(raster::mosaic, raster.list)
    
    #set crs of output
    crs(mos) <- crs(x = raster(rasList[1]))
    return(mos)
  }
  
  if(length(raster_files) > 1) {
    world_layer <- mosaicList(raster_files)
    
  } else {
    world_layer <- raster(raster_files)
    world_layer <- projectRaster(world_layer, snap, method = "bilinear")
  }
  
  dir.create(path = speciesname)
  output_dir <- paste0(speciesname, "/")
  writeRaster(world_layer, paste0(output_dir, "S1bin", speciesname, ".tif"), NAflag=-9999)
  
  # world <- getMap(resolution = "low")
  # png(world_layer, file = paste0(output_dir, "S2bin", speciesname, ".png"), width = 11, height = 8.5,res=110, units="in")
  # plot <- plot(world_layer)
  # plot <- plot(world, add = TRUE)
  # print(plot)
  # dev.off()
  
}

# In parallel
library(parallel)
library(foreach)

SpeciesList <- read.csv("S1.SpeciesList.csv", header = T)

species <- unique(SpeciesList$species)

tifs <- list.files(path = "Rasters_SeasonalEcoregional", pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)

tifs <- tifs[stringr::str_detect(tifs, "S4")]
tifs <- tifs[stringr::str_detect(tifs, "binary")]

# Define variables
n_threads <- 2

# Setup parallel cluster
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
  library(parallel)
  library(foreach)
})
clusterExport(cl, c("tifs")) # send to cluster

# Main processing
result <- try(parLapply(cl, species, function(i) {
      print(i)
    speciesname <- gsub(" ", "_", i)
    r <- tifs[stringr::str_detect(tifs, as.character(i))]
    head(r)
    raster_files <- r[stringr::str_detect(r, "binary")]
    
    snap <- raster(resolution = c(0.04166667, 0.04166667), 
                   xmn = -180, xmx = 180, ymn = -90, ymx = 90, 
                   crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 
    
    mosaicList <- function(rasList){
      
      #Internal function to make a list of raster objects from list of files.
      ListRasters <- function(list_names) {
        raster_list <- list() # initialise the list of rasters
        
        for (j in 1:(length(list_names))){ 
          grd_name <- list_names[j] # list_names contains all the names of the images in .grd format
          raster_file <- raster::raster(grd_name)
          
          raster_file <- projectRaster(raster_file, snap, method = "bilinear")
          
        }
        raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
        
        
      }
      
      #convert every raster path to a raster object and create list of the results
      raster.list <- sapply(rasList, FUN = ListRasters)
      
      # edit settings of the raster list for use in do.call and mosaic
      names(raster.list) <- NULL
      #####This function deals with overlapping areas
      raster.list$fun <- sum
      # raster.list$tolerance <- 0.1
      
      #run do call to implement mosaic over the list of raster objects.
      mos <- do.call(raster::mosaic, raster.list)
      
      #set crs of output
      crs(mos) <- crs(x = raster(rasList[1]))
      return(mos)
    }
    
    if(length(raster_files) > 1) {
      world_layer <- mosaicList(raster_files)
      
    } else {
      world_layer <- raster(raster_files)
      world_layer <- projectRaster(world_layer, snap, method = "bilinear")
    }
    
    dir.create(path = speciesname)
    output_dir <- paste0(speciesname, "/")
    writeRaster(world_layer, paste0(output_dir, "S4bin", speciesname, ".tif"), NAflag=-9999, overwrite = TRUE)
  
}), silent = FALSE)

# Stop cluster
cl <- stopCluster(cl)

