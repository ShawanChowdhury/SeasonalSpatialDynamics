options(java.parameters = "-Xmx6g")

library(sp)
library(raster)
library(dismo)
library(rJava)
library(dplyr)
library(sf)
library(rgeos)
library(spocc)
library(ENMeval)
library(data.table)
library(stringr)
library(rgdal)
library(data.table)

# This script is for S1 Palaearctic, we repeated the same procedure for 11 realms and four seasons.

# Reading csv file
dc1 <- read.csv("S1.Palaearctic.csv", header = T)

dc_cl <- dc1 %>%
  dplyr::select("species", "decimalLon", "decimalLat")


# Remove records without coordinates
dc_cl <- dc_cl%>%
  filter(!is.na(decimalLon))%>%
  filter(!is.na(decimalLat))%>%
  filter(!is.na(species))

# Remove species with low occurrence records
dc_cl <- dc_cl %>%
  group_by(species) %>%
  filter(n() > 3) %>%
  ungroup()

# Remove duplicated records
dc_cl <- dc_cl[!duplicated(dc_cl),]

#attaching climatic layers
vars<-c("bio1.tif","bio2.tif","bio3.tif","bio4.tif","bio5.tif","bio6.tif","bio7.tif","bio8.tif","bio9.tif","bio10.tif","bio11.tif","bio12.tif","bio13.tif","bio14.tif","bio15.tif","bio16.tif","bio17.tif","bio18.tif","bio19.tif")
clim.stack <- stack(paste(getwd(),"/clim/", vars, sep=""))

species <- unique(dc_cl$species)

for (currspecies in species) try({
  print(currspecies)
  speciesname <- gsub(" ", "_", currspecies)
  dc <- dc_cl %>% filter(species==currspecies)
  d.c<- dc[,2:3]
  d.c <- as.data.frame(d.c)
  d.c.sp <- SpatialPoints(d.c)
  bb <- bbox(d.c.sp)
  bb.buf <- extent(bb[1]-10, bb[3]+10, bb[2]-10, bb[4]+10)
  #clim.stack <- crop(clim.stack, bb.buf)
  bg<-randomPoints(clim.stack[[19]], n=10000)
  bg <- as.data.frame(bg)
  pred.mod.allyr <- ENMevaluate(d.c, clim.stack, bg, method='checkerboard2', 
                                RMvalues=c(1,2), fc=c('L','LQ','LQP'), parallel=TRUE, algorithm='maxent.jar')
  
  aic.opt <- pred.mod.allyr@models[[which.max((pred.mod.allyr@results$avg.test.AUC))]]
  dir.create(path = speciesname)
  
  output_dir <- paste0(speciesname, "/")
  
  ModelOutput <- 
    as.data.frame(aic.opt@results) %>% 
    rownames_to_column(var = "param") %>%
    spread(key = param, value = V1) %>%
    mutate(species = speciesname) %>%
    dplyr::select(species, everything())
  
  
  write.csv(ModelOutput, file = paste0(output_dir, "S1PalaearcticModelOutput_", speciesname, ".csv"))
  
  VariableContribution <- var.importance((aic.opt))
  write.csv(VariableContribution, file = paste0(output_dir, "S1PalaearcticVariableContribution_", speciesname, ".csv"))
  
  r <- dismo::predict(aic.opt, clim.stack)
  writeRaster(r, paste0(output_dir, "S1PalaearcticMap_", speciesname, ".tif"), NAflag=-9999, overwrite = TRUE)
  
  r_bin <- r >= ModelOutput$Maximum.training.sensitivity.plus.specificity.Cloglog.threshold
  writeRaster(r_bin, paste0(output_dir, "S1PalaearcticMapbinary", speciesname, ".tif"), NAflag=-9999, overwrite = TRUE)
  
  dc_cl <- dc_cl[!(dc_cl$species %in% c(currspecies)), ]
  
  rm(aic.opt, bb, bb.buf, bg, d.c, d.c.sp, ModelOutput, pred.mod.allyr, r, r_bin, 
     VariableContribution, dc, currspecies, output_dir)
  
}, silent = FALSE)
