options(java.parameters = "-Xmx6g")

library(dplyr)
library(countrycode)
library(CoordinateCleaner)
library(data.table)

# Data DOI: https://doi.org/10.26197/ala.d21a654d-fcb4-4193-9879-798baf7d4a08

# We downloaded all the butterfly records from the ALA and manually selected
# 25 species based on opinions from experts (Michael Braby and Myron Zalucki)

df <- read.csv("NonMigrants/NonMigrants.csv")
head(df)

# Remove records without coordinates
df <- df%>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# Remove blank cells
df <- df[!(is.na(df$species) | df$species == ""),]
df <- df[!(is.na(df$decimalLongitude) | df$decimalLongitude == ""),]
df <- df[!(is.na(df$decimalLatitude) | df$decimalLatitude == ""),]
df <- df[!(is.na(df$month) | df$month == ""),]

# Remove duplicated records
df1 <- df[!duplicated(df),]

rm(df)

# Remove species with low occurrence records
df2 <- df1 %>%
  group_by(species, season) %>%
  filter(n() > 3) %>%
  ungroup()

write.csv(df2, "NonMigrants/CleanedRecords_NonMigrants.csv")

# Select groups
dc_cl <- read.csv("NonMigrants/CleanedRecords_NonMigrants.csv")
head(dc_cl)

dc_cl <- dc_cl %>%
  dplyr::select("species", "decimalLongitude", "decimalLatitude", "month")

# Group by seasons
S1 <- dc_cl %>%
  filter(month %in% c(12, 1, 2))

S2 <- dc_cl %>%
  filter(month %in% c(3, 4, 5))

S3 <- dc_cl %>%
  filter(month %in% c(6, 7, 8))

S4 <- dc_cl %>%
  filter(month %in% c(9, 10, 11))

write.csv(S1, "NonMigrants/S1.csv")
write.csv(S2, "NonMigrants/S2.csv")
write.csv(S3, "NonMigrants/S3.csv")
write.csv(S4, "NonMigrants/S4.csv")

##################################################
rm(list = ls())

# Running habitat suitability maps

library(sp)
library(raster)
library(dismo)
library(dplyr)
library(sf)
library(rgeos)
library(spocc)
library(ENMeval)
library(data.table)
library(stringr)
library(rgdal)


#attaching climatic layers
df <- read.csv("NonMigrants/SpeciesList.csv")
dc_cl <- read.csv("NonMigrants/S1.csv")
dc_cl <- dc_cl %>%
  dplyr::select("species", "decimalLongitude", "decimalLatitude")

vars<-c("bio1.tif","bio2.tif","bio3.tif","bio4.tif","bio5.tif","bio6.tif","bio7.tif","bio8.tif","bio9.tif","bio10.tif","bio11.tif","bio12.tif","bio13.tif","bio14.tif","bio15.tif","bio16.tif","bio17.tif","bio18.tif","bio19.tif")
clim.stack <- stack(paste(getwd(),"/NonMigrants/clim/", vars, sep=""))

species <- unique(dc_cl$species)

for (currspecies in species) try( {
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
  
  
  # write.csv(ModelOutput, file = paste0(output_dir, "S1AustralianModelOutput_", speciesname, ".csv"))
  
  VariableContribution <- var.importance((aic.opt))
  # write.csv(VariableContribution, file = paste0(output_dir, "S1AustralianVariableContribution_", speciesname, ".csv"))
  
  r <- dismo::predict(aic.opt, clim.stack)
  # writeRaster(r, paste0(output_dir, "S1AustralianMap_", speciesname, ".tif"), NAflag=-9999)
  
  r_bin <- r >= ModelOutput$Maximum.training.sensitivity.plus.specificity.Cloglog.threshold
  writeRaster(r_bin, paste0(output_dir, "S1AustralianMapbinary", speciesname, ".tif"), NAflag=-9999)
  
  dc_cl <- dc_cl[!(dc_cl$species %in% c(currspecies)), ]
  
  rm(aic.opt, bb, bb.buf, bg, d.c, d.c.sp, ModelOutput, pred.mod.allyr, r, r_bin, 
     VariableContribution, dc, dc1, AOO, output_dir)
  
}, silent = FALSE)


########################
# Migratory percentages

tifs <- list.files(path = "NonMigrants/SDM/",pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
tifs <- tifs[stringr::str_detect(tifs, "binary")]


SpeciesList <- read.csv("NonMigrants/SpeciesList.csv", header = T)

sp <- unique(SpeciesList$species)

for (i in sp) try( {
  
  print(i)
  speciesname <- gsub(" ", "_", i)
  r <- tifs[stringr::str_detect(tifs, as.character(i))]
  
  raster_files <- r[stringr::str_detect(r, "binary")]
  
  if(length(raster_files) > 1) {
    #convert every raster path to a raster object and create list of the results
    r_list <- sapply(raster_files, FUN = raster::raster)
    
    # stack seasonal raster layers
    r_stack <- raster::stack(x = r_list)
    
    # sum each pixel (range should be 0-4)
    r_sum <- sum(r_stack)
    
    r_migr_bin <- r_sum
    
    n_seasons <- max(r_migr_bin[], na.rm = T)
    
    r_migr_bin[r_migr_bin == 0] <- -1
    r_migr_bin[r_migr_bin %in% 1:(n_seasons - 1)] <- 1
    r_migr_bin[r_migr_bin == n_seasons] <- 0
    
    
    writeRaster(r_migr_bin, paste0("NonMigrants/MPR/", speciesname, "_Migr.tif"), NAflag=-9999, overwrite = TRUE)
    
    
    results <- table(getValues(r_migr_bin))
    
    ModelOutput <-
      as.data.frame(results) %>%
      mutate(species = speciesname) %>%
      dplyr::select(species, everything())
    
    write.csv(ModelOutput, file = paste0("NonMigrants/MPR/MP_", speciesname, ".csv"))
    
  }
  
}, silent = FALSE)



# Merge csvs
input_folder <- "NonMigrants/MPR"
range <- dir(input_folder, "^.*\\.csv$", full.names = TRUE)
migr.per <- plyr::ldply(range, readr::read_csv)
write.csv(migr.per, "NonMigrants/MPR_NM.csv")

# Calculation
mp <- read.csv("NonMigrants/MPR_NM.csv", header = T)

migr.per <- mp %>%
  group_by(species) %>%
  summarise(TotalCells = sum(Freq), migratory = sum(migr), mpr = (migratory/TotalCells)*100) %>%
  ungroup()

write.csv(migr.per, "NonMigrants/MPR_calc.csv")


# Boxplot
df1 <- read.csv("Migr&NM_MPR.csv")

ggplot(df1, aes(x = group, y = mpr)) +
  geom_boxplot() + theme_classic() + ylab("Seasonal switching") + xlab("")