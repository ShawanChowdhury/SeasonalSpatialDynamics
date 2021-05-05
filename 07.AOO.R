options(java.parameters = "-Xmx6g")

library(sp)
library(dplyr)
library(sf)
library(raster)
library(dismo)
library(maptools)
library(rgeos)
library(data.table)
library(stringr)
library(rgdal)
library(rworldmap)
library(tmap)
library(cowplot)
library(parallel)
library(foreach)

# Estimating AOO
SpeciesList <- read.csv("SpeciesList_updated.csv", header = T)
sp <- unique(SpeciesList$species)

tifs <- list.files(path = "Rasters_Global",pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
tifs <- tifs[stringr::str_detect(tifs, "S1bin")]


n_threads <- 16
cl <- makeCluster(n_threads, "PSOCK") # create workers
clusterEvalQ(cl, { # load packages into workers
  library(sp)
  library(dplyr)
  library(sf)
  library(raster)
  library(dismo)
  library(maptools)
  library(rgeos)
  library(data.table)
  library(stringr)
  library(rgdal)
  library(rworldmap)
  library(tmap)
  library(cowplot)
  library(parallel)
  library(foreach)
})
clusterExport(cl, c("tifs")) # send to cluster

# Main processing
result <- parLapply(cl, sp, function(i) {
  print(i)
  speciesname <- gsub(" ", "_", i)
  r <- tifs[stringr::str_detect(tifs, as.character(i))]
  raster_files <- r[(stringr::str_detect(r, "S1bin"))]
  
  if(length(raster_files) > 0) {
    world_layer <- raster(raster_files)
    AOO <- cellStats(world_layer, "sum")
    write.csv(AOO, file = paste0("S1AOO_", speciesname, ".csv"))
    
  }
})



# Stop cluster
cl <- stopCluster(cl)



# Merge csvs
input_folder <- "AOO"
range <- dir(input_folder, "^.*\\.csv$", full.names = TRUE)
migr.per <- plyr::ldply(range, readr::read_csv)
migr.per <- 
  migr.per %>% 
  mutate(species = range) %>%
  dplyr::select(species, everything())

write.csv(migr.per, "AOO.csv")


# Calcuating mean, range and fluctuation

aoo <- read.csv("AOO.csv", header = T)

aooR <- aoo %>%
  group_by(species) %>%
  summarise(seasons = n(),mean_AOO = mean(aoo), range = (max(aoo) - min(aoo)), fluctuation = max(aoo)/min(aoo)) %>%
  ungroup()

write.csv(aooR, "AOO_Calc.csv")



# Figure
aoo <- read.csv("AOO_Calc.csv", header = T)

# Conditional colouring
aoo %>% mutate(Color = ifelse(fluctuation > 10, "red", "blue")) %>%
  ggplot(aes(mean_AOO_m, fluctuation, color = Color)) +
  geom_point(size = 2) + theme_classic() +
  xlab("Mean AOO (million km2)") + ylab("Magnitude of fluctuation") +
  scale_color_identity() + scale_y_log10()  + 
  geom_hline(yintercept = 10, linetype = "dashed")


# Latitudinal and elevational plots
RasterData <- read.delim(file = "RasterData.txt", sep = "", header = T)

MP <- raster("Layers/Ov_migr_final.tif")
elev <- raster("GlobalElevationData.tif")

s <- stack(MP, elev)
v <- data.frame(na.omit(values(s)))

# Hexagon

g1 <- ggplot(RasterData, aes(y, vals)) + 
  geom_hex(bins = 100) +
  theme_classic() + 
  labs(x = "Latitude", y = "% of species showing seasonal switching") + ylim(c(0, 50)) +
  scale_fill_gradient(low = "blue", high = "red", trans = "log10") +
  stat_summary(fun.y = median, geom = "point", shape = 18, size = 4, color = "black")



g2 <- ggplot(v, aes(GlobalElevationData, Ov_migr_final)) +
  geom_hex(bins = 40) +
  theme_classic() + 
  labs(x = "Elevation (m)", y = "% of species showing seasonal switching") + ylim(c(0, 50)) +
  scale_fill_gradient(low = "blue", high = "red", trans = "log10")


cowplot::plot_grid(g1, g2, ncol = 2)  

g1 <- ggplot(RasterData, aes(y, vals)) +
  geom_hex(bins = 100) + theme_classic() + 
  labs(x = "Latitude", y = "Migratory species (%)") +
  scale_fill_gradient(low = "lightblue1", high = "darkblue", trans = "log10") + 
  ylim(c(0, 50))



g2 <- ggplot(v, aes(GlobalElevationData, Ov_Pre)) +
  geom_hex(bins = 100) + theme_classic() + 
  labs(x = "Elevation (m)", y = "Migratory species (%)") +
  scale_fill_gradient(low = "lightblue1", high = "darkblue", trans = "log10") + 
  ylim(c(0, 50))

cowplot::plot_grid(g1,g2)