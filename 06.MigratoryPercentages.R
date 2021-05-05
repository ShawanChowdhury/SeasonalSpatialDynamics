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

tifs <- list.files(path = "Rasters_Global",pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
tifs <- tifs[stringr::str_detect(tifs, "bin")]


SpeciesList <- read.csv("SpeciesList_updated.csv", header = T)

sp <- unique(SpeciesList$species)

for (i in sp) try( {
  
  print(i)
  speciesname <- gsub(" ", "_", i)
  r <- tifs[stringr::str_detect(tifs, as.character(i))]
  
  raster_files <- r[stringr::str_detect(r, "bin")]
  
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
   
    
    writeRaster(r_migr_bin, paste0("", speciesname, "_Migr.tif"), NAflag=-9999, overwrite = TRUE)
    
    
    results <- table(getValues(r_migr_bin))
    
    ModelOutput <-
      as.data.frame(results) %>%
      mutate(species = speciesname) %>%
      dplyr::select(species, everything())
    
    write.csv(ModelOutput, file = paste0("MP_", speciesname, ".csv"))
    
  }
  
}, silent = FALSE)



# Merge csvs
input_folder <- "MPR (3groups)"
range <- dir(input_folder, "^.*\\.csv$", full.names = TRUE)
migr.per <- plyr::ldply(range, readr::read_csv)
write.csv(migr.per, "MPR.csv")

# Calculation
mp <- read.csv("MPR.csv", header = T)

migr.per <- mp %>%
  group_by(species) %>%
  summarise(TotalCells = sum(Freq), migratory = sum(migr), mpr = (migratory/TotalCells)*100) %>%
  ungroup()

write.csv(migr.per, "MPR_calc.csv")


# Migratory Cell Figure
MigrCell <- read.csv("MPR_calc.csv", header = T)

library(ggridges)
library(hrbrthemes)

# Frequency distribution

MigrCell$family <- factor(MigrCell$family, levels=c("Hesperiidae", "Lycaenidae",
                                                    "Nymphalidae", "Papilionidae", 
                                                    "Pieridae", "Overall"))
ggplot(MigrCell, aes(mpr, fill = family)) +
  geom_histogram(binwidth = 5, colour = "black") +
  facet_wrap(~family) + theme_classic() + 
  theme(legend.position = "none") +
  xlab("Seasonal switching (%)") +
  ylab("No. of species") + geom_density(alpha=0.6)



# Density Plot
MigrCell <- read.csv("MPR_calc.csv", header = T)

ggplot(data = MigrCell, aes(x=migr.per, group=family, fill=family)) +
  geom_density(adjust=1.5) +
  theme_classic() +facet_wrap(~family) +
  theme(
    legend.position = "none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x = element_blank()
  ) + xlab("Migratory percentages") + ylab("Density")


ggplot(MigrCell, aes(x = migr.per, y = family, fill = ..x..)) +
  geom_point() +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01) +
  scale_fill_viridis() +
  labs(title = 'Migratory Cell (%)') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.2, "lines"),
    strip.text.x = element_text(size = 8)
  )


# MPR (2 groups) [positive] [repeat this process for the negatives]
tifs <- list.files(path = "MPR (3groups)",pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
tifs <- tifs[stringr::str_detect(tifs, "Migr")]


SpeciesList <- read.csv("SpeciesList_updated.csv", header = T)

sp <- unique(SpeciesList$species)

for (i in sp) try({
  print(i)
  speciesname <- gsub(" ", "_", i)
  r_migr_bin <- tifs[stringr::str_detect(tifs, as.character(i))]
  
  r_migr_bin <- raster(r_migr_bin)
  
  #n_seasons <- max(r_migr_bin[], na.rm = T)
  
  r_migr_bin[r_migr_bin == 0] <- 0
  r_migr_bin[r_migr_bin > 0] <- 1 # 0 for the negatives
  r_migr_bin[r_migr_bin < 0] <- 0 # -1 for the negatives
  
  #r_migr_bin[r_migr_bin %in% 1:(n_seasons - 1)] <- 1
  #r_migr_bin[r_migr_bin == n_seasons] <- 0
  
  writeRaster(r_migr_bin, paste0("", speciesname, "_Migr_Pr.tif"), NAflag=-9999)
  
}, silent = FALSE)


# In parallel
n_threads <- 8

tifs <- list.files(path = "MPR (3groups)",pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
tifs <- tifs[stringr::str_detect(tifs, "Migr")]


SpeciesList <- read.csv("SpeciesList_updated.csv", header = T)

head(SpeciesList)
sp <- unique(SpeciesList$species)

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
clusterExport(cl, c("tifs"))

# Main processing
result <- parLapply(cl, sp, function(i) {
  print(i)
  speciesname <- gsub(" ", "_", i)
  r_migr_bin <- tifs[stringr::str_detect(tifs, as.character(i))]
  
  if(length(r_migr_bin) > 0) {
    r_migr_bin <- raster(r_migr_bin)
    
    #n_seasons <- max(r_migr_bin[], na.rm = T)
    
    r_migr_bin[r_migr_bin == 0] <- 0
    r_migr_bin[r_migr_bin > 0] <- 0 # 1 for the positives
    r_migr_bin[r_migr_bin < 0] <- -1 # 0 for the positives
    
    writeRaster(r_migr_bin, paste0("", speciesname, "_Migr_Ab.tif"), 
                NAflag=-9999, overwrite = TRUE)
    
  }
  
  
    
  })

# Stop cluster
cl <- stopCluster(cl)