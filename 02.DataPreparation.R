options(java.parameters = "-Xmx6g")

library(sp)
library(dplyr)
library(sf)
library(data.table)

# Converted other shapefiles in ArcGIS to save some time
# Data to shapefile
wdpa_crs <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"

dc1 <- read.csv("S3.csv", header = T)

class(dc1)
coordinates(dc1) <- ~decimalLongitude + decimalLatitude
projection(dc1) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
shapefile(dc1, "S3.shp")

dc1 <- read.csv("S4.csv", header = T)

class(dc1)
coordinates(dc1) <- ~decimalLongitude + decimalLatitude
projection(dc1) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"
shapefile(dc1, "S4.shp")


# Intersection
s4 <- st_read("Layers/S4.shp")
Nearctic <- st_read("ZoogeographicRegions/Nearctic.shp")

S4.Nearctic <- intersect(s4, Nearctic)
st_write(S4.Nearctic, "S4.Nearctic.shp")

# Shapefile to csvs
s3.pal <- st_read("S3.Palaearctic.shp")

s1.afro1 <- s1.afro %>%
  dplyr::select("species", "decimalLon", "decimalLat")

write.csv(s1.afro1, "S1.afrotropical.csv")

dc <- read.csv("EcoregionsList.csv", header = T)

species <- unique(dc$species)

for (i in species) try( {
  print(i)
  shp <- st_read(i)
  
  shp <- shp %>%
    dplyr::select("species", "decimalLon", "decimalLat")
  
  write.csv(shp, file = paste0("", i, ".csv"))
}, silent = FALSE)


# NonMigrants
df <- read.csv("NonMigrants/NonMigrants_Australia_table.csv")

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

write.csv(df1, "NonMigrants/NonMigrants.csv")
