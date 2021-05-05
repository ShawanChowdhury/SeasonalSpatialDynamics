options(java.parameters = "-Xmx6g")

library(sp)
library(dplyr)
library(sf)
library(countrycode)
library(CoordinateCleaner)
library(maptools)
library(rgeos)
library(data.table)
library(rworldmap)

## GBIF DOI
# https://doi.org/10.15468/dl.3b9m6s; https://doi.org/10.15468/dl.vxvfvj;
# https://doi.org/10.15468/dl.q8qjw8; https://doi.org/10.15468/dl.rmb27u;
# https://doi.org/10.15468/dl.83q9a6; https://doi.org/10.15468/dl.xt74fw;
# https://doi.org/10.15468/dl.spkfhu

df1 <- read_delim(file = 'occ/1/occurrence.txt', delim = "\t")
df1 <- df1 %>%
  dplyr::select(species, scientificName, verbatimScientificName, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year, month,
                basisOfRecord, institutionCode, datasetName) 
df2 <- read_delim(file = 'occ/2/occurrence.txt', delim = "\t")
df2 <- df1 %>%
  dplyr::select(species, scientificName, verbatimScientificName, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year, month,
                basisOfRecord, institutionCode, datasetName)

df3 <- read_delim(file = 'occ/3/occurrence.txt', delim = "\t")
df3 <- df3 %>%
  dplyr::select(species, scientificName, verbatimScientificName, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year, month,
                basisOfRecord, institutionCode, datasetName)

df4 <- read_delim(file = 'occ/4/occurrence.txt', delim = "\t")
df4 <- df4 %>%
  dplyr::select(species, scientificName, verbatimScientificName, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year, month,
                basisOfRecord, institutionCode, datasetName)


df5 <- read_delim(file = 'occ/5/occurrence.txt', delim = "\t")
df5 <- df5 %>%
  dplyr::select(species, scientificName, verbatimScientificName, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year, month,
                basisOfRecord, institutionCode, datasetName)

df6 <- read_delim(file = 'occ/6/occurrence.txt', delim = "\t")
df6 <- df6 %>%
  dplyr::select(species, scientificName, verbatimScientificName, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year, month,
                basisOfRecord, institutionCode, datasetName)

df7 <- read_delim(file = 'occ/7/occurrence.txt', delim = "\t")
df7 <- df7 %>%
  dplyr::select(species, scientificName, verbatimScientificName, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year, month,
                basisOfRecord, institutionCode, datasetName)

df <- rbind(df1, df2, df3, df4, df5, df6, df7)
write.csv(df, "SecondChapterCompleteData.11thMay2020.csv")

rm(df1, df2, df3, df4, df5, df6, df7)


# Removing records without coordinates
df <- df%>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# Removing blank cells
df <- df[!(is.na(df$species) | df$species == ""),]
df <- df[!(is.na(df$decimalLongitude) | df$decimalLongitude == ""),]
df <- df[!(is.na(df$decimalLatitude) | df$decimalLatitude == ""),]
df <- df[!(is.na(df$month) | df$month == ""),]

# Removing duplicated records
df1 <- df[!duplicated(df),]

rm(df)

# Removing species with low occurrence records
df1 <- df1 %>%
  group_by(species) %>%
  filter(n() > 3) %>%
  ungroup()

# Removing invalid and imprecise coordinates
# Convert country code from ISO2c to ISO3c
df1$countryCode <-  countrycode(df1$countryCode, origin =  'iso2c', destination = 'iso3c')

# Flag problems
dc1 <- data.frame(df1)
flags <- clean_coordinates(x = dc1, lon = "decimalLongitude", lat = "decimalLatitude",
                           countries = "countryCode", 
                           species = "species",
                           tests = c("capitals", "centroids", "gbif", "institutions",
                                     "zeros", "countries")) # most test are on by default

summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

# Exclude problematic records
dc_cl <- dc1[flags$.summary,]


# Select columns of interest
dc_cl <- dc_cl1 %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, month)

# Seasonal distribution
dc_cl <- read.csv("CleanedRecords.csv", header = T)

s1 <- dc_cl %>%
  filter(month == c(12, 1, 2))
write.csv(s1, "S1.csv")

s2 <- dc_cl %>%
  filter(month == c(3, 4, 5))
write.csv(s2, "S2.csv")

s3 <- dc_cl %>%
  filter(month == c(6, 7, 8))
write.csv(s3, "S3.csv")

s4 <- dc_cl %>%
  filter(month == c(9, 10, 11))
write.csv(s4, "S4.csv")