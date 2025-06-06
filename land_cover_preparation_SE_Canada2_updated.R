
##Script name: land cover data preparation for southeastern canada
##
## Purpose of script: The primary goal of this script to generate explanatory variables on land cover data for southeastern Canada in terms of proportion of habitat per site (at a 1km X 1km spatial scale).
##
## Author: Yisa Ginath Yuh, adapted from a pipeline developed by James Paterson
##
## 
## Email: gyisa@uottawa.ca
##
## ---------------------------

##### 1. Load libraries and data ----------------------------------------------

# Load libraries
library(dplyr)
library(sf)
library(terra)

# Load the study area for masking land cover data
setwd("F:/Uottawa_data/study_area/AOI_50kmBuffer/AOI_50kmBuffer")# Set working directory to location of study area
study_area <- sf::st_read(dsn = "AOI_Buffer.shp")

# Load the land cover datasets and visualize
#SECanadaLandCover <- rast("E:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/SECanada_landcover_resampled_30m.tif")
SECanadaLandCover <- rast("F:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/LC2020_Final.tif")
plot(SECanadaLandCover)
SECanadaLandCover

# Create a base raster with some land use values for each pixel
# raster x (blank but with 1000m cells) can be used to "resample" rasters of different resolution or calculate variables
# Use CRS of the land cover data, and insert list of values based on the number of cells defined within the base raster
x <- terra::rast(resolution = 1000,
                 extent = ext(SECanadaLandCover),
                 crs = terra::crs(SECanadaLandCover),
                 vals = 1:43941744)
names(x) <- "cell_id"

ncell(x)# check number of raster cells to be inserted within the vals range of the base raster and modify the base raster values


##### 3. Use land cover raster to measure proportion of habitat per site -----------------

# Start by matching the CRS of the study area to that of the land cover data, and convert to a vector object
SECanadaVectTemp <- study_area %>%
  st_transform(., terra::crs(SECanadaLandCover)) %>%
  vect(.)

# Mask to the land cover data to the new stucy area (vector object) so that aggregate doesn't have to do so many cells
# use the terra crop and mask functions in a 2 step process
SECanadaLandCoverCropped <- terra::crop(x = SECanadaLandCover, 
                                        y = SECanadaVectTemp)
plot(SECanadaLandCoverCropped, main = "cropped")

SECanadaLandCover <- terra::mask(x = SECanadaLandCoverCropped, 
                                 mask = SECanadaVectTemp)
plot(SECanadaLandCover, main = "masked")



# Mask (if landcover is incomplete or to adjust to border of the study area)
x_lc <- x %>%
  terra::project(., terra::crs(SECanadaLandCover)) %>%
  terra::mask(., 
              mask = SECanadaVectTemp)
plot(x_lc)

# List the legend of the land cover data to ease selection of classes to be aggregated and calculate proportion of habitat per site
# land cover Legend 
# 1 = Lake
# 2 = Freshwater Pond
# 3 = Freshwater Emergent Wetland
# 4 = Freshwater Forested/Shrub Wetland
# 5 =  Riverine
# 6 = Other
# 7 = Estuarine and Marine Wetland
# 8 = Estuarine and Marine Deepwater
# 21 = Settlement
# 22 = High Reflectance Settlement
# 24 = Settlement Forest
# 25 = Roads
# 28 = Vegetated Settlement
# 29 = Very High Reflectance Settlement
# 31 = Water 
# 41 = Forest
# 42 = Forest Wetland
# 43 = Forest Regenerating after harvest < 20years
# 44 = Forest wetland regenerating after harvest <20 years
# 47 = Forest Regenerating after Harvest 20-29 years
# 48 = Forest Wetland Regenerating after Harvest 20-29 years
# 49 = Forest regenerating after fire <20 years
# 51 = Perennial Cropland (post 2015)
# 52 = Annual Cropland (post 2015)
# 55 = Land Converted to Cropland
# 56 = Land Converted to Annual Cropland
# 61 = Grassland Managed
# 62 = Grassland unmanaged
# 71 = SDLU Wetland
# 81 = Newly-Detected Settlement <10 years
# 82 = Newly Detected High Reflectance Settlement <10 years
# 84 = Newly Detected Settlement Forest <10 years
# 88 = Newly Detected Vegetated Settlement  <10 years
# 89 = Newly Detected Very High Reflectance Settlement <10 years
# 91 = Other land (rock, beaches, ice, barren land)
# 92 = Snow and Ice
# 100 = Emergent Herbaceous Wetlands

# Calculate ratio for aggregating data from smaller cells into bigger cells
Nratio <- as.integer(terra::res(x_lc)/terra::res(SECanadaLandCover))[1] 
Nratio

# -------------------proportion of habitat per site calculations--------------

# A. Native grass

# Lets start with calculating habitat per site for native grass
# Native grass has the below codes
# 61 = Grassland Managed
# 62 = Grassland unmanaged
# sum with perenial croplands (51, 55)

# Now, calculate the proportion of habitat per site for native grass
system.time(
  SENativegrass_temp <- aggregate(SECanadaLandCover,
                                  Nratio, 
                                  # Function sums the number of grassland cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                                  fun = function(x, na.rm=T) {(sum(x==61, na.rm = na.rm) + 
                                                                 sum(x==62, na.rm = na.rm) + 
                                                                 sum(x==51, na.rm = na.rm) +
                                                                 sum(x==55, na.rm = na.rm) )*10*10/(1000*1000)})
)

plot(SENativegrass_temp)

# Resample data to correct resolution and extent.
SENativegrass_resampled <- terra::resample(SENativegrass_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SENativegrass_resampled) <- "nativegrass"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SENativegrass_resampled[is.na(SENativegrass_resampled)] <- 0
SENativegrass <- terra::mask(x = SENativegrass_resampled, 
                             mask = study_area %>%
                               st_transform(., terra::crs(SENativegrass_resampled)) %>%
                               vect(.))
# Plot layer
plot(SENativegrass)

# Save results
output_dir <- "F:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/updated1/"
writeRaster(SENativegrass, 
            filename = paste0(output_dir, "native_grass_proportion_1km.tif"), 
            overwrite = TRUE)


# use the below function to clear memory if too full
#terra::tmpFiles(remove=TRUE)  # Remove temporary raster files
#rm(list=ls())  # Remove all objects from the work space
#gc()

# ---------------------------------------------------------------------------

# B. Forest

# codes = 41, 42, 43, 44, 47, 48, 49
# Lets calculate for forest using forest codes
# 41 = Forest
# 42 = Forest wetland
# 43 = Forest regenerating after harvest < 20years
# 44 = Forest wetland regenerating after harvest <20 years
# 47 = Forest Regenerating after Harvest 20-29 years
# 48 = Forest Wetland Regenerating after Harvest 20-29 years
# 49 = Forest regenerating after fire <20 years

# Calculate new categories for forest by excluding wetland categories 
# 41 = Forest
# 43 = Forest regenerating after harvest < 20years
# 47 = Forest Regenerating after Harvest 20-29 years
# 49 = Forest regenerating after fire <20 years


# Now, calculate the proportion of habitat per site for forest
#system.time(
#SEforest_temp <- aggregate(SECanadaLandCover,
#Nratio, 
                             # Function sums the number of forest cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
#fun = function(x, na.rm=T) {(sum(x==41, na.rm = na.rm) + 
#sum(x==42, na.rm = na.rm)+
#sum(x==43, na.rm = na.rm)+
#sum(x==44, na.rm = na.rm)+
#sum(x==47, na.rm = na.rm)+
#sum(x==48, na.rm = na.rm)+
#sum(x==49, na.rm = na.rm))*10*10/(1000*1000)})
#)


system.time(
  SEforest_temp <- aggregate(SECanadaLandCover,
                             Nratio, 
                             # Function sums the number of forest cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                             fun = function(x, na.rm=T) {(sum(x==41, na.rm = na.rm) + 
                                                            
                                                            sum(x==43, na.rm = na.rm)+
                                                           
                                                            sum(x==47, na.rm = na.rm)+
                                                            
                                                            sum(x==49, na.rm = na.rm))*10*10/(1000*1000)})
)

plot(SEforest_temp)

# Resample data to correct resolution and extent.
SEforest_resampled <- terra::resample(SEforest_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SEforest_resampled) <- "forest"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SEforest_resampled[is.na(SEforest_resampled)] <- 0
SEforest <- terra::mask(x = SEforest_resampled, 
                        mask = study_area %>%
                          st_transform(., terra::crs(SEforest_resampled)) %>%
                          vect(.))
# Plot layer
plot(SEforest)

# Save results
output_dir <- "F:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/updated1/"
writeRaster(SEforest, 
            filename = paste0(output_dir, "forest_proportion_1km.tif"), 
            overwrite = TRUE)

# ---------------------------------------------------------------------------

# C. Wetlands
# categories are 2, 3, 4, 7, 71, 100
# wetland legends
# 2 = Freshwater Pond
# 3 = Freshwater Emergent Wetland
# 4 = Freshwater Forested/Shrub Wetland
# 7 = Estuarine and Marine Wetland
# 71 = SDLU Wetland
# 100 = Emergent Herbaceous Wetlands


# Third, calculate the proportion of habitat per site for wetlands
#system.time(
#SEwetlands_temp <- aggregate(SECanadaLandCover,
#Nratio, 
                               # Function sums the number of wetland cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
#fun = function(x, na.rm=T) {(sum(x==2, na.rm = na.rm) + 
#sum(x==3, na.rm = na.rm)+
#sum(x==4, na.rm = na.rm)+
#sum(x==7, na.rm = na.rm)+
#sum(x==71, na.rm = na.rm)+
#sum(x==100, na.rm = na.rm))*10*10/(1000*1000)})
#)


# C. Wetlands updated
# categories are  71, 100, 101
# wetland legends
# 71 = SDLU Wetland
# 100 = CWI Freshwater Wetland (omitting CWI Shallow Open Water category)
# 101 = CWI Saltwater/Brackish Wetlands
# Add forested wetland categories (42, 44, 48)

# Third, calculate the proportion of habitat per site for wetlands
system.time(
  SEwetlands_temp <- aggregate(SECanadaLandCover,
                               Nratio, 
                               # Function sums the number of wetland cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                               fun = function(x, na.rm=T) {(sum(x==71, na.rm = na.rm) + 
                                                              sum(x==100, na.rm = na.rm)+
                                                              sum(x==101, na.rm = na.rm) +
                                                              sum(x==42, na.rm = na.rm) +
                                                              sum(x==44, na.rm = na.rm) +
                                                              sum(x==48, na.rm = na.rm))*10*10/(1000*1000)})
)


plot(SEwetlands_temp)

# Resample data to correct resolution and extent.
SEwetlands_resampled <- terra::resample(SEwetlands_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SEwetlands_resampled) <- "wetlands"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SEwetlands_resampled[is.na(SEwetlands_resampled)] <- 0
SEwetlands <- terra::mask(x = SEwetlands_resampled, 
                          mask = study_area %>%
                            st_transform(., terra::crs(SEwetlands_resampled)) %>%
                            vect(.))
# Plot layer
plot(SEwetlands)

# Save results
output_dir <- "F:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/updated1/"
writeRaster(SEwetlands, 
            filename = paste0(output_dir, "wetlands_proportion_1km.tif"), 
            overwrite = TRUE)

# ---------------------------------------------------------------------------

# D. Water
# Water categories = 1, 5, 8, 31
# 1 = Lake
# 5 =  Riverine
# 8 = Estuarine and Marine Deepwater
# 31 = Water 

# Third, calculate the proportion of habitat per site for wetlands
system.time(
  SEwater_temp <- aggregate(SECanadaLandCover,
                            Nratio, 
                            # Function sums the number of water body cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                            fun = function(x, na.rm=T) {(sum(x==1, na.rm = na.rm) + 
                                                           sum(x==5, na.rm = na.rm)+
                                                           sum(x==8, na.rm = na.rm)+
                                                           sum(x==31, na.rm = na.rm))*10*10/(1000*1000)})
)


# D. Fresh Water
# Water categories = 1
# 1 = Fresh water

# Third, calculate the proportion of habitat per site for wetlands
system.time(
  SEwater_temp <- aggregate(SECanadaLandCover,
                            Nratio, 
                            # Function sums the number of water body cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                            fun = function(x, na.rm=T) {(sum(x==1, na.rm = na.rm))*10*10/(1000*1000)})
)

plot(SEwater_temp)

# Resample data to correct resolution and extent.
SEwater_resampled <- terra::resample(SEwater_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SEwater_resampled) <- "water"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SEwater_resampled[is.na(SEwater_resampled)] <- 0
SEwater <- terra::mask(x = SEwater_resampled, 
                       mask = study_area %>%
                         st_transform(., terra::crs(SEwater_resampled)) %>%
                         vect(.))
# Plot layer
plot(SEwater)

# Save results
output_dir <- "E:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/"
writeRaster(SEwater, 
            filename = paste0(output_dir, "freshwater_proportion_1km.tif"), 
            overwrite = TRUE)



# D. Ocean Water
# Water categories = 103
# 103 = Ocean water

# Third, calculate the proportion of habitat per site for wetlands
system.time(
  SEwater_temp <- aggregate(SECanadaLandCover,
                            Nratio, 
                            # Function sums the number of water body cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                            fun = function(x, na.rm=T) {(sum(x==103, na.rm = na.rm))*10*10/(1000*1000)})
)

plot(SEwater_temp)

# Resample data to correct resolution and extent.
SEwater_resampled <- terra::resample(SEwater_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SEwater_resampled) <- "water"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SEwater_resampled[is.na(SEwater_resampled)] <- 0
SEwater <- terra::mask(x = SEwater_resampled, 
                       mask = study_area %>%
                         st_transform(., terra::crs(SEwater_resampled)) %>%
                         vect(.))
# Plot layer
plot(SEwater)

# Save results
output_dir <- "E:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/"
writeRaster(SEwater, 
            filename = paste0(output_dir, "oceanwater_proportion_1km.tif"), 
            overwrite = TRUE)

# ---------------------------------------------------------------------------

# E. Settlements
# settlements codes = 21, 22, 24, 28, 29, 81, 82, 84, 88, 89

# 21 = Settlement
# 22 = High Reflectance Settlement
# 24 = Settlement Forest
# 28 = Vegetated Settlement
# 29 = Very High Reflectance Settlement
# 81 = Newly-Detected Settlement <10 years
# 82 = Newly Detected High Reflectance Settlement <10 years
# 84 = Newly Detected Settlement Forest <10 years
# 88 = Newly Detected Vegetated Settlement  <10 years
# 89 = Newly Detected Very High Reflectance Settlement <10 years

# Now, lets calculate the proportion of habitat per site for settlements
#system.time(
#SEsettlements_temp <- aggregate(SECanadaLandCover,
# Nratio, 
                                  # Function sums the number of forest cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
#fun = function(x, na.rm=T) {(sum(x==21, na.rm = na.rm) + 
#  sum(x==22, na.rm = na.rm)+
# sum(x==24, na.rm = na.rm)+
# sum(x==28, na.rm = na.rm)+
# sum(x==29, na.rm = na.rm)+
#sum(x==81, na.rm = na.rm)+
# sum(x==82, na.rm = na.rm)+
#sum(x==84, na.rm = na.rm)+
# sum(x==88, na.rm = na.rm)+
#sum(x==89, na.rm = na.rm))*10*10/(1000*1000)})
#)

#plot(SEsettlements_temp)



#calculate the proportion of habitat per site for settlements without vegetation
# (21, 22, 29, 81, 82, 89)

system.time(
  SEsettlements_temp <- aggregate(SECanadaLandCover,
                                  Nratio, 
                                  # Function sums the number of forest cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                                  fun = function(x, na.rm=T) {(sum(x==21, na.rm = na.rm) + 
                                                                 sum(x==22, na.rm = na.rm)+
                                                                 
                                                                 sum(x==29, na.rm = na.rm)+
                                                                 sum(x==81, na.rm = na.rm)+
                                                                 sum(x==82, na.rm = na.rm)+
                                                                
                                                                 sum(x==89, na.rm = na.rm))*10*10/(1000*1000)})
)

plot(SEsettlements_temp)


# Resample data to correct resolution and extent.
SEsettlements_resampled <- terra::resample(SEsettlements_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SEsettlements_resampled) <- "settlements"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SEsettlements_resampled[is.na(SEsettlements_resampled)] <- 0
SEsettlements <- terra::mask(x = SEsettlements_resampled, 
                             mask = study_area %>%
                               st_transform(., terra::crs(SEsettlements_resampled)) %>%
                               vect(.))
# Plot layer
plot(SEsettlements)

# Save results
#output_dir <- "E:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/"
writeRaster(SEsettlements, 
            filename = paste0(output_dir, "settlements_proportion_1km.tif"), 
            overwrite = TRUE)

# Now, do for vegetated settlements (24, 28, 84, 88)

system.time(
  veg_SEsettlements_temp <- aggregate(SECanadaLandCover,
                                  Nratio, 
                                  # Function sums the number of forest cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                                  fun = function(x, na.rm=T) {(sum(x==24, na.rm = na.rm) + 
                                                                 sum(x==28, na.rm = na.rm)+
                                                                 
                                                                 sum(x==84, na.rm = na.rm)+
                                                                 
                                                                 sum(x==88, na.rm = na.rm))*10*10/(1000*1000)})
)

plot(veg_SEsettlements_temp)


# Resample data to correct resolution and extent.
veg_SEsettlements_resampled <- terra::resample(veg_SEsettlements_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(veg_SEsettlements_resampled) <- "vegetated_settlements"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
veg_SEsettlements_resampled[is.na(veg_SEsettlements_resampled)] <- 0
veg_SEsettlements <- terra::mask(x = veg_SEsettlements_resampled, 
                             mask = study_area %>%
                               st_transform(., terra::crs(veg_SEsettlements_resampled)) %>%
                               vect(.))
# Plot layer
plot(veg_SEsettlements)

# Save results
#output_dir <- "E:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/"
writeRaster(veg_SEsettlements, 
            filename = paste0(output_dir, "vegetated_settlements_proportion_1km.tif"), 
            overwrite = TRUE)

# ---------------------------------------------------------------------------

# F. Annual croplands
# cropland codes = 52, 55, 56
# 52 = Annual Cropland (post 2015)
# 55 = Land Converted to Cropland
# 56 = Land Converted to Annual Cropland

# New data = 52, 56
# Calculate the proportion of habitat per site for annual croplands
system.time(
  SEannual_croplands_temp <- aggregate(SECanadaLandCover,
                                       Nratio, 
                                       # Function sums the number of annual cropland cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                                       fun = function(x, na.rm=T) {(sum(x==52, na.rm = na.rm) + 
                                                                      
                                                                      sum(x==56, na.rm = na.rm))*10*10/(1000*1000)})
)

plot(SEannual_croplands_temp)

# Resample data to correct resolution and extent.
SEannual_croplands_resampled <- terra::resample(SEannual_croplands_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SEannual_croplands_resampled) <- "annual croplands"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SEannual_croplands_resampled[is.na(SEannual_croplands_resampled)] <- 0
SEannual_croplands <- terra::mask(x = SEannual_croplands_resampled, 
                                  mask = study_area %>%
                                    st_transform(., terra::crs(SEannual_croplands_resampled)) %>%
                                    vect(.))
# Plot layer
plot(SEannual_croplands)

# Save results
#output_dir <- "E:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/"
writeRaster(SEannual_croplands, 
            filename = paste0(output_dir, "annual_croplands_proportion_1km.tif"), 
            overwrite = TRUE)

# ---------------------------------------------------------------------------

# G. perennial croplands
# Perennial cropland (category 51)
# 51 = Perennial Cropland (post 2015)

# Calculate the proportion of habitat per site for perennial croplands
system.time(
  SEperennial_croplands_temp <- aggregate(SECanadaLandCover,
                                          Nratio, 
                                          # Function sums the number of annual cropland cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                                          fun = function(x, na.rm=T) {(sum(x==51, na.rm = na.rm))*10*10/(1000*1000)})
)

plot(SEperennial_croplands_temp)

# Resample data to correct resolution and extent.
SEperennial_croplands_resampled <- terra::resample(SEperennial_croplands_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SEperennial_croplands_resampled) <- "perennial croplands"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SEperennial_croplands_resampled[is.na(SEperennial_croplands_resampled)] <- 0
SEperennial_croplands <- terra::mask(x = SEperennial_croplands_resampled, 
                                     mask = study_area %>%
                                       st_transform(., terra::crs(SEperennial_croplands_resampled)) %>%
                                       vect(.))
# Plot layer
plot(SEperennial_croplands)

# Save results
output_dir <- "E:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/"
writeRaster(SEperennial_croplands, 
            filename = paste0(output_dir, "perennial_croplands_proportion_1km.tif"), 
            overwrite = TRUE)

# ---------------------------------------------------------------------------

# H. Barren areas
# Barren land (category 91)
# 91 = Other land (rock, beaches, ice, barren land)

# Calculate the proportion of habitat per site for barren lands
system.time(
  SEbarren_lands_temp <- aggregate(SECanadaLandCover,
                                   Nratio, 
                                   # Function sums the number of barren land cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                                   fun = function(x, na.rm=T) {(sum(x==91, na.rm = na.rm))*10*10/(1000*1000)})
)

plot(SEbarren_lands_temp)

# Resample data to correct resolution and extent.
SEbarren_lands_resampled <- terra::resample(SEbarren_lands_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SEbarren_lands_resampled) <- "barren lands"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SEbarren_lands_resampled[is.na(SEbarren_lands_resampled)] <- 0
SEbarren_lands <- terra::mask(x = SEbarren_lands_resampled, 
                              mask = study_area %>%
                                st_transform(., terra::crs(SEbarren_lands_resampled)) %>%
                                vect(.))
# Plot layer
plot(SEbarren_lands)

# Save results
output_dir <- "E:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/"
writeRaster(SEbarren_lands, 
            filename = paste0(output_dir, "barren_lands_proportion_1km.tif"), 
            overwrite = TRUE)

# ---------------------------------------------------------------------------

# I. Roads
# roads (category 25)

# Calculate the proportion of habitat per site for roads
system.time(
  SEroads_temp <- aggregate(SECanadaLandCover,
                            Nratio, 
                            # Function sums the number of road cells x cell area (10 x 10) / area of environmental layer raster cells (1000 x 1000)
                            fun = function(x, na.rm=T) {(sum(x==25, na.rm = na.rm))*10*10/(1000*1000)})
)

plot(SEroads_temp)

# Resample data to correct resolution and extent.
SEroads_resampled <- terra::resample(SEroads_temp, x_lc, method = "bilinear") 

# Before saving, make the name of the layer correspond to the variable
names(SEroads_resampled) <- "roads"

# Need to mask landcover resampled rasters (so that 0 within mask, NA outside)
# Make NA values 0, then re-mask by study area outline
SEroads_resampled[is.na(SEroads_resampled)] <- 0
SEroads <- terra::mask(x = SEroads_resampled, 
                       mask = study_area %>%
                         st_transform(., terra::crs(SEroads_resampled)) %>%
                         vect(.))
# Plot layer
plot(SEroads)

# Save results
output_dir <- "E:/Uottawa_data/study_area/land_cover_data/Eastern Expansion/resampled_data2/"
writeRaster(SEroads, 
            filename = paste0(output_dir, "roads_proportion_1km.tif"), 
            overwrite = TRUE)

# --------------------------------End---------------------------------------

