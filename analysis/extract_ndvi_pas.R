## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Extract annual NDVI data over WDPA PAs
## Amina Ly, Jan 2023
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Script return annual mean temperature of PAs based on gridded CMIP6
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ifelse(dir.exists("~/Box Sync/biodiversity_cmip6"),
       setwd("~/Box Sync/biodiversity_cmip6"),
       setwd("/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6"))

source(paste0(getwd(), "/analysis/util.R"))

library(ncdf4)
library(dplyr)
library(stringr)
library(sf)
library(raster)
library(lfe)
library(lubridate)
library(exactextractr)


## test run? If so, just perform this on KBAs in Brazil (or whatever)
TEST <- FALSE

#load in PAs, subset if necessary, and clean up
ifelse(file.exists("processed_data/wdpa/clean_wdpa_terrestrial.shp"),  
       pas <- st_read(dsn = "processed_data/wdpa/clean_wdpa_terrestrial.shp", stringsAsFactors = F, crs = 4326), 
       pas <- clean_pas("raw_data/WDPA"))

if(TEST) pas <- pas %>% filter(grepl("BRA", ISO3)) 

all_data <- c()

#### extract NDVI  Data ----

## Run through CMIP6 temperature brick and extract over the buffers
print("starting NDVI")
for(y in 1981:2022) {
  
  temp <- pas %>% dplyr::select(WDPAID, NAME, DESIG, DESIG_TYPE, STATUS_YR, 
                                GOV_TYPE, ISO3) %>% 
    mutate(year = y) %>% st_drop_geometry()

  #extract and add avg mean
  file_mean <- brick(paste0(getwd(), "/processed_data/ndvi/", y, "_mean.nc"))
  ev <- exact_extract(file_mean, pas, "mean")
  temp <- cbind(temp, mean_ndvi = ev)
  
  temp2 <- pas %>% dplyr::select(WDPAID) %>% st_drop_geometry()
  
  #extract and add avg max
  file_max <- brick(paste0("processed_data/ndvi/", y, "_max.nc"))
  ev <- exact_extract(file_max, pas, "mean")
  temp2 <- cbind(temp2, max_ndvi = ev)
  
  temp <- left_join(temp, temp2, by = "WDPAID")
  
  all_data <- rbind(all_data, temp)
  
}

print("finished ndvi")

#save this out to make my life easier
file_name <- paste0("./processed_data/ndvi/ndvi_pa_ovl.csv")
write.csv(all_data, file_name)

