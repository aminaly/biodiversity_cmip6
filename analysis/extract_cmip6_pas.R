## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Extract annual CMIP6  data from WDPA PAs
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

## pick up args from commandline/sbatch
args <- commandArgs(trailingOnly = TRUE)
rep <- as.numeric(args[1])

## test run? If so, just perform this on KBAs in USA
TEST <- FALSE

#get list of allCMIP and select the one for this task
cmip_files <- list.files("raw_data/CMIP6_for_Amina", pattern = "*.nc", full.names = T)
i <- cmip_files[rep]
file <- brick(i)
file_source <- str_extract(filename(file), "[^/]*$")
file_source <- gsub("\\.[^.]*$", "", file_source)

#### extract over all PAs ----
#load in PAs, subset if necessary, and clean up
ifelse(file.exists("processed_data/wdpa/clean_wdpa_terrestrial.shp"),  
       pas <- st_read(dsn = "processed_data/wdpa/clean_wdpa_terrestrial.shp", stringsAsFactors = F, crs = 4326), 
       pas <- clean_pas("raw_data/WDPA"))

if(TEST) pas <- pas %>% filter(grepl("BRA", ISO3)) 

## make sure the extents match up (CMIP 6 is weird)
extent(file) <- extent(-180, 180, -79.21161, 81.36814)

## Run through CMIP6 temperature brick and extract over the buffers
all_data <- c()
print("starting cmip6 extract")
for(j in 1:length(names(file))) {
  temp <- pas %>% dplyr::select(WDPAID) %>% 
    mutate(year = getZ(file[[j]]), source = file_source) %>% st_drop_geometry()
  
  ## make sure the layers align
  extent(file[[j]]) <- extent(-180, 180, -79.21161, 81.36814)
  
  ev <- exact_extract(file[[j]], pas, c("mean", "min", "max"))
  
  temp <- cbind(temp, ev)
  
  all_data <- rbind(all_data, temp)
  
}
print("finished")

#rename columns for later
all_data <- all_data %>% rename(mean_temp = mean, min_temp = min, max_temp = max)

#save this out to make my life easier
file_name <- paste0(getwd(), "/processed_data/cmip6_pa_ovl/", file_source, ".csv")
write.csv(all_data, file_name)


