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
TEST <- TRUE

#get list of allCMIP and select the one for this task
cmip_files <- list.files("raw_data/CMIP6_for_Amina", pattern = "*.nc", full.names = T)
i <- cmip_files[rep]
file <- brick(i)
file_source <- str_extract(filename(file), "[^/]*$")

#### extract over all PAs ----
#load in PAs, subset if necessary, and clean up
pas <- st_read(dsn = paste0(getwd(), "/raw_data/WDPA/WDPA_Nov2020_Public_shp/WDPA_poly_Nov2020_filtered.gdb"))
if(TEST) pas <- pas %>% filter(grepl("BRA", ISO3)) 
pas <- clean_pas(pas)

## make sure the extents match up (CMIP 6 is weird)
extent(file) <- extent(pas)

## Run through CMIP6 temperature brick and extract over the buffers
all_data <- c()
print("starting kbas")
for(j in 1:length(names(file))) {
  temp <- pas %>% dplyr::select(WDPAID, Name, DESIG, DESIG_TYPE, STATUS_YR, 
                         GOV_TYPE, ISO3) %>% 
    mutate(year = getZ(file[[j]]), source = file_source) %>% st_drop_geometry()
  
  ## make sure the layers align
  extent(file[[j]]) <- extent(pas)
  
  ev <- exact_extract(file[[j]], pas, c("mean", "min", "max"))
  
  temp <- cbind(temp, ev)
  
  all_data <- rbind(all_data, temp)
  
}
print("finished")


## now we do this with NDVI


## new we combine NDVI and temp

## now we save out 


#save this out to make my life easier
file_name <- paste0("./processed_data/pa_cmip_ovl/", file_source, ".rds")
saveRDS(all_data, file_name)

## unload all KBAs and all_data to save memory 
rm(all_data)
