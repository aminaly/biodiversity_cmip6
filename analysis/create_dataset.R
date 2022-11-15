ifelse(dir.exists("~/Box Sync/biodiversity_cmip6"),
       setwd("~/Box Sync/biodiversity_cmip6"),
       setwd("/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6"))

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
TEST = TRUE

#get list of allCMIP
cmip_files <- list.files("data/CMIP6_for_Amina", pattern = "*.nc", full.names = T)

#load in KBA
kbas <- st_read(dsn = paste0(getwd(), "/data/KBA/KBA2020/KBAsGlobal_2020_September_02_POL.shp"), stringsAsFactors = F, crs = 4326) 
if(TEST) kbas <- kbas %>% filter(ISO3 == "USA") %>% slice_head(n = 25)

# Run through temperature brick and extract over the buffers
all_data <- c()

i <- cmip_files[rep]

print(i)
file <- brick(i)
crs(file) <- CRS("+init=epsg:4326")

for(j in 1:length(names(file))) {
  temp <- c()
  nms <- as.numeric(substring(as.character(names(file[[j]])),2))
  temp$date <- rep(as.Date(nms, origin= "1900-01-01"), nrow(block_group))
  temp <- as.data.frame(temp)
  temp$county <- block_group$COUNTYFP
  temp$census_block_group <- block_group$GEOID
  temp$fips <- block_group$fips
  temp$measure <- rep(substring(i, 26, 29), nrow(block_group))
  temp$year <- year(temp$date)
  
  
  extracted_vals <-  exact_extract(file[[j]], block_group)
  
  temp$mean_measure <- lapply(extracted_vals, function(x){mean(as.numeric(x$value), na.rm = T)}) %>% unlist()  
  temp$max_measure <- lapply(extracted_vals, function(x){max(as.numeric(x$value), na.rm = T)}) %>% unlist()
  
  all_data <- bind_rows(all_data, temp)
  
}

#save this out to make my life easier
saveRDS(all_data, file_name)

