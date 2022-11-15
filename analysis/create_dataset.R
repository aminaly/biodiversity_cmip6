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

#get list of allCMIP and select the one for this task
cmip_files <- list.files("data/CMIP6_for_Amina", pattern = "*.nc", full.names = T)
i <- cmip_files[rep]
file <- brick(i)
file_source <- str_extract(filename(file), "[^/]*$")

#### extract over all KBAs ----
#load in KBA
kbas <- st_read(dsn = paste0(getwd(), "/data/KBA/KBA2020/KBAsGlobal_2020_September_02_POL.shp"), stringsAsFactors = F, crs = 4326) 
if(TEST) kbas <- kbas %>% filter(ISO3 == "USA") %>% slice_head(n = 25)

# Run through temperature brick and extract over the buffers
all_data <- c()

for(j in 1:length(names(file))) {
  
  kbas %>% dplyr::select(SitRecID, Country, ISO3, NatName, IntName, 
                         SitArea, AddedDate) %>% 
    mutate(year = getZ(file[[j]]), source = file_source)

  extracted_vals <-  exact_extract(file[[j]], kbas)
  
  temp$mean_temp <- lapply(extracted_vals, function(x){mean(as.numeric(x$value), na.rm = T)}) %>% unlist()  
  temp$max_temp <- lapply(extracted_vals, function(x){max(as.numeric(x$value), na.rm = T)}) %>% unlist()
  temp$min_temp <- lapply(extracted_vals, function(x){min(as.numeric(x$value), na.rm = T)}) %>% unlist()
  
  all_data <- bind_rows(all_data, temp)
  
}

#save this out to make my life easier
file_name <- paste0("./biodiversity_cmip6/processed_data/originalkba_cmip_ovl/", file_source, ".rds")
saveRDS(all_data, file_name)

## unload all KBAs and all_data to save memory 
rm(kbas,all_data)

#### extract over mountainous KBAs (source file to create it pulled from SDG Calculator https://github.com/GMBA-biodiversity/SDG15.4.1_Calculator) ----

#### load in cleaned KBA, and if not source file to make it ----
kba_loc <- paste0(folder,"/data/KBA/KBA2020/KBAsGlobal_2020_September_02_POL_noOverlaps.shp")

if(file.exists(kba_loc)) {
  kbas <- st_read(dsn = kba_loc, stringsAsFactors = F, crs = 4326) 
} else {
  source(paste0(folder, "/analysis/kba_cleaning.R"))
  kbas <- st_read(dsn = kba_loc, stringsAsFactors = F, crs = 4326) 
}

## clean up column names and shorten if testing
coln <- c("SitRecID", "Country", "ISO3", "NatName", "IntName", "SitArea", "IbaStatus",
          "KBAStatus", "AzeStatus", "AddedDate", "ChangeDate", "Source", "DelTxt",
          "DelGeom", "Shape_Leng", "Shape_Area", "original_area", "kba_notes", 
          "akba", "geometry")
colnames(kbas) <- coln
if(TEST) kbas <- kbas %>% filter(ISO3 == "USA") %>% slice_head(n = 25)

#### now extract over these cleaned up KBAs and save ----
all_data <- c()

for(j in 1:length(names(file))) {
  
  kbas %>% dplyr::select(SitRecID, Country, ISO3, NatName, IntName, 
                         SitArea, AddedDate) %>% 
    mutate(year = getZ(file[[j]]), source = file_source)
  
  extracted_vals <-  exact_extract(file[[j]], kbas)
  
  temp$mean_temp <- lapply(extracted_vals, function(x){mean(as.numeric(x$value), na.rm = T)}) %>% unlist()  
  temp$max_temp <- lapply(extracted_vals, function(x){max(as.numeric(x$value), na.rm = T)}) %>% unlist()
  temp$min_temp <- lapply(extracted_vals, function(x){min(as.numeric(x$value), na.rm = T)}) %>% unlist()
  
  all_data <- bind_rows(all_data, temp)
  
}

#save this out to make my life easier
file_name <- paste0("./biodiversity_cmip6/processed_data/cleankba_cmip_ovl/",file_source, ".rds")
saveRDS(all_data, file_name)

