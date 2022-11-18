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
TEST <- FALSE

#get list of allCMIP and select the one for this task
cmip_files <- list.files("raw_data/CMIP6_for_Amina", pattern = "*.nc", full.names = T)
i <- cmip_files[rep]
file <- brick(i)
file_source <- str_extract(filename(file), "[^/]*$")

#### extract over all KBAs ----
#load in KBA
kbas <- st_read(dsn = paste0(getwd(), "/raw_data/KBA2020/KBAsGlobal_2020_September_02_POL_valid.shp"), stringsAsFactors = F) 
if(TEST) kbas <- kbas %>% filter(ISO3 == "USA") %>% slice_head(n = 25)
if(sum(st_is_valid(kbas)) < nrow(kbas)) kbas <- st_make_valid(kbas)

# Run through temperature brick and extract over the buffers
all_data <- c()
print("starting kbas")
for(j in 1:length(names(file))) {
  temp <- kbas %>% dplyr::select(SitRecID, Country, ISO3, NatName, IntName, 
                         SitArea, AddedDate) %>% 
    mutate(year = getZ(file[[j]]), source = file_source)
  
  temp <- temp %>%
    mutate(mean_temp = raster::extract(file[[j]], temp, fun = mean)) %>%
    mutate(max_temp = raster::extract(file[[j]], temp, fun = max)) %>%
    mutate(min_temp = raster::extract(file[[j]], temp, fun = min))
  
  all_data <- rbind(all_data, temp)
  
}
print("finished")
#save this out to make my life easier
file_name <- paste0("./biodiversity_cmip6/processed_data/originalkba_cmip_ovl/", file_source, ".rds")
saveRDS(all_data, file_name)

## unload all KBAs and all_data to save memory 
rm(kbas,all_data)

#### extract over mountainous KBAs (source file to create it pulled from SDG Calculator https://github.com/GMBA-biodiversity/SDG15.4.1_Calculator) ----

#### load in cleaned KBA, and if not source file to make it ----
kba_loc <- paste0(getwd(),"/raw_data/KBA2020/KBAsGlobal_2020_September_02_POL_noOverlaps.shp")

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
print("starting fied kbas")
for(j in 1:length(names(file))) {
  temp <- kbas %>% dplyr::select(SitRecID, Country, ISO3, NatName, IntName, 
                                 SitArea, AddedDate) %>% 
    mutate(year = getZ(file[[j]]), source = file_source)
  
  temp <- temp %>%
    mutate(mean_temp = raster::extract(file[[j]], temp, fun = mean)) %>%
    mutate(max_temp = raster::extract(file[[j]], temp, fun = max)) %>%
    mutate(min_temp = raster::extract(file[[j]], temp, fun = min))
  
  all_data <- rbind(all_data, temp)
  
}
print("finished")

#save this out to make my life easier
file_name <- paste0("./biodiversity_cmip6/processed_data/cleankba_cmip_ovl/",file_source, ".rds")
saveRDS(all_data, file_name)

