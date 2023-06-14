## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Extract annual CMIP6 data over KBAs 
## Amina Ly, Jan 2023
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Script return annual mean temperature of PAs based on gridded CMIP6

## NOTE if you are running this with a new KBA file, ensure file paths here AND 
## those in analysis/kba_cleaning.R are correct!!! 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

## pick up args from command line/sbatch
args <- commandArgs(trailingOnly = TRUE)
rep <- as.numeric(args[1])

## test run? If so, just perform this on KBAs in subset of locations
TEST <- TRUE
locs <- c("ZAF")
models <- c("cwd", "cdd", "txx")

#get list of allCMIP and select the one for this task
cmip_files <- list.files("raw_data/extreme-indexes-cmip6-4", pattern = paste(models, collapse = "|"), full.names = T)
i <- cmip_files[rep]
file <- brick(i)
file_source <- str_extract(filename(file), "[^/]*$")
file_source <- gsub("\\.[^.]*$", "", file_source)

#### load in cleaned KBA, and if not source file to make it ----
kba_loc <- paste0(getwd(),"/processed_data/kba/KBAsGlobal_2022_September_02_POL_noOverlaps.shp")

if(file.exists(kba_loc)) {
  kbas <- st_read(dsn = kba_loc, stringsAsFactors = F, crs = 4326) 
} else {
  source(paste0(getwd(), "/analysis/kba_cleaning.R"))
  kbas <- st_read(dsn = kba_loc, stringsAsFactors = F, crs = 4326) 
}

## clean up column names, select terrestrial, and shorten if testing
coln <- c("SitRecID", "Region", "Country", "ISO3", "NatName", "IntName", "FinCode", 
          "SitLat", "SitLong", "GISArea", "IbaStatus", "KBAStatus", 
          "AzeStatus", "AddedDate", "ChangeDate", "Source", "DelTxt",
          "DelGeom", "KBA_Quality", "Shape_Long", "Shape_Area", "LegacyKBA", "Criteria",
          "original_area", "kba_notes", "akba", "class", "geometry")
names(kbas) <- coln
if(TEST) kbas <- kbas %>% filter(ISO3 %in% locs) 

kbas <- kbas %>% filter(class == "Terrestrial")

extent(file) <- extent(kbas)

#### now extract over these cleaned up KBAs and save ----

gcm_name <- str_split(file_source, "_")[[1]][3]
msure <- str_split(file_source, "_")[[1]][1]

all_data <- c()
print("starting clean kbas")
for(j in 1:length(names(file))) {
  temp <- kbas %>% dplyr::select(SitRecID) %>% 
    mutate(date = getZ(file[[j]]), source = file_source,
           gcm  = gcm_name,
           measure = msure) %>% st_drop_geometry()
  
  ## make sure the layers align
  extent(file[[j]]) <- extent(kbas)
  
  ev <- exact_extract(file[[j]], kbas, c("mean", "min", "max"))
  
  temp <- cbind(temp, ev)
  
  all_data <- rbind(all_data, temp)
  
}
print("finished")

#save this out to make my life easier
file_name <- paste0("./processed_data/extremes_kbas/",file_source, ".csv")
write.csv(all_data, file_name)


# CODE FOR EXTRACTING ORIGINAL KBAs w/o cleaning

# #### extract over all KBAs ----
# #load in KBAs
# kba_loc <- paste0(getwd(), "/raw_data/KBA2022/KBAsGlobal_2022_September_02_POL_valid.shp") 
# 
# # if we haven't made these valid yet, do that and save it 
# if(file.exists(kba_loc)) {
#   kbas <- st_read(dsn = kba_loc, stringsAsFactors = F, crs = 4326) 
# } else {
#   kbas <- st_read(dsn = paste0(getwd(), "/raw_data/KBA2022/KBAsGlobal_2022_September_02_POL.shp") , stringsAsFactors = F, crs = 4326) 
#   if(sum(st_is_valid(kbas)) < nrow(kbas)) kbas <- st_make_valid(kbas)
#   st_write(kbas, dsn = kba_loc)
# }
# 
# if(TEST) kbas <- kbas %>% filter(ISO3 %in% c("BTN", "CHE")) 
# 
# extent(file) <- extent(kbas)
# ## Run through temperature brick and extract over the buffers
# all_data <- c()
# print("starting kbas")
# for(j in 1:length(names(file))) {
#   temp <- kbas %>% dplyr::select(SitRecID, Country, ISO3, NatName, IntName, 
#                                  SitArea, AddedDate) %>% 
#     mutate(year = getZ(file[[j]]), source = file_source)
#   
#   ## make sure the layers align
#   extent(file[[j]]) <- extent(kbas)
#   
#   ev <- exact_extract(file[[j]], kbas, c("mean", "min", "max"))
#   
#   temp <- cbind(temp, ev)
#   
#   all_data <- rbind(all_data, temp)
#   
# }
# print("finished")
# #save this out to make my life easier
# file_name <- paste0("./processed_data/originalkba_cmip_ovl/", file_source, ".rds")
# saveRDS(all_data, file_name)
# 
# ## unload all KBAs and all_data to save memory 
# rm(kbas,all_data)
# 
# #### extract over mountainous KBAs (source file to create it pulled from SDG Calculator https://github.com/GMBA-biodiversity/SDG15.4.1_Calculator) ----
# 
# 
