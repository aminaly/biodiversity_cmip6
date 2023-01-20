## Functions used throughout analysis

#libraries required
library(raster)
library(sp)
library(dplyr)

## Clean PAs -- strategy imitates official methodology used to clean PAs
## for SDG 15.4.1 (without ISO3 edits)
clean_pas <- function(location) {
  
  ## read in 3 WDPA files
  p1 <- st_read(dsn = paste0(location, "/WDPA_Jan2023_Public_shp/WDPA_Jan2023_Public_shp_0/WDPA_Jan2023_Public_shp-polygons.shp"), stringsAsFactors = F, crs = 4326)
  p2 <- st_read(dsn = paste0(location, "/WDPA_Jan2023_Public_shp/WDPA_Jan2023_Public_shp_1/WDPA_Jan2023_Public_shp-polygons.shp"), stringsAsFactors = F, crs = 4326)
  p3 <- st_read(dsn = paste0(location, "/WDPA_Jan2023_Public_shp/WDPA_Jan2023_Public_shp_2/WDPA_Jan2023_Public_shp-polygons.shp"), stringsAsFactors = F, crs = 4326)
  
  pas <- dplyr::bind_rows(p1, p2, p3)
  
  ## select only terrestrial PAs for our specific analysis
  pas <- pas %>% filter(MARINE == 0)
  
  ## remove sites not typically used (proposed, no designation year, and no UNESCO Man and Biosphere Reserves)
  pas %>% filter(!STATUS %in% c("Proposed", "Not Reported"), STATUS_YR == 0,
                 !grepl("UNESCO-MAB", DESIG))
  
  #clean up geometries
  pas <- pas[!is.na(st_dimension(pas)),]
  as.character(unique(st_geometry_type(st_geometry(pas)))) ## what geometries are in the dataset
  
  #check for and repair any geometry issues
  if(sum(st_is_valid(pas)) < nrow(pas)) pas <- st_make_valid(pas)
  
  ## convert factors to characters in the dataframes
  ## PAs dataframe
  pas$ISO3 <- as.character(pas$ISO3)
  pas$PARENT_ISO <- as.character(pas$PARENT_ISO)
  str(pas)
  
  ## make sure all unassigned PAs are the same (convert to NAs)
  pas <- pas %>% mutate(ISO3 = ifelse((ISO3 == " " | ISO3 == '---'), NA, ISO3))
  
  ## save out for the future
  st_write(pas, "processed_data/WDPA/clean_wdpa_terrestrial.shp")

  return(pas)
  
}

## creates an input file to be used in the ML model on python (or here idk yet) 
## note model name should have .nc but not .csv
get_input_hist <- function(model, processed_data_loc = NA) {
  
  final_loc <- paste0(getwd(), "/processed_data/input_data/", model, ".csv")
  if(is.na(processed_data_loc)) processed_data_loc <- paste0(getwd(), "/processed_data/") 
  
  ## check to see if this has been done before. If so, pick it up
  if(file.exists(final_loc)) {
    return(read.csv(final_loc))
  }
  
  ##otherwise, lets make the dataset, save it, and return it
  #pickup wd for reset later
  oldwd <- getwd()

  #read in data needed, combine, and filter
  cmip6 <- read.csv(paste0(processed_data_loc, "cmip6_pa_ovl/", model, ".csv"))
  ndvi <- read_csv(paste0(processed_data_loc, "ndvi/ndvi_pa_ovl.csv"))
  
  #combine
  data <- left_join(ndvi[,-1], cmip6[,-1], by = c("WDPAID", "year"))
  
  #remove any years where there is a mismatch
  data <- data %>% filter(!is.na(source))
  write.csv(data, final_loc, row.names =  F)
  
  #reset wd
  setwd(oldwd)
  
  return(data)
}

















