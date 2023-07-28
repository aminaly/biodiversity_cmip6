## Functions used throughout analysis

#libraries required
library(raster)
library(sp)
library(dplyr)
library(sf)
library(dplyr)
library(tidyverse)
library(lwgeom)
library(ggplot2)
#library(transformr)
library(RColorBrewer)

#### Clean PAs -- strategy imitates official methodology used to clean PAs ----
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

#### Get cummulative KBA Coverage given output of kba_pa_ovl ----

cummulative_kba <- function(data, years = c(2022), level = "kba") {
  ## set working directory
  ifelse(dir.exists("~/Box Sync/biodiversity_cmip6"),
         setwd("~/Box Sync/biodiversity_cmip6"),
         setwd("/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6"))
  
  td <- format(Sys.Date(), "%m_%Y")
  
  fin <- paste0("./processed_data/pa_kba_ovl/cummulative_coverage_", level, min(years), max(years), ".csv")
  if(file.exists(fin)) return(read.csv(fin))
    
  ## Set NA ovl to 0
  results <- data %>% mutate(ovl = ifelse(is.na(ovl), 0, ovl))
  
  ## create dataframe of all years with all ID vars
  uniqids <- unique(results %>% dplyr::select(SitRecID, kba, Country))
  results_all_years <- merge(uniqids , c(min(results$year, na.rm = T):2022)) 
  results_all_years <- results_all_years %>% rename(year = y)
  
  ## join with full run, get rid of NA years (those with issues calculating) , 
  #and set NA ovl to 0 (no additional coverage that year)
  results_all_years <- left_join(results_all_years, results, 
                                 by = c(names(uniqids), "year"))
  results_all_years <- results_all_years %>% 
    mutate(ovl = ifelse(is.na(ovl), 0, ovl)) %>% 
    dplyr::select(-c(percPA))
  
  ## calculate cumulative coverage 
  results_all_years <- results_all_years  %>% 
    mutate(percPA = (ovl/kba) * 100) %>% 
    group_by(SitRecID, Country) %>% 
    mutate(cum_overlap = cumsum(ovl)) %>%
    mutate(cum_percPA = (cum_overlap/kba) * 100)
  
  ## select only the years necessary. If no year given, 2022
  return_data <- results_all_years %>% filter(year %in% years)
  
  ## this will be the data that creates all the other aggregations
  if(level == "kba") {
    write.csv(return_data, fin)
    return(return_data)
  } else {
    write.csv(return_data, 
              paste0("./processed_data/pa_kba_ovl/cummulative_coverage_kba", 
                     min(years), max(years), ".csv")) 
  }
  
  ## calculate cummulative country 
  results_all_years_country <- results_all_years %>% 
    group_by(year, Country) %>%
    summarize(kba_area = sum(unique(kba), na.rm = T),
              cum_overlap = sum(cum_overlap, na.rm = T)) %>%
    mutate(cum_percPA = (cum_overlap/kba_area) * 100) %>%
    filter(cum_year %in% years)
  
  write.csv(results_all_years_country, fin)
  return(results_all_years_country)
}

#### Create an input file to be used in the ML model on python (or here idk yet)  ----
## note model name should include it's extension
get_input_hist <- function(model, processed_data_loc = NA) {
  
  final_loc <- paste0(getwd(), "/processed_data/input_data/", model)
  if(is.na(processed_data_loc)) processed_data_loc <- paste0(getwd(), "/processed_data/") 
  
  ## check to see if this has been done before. If so, pick it up
  if(file.exists(final_loc)) {
    return(read.csv(final_loc))
  }
  
  ##otherwise, lets make the dataset, save it, and return it
  #pickup wd for reset later
  oldwd <- getwd()

  #read in data needed, combine, and filter
  cmip6 <- read.csv(paste0(processed_data_loc, "cmip6_pa_ovl/", model))
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

#### Gets the WDPAID. Used in an apply function in Figures.R ----
get_wdpaid <- function(index, pas) {
  pas[index,]
}



















#### Model Agreement ----
mod_agreement <- function(extreme_data, measures, sites, reps) {
  model_agreement <- c()
  ## bootstrap each site
  for(m in 1:length(measures)) {
    me <- measures[m]
    print(me)
    
    ## pull this measure
    hist <- extreme_data %>% 
      filter(measure == me,
             scenario == "historical",
             year %in% c(1995:2014))  
    comp1 <- extreme_data %>% 
      filter(measure == me,
             year %in% c(2015:2025)) %>% 
      filter(if(ERA != "") scenario == ERA)
    comp2 <- extreme_data %>% 
      filter(measure == me,
             year %in% c(2026:2036)) %>% 
      filter(if(ERA != "") scenario == ERA)
    
    ## pull this site, sample 1000 times and exit loop with 95% CI, % >0, and %<0 
    for(s in 1:length(sites)) {
      site <- sites[s]
      print(site)
      
      hist_m <- hist %>% filter(SitRecID == site)
      comp1_m <- comp1 %>% filter(SitRecID == site)
      comp2_m <- comp2 %>% filter(SitRecID == site)
      
      boots <- c()
      
      for(boot in 1:reps) {
        
        h <- sample_n(hist_m, nrow(hist_m), replace = T) %>% 
          summarize(mean_index = mean(mean, na.rm = T)) %>% pull(mean_index)
        c1 <- sample_n(comp1_m, nrow(comp1_m), replace = T) %>%
          summarize(mean_index = mean(mean, na.rm = T)) %>% pull(mean_index)
        c2 <- sample_n(comp2_m, nrow(comp2_m), replace = T) %>%
          summarize(mean_index = mean(mean, na.rm = T)) %>% pull(mean_index)
        
        if(measure == "txxETCCDI") {
          boots <- rbind(boots, cbind(firstdecade = (c1 - h),
                                      seconddecade = (c2 - h)))
        } else{
          boots <- rbind(boots, cbind(firstdecade = (c1 - h)/h,
                                      seconddecade = (c2 - h)/h))          
        } 
      }
      
      ## using the 1000 boots of all 4, return summary row for this site and add to model agreement
      low <- apply(boots, 2, quantile, c(0.025), na.rm = T)

      high <- apply(boots, 2, quantile, c(0.975), na.rm = T)
      
      mid <- apply(boots, 2, quantile, c(.5), na.rm = T)
      
      over0 <- apply(boots, 2, function(x) {sum(x > 0)/ reps})
      
      under0 <- apply(boots, 2, function(x) {sum(x < 0)/ reps})
      
      over100 <- apply(boots, 2, function(x) {sum(abs(x) > 1)/ reps})
      
      model_agreement <- rbind(model_agreement,
                               cbind(SitRecID = site, measure = me, low, high,
                                     over0, under0, over100))
      
      
    }
    
    
  }
  
  write.csv(model_agreement, paste0("./processed_data/model_agreement.csv"))

  return(model_agreement)
}