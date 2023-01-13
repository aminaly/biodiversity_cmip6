## Functions used throughout analysis


## Clean PAs -- strategy imitates official methodology used to clean PAs
## for SDG 15.4.1 (without ISO3 edits)
clean_pas <- function(pas) {
  
  # select only terrestrial PAs 
  pas <- pas %>% filter(MARINE == 0)
  
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
  
  return(pas)
  
}
