
### Auxillary Code to Get Rid of KBA overlaps. Sourced in extract_over_kbas

#### Part 1.1 load packages ----
# if you do not have any of the packages, please install them before running this code

library(sf)
library(dplyr)
library(tidyverse)
library(lwgeom)
library(rgdal)
sf::sf_use_s2(FALSE) ## to deal with some issues not fixable with st_make_valid

#### Part 1.2 Working Directory & Files ----

## set the working directory
ifelse(dir.exists("~/Box Sync/biodiversity_cmip6"),
       setwd("~/Box Sync/biodiversity_cmip6"),
       setwd("/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6"))
finfile <- paste0(getwd(), "/processed_data/kba/KBAsGlobal_2022_September_02_POL_noOverlaps.shp") #folder where the files per country will be saved

#read in KBAs
kbas <- st_read(dsn = paste0(getwd(), "/raw_data/KBA2022/KBAsGlobal_2022_September_02_POL.shp"), stringsAsFactors = F, crs = 4326) 
if(sum(st_is_valid(kbas)) < nrow(kbas)) kbas <- st_make_valid(kbas)

## get all KBA areas
kbas$akba <- as.numeric(suppressWarnings(tryCatch({st_area(kbas$geometry, byid = FALSE)}, error=function(e){})))
kbas$kba_notes <- ""
new_kbas <- c()

intersecs_all <- st_intersects(kbas, sparse = F)

## run through intersections
for(k in 1:nrow(intersecs_all)) {
  
  print("KBA K # is")
  print(k)
  #get this KBA
  kba <- kbas[k,]

  ## for this KBA, check to see if it intersects with any other KBA, select them all and loop through 
  inter_kbas <- kbas[which(intersecs_all[k,]==T),]
  
  ## start loop for KBA's that intersect with this one 
  for(i in 1:nrow(inter_kbas)) {
    print(i)
    intersec <- inter_kbas[i,]
    #is this the same polygon? if so, skip
    if(intersec$SitRecID == kba$SitRecID) next
    
    #intersect the two
    overlap <- st_intersection(kba, intersec)
    overlap_area <- as.numeric(st_area(overlap$geometry))
    
    #if this KBA has been cut already so now the overlap is no longer there, next
    if(nrow(overlap) == 0) next

    ## is the overlapping area > 2% of this KBA's area?
    if(0.02 < (overlap_area/ kba$akba)) {
      
      ## if the overlapping area is within 2% of the other
      if(0.98 <= (overlap_area/ kba$akba) & 1.02 >= (overlap_area/ kba$akba)) {
      
        ## if this KBA was created after it's duplicate, mark it to be removed, otherwise note it's partner
        kba$kba_notes <- ifelse(kba$AddedDate > intersec$AddedDate,
                                "remove duplicate",
                                paste(kba$kba_notes, "removed duplicate", intersec$SitRecID))
        
        } else if(kba$akba < intersec$akba) {
          
          kba_diff <- st_difference(kba, intersec)
          ## if this overlap ended up getting rid of this piece entirely, mark to be dropped
          if(nrow(kba_diff) == 0) {
            kba$kba_notes <- paste("remove -- fully overlapped")
          } else { ## else crop it
            kba <- kba_diff
            if(!st_is_valid(kba)) kba <- st_make_valid(kba) ## make sure the difference is still valid
            print("KBA difference")
            print(kba)
            kba$kba_notes <- paste(kba$kba_notes, "clipped by:", intersec$SitRecID, ";")
          }
        }
      }
    #get rid of the info from the second kba/select columns 
    kba <- kba %>% select(SitRecID, Region, Country, ISO3, NatName, IntName, FinCode,
                          SitLat, SitLong, GISArea, IbaStatus, KbaStatus, AzeStatus, 
                          AddedDate, ChangeDate, Source, DelTxt, DelGeom, KBA_Qual, 
                          Shape_Leng, Shape_Area, LegacyKBA, Criteria, geometry, 
                          akba, kba_notes)
  }
  # now we've done all the kba adjustments, add it in
  new_kbas <- rbind(kba, new_kbas)
  
}

new_kbas <- new_kbas %>% rename(original_area = akba) %>% 
  filter(!kba_notes == "remove duplicate") %>%
  filter(!kba_notes == "remove -- fully overlapped")

saveRDS(new_kbas, finfile)
st_write(new_kbas, dsn = finfile)

new_kbas$akba <- as.numeric(suppressWarnings(tryCatch({st_area(new_kbas$geometry, byid = FALSE)}, error=function(e){})))

st_write(new_kbas, dsn = finfile)

