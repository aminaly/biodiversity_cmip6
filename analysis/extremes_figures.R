#### Set WD, source utils, and load in libraries ----
ifelse(dir.exists("~/Box Sync/biodiversity_cmip6"),
       setwd("~/Box Sync/biodiversity_cmip6"),
       setwd("/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6"))

source("./analysis/util.R")

library(dplyr)
library(lfe)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(reshape2)
library(stringr)
library(broom)
library(gridExtra)
library(interactions)
library(quantreg)
library(jtools)
library(sf)
#library(gridExtra)
#library(report)

#### Set Global Vars ----
CLEAN = FALSE ##switch to false if you want to use original KBAs
td <- format(Sys.Date(), "%m_%d_%Y")
ERA <- "ssp370"
COUNTRY <- "BRA"
TYPE <- "Terrestrial"
models <- c("ACCESS-CM2", "CNRM-CM6-1", "CanESM5", "GFDL-ESM4") ## list of the models we will use in this case
sf::sf_use_s2(FALSE) ## to deal with some issues not fixable with st_make_valid

#### Get data ----
world <- st_read(dsn = "./raw_data/WB_countries_Admin0_10m/WB_countries_Admin0_10m.shp", stringsAsFactors = F, crs = 4326) 
kba_class <- read_csv("./raw_data/KBA2022/kba_class_2022_clean.csv")[,-1]
kba_protected_area <- read.csv("./processed_data/pa_kba_ovl/all_countries_2022.csv")
kba_geometry <- st_read(dsn = "./processed_data/kba/KBAsGlobal_2022_September_02_POL_noOverlaps.shp", stringsAsFactors = F, crs = 4326) 
#climate_zones  <- st_read(dsn = "./raw_data/climate_zones/other_climate_2007_koppen_geiger.shp")
pas <- st_read(dsn = "./processed_data/wdpa/clean_wdpa_terrestrial.shp")
ndvi <- read.csv("./processed_data/ndvi/ndvi_pa_ovl.csv")

## fix column names 
coln <- c("SitRecID", "Region", "Country", "ISO3", "NatName", "IntName", "FinCode", 
          "SitLat", "SitLong", "GISArea", "IbaStatus", "KBAStatus", 
          "AzeStatus", "AddedDate", "ChangeDate", "Source", "DelTxt",
          "DelGeom", "KBA_Quality", "Shape_Long", "Shape_Area", "LegacyKBA", "Criteria",
          "original_area", "kba_notes", "akba", "class", "geometry")
names(kba_geometry) <- coln

#### Clean and combine datasets ----

## subset all read in data based on global vars
kba_class <- kba_class %>% filter(ISO == COUNTRY & Type == TYPE)
kba_protected_area <- kba_protected_area %>% filter(SitRecID %in% unique(kba_class$SitRecID))
kba_geometry <- kba_geometry %>% filter(SitRecID %in% unique(kba_class$SitRecID))

kba_2022 <- cummulative_kba(kba_protected_area, years = c(2022), level ="kba")
kba_2022 <- kba_2022 %>% dplyr::select(SitRecID, kba, Country, cum_year = year,
                                       cum_overlap, cum_percPA) %>%
  mutate(protected = ifelse(cum_percPA < 2, "NP", ifelse(cum_percPA >= 98, "FP", "P")))
kba_country_2022 <- cummulative_kba(kba_protected_area, years = c(2022), level ="country")
kba_country_2022 <- kba_country_2022 %>% dplyr::select(kba, Country, cum_year,
                                                       cum_overlap, cum_percPA)  %>%
  mutate(protected = ifelse(cum_percPA < 2, "NP", ifelse(cum_percPA >= 98, "FP", "P")))

## get a list of WDPAIDs that match with KBAs
intersections <- st_intersects(kba_geometry, pas)
kba_wdpa <- c()
for(i in 1:nrow(intersections)) {
  indeces <- pas %>% slice(intersections[[i]])
  ifelse(nrow(indeces) == 0, WDPA <- NA, WDPA <- indeces %>% pull(WDPAID))
  kba_wdpa <- rbind(kba_wdpa,
                    cbind(SitRecID = kba_geometry[i,] %>% st_drop_geometry() %>% pull(SitRecID), 
                          WDPAID = WDPA))
}
kba_wdpa <- as.data.frame(kba_wdpa)

kba_wdpa <- left_join(kba_2022, kba_wdpa, by = "SitRecID")
kba_wdpa <- left_join(kba_wdpa, ndvi, by = c("WDPAID"))

#### Start Plots ----
pdf(paste0("./visuals/view_", td, ".pdf")) ## start pdf up here

#### Extremes