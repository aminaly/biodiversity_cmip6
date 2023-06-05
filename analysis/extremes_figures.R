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
COUNTRY <- "ZAF"
TYPE <- "Terrestrial"
indexes <- c("*") ## extreme indexes for pattern matching (* if all)
sf::sf_use_s2(FALSE) ## to deal with some issues not fixable with st_make_valid
ids <- c(7076, 7086, 7090, 7092, 7095, 7097, 7100, 7101, 7102, 7103, 7107,
         7116, 7117, 7124, 7131, 7132, 7136, 7139, 7142, 7147, 7149, 7152, 
         7153, 7154, 7155, 7157, 7158, 7159, 7161, 7162, 7163, 7164, 7165,
         7166, 7168, 7170, 7171, 7172, 7174, 7175, 32050, 32058, 44661, 44671)

#### Get data ----
world <- st_read(dsn = "./raw_data/WB_countries_Admin0_10m/WB_countries_Admin0_10m.shp", stringsAsFactors = F, crs = 4326) 
world <- world %>% filter(ISO_A3 == COUNTRY)
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

## KBA dataset
kbas <- left_join(kba_class %>% filter(ISO == COUNTRY, Type == TYPE),
          kba_protected_area, by = c("Country", "SitRecID"))
kbas <- kbas %>% mutate(climate_threat = ifelse(SitRecID %in% ids, T, F))
kba_geometry <- kba_geometry %>% filter(ISO3 == COUNTRY)

plot_data <- c()
for(index in indexes) {
  if(index == "*") plot_data <- read.csv("./processed_data/extremes_kbas/all_plot_data.csv") break
  data_list <- list.files("./processed_data/extremes_kbas/", full.names = T, pattern = index)
  for(file in data_list) {
    data <- read_csv(file)[,-1]
    data <- data %>% filter(SitRecID %in% kbas$SitRecID) %>% 
      mutate(year = year(date), year_group = cut(year, 8, labels = F)) %>% select(-source)
    plot_data <- rbind(plot_data, data)
  }
}

#### Start Plots ----
pdf(paste0("./visuals/view_", td, ".pdf")) ## start pdf up here

#### Loop through and make the same figures for each measure
next_10 <- plot_data %>% group_by(SitRecID, year_group, measure) %>%
  summarize(mean_sd = sd(mean), mean_index = mean(mean)) %>% filter(year_group < 3) %>% 
  pivot_wider(id_cols = c(SitRecID, measure), names_from = year_group, values_from = c("mean_sd", "mean_index")) %>% 
  mutate(diff_mean = mean_index_2 - mean_index_1, 
         climate_threat = ifelse(SitRecID %in% ids, "Y", "N"))

for(index in unique(next_10$measure)) {
  
  ### Average of index for every 10 years (about)
  next_10_i <- next_10 %>% filter(measure == index)
  next_10_i <- left_join(next_10_i, kba_geometry) %>% st_set_geometry("geometry")
  
  paste(ggplot(data = next_10_i) +
    ggtitle("Change in Index (avg 2036-2026 minus 2025-2015)") +
    geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
    geom_sf(data = next_10_i, size = 0.0002, aes(fill = diff)) +
    coord_sf(ylim = c(-22, -35)) +
    scale_fill_continuous(na.value = "grey") +
    labs(fill = "Change in Index") +
    theme_bw())
  
  paste(ggplot(data = next_10_i) +
    ggtitle("SD 2025-2015") +
    geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
    geom_sf(data = next_10_i, size = 0.0002, aes(fill = mean_sd_1)) +
    coord_sf(ylim = c(-22, -35)) +
    scale_fill_continuous(na.value = "grey") +
    labs(fill = "Standard Deviation") +
    theme_bw())
  
  paste(ggplot(data = next_10_i) +
    ggtitle("SD 2026-2036") +
    geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
    geom_sf(data = next_10_i, size = 0.0002, aes(fill = mean_sd_2)) +
    coord_sf(ylim = c(-22, -35)) +
    scale_fill_continuous(na.value = "grey") +
    labs(fill = "Standard Deviation") +
    theme_bw())
  
  paste(ggplot(data = next_10_i, aes(x = mean_index_1, y = mean_index_2)) +
    geom_point(aes(color = climate_threat)) +
    facet_wrap(~ measure) +
    xlab("Average Index 2015-2025") + ylab("Average Index 2026-2036") +
    theme_bw())
  
}

dev.off()



