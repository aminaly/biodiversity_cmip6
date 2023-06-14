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
COUNTRY <- "ZAF"
TYPE <- "Terrestrial"
ERA <- "ssp370"
indexes <- c("r99p", "r95p", "tx90p", "tn90p", "wsdi")  ## extreme indexes for pattern matching (* if all)
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

#### Create plot datasets ----
## get KBAs in this country and current protections
kbas <- cummulative_kba(kba_protected_area, years = c(2022), level ="kba")
kbas <- kbas %>% dplyr::select(SitRecID, kba, Country, cum_year = year,
                                       cum_overlap, cum_percPA) %>%
  mutate(protected = ifelse(cum_percPA < 2, "NP", ifelse(cum_percPA >= 98, "FP", "P")))
kbas <- left_join(kba_class %>% filter(ISO == COUNTRY, Type == TYPE),
                  kbas, by = c("Country", "SitRecID"))
kbas <- kbas %>% mutate(climate_threat = ifelse(SitRecID %in% ids, T, F))
kba_geometry <- kba_geometry %>% filter(ISO3 == COUNTRY)

## create dataset of extreme indexes (by kba)
extreme_data <- c()
if(file.exists("./processed_data/extremes_kbas/all_plot_data.csv")) {
  extreme_data <- read.csv("./processed_data/extremes_kbas/all_plot_data.csv")  
} else {
  data_list <- list.files("./processed_data/extremes_kbas/", full.names = T, pattern = paste(indexes, collapse = "|"))
  for(file in data_list) {
    data <- read_csv(file)[,-1]
    data <- data %>% filter(SitRecID %in% kbas$SitRecID) %>% 
      mutate(year = year(date), #year_group = cut(year, 8, labels = F),
             scenario = str_split(source, "_")[[1]][4]) %>% select(-source)
    extreme_data <- rbind(extreme_data, data)
  }
}

extreme_data <- extreme_data %>% mutate(climate_threat = ifelse(SitRecID %in% ids, "Y", "N"))

hist <- extreme_data %>% 
  filter(scenario == "historical",
         year %in% c(2001:2011)) %>% 
  group_by(SitRecID, scenario, measure, climate_threat) %>%
  summarize(standard_dev = sd(mean), mean_index = mean(mean)) %>%
  mutate(year_group = "hist")

comp <- extreme_data %>% 
  filter(year %in% c(2015:2036),
         measure == paste0(index, "ETCCDI")) %>% 
  mutate(year_group = cut(year, 2, labels = c("first", "second"))) %>%
  group_by(SitRecID, scenario, year_group, measure) %>%
  summarize(standard_dev = sd(mean), mean_index = mean(mean)) %>%
  pivot_wider(id_cols = c(SitRecID, scenario), 
              names_from = year_group, values_from = c("standard_dev", "mean_index")) %>%
  filter(if(ERA != "") scenario == ERA)

extreme_comp_data <- left_join(hist %>% ungroup %>% select(-scenario), comp, by = "SitRecID") %>%
  mutate(diff_first = mean_index_first - mean_index,
         diff_second = mean_index_second - mean_index)

if(ERA != "") extreme_comp_data <- extreme_comp_data %>% filter(scenario == ERA)

## create dataset of goverance types
intersections <- st_intersects(kba_geometry, pas)
governance <- c()
for(i in 1:nrow(intersections)) {
  indeces <- pas %>% slice(intersections[[i]])
  ifelse(nrow(indeces) == 0, WDPA <- NA, WDPA <- indeces %>% pull(WDPAID))
  governance <- rbind(governance,
                    cbind(SitRecID = kba_geometry[i,] %>% st_drop_geometry() %>% pull(SitRecID), 
                          WDPAID = WDPA))
}
governance <- as.data.frame(governance)

governance <- left_join(kbas, governance, by = "SitRecID")
governance <- left_join(governance, ndvi %>%
                          filter(year == 2022, ISO3 == COUNTRY) %>% 
                          select(WDPAID, GOV_TYPE, DESIG_TYPE), 
                        by = c("WDPAID"))

#### Start Plots ----
#### Plot extremes ----
pdf(paste0("./visuals/extremes_", td, ".pdf"), width=7,height=7,pointsize=12) ## start pdf up here
for(index in indexes) {
  
  m <- paste0(index, "ETCCDI")
  data <- left_join(extreme_comp_data %>% filter(measure == m), 
                    kbas %>% select(SitRecID, percent_protected = cum_percPA, protected)) %>% 
    filter(measure == measure)
  data_i <- left_join(data, 
                      kba_geometry %>% 
                        filter(SitRecID %in% unique(data$SitRecID)) %>% 
                        select(SitRecID, geometry),
                      by = "SitRecID") %>% st_set_geometry("geometry") 
  
  ## historic extremes
  lims <- range(data_i$mean_index, data_i$mean_index_first, data_i$mean_index_second)
  a <- (ggplot(data = data_i) +
          ggtitle(paste("Historic ('01-11)", index)) +
          geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index)) +
          coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", 
                              na.value = "grey", limits = lims) +
          theme_bw())
  
  b <- (ggplot(data = data_i) +
          ggtitle(paste("'15-25", index)) +
          geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index_first)) +
          coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", 
                              na.value = "grey", limits = lims) +
          theme_bw())
  
  c <- (ggplot(data = data_i) +
          ggtitle(paste("26-36", index)) +
          geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index_second)) +
          coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", 
                              na.value = "grey", limits = lims) +
          theme_bw())
  
  print(grid.arrange(a, b, c))
  
  ## historic vs next 10 
  lims <- range(data_i$mean_index_first, data_i$mean_index_second)
  first <- ggplot(data = data_i, aes(x = mean_index, y = mean_index_first)) +
    ggtitle(paste("Expected '15-25 + climate threat ", m)) +
    geom_point(aes(color = climate_threat)) + ylim(lims) +
    xlab("Average Index `01-11") + ylab("Average Index 2015-2025") +
    theme_bw()
  second <- ggplot(data = data_i, aes(x = mean_index, y = mean_index_second)) +
    ggtitle(paste("Expected '26-36 + climate threat ", m)) +
    geom_point(aes(color = climate_threat)) + ylim(lims) +
    xlab("Average Index `01-11") + ylab("Average Index 2015-2025") +
    theme_bw()
  
  print(grid.arrange(first, second))
  
  #historic vs next 10, but protected
  lims <- range(data_i$mean_index_first, data_i$mean_index_second)
  first <- ggplot(data = data_i, aes(x = mean_index, y = mean_index_first)) +
    ggtitle(paste("Expected '15-25 + protection ", m)) +
    geom_point(aes(color = percent_protected)) + ylim(lims) +
    xlab("Average Index `01-11") + ylab("Average Index 2015-2025") +
    theme_bw()
  second <- ggplot(data = data_i, aes(x = mean_index, y = mean_index_second)) +
    ggtitle(paste("Expected '26-36 + protection ", m)) +
    geom_point(aes(color = percent_protected)) + ylim(lims) +
    xlab("Average Index `01-11") + ylab("Average Index 2015-2025") +
    theme_bw()
  
  print(grid.arrange(first, second))
  
  ## difference in each decade to historical
  lims <- range(data_i$diff_first, data_i$diff_second)
  first <- ggplot(data = data, aes(x = mean_index, y = diff_first)) +
    ggtitle(paste("Change historical to '15-25 + climate threat ", m)) +
    geom_point(aes(color = climate_threat)) + ylim(lims) +
    xlab("Average Index `01-11") + ylab("Change in Average Index (10yr)") +
    theme_bw()
  second <- ggplot(data = data, aes(x = mean_index, y = diff_second)) +
    ggtitle(paste("Change historical to '26-36 + climate threat ", m)) +
    geom_point(aes(color = climate_threat)) + ylim(lims) +
    xlab("Average Index `01-11") + ylab("Change in Average Index (20 yr)") +
    theme_bw()
  
  print(grid.arrange(first, second))
  
  
}

dev.off()

#### Plot protection  ----
pdf(paste0("./visuals/protection_", td, ".pdf")) ## start pdf up here

data_i <- left_join(kbas, 
                    kba_geometry %>% 
                      filter(SitRecID %in% unique(kbas$SitRecID)) %>% 
                      select(SitRecID, geometry),
                    by = "SitRecID") %>% st_set_geometry("geometry") %>% 
  filter(!is.na(cum_percPA)) %>% rename(percent_protected = cum_percPA) 

## current protection status
print(ggplot(data = data_i) +
        ggtitle(paste("Current Protection Status")) +
        geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
        geom_sf(data = data_i, size = 0.0002, aes(fill = percent_protected)) +
        coord_sf(ylim = c(-22, -35)) +
        labs(fill = "Percent of KBA covered by Protected Areas") +
        scale_fill_gradient(low = "#edf8fb", high = "#006d2c", na.value = "grey") +
        theme_bw())

data_i <- data_i %>% mutate(protected = fct_relevel(protected, c("FP", "P", "NP")))
print(ggplot(data = data_i) +
        ggtitle(paste("Current Coverage Group", index)) +
        geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
        geom_sf(data = data_i, size = 0.0002, aes(fill = protected)) +
        coord_sf(ylim = c(-22, -35)) +
        labs(fill = "Change in Index") +
        scale_fill_manual(values = c("#4d9221", "#66c2a4", "#c51b7d")) +
        theme_bw())

ggplot(data = data_i, aes(protected)) +
  geom_bar(aes(fill = factor(protected, levels = c("FP", "P", "NP")))) +
  labs(xlab = "°C Change from 2015 to 2100", ylab = "Number of KBAs", fill = "KBA Coverage") +
  scale_fill_manual(values = c("#4d9221", "#66c2a4", "#c51b7d")) +
  theme_bw()

## and who is in charge (governance type)
pas_i <- left_join(governance, pas %>% select(WDPAID, geometry), by = "WDPAID") %>% 
  st_set_geometry("geometry") %>%
  filter(!is.na(GOV_TYPE))
ggplot(data = pas_i) +
  ggtitle(paste("Governance Type (with KBA boundaries")) +
  geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
  geom_sf(data = pas_i, size = 0.0002, aes(fill = GOV_TYPE)) +
  geom_sf(data = data_i, size = 0.002, fill = "transparent") +
  coord_sf(ylim = c(-22, -35)) +
  scale_fill_discrete(na.value = "grey") +
  labs(fill = "Gov Type") +
  theme_bw()

ggplot(data = pas_i) +
  ggtitle(paste("Governance Type (w/o KBA boundaries)")) +
  geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
  geom_sf(data = pas_i, size = 0.0002, aes(fill = GOV_TYPE)) +
  #geom_sf(data = data_i, size = 0.002, fill = "transparent") +
  coord_sf(ylim = c(-22, -35)) +
  scale_fill_discrete(na.value = "grey") +
  labs(fill = "Gov Type") +
  theme_bw()

dev.off()
#### Loop through and plot historic + 2015-2025 for each SSP for each index ----
for(index in indexes) {
  
  # select correct measure
  data <- plot_data %>% filter(measure == paste0(index, "ETCCDI"))
  data_i <- left_join(data, 
                      kba_geometry %>% filter(SitRecID %in% unique(data$SitRecID)),
                      by = "SitRecID") %>% st_set_geometry("geometry")
  
  
  print(ggplot(data = data_i) +
          ggtitle(paste("Average 2015-2025", index)) +
          geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index_first)) +
          coord_sf(ylim = c(-22, -35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", na.value = "grey") +
          facet_wrap(~ scenario.y, nrow = 3) +
          theme_bw())
  
  print(ggplot(data = data_i) +
          ggtitle(paste("Average 2026-2036", index)) +
          geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index_second)) +
          coord_sf(ylim = c(-22, -35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", na.value = "grey") +
          facet_wrap(~ scenario.y, nrow = 3) +
          theme_bw())
  
  print(ggplot(data = data_i) +
          ggtitle(paste("Diff Historical + 2015-2025", index)) +
          geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = diff_first)) +
          coord_sf(ylim = c(-22, -35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", na.value = "grey") +
          facet_wrap(~ scenario.y, nrow = 3) +
          theme_bw())
  
  print(ggplot(data = data_i) +
          ggtitle(paste("Diff Historical + 2026-2036", index)) +
          geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = diff_second)) +
          coord_sf(ylim = c(-22, -35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", na.value = "grey") +
          facet_wrap(~ scenario.y, nrow = 3) +
          theme_bw())
  
  reg <- lm(mean_index_first ~ mean_index, data = data)
  print(ggplot(data = data, aes(x = mean_index, y = mean_index_first)) +
          ggtitle(paste("Change historical to '15-25 + climate threat ", index)) +
          geom_point(aes(color = climate_threat)) +
          geom_abline(slope=reg$coefficients[2], intercept = reg$coefficients[1]) +
          facet_wrap(~ scenario.y, nrow = 3) +
          xlab("Average Index Hist") + ylab("Average Index 2015-2025") +
          theme_bw())
  
  
  reg <- lm(mean_index_second ~ mean_index, data = data)
  print(ggplot(data = data, aes(x = mean_index, y = mean_index_second)) +
          ggtitle(paste("Change historical to '26-36 + climate threat ", index)) +
          geom_point(aes(color = climate_threat)) +
          geom_abline(slope=reg$coefficients[2], intercept = reg$coefficients[1]) +
          facet_wrap(~ scenario.y, nrow = 3) +
          xlab("Average Index Hist") + ylab("Average Index 2026-2036") +
          theme_bw())
  
}

dev.off()

#### Do some ANOVAs that will plot out differences in models ----

for(s in unique(plot_data$scenario)) {
  if(scenario == "historic") next
  for(m in unique(plot_data$measure)){
    d <- plot_data %>% filter(scenario == s & measure == m, 
                              year %in% c(2015:2025))
    aov(mean ~ year * gcm, data = d)
    
  }
}


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
  
  print(ggplot(data = next_10_i) +
    ggtitle(paste("Change in Index (avg 2036-2026 minus 2025-2015) ", index)) +
    geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
    geom_sf(data = next_10_i, size = 0.0002, aes(fill = diff_mean)) +
    coord_sf(ylim = c(-22, -35)) +
    scale_fill_continuous(na.value = "grey") +
    labs(fill = "Change in Index") +
    theme_bw())
  
  print(ggplot(data = next_10_i) +
    ggtitle(paste("SD 2025-2015 ", index)) +
    geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
    geom_sf(data = next_10_i, size = 0.0002, aes(fill = mean_sd_1)) +
    coord_sf(ylim = c(-22, -35)) +
    scale_fill_continuous(na.value = "grey") +
    labs(fill = "Standard Deviation") +
    theme_bw())
  
  print(ggplot(data = next_10_i) +
    ggtitle(paste("SD 2026-2036 ", index)) +
    geom_sf(data = world, size = 0.002, fill = "#d9f0a3") +
    geom_sf(data = next_10_i, size = 0.0002, aes(fill = mean_sd_2)) +
    coord_sf(ylim = c(-22, -35)) +
    scale_fill_continuous(na.value = "grey") +
    labs(fill = "Standard Deviation") +
    theme_bw())
  
}

dev.off()



