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
library(gridExtra)
library(report)

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

#### Plots dealing with CMIP6 models ----
#### loop through selected models, make plots for individuals and aggregate data ----
## run through all CMIP models overlapped with KBA and aggregate
all_cmips <- list.files(paste0(getwd(),"/processed_data/cleankba_cmip_ovl"), full.names = T, pattern = ERA)
kba_2022_cmip <- c()
kba_country_2022_cmip <- c()

for(file in all_cmips) {
  name <- str_split(file[1], "_")[[1]][8]
  if(sum(grepl(name, models)) == 0) next()
  kba_temps <- read_csv(file)
  kba_temps <- kba_temps %>% filter(SitRecID %in% unique(kba_class$SitRecID))
  
  kba_temps <- kba_temps %>% 
    mutate(yr_grp = cut(year,c(2015, 2030, 2045, 2060, 2075, 2090, 2100), 
                        labels=c("2015-2030", "2030-2045", 
                                 "2045-2060", "2060-2075", "2075-2090", 
                                 "2090-2100"), include.lowest = T)) %>%
    group_by(SitRecID, yr_grp) %>% mutate(model = name) %>% 
    filter(SitRecID %in% kba_class$SitRecID) %>%
    dplyr::select(SitRecID, year, mean, yr_grp, model)
  
  kc <- left_join(kba_2022, kba_temps, by = c("SitRecID"))
  kba_2022_cmip <- rbind(kba_2022_cmip, kc)
  
  ### plot this model difference 2015 to 2100 ----
  a <- kc %>% filter(year==2015) %>% ungroup() %>% dplyr::select(SitRecID, mean, protected) 
  b <- kc %>% filter(year==2100) %>% ungroup() %>% dplyr::select(SitRecID, mean, protected)
  kt <- left_join(a, b, by = c("SitRecID", "protected"))
  kt <- kt %>% mutate(dif = mean.y - mean.x) %>% mutate(group = cut(dif, 4))
  kt_g <- left_join(kt, kba_geometry, by = "SitRecID")
  kt_g <- st_set_geometry(kt_g, "geometry")
  
  country <- world %>% filter(ISO_A3 == COUNTRY)
  
  #### expected change in temp from 2015 to 2100 
  paste(ggplot(data = kt_g) +
    ggtitle(paste(name, "Expected Change in Temperature 2015-2100")) +
    geom_sf(data = country, size = 0.002, fill = "#d9f0a3") +
    geom_sf(data = kt_g, size = 0.0002, aes(fill = dif)) +
    scale_fill_distiller(palette = "RdPu", direction = 1, na.value = "grey") +
    labs(colour="2100 - 2015 Annual Mean Temperature", fill = "Temperature \n Difference °C") +
    theme_bw())
  
  #### protection status by expected temperature change
  kt <- kt %>% mutate(protected = fct_relevel(protected, c("FP", "P", "NP")))
  paste(ggplot(data = kt, aes(group, group = protected)) +
    geom_bar(aes(fill = factor(protected))) +
    labs(title = paste(name, "KBA count by expected temperature change"), fill = "KBA Coverage") +
    xlab("°C Change from 2015 to 2100") +
    ylab("Number of KBAs") +
    scale_fill_manual(values = c("#4d9221", "#e6f5d0", "#c51b7d")) +
    theme_bw())
  
  #### avg temp for each 15 years 
  kt <- kc %>% group_by(yr_grp, SitRecID, protected) %>% summarize(mean = mean(mean, na.rm = T))
  kt_g <- left_join(kt, kba_geometry, by = "SitRecID")
  kt_g <- st_set_geometry(kt_g, "geometry")
  
  paste(ggplot(data = kt_g) +
    ggtitle(paste(name, "Avg Temperature every 15 years")) +
    geom_sf(data = country, size = 0.002, fill = "#d9f0a3") +
    geom_sf(data = kt_g, size = 0.0002, aes(fill = mean)) +
    scale_fill_distiller(palette = "RdPu", direction = 1, na.value = "grey") +
    labs(colour="2100 - 2015 Annual Mean Temperature", fill = "Average \n Temperature °C") +
    facet_wrap(~ yr_grp, nrow = 2) + 
    theme_bw())
  
  #### protection status as of 2022
  paste(ggplot(data = kt_g %>% filter(yr_grp == "2015-2030")) +
    ggtitle(paste(name, "Avg Temperature every 15 years")) +
    geom_sf(data = country, size = 0.002, fill = "#d9f0a3") +
    geom_sf(data = kt_g, size = 0.0002, aes(fill = protected)) +
    #scale_fill_discrete(palette = "RdPu", direction = 1, na.value = "grey") +
    labs(colour="2100 - 2015 Annual Mean Temperature", fill = "Temperature \n Difference °C") +
    theme_bw())
  
}

#### 4 models box and whisker ----
ggplot(kba_2022_cmip, aes(model, mean, group = model)) +
  geom_boxplot(aes(fill = model), outlier.alpha = .1) +
  #scale_fill_brewer(palette = "Oranges") +
  facet_wrap(~ yr_grp, nrow = 3)
  labs(title = "Annual Mean Temperatures, every 15 years")

a <- kba_2022_cmip %>% filter(year==2015) %>% ungroup() %>% dplyr::select(SitRecID, mean, model) 
b <- kba_2022_cmip %>% filter(year==2100) %>% ungroup() %>% dplyr::select(SitRecID, mean, model)
kt <- left_join(a, b, by = c("SitRecID", "model"))
kt <- kt %>% mutate(dif = mean.y - mean.x)
kt <- left_join(kt, kba_2022 %>% select(SitRecID, protected), by = "SitRecID")

ggplot(kt, aes(model, dif)) +
  geom_boxplot(outlier.alpha = .1, aes(fill = model)) +
  #scale_fill_brewer(palette = "Oranges") +
  labs(title = "Expected Temperature Difference 2015-2100")


#### Plots having to do with protections & governance ----
#### protection proportions ----
kt <- kt %>% mutate(protected = fct_relevel(protected, c("FP", "P", "NP")))
ggplot(data = kt, aes(protected)) +
  geom_bar(aes(fill = factor(protected, levels = c("FP", "P", "NP")))) +
  labs(xlab = "°C Change from 2015 to 2100", ylab = "Number of KBAs", fill = "KBA Coverage") +
  scale_fill_manual(values = c("#4d9221", "#e6f5d0", "#c51b7d")) +
  theme_bw()

ggplot(data = kt %>% mutate(temp_grp = cut(dif, 4)), aes(temp_grp, group = protected)) +
  geom_bar(aes(fill = factor(protected, levels = c("FP", "P", "NP")))) +
  labs(xlab = "°C Change from 2015 to 2100", ylab = "Number of KBAs", fill = "KBA Coverage") +
  scale_fill_manual(values = c("#4d9221", "#e6f5d0", "#c51b7d")) +
  facet_wrap(~model)
  theme_bw()

#### governance ----
kba_wdpa_g <- left_join(kba_wdpa, pas %>% select(WDPAID, geometry))
kba_wdpa_g <- kba_wdpa_g %>% st_set_geometry("geometry") %>% filter(year = 2022)

ggplot(data = kba_wdpa_g) +
  ggtitle(paste(name, "Governance Type")) +
  geom_sf(data = country, size = 0.002, fill = "#d9f0a3") +
  geom_sf(data = kba_wdpa_g, size = 0.0002, aes(fill = GOV_TYPE)) +
  geom_sf(data = kt_g, size = 0.0002, fill = "transparent") +
  scale_fill_discrete(na.value = "grey") +
  labs(colour="2100 - 2015 Annual Mean Temperature", fill = "Temperature \n Difference °C") +
  theme_bw()

#### plots having to do with NDVI ----
ndvi_vars <- kba_wdpa %>% select(GOV_TYPE, max_ndvi, mean_ndvi)

ggplot(ndvi_vars, aes(x = GOV_TYPE, y = max_ndvi)) + 
  geom_boxplot(aes(color = GOV_TYPE)) + 
  labs(title = "Max NDVI") +
  xlab("Governance Type") +
  ylab("Maximum NDVI") +
  theme_bw()

ggplot(ndvi_vars, aes(x = GOV_TYPE, y = mean_ndvi)) + 
  geom_boxplot(aes(color = GOV_TYPE)) + 
  labs(title = "Mean NDVI") +
  xlab("Governance Type") +
  ylab("Mean NDVI") +
  theme_bw()

aov_ndvi <- aov(mean_ndvi ~ GOV_TYPE, data = ndvi_vars)
report(aov_ndvi)
plot(TukeyHSD(aov_ndvi), las=1,cex.axis=0.4, sub = "mean_ndvi")

aov_ndvi <- aov(max_ndvi ~ GOV_TYPE, data = ndvi_vars)
report(aov_ndvi)
plot(TukeyHSD(aov_ndvi), las=1,cex.axis=0.4, sub = "max ndvi")

dev.off()

#### general plots ----
kba_wdpa_g <- left_join(kba_wdpa, pas %>% select(WDPAID, geometry))
kba_wdpa_g <- kba_wdpa_g %>% group_by(SitRecID) %>% summarize(mn = mean(mean_ndvi, na.rm = T))

ggplot(data = kba_wdpa_g) +
  ggtitle(paste(name, "Governance Type")) +
  geom_sf(data = country, size = 0.002, fill = "#d9f0a3") +
  geom_sf(data = kba_wdpa_g, size = 0.0002, aes(fill = mn)) +
  scale_fill_continuous(na.value = "grey", color = "greens") +
  labs(colour="Average NDVI Values (protected areas only)", 
       fill = "Average NDVI") +
  theme_bw()

dev.off()
