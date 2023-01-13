ifelse(dir.exists("~/Box Sync/biodiversity_cmip6"),
       setwd("~/Box Sync/biodiversity_cmip6"),
       setwd("/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6"))

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

#### Set Global Vars ----
CLEAN = FALSE ##switch to false if you want to use original KBAs
td <- format(Sys.Date(), "%m_%d_%Y")
COUNTRY <- "CHE"

pdf(paste0("./visuals/view_", td, ".pdf")) ## start pdf up here

#### Download global data ----
world <- st_read(dsn = "./raw_data/WB_countries_Admin0_10m/WB_countries_Admin0_10m.shp", stringsAsFactors = F, crs = 4326) 
kba_class <- read.csv("./raw_data/kba_class_2020.csv")

#### Just testing for now ----

## move to folder
loc <- ifelse(CLEAN, "/processed_data/cleankba_cmip_ovl/",
       "/processed_data/originalkba_cmip_ovl/")

all <- list.files(paste0(getwd(), loc), full.names = T, pattern = "ssp370")

kba_temps <- readRDS(all[1])
name <- str_split(all[1], "_")[[1]][8]

#### plot difference from 2015 to 2100 ----
a <- kba_temps %>% filter(year==2015 & ISO3 == COUNTRY) %>% dplyr::select(SitRecID, mean)
b <- kba_temps %>% filter(year==2100 & ISO3 == COUNTRY) %>% dplyr::select(SitRecID, mean) %>% st_drop_geometry()
kt <- left_join(a, b, by = "SitRecID")
country <- world %>% filter(ISO_A3 == COUNTRY)
kt <- left_join(kt, kba_class, by = "SitRecID")
kt <- kt %>% filter(terrestrial == 1 & marine == 0) %>% rename(cov = PA.Coverage) %>% mutate(dif = mean.y - mean.x)

pdf(paste0("./visuals/view_", COUNTRY, td, ".pdf")) ## start pdf up here

ggplot(data = kt) +
  ggtitle(paste(name, "Expected Change in Temperature 2015-2100")) +
  geom_sf(data = country, size = 0.002, fill = "#d9f0a3") +
  geom_sf(data = kt, size = 0.0002, aes(fill = dif)) +
  scale_fill_distiller(palette = "RdPu", direction = 1, na.value = "grey") +
  labs(colour="2100 - 2015 Annual Mean Temperature", fill = "Temperature \n Difference °C") +
  theme_bw()

### plot temperature trend
kt_nogeo <- kt %>% mutate(group = cut(dif, 5)) %>% st_drop_geometry 

ggplot(data = kt_nogeo, aes(group)) +
  geom_bar(aes(fill = factor(cov, levels = c("complete", "partial", "none")))) +
  labs(xlab = "°C Change from 2015 to 2100", ylab = "Number of KBAs", fill = "KBA Coverage") +
  theme_bw()


dev.off()
