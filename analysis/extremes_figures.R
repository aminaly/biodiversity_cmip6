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
library(forcats)
library(nngeo)
#library(ggpattern)
#library(ggridges)
#library(gridExtra)
#library(report)

#### Set Global Vars ----

## pick up args from command line/sbatch
args <- commandArgs(trailingOnly = TRUE)
measurerep <- as.numeric(args[1])

CLEAN = FALSE ##switch to false if you want to use original KBAs
td <- format(Sys.Date(), "%m_%d_%Y")
COUNTRY <- "ZAF"
TYPE <- "Terrestrial"
ERA <- "ssp370"
indexes <- c("cddETCCDI", "r95pETCCDI", "txxETCCDI", "wsdiETCCDI")  ## extreme indexes for pattern matching (* if all)
sf::sf_use_s2(FALSE) ## to deal with some issues not fixable with st_make_valid
ids <- c(7076, 7086, 7090, 7092, 7095, 7097, 7100, 7101, 7102, 7103, 7107,
         7116, 7117, 7124, 7131, 7132, 7136, 7139, 7142, 7147, 7149, 7152, 
         7153, 7154, 7155, 7157, 7158, 7159, 7161, 7162, 7163, 7164, 7165,
         7166, 7168, 7170, 7171, 7172, 7174, 7175, 32050, 32058, 44661, 44671)
num_k <- 5 # set number of locations to consider in vulnerability

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
  #### get KBAs in this country and current protections ----
  kbas <- cummulative_kba(kba_protected_area, years = c(2022), level ="kba")
  kbas <- kbas %>% dplyr::select(SitRecID, kba, Country, cum_year = year,
                                         cum_overlap, cum_percPA) %>%
    mutate(protected = ifelse(cum_percPA < 2, "NP", ifelse(cum_percPA >= 98, "FP", "P")))
  kbas <- left_join(kba_class %>% filter(ISO == COUNTRY, Type == TYPE),
                    kbas, by = c("Country", "SitRecID"))
  kbas <- kbas %>% mutate(climate_threat = ifelse(SitRecID %in% ids, T, F))
  kba_geometry <- kba_geometry %>% filter(ISO3 == COUNTRY)
  
  #### create dataset of extreme indexes (by kba) ----
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
    write.csv(extreme_data, "./processed_data/extremes_kbas/all_plot_data.csv")
  }
  
  extreme_data <- extreme_data %>% mutate(climate_threat = ifelse(SitRecID %in% ids, "Y", "N"))
  
  hist <- extreme_data %>% 
    filter(scenario == "historical",
           year %in% c(1995:2014)) %>% 
    group_by(SitRecID, scenario, measure, climate_threat) %>%
    summarize(standard_dev = sd(mean, na.rm = T), mean_index = median(mean, na.rm = T)) %>%
    mutate(year_group = "hist")
  
  comp <- extreme_data %>% 
    filter(year %in% c(2015:2036)) %>% 
    mutate(year_group = cut(year, 2, labels = c("first", "second"))) %>%
    group_by(SitRecID, scenario, year_group, measure) %>%
    summarize(standard_dev = sd(mean, na.rm = T), mean_index = median(mean, na.rm = T)) %>%
    pivot_wider(id_cols = c(SitRecID, scenario, measure), 
                names_from = year_group, values_from = c("standard_dev", "mean_index")) %>%
    filter(if(ERA != "") scenario == ERA) %>% ungroup %>% select(-scenario)
  
  extreme_comp_data <- left_join(hist %>% ungroup %>% dplyr::select(-scenario), 
                                 comp, by = c("SitRecID", "measure")) %>%
    mutate(diff_first = ifelse(measure %in% c("tn90pETCCDI", "tx90pETCCDI", "txxETCCDI"), 
                               mean_index_first - mean_index,
                               ((mean_index_first - mean_index)/mean_index)*100),
           diff_abs_first = mean_index_first - mean_index,
           diff_abs_second = mean_index_second - mean_index,
           diff_second = ifelse(measure %in% c("tn90pETCCDI", "tx90pETCCDI", "txxETCCDI"), 
                                mean_index_second - mean_index,
                                ((mean_index_second - mean_index)/mean_index)*100))
  
  categories <- as.data.frame(cbind(measure = unique(extreme_comp_data$measure),
                                    category = c("precip", "precip", "temp", "temp")))
  
  extreme_comp_data <- left_join(extreme_comp_data, categories, by = "measure")
  
  
  #### create dataset of goverance types ----
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
  
  governance <- left_join(kbas %>% filter(!is.na(kba)), governance, by = "SitRecID")
  governance <- left_join(governance, ndvi %>%
                            filter(year == 2022, ISO3 == COUNTRY) %>%
                            select(WDPAID, GOV_TYPE, DESIG_TYPE),
                          by = c("WDPAID"))
  
  
  #### Model Agreement ----
  if(file.exists("./processed_data/model_agreement/model_agreement.csv")) {
    model_agreement <- read.csv("./processed_data/model_agreement/model_agreement.csv")
  } else {
    ## as in Singh et al 2014 and Horton et al. 2014 which use IPCC thresholds of 66% agreement 
    extreme_data <- extreme_data %>% mutate(climate_threat = ifelse(SitRecID %in% ids, "Y", "N"))
    measures <- unique(extreme_data$measure)[measurerep]
    sites <- unique(extreme_data$SitRecID)
    reps <- 1000
    model_agreement <- mod_agreement(extreme_data, measures, sites, reps)
    
    
    model_agreement <- model_agreement %>%
      mutate(confidence = ifelse(over0 >= .99 | under0 >= .99, "virtually certain",
                                 ifelse(over0 >= .9 | under0 >= .9, "very likely",
                                        ifelse(over0 >= .66 | under0 >= .66, "likely",
                                               ifelse(over0 >= .33 | under0 >= .33, "neither",
                                                      "unlikely"))))) %>%
      mutate(confidence = fct_relevel(confidence, c("unlikely", "neither", "likely", "very likely", "virtually certain")))
  }
  
  #### measures of change in climate hazard ----
  
  #get zscore and normalized value (divided by largest in group)
  climate_hazard <- extreme_comp_data %>% group_by(measure) %>%
    mutate(z_score_hist = (mean_index - mean(mean_index)) / sd(mean_index), 
           z_score_first = (mean_index_first - mean(mean_index_first)) / sd(mean_index_first),
           z_score_second = (mean_index_second - mean(mean_index_second)) / sd(mean_index_second)) %>% ##zscore of this value within the distribution of the measure in the respective decade
    mutate(normalized_hist = mean_index/max(mean_index),
           normalized_first = mean_index_first/max(mean_index_first),
           normalized_second = mean_index_second/max(mean_index_second)) %>% ## normalized data by dividing by the max value
    ungroup %>% select(SitRecID, measure, z_score_hist, z_score_first, z_score_second,
                       normalized_hist, normalized_first, normalized_second)
  
  ## for each site, get the euclidian distance between the zscore and normalized values from the first decade to the second decade
  climate_hazard_zscore <-  climate_hazard %>% 
    group_by(SitRecID) %>% summarize(zscore_euclid_climate_first = dist(rbind(z_score_hist,z_score_first))[1],
                                     zscore_euclid_climate_second = dist(rbind(z_score_hist,z_score_second))[1])
  
  climate_hazard_normalized <-  climate_hazard %>% 
    group_by(SitRecID) %>% summarize(normalized_euclid_climate_first = dist(rbind(normalized_hist,normalized_first))[1],
                                     normalized_euclid_climat_second = dist(rbind(normalized_hist,normalized_second))[1])

#### measures of vulnerability (distance to other KBAs and protection) ----

climate_vulnerability <- st_nn(kba_geometry, pas, k = num_k, parallel = 5, returnDist = T)$dist
cv <- unlist(lapply(climate_vulnerability, FUN = mean))
cv <- cbind(kba_geometry$SitRecID %>% st_drop_geometry(), cv)
saveRDS(cv, "./processed_data/kba_dist_pas.rds")

climate_vulnerability <- st_nn(kba_geometry, kba_geometry, k = num_k + 1 , parallel = 5, returnDist = T)$dist
cv <- unlist(lapply(climate_vulnerability, FUN = mean))
cv <- cbind(kba_geometry$SitRecID %>% st_drop_geometry(), cv)
saveRDS(cv, "./processed_data/kba_dist_kbas.rds")

#### rank climate hazard based on percentile changes ----
model_agreement_rank_first <- left_join(model_agreement %>% filter(X == "firstdecade"), 
                                  kbas %>% select(SitRecID, 
                                                  percent_protected = cum_percPA),
                                  by = "SitRecID") %>%
  pivot_wider(id_cols = c(SitRecID, percent_protected),
              names_from = measure,
              values_from = c("mid")) %>%
  mutate(cdd_rank = percent_rank(abs(cddETCCDI)), txx_rank = percent_rank(abs(txxETCCDI)),
         r95p_rank = percent_rank(abs(r95pETCCDI)), wsdi_rank = percent_rank(abs(wsdiETCCDI)),
         protected_rank = percent_rank(desc(percent_protected)),
         total_rank = cdd_rank + txx_rank + r95p_rank + wsdi_rank) %>% 
  arrange(desc(total_rank)) %>% mutate(risk_group = cut(total_rank, 6, labels = F))
  
model_agreement_rank_second <- left_join(model_agreement %>% filter(X == "seconddecade"), 
                                         kbas %>% select(SitRecID, 
                                                         percent_protected = cum_percPA),
                                         by = "SitRecID") %>%
  pivot_wider(id_cols = c(SitRecID, percent_protected),
              names_from = measure,
              values_from = c("mid")) %>%
  mutate(cdd_rank = percent_rank(abs(cddETCCDI)), txx_rank = percent_rank(abs(txxETCCDI)),
         r95p_rank = percent_rank(abs(r95pETCCDI)), wsdi_rank = percent_rank(abs(wsdiETCCDI)),
         protected_rank = percent_rank(desc(percent_protected)),
         total_rank = cdd_rank + txx_rank + r95p_rank + wsdi_rank) %>% 
  arrange(desc(total_rank)) %>% mutate(risk_group = cut(total_rank, 6, labels = F))
  
  
#### Start Plots ----
pdf(paste0("./visuals/protection_", td, ".pdf")) ## start pdf up here

#### Figure plot protection  ----

data_i <- left_join(kbas,
                    kba_geometry %>%
                      filter(SitRecID %in% unique(kbas$SitRecID)) %>%
                      select(SitRecID, geometry),
                    by = "SitRecID") %>% st_set_geometry("geometry") %>%
  filter(!is.na(cum_percPA)) %>% rename(percent_protected = cum_percPA)

## current protection status by percent
print(ggplot(data = data_i) +
        ggtitle(paste("Current Protection Status")) +
        geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
        geom_sf(data = data_i, size = 0.0002, aes(fill = percent_protected)) +
        coord_sf(ylim = c(-22, -35)) +
        labs(fill = "Percent of KBA covered by Protected Areas") +
        scale_fill_gradient(low = "#edf8fb", high = "#006d2c", na.value = "grey") +
        theme_bw())

## curernt protection status by group
data_i <- data_i %>% mutate(protected = fct_relevel(protected, c("FP", "P", "NP")))
print(ggplot(data = data_i) +
        ggtitle(paste("Current Coverage Group")) +
        geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
        geom_sf(data = data_i, size = 0.0002, aes(fill = protected)) +
        coord_sf(ylim = c(-22, -35)) +
        labs(fill = "Protection Progress") +
        scale_fill_manual(values = c("#4d9221", "#66c2a4", "#c51b7d")) +
        theme_bw())

## map of climate threats
print(ggplot(data = data_i) +
        ggtitle(paste("Climate Threat")) +
        geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
        geom_sf(data = data_i, size = 0.0002, aes(fill = climate_threat)) +
        coord_sf(ylim = c(-22, -35)) +
        labs(fill = "IUCN Climate Threatened") +
        scale_fill_manual(values = c("#01665e", "#8c510a")) +
        theme_bw())

## bar of protections
ggplot(data = data, aes(protected)) +
  geom_bar(aes(fill = factor(protected, levels = c("FP", "P", "NP")))) +
  labs(xlab = "protected", ylab = "Count of KBAs", fill = "KBA Coverage") +
  scale_fill_manual(values = c("#4d9221", "#66c2a4", "#c51b7d")) +
  theme_bw()

## bar of protections and climate threats
ggplot(data = data, aes(x = protected, alpha = climate_threat)) +
  geom_bar(position = "dodge", aes(fill = factor(protected, levels = c("FP", "P", "NP")))) +
  labs(xlab = "Protection Status", ylab = "Count of KBAs", fill = "KBA Coverage") +
  scale_fill_manual(values = c("#4d9221", "#66c2a4", "#c51b7d")) +
  theme_bw()


## and who is in charge (governance type)
pas_i <- left_join(governance, pas %>% select(WDPAID, geometry), by = "WDPAID") %>%
  st_set_geometry("geometry") %>%
  filter(!is.na(GOV_TYPE))
ggplot(data = pas_i) +
  ggtitle(paste("Governance Type (with KBA boundaries)")) +
  geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
  geom_sf(data = pas_i, size = 0.0002, aes(fill = GOV_TYPE)) +
  geom_sf(data = data_i, size = 0.002, fill = "transparent") +
  coord_sf(ylim = c(-22, -35)) +
  scale_fill_discrete(na.value = "grey") +
  labs(fill = "Gov Type") +
  theme_bw()

ggplot(data = pas_i) +
  ggtitle(paste("Governance Type (w/o KBA boundaries)")) +
  geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
  geom_sf(data = pas_i, size = 0.0002, aes(fill = GOV_TYPE)) +
  #geom_sf(data = data_i, size = 0.002, fill = "transparent") +
  coord_sf(ylim = c(-22, -35)) +
  scale_fill_discrete(na.value = "grey") +
  labs(fill = "Gov Type") +
  theme_bw()

dev.off()
#### Figure plot 4 climate vars and model agreement ---- 
pdf(paste0("./visuals/model_agreement_", td, ".pdf")) ## start pdf up here

## plot map of all KBAs
data_i <- left_join(model_agreement,
                    kba_geometry %>%
                      filter(SitRecID %in% unique(data$SitRecID)) %>%
                      select(SitRecID, geometry),
                    by = "SitRecID") %>% st_set_geometry("geometry")

world_crop <- st_crop(world, xmin = 25, xmax = 33, ymin = -34, ymax = -25)

for(dec in c("firstdecade", "seconddecade")) {
  for(ind in 1:length(indexes)) {
    i <- indexes[ind]
    data_ii <- data_i %>% filter(measure == i, X == dec)
    #data_ii_crop <- st_crop(data_ii, xmin = 25, xmax = 33, ymin = -34, ymax = -25)
    
    lims <- range(data_i %>% filter(measure == i) %>% pull(mid))
    
    # ### plot this index's change in the first decade with inset
    print(ggplot(data = data_ii) +
      geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
      geom_sf(data = data_ii, color = NA, aes(fill = mid)) +
      scale_fill_gradientn(colors = c("#f1b6da", "#de77ae", "#8e0152"), na.value = "grey", limits = lims) +
      geom_sf(data = world, size = 0.002, fill = NA) +
      ggtitle(paste(dec, i)) +
      coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
      labs(fill = "Average Change") +
      theme_bw())
    # 
    # print(ggplot(data = data_ii_crop) +
    #   geom_sf(data = world_crop, size = 0.002, fill = "#e1e7ce") +
    #   geom_sf(data = data_ii_crop, color = NA, aes(fill = mid*100)) +
    #   scale_fill_continuous(high = "#ef8a62", na.value = "grey") +
    #   ggtitle(paste(dec, i)) +
    #   coord_sf(ylim = c(-25, -34), xlim = c(25, 33)) +
    #   labs(fill = "Average Change") +
    #   theme_bw())
    
    ## plot the model certainty with inset
    # print(ggplot(data = data_ii) +
    #   geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
    #   geom_sf(data = data_ii, color = NA, aes(fill = confidence)) +
    #   #scale_fill_manual(values = c("#d0d1e6", "#67a9cf", "#3690c0", "02818a", "016c59")) +
    #   ggtitle(paste("Model Agreement", i, dec)) +
    #   coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
    #   #geom_rect(aes(xmin = 25, xmax = 33, ymin = -34, ymax = -25), 
    #   #          color = "purple", linewidth = 0.2, inherit.aes = FALSE, fill = NA) +
    #   labs(fill = "Model Confidence") +
    #   theme_bw())
    
    # print(ggplot(data = data_ii_crop) +
    #   geom_sf(data = world_crop, size = 0.002, fill = "#e1e7ce") +
    #   geom_sf(data = data_ii_crop, color = NA, aes(fill = confidence)) +
    #   #scale_fill_manual(values = c("#01665e", "#02818a", "#3690c0", "67a9cf", "d0d1e6")) +
    #   ggtitle(paste("Model Agreement", i, dec)) +
    #   coord_sf(ylim = c(-25, -34), xlim = c(25, 33), expand = F) +
    #   labs(fill = "Model Confidence") +
    #   theme_bw())
    
    # print(ggplot(data = data_ii %>% arrange(mid),
    #        aes(x = measure, y = mid, ymin = low, ymax = high)) +
    #   geom_point(position = position_dodge2(1)) +
    #   geom_errorbar(width = 1, position = position_dodge2(1)) +
    #   labs(x = "Index", y = "Estimated Range of Change Across All Models and KBAs",
    #        title = paste0("Average 95% CI of Change for all KBAs")) +
    #   theme_bw())

    print(ggplot(data = data_ii %>% st_drop_geometry(), aes(x = measure, group = fct_rev(confidence))) +
            geom_bar(position = "dodge", aes(fill = factor(confidence))) +
            labs(xlab = "Model Agreement", ylab = "Count of KBAs", fill = "Confidence") +
            geom_text(stat='count', aes(label=..count..), vjust= -1, position = position_dodge(width = .9),
                      hjust = 0) +
            scale_fill_manual(values = c("#01665e", "#02818a", "#3690c0", "67a9cf", "d0d1e6")) +
            facet_grid(~.dir) +
            theme_bw())

  }
}

dev.off()

pdf(paste0("./visuals/model_agreement_CIs", td, ".pdf")) ## start pdf up here


for(ind in 1:length(indexes)) {
  i <- indexes[ind]
  data_ii <- data_i %>% filter(measure == i) %>% mutate(ids = as.character(SitRecID)) %>% arrange(mid)
  
  print(ggplot(data = data_ii,
               aes(x = reorder(ids, mid, FUN=mean), y = mid, ymin = low, ymax = high, group = X)) +
          #geom_point(position = position_dodge2(1), aes(color = X)) +
          geom_errorbar(width = 1, position = position_dodge2(1), aes(color = X)) +
          scale_x_discrete(breaks=1:8) +
          labs(x = "Index", y = "Estimated Range of Change Across All Models and KBAs",
               title = paste0("Average 95% CI of Change for all KBAs ", i)) +
          theme_bw())
  
}

dev.off()

#### Figure plot high heat (wdsi + txx + climate threat scatterplots) ----
pdf(paste0("./visuals/f3_scatterplots_", td, ".pdf")) ## start pdf up here

data <- left_join(model_agreement,
                  kbas %>% 
                    select(SitRecID, percent_protected = cum_percPA, protected, climate_threat), 
                  by = "SitRecID") %>%
  mutate(protected = fct_relevel(protected, c("FP", "P", "NP")),
         protected_group = cut(percent_protected, 7, labels = F ))

## plot heat
xlims <- range(data %>% filter(measure == "txxETCCDI") %>% pull(mid))
ylims <- range(data %>% filter(measure == "wsdiETCCDI") %>% pull(mid)) 

for(dec in c("firstdecade", "seconddecade")) {
  data_piv <- data %>% filter(X == dec) %>%
    pivot_wider(id_cols = c(SitRecID, climate_threat, protected, percent_protected),
                                   names_from = measure,
                                   values_from = c("mid"))
  
  m <- lm(wsdiETCCDI~ txxETCCDI, data = data_piv)
  ## scatterplot with climate threat
  print(ggplot(data = data_piv %>% filter(confidence != "unlikely"), 
         aes(x = txxETCCDI, y = wsdiETCCDI)) +
    ggtitle(paste("Average % Change txx + wsdi", dec, "r-squared", round(summary(m)$r.squared, 2))) +
    geom_point(aes(color = climate_threat), size = 2.5) +
    scale_colour_manual(values=c("#f6e8c3", "#01665e")) +
    geom_rug(col=rgb(.5,0,0,alpha=.2)) +
    geom_abline(intercept = m$coefficients[1], slope = m$coefficients[2], color = "darkblue") +
    xlim(xlims) + ylim(ylims) +
    xlab("% Change in txxETCCDI from `95-14") +
    ylab("% Change in wsdiETCCDI \n from `95-14") +
    theme_bw())
  
  ## grid data 
  print(ggplot(data_piv, aes(x = txxETCCDI, y = wsdiETCCDI)) +
    stat_bin2d(aes(fill = after_stat(count), alpha = after_stat(count)), bins = 4, na.rm = T) +
    scale_fill_gradient(high = "#01665e") +
    theme_bw())
}

dev.off()
#### Figure plot wet /dry (cdd, txx + climate threat scatterplots ) ----
pdf(paste0("./visuals/f4_scatterplots_", td, ".pdf")) ## start pdf up here

xlims <- range(data %>% filter(measure == "cddETCCDI") %>% pull(mid))
ylims <- range(data %>% filter(measure == "r95pETCCDI") %>% pull(mid)) 

for(dec in c("firstdecade", "seconddecade")) {
  data_piv <- data %>% filter(X == dec) %>%
    pivot_wider(id_cols = c(SitRecID, climate_threat, protected, percent_protected),
                names_from = measure,
                values_from = c("mid"))
  m <- lm(wsdiETCCDI~ txxETCCDI, data = data_piv)
  
   ## scatterplot with climate threat
  print(ggplot(data = data_piv, 
         aes(x = cddETCCDI, y = r95pETCCDI)) +
    ggtitle(paste("Average % Change cdd + r95p", dec, "r-squared", round(summary(m)$r.squared, 2))) +
    geom_point(aes(color = climate_threat)) +
    scale_colour_manual(values=c("#f6e8c3", "#01665e")) +
    geom_rug(col=rgb(.5,0,0,alpha=.2)) +
    geom_abline(intercept = m$coefficients[1], slope = m$coefficients[2], color = "darkblue") +
    xlim(xlims) + ylim(ylims) +
    xlab("% Change in txxETCCDI from `95-14") +
    ylab("% Change in wsdiETCCDI \n from `95-14") +
    theme_bw())
  
  ## grid data 
  print(ggplot(data_piv, aes(x = cddETCCDI, y = r95pETCCDI)) +
    stat_bin2d(aes(fill = after_stat(count), alpha = after_stat(count)), bins = 4, na.rm = T) +
    scale_fill_gradient(high = "#01665e") +
    theme_bw())
}

dev.off()


#### Figure plot risk rank for KBAs ----
pdf(paste0("./visuals/risk_ranking_", td, ".pdf")) ## start pdf up here

data_i <- left_join(model_agreement_rank_first,
                    kba_geometry %>%
                      filter(SitRecID %in% unique(model_agreement_rank_first$SitRecID)) %>%
                      select(SitRecID, geometry),
                    by = "SitRecID") %>% st_set_geometry("geometry")
top10 <- model_agreement_rank_first[120:160,] %>% pull(SitRecID)

# ### plot this index's change in the first decade with inset
ggplot(data = data_i) +
  geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
  geom_sf(data = data_i, color = NA, aes(fill = risk_group)) +
  scale_fill_gradientn(colors = c("#d73027", "#8856a7")) +
  geom_sf(data = data_i %>% filter(SitRecID %in% top10), fill = "orange", size = 0.2) +
  geom_sf(data = world, size = 0.002, fill = NA) +
  ggtitle("Combined Ranked Risk 2015-2025") +
  coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
  labs(fill = "Risk Group") +
  theme_bw()

write.csv(model_agreement_rank_first[1:10,], "./processed_data/model_agreement/ranking_2015_2025.csv")

data_i <- left_join(model_agreement_rank_second,
                    kba_geometry %>%
                      filter(SitRecID %in% unique(model_agreement_rank_first$SitRecID)) %>%
                      select(SitRecID, geometry),
                    by = "SitRecID") %>% st_set_geometry("geometry")
top10 <- model_agreement_rank_second[120:160,] %>% pull(SitRecID)

# ### plot this index's change in the first decade with inset
ggplot(data = data_i) +
  geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
  geom_sf(data = data_i, color = NA, aes(fill = risk_group)) +
  scale_fill_gradientn(colors = c("#d73027", "#8856a7")) +
  geom_sf(data = data_i %>% filter(SitRecID %in% top10), fill = "orange", size = 0.2) +
  geom_sf(data = world, size = 0.002, fill = NA) +
  ggtitle("Combined Ranked Risk 2026-2036") +
  coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
  labs(fill = "Risk Group") +
  theme_bw()

write.csv(model_agreement_rank_second[1:10,], "./processed_data/model_agreement/ranking_2026_2036.csv")


dev.off()


#### Supp Figure - Spread of % change (old fig 2)----
pdf(paste0("./visuals/f2_boxplot_", td, ".pdf")) ## start pdf up here

ggplot(data = extreme_comp_data,
       aes(x=measure, y = diff_first, color = measure)) +
  geom_boxplot() +
  labs(x = "Date", y = "% Change from Historical", title = "% Change in Index Historical to 2015-25") +
  facet_wrap(~category, scales = "free", nrow = 2) +
  theme_bw()
ggplot(data = extreme_comp_data, aes(x=measure, y = diff_second, color = measure)) +
  geom_boxplot() +
  labs(x = "Date", y = "% Change from Historical", title = "% Change in Index Historical to 2026-36") +
  facet_wrap(~category, scales = "free", nrow = 2) +
  theme_bw()

ggplot(data = extreme_comp_data,
       aes(x=measure, y = diff_abs_first, color = measure)) +
  geom_boxplot() +
  labs(x = "Date", y = "Absolute Change from Historical", title = "Absolute Change in Index Historical to 2015-25") +
  facet_wrap(~category, scales = "free", nrow = 2) +
  theme_bw()
ggplot(data = extreme_comp_data, aes(x=measure, y = diff_abs_second, color = measure)) +
  geom_boxplot() +
  labs(x = "Date", y = "Absolute Change from Historical", title = "Absolute Change in Index Historical to 2026-36") +
  facet_wrap(~category, scales = "free", nrow = 2) +
  theme_bw()


dev.off()
#### Loop through and plot 2 decades of % change ----
for(index in indexes) {
  
  # select correct measure
  data <- extreme_comp_data %>% filter(measure == paste0(index, "ETCCDI"))
  model_agreement %>% pivot_wider(id_cols = c(SitRecID, measure), 
                                  names_from = X,
                                  values_from = c(low, high, over0, under0, over100))
  data <- left_join(data, model_agreement, by = c("SitRecID", "measure"))
  data_i <- left_join(data,
                      kba_geometry %>% filter(SitRecID %in% unique(data$SitRecID)),
                      by = "SitRecID") %>% st_set_geometry("geometry")
  
  
  print(ggplot(data = data_i) +
          ggtitle(paste("Average 2015-2025", index)) +
          geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index_first)) +
          coord_sf(ylim = c(-22, -35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", na.value = "grey") +
          facet_wrap(~ scenario.y, nrow = 3) +
          theme_bw())
  
  print(ggplot(data = data_i) +
          ggtitle(paste("Average 2026-2036", index)) +
          geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index_second)) +
          coord_sf(ylim = c(-22, -35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", na.value = "grey") +
          facet_wrap(~ scenario.y, nrow = 3) +
          theme_bw())
  
  print(ggplot(data = data_i) +
          ggtitle(paste("Diff Historical + 2015-2025", index)) +
          geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
          geom_sf(data = data_i, size = 0.0002, aes(fill = diff_first)) +
          coord_sf(ylim = c(-22, -35)) +
          labs(fill = "Change in Index") +
          scale_fill_gradient(low = "#998ec3", high = "#f1a340", na.value = "grey") +
          facet_wrap(~ scenario.y, nrow = 3) +
          theme_bw())
  
  print(ggplot(data = data_i) +
          ggtitle(paste("Diff Historical + 2026-2036", index)) +
          geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
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
#### Scratch ----
m <- model_agreement_rank_first[,7:11]
model_agreement_rank_first$highest_rank <- colnames(m)[apply(m,1,which.max)]
m <- model_agreement_rank_second[,7:11]
model_agreement_rank_second$highest_rank <- colnames(m)[apply(m,1,which.max)]
pdf("visuals/highest_rank.pdf")

data_i <- left_join(model_agreement_rank_first,
                    kba_geometry %>%
                      filter(SitRecID %in% unique(model_agreement_rank_first$SitRecID)) %>%
                      select(SitRecID, geometry),
                    by = "SitRecID") %>% st_set_geometry("geometry")

# ### plot this index's change in the first decade with inset
ggplot(data = data_i) +
  geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
  geom_sf(data = data_i, color = NA, aes(fill = highest_rank)) +
  #scale_fill_gradientn(colors = c("#e0ecf4", "#8856a7")) +
  geom_sf(data = world, size = 0.002, fill = NA) +
  ggtitle("Highest Ranked Risk 2015-2025") +
  coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
  labs(fill = "Risk Group") +
  theme_bw()

data_i <- left_join(model_agreement_rank_second,
                    kba_geometry %>%
                      filter(SitRecID %in% unique(model_agreement_rank_second$SitRecID)) %>%
                      select(SitRecID, geometry),
                    by = "SitRecID") %>% st_set_geometry("geometry")

# ### plot this index's change in the first decade with inset
ggplot(data = data_i) +
  geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
  geom_sf(data = data_i, color = NA, aes(fill = highest_rank)) +
  #scale_fill_gradientn(colors = c("#e0ecf4", "#8856a7")) +
  geom_sf(data = world, size = 0.002, fill = NA) +
  ggtitle("Combined Ranked Risk 2026-2036") +
  coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
  labs(fill = "Risk Group") +
  theme_bw()

dev.off()
# #### Plot maps  ----
# pdf(paste0("./visuals/maps_", td, ".pdf")) ## start pdf up here
# 
# ids <- c("r95p", "tx90p", "wsdi", "cwd", "cdd", "txx")
# for(index in ids) {
# 
# #### data setup ----
#   m <- paste0(index, "ETCCDI")
#   data <- left_join(extreme_comp_data %>% filter(measure == m),
#                     kbas %>% select(SitRecID, percent_protected = cum_percPA, protected)) %>%
#     filter(measure == measure) %>%
#     mutate(protected = fct_relevel(protected, c("FP", "P", "NP")),
#            protected_group = cut(percent_protected, 7, labels = F ))
# 
#   data_i <- left_join(data,
#                       kba_geometry %>%
#                         filter(SitRecID %in% unique(data$SitRecID)) %>%
#                         select(SitRecID, geometry),
#                       by = "SitRecID") %>% st_set_geometry("geometry") %>%
#     mutate(climate_threat_fill = ifelse(climate_threat, TRUE, NA)) %>%
#     mutate(future_fill = ifelse(diff_first >= .5, "first",
#                                 ifelse(diff_second >= .5, "second", "none"))) %>%
#     mutate(future_fill = fct_relevel(future_fill, c("none", "second", "first")))
# 
#   ggplot(data = data_i) +
#     ggtitle(paste("KBAs at risk for >50% Change in ", index)) +
#     geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
#     geom_sf(data = data_i, size = 0.0002, aes(fill = future_fill)) +
#     coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
#     labs(fill = "Era for Risk") +
#     scale_fill_manual(values = c("#9ebcda", "#8856a7", "#810f7c"),
#                         na.value = "grey") +
#     theme_bw()
# }
# 
# dev.off()

# ####  extremes maps ----
# lims <- range(data_i$mean_index, data_i$mean_index_first, data_i$mean_index_second)
#
# a <- (ggplot(data = data_i) +
#         ggtitle(paste("Historic ('01-11)", index)) +
#         geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
#         geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index)) +
#         coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
#         labs(fill = "Index Value") +
#         scale_fill_gradient(low = "#998ec3", high = "#f1a340",
#                             na.value = "grey", limits = lims) +
#         theme_bw())
#
# b <- (ggplot(data = data_i) +
#         ggtitle(paste("'15-25", index)) +
#         geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
#         geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index_first)) +
#         coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
#         labs(fill = "Index Value") +
#         scale_fill_gradient(low = "#998ec3", high = "#f1a340",
#                             na.value = "grey", limits = lims) +
#         theme_bw())
#
# c <- (ggplot(data = data_i) +
#         ggtitle(paste("26-36", index)) +
#         geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
#         geom_sf(data = data_i, size = 0.0002, aes(fill = mean_index_second)) +
#         coord_sf(ylim = c(-22, -35), xlim = c(15, 35)) +
#         labs(fill = "Index Value") +
#         scale_fill_gradient(low = "#998ec3", high = "#f1a340",
#                             na.value = "grey", limits = lims) +
#         theme_bw())
#
# print(grid.arrange(a, b, c))
#
# #### scatter plots ----
# ## historic vs next 10
# lims <- range(data$diff_first, data$diff_second)
# first <- ggplot(data = data, aes(x = diff_first, y = percent_protected)) +
#   ggtitle(paste("% Change '15-25", m)) +
#   geom_point(aes(color = climate_threat)) + xlim(lims) +
#   geom_rug(col=rgb(.5,0,0,alpha=.2)) +
#   xlab("% Change in Index from `01-11") + ylab("% Area protected") +
#   theme_bw()
# second <- ggplot(data = data, aes(x = diff_second, y = percent_protected)) +
#   ggtitle(paste("% Change'26-36", m)) +
#   geom_point(aes(color = climate_threat)) + xlim(lims) +
#   geom_rug(col=rgb(.5,0,0,alpha=.2)) +
#   xlab("% Change in Index from `01-11") + ylab("% Area protected") +
#   theme_bw()
#
# print(grid.arrange(first, second))
#
# #### density plots ----
# #historic vs next 10, but protected
# lims <- range(0, data$diff_first, data$diff_second)
# first <- ggplot(data = data, aes(x = diff_first, y = protected_group, group = protected_group)) +
#   geom_density_ridges(aes(fill = protected_group)) +
#   scale_fill_gradient(low = "#762a83", high = "#1b7837", na.value = "grey", name = "Protected Group \n (7 High)") +
#   xlim(lims) +
#   ggtitle(paste("% Change '15-25", m)) +
#   xlab("% Change in Index from `01-11") + ylab("% Area protected") +
#   theme_ridges()
# second <- ggplot(data = data, aes(x = diff_second, y = protected_group, group = protected_group)) +
#   geom_density_ridges(aes(fill = protected_group)) +
#   scale_fill_gradient(low = "#762a83", high = "#1b7837", na.value = "grey", name = "Protected Group \n (7 High)") +
#   xlim(lims) +
#   ggtitle(paste("% Change '26-36", m)) +
#   xlab("% Change in Index from `01-11") + ylab("% Area protected") +
#   theme_ridges()
#
# print(grid.arrange(first, second))
#




dev.off()

# #### Do some ANOVAs that will plot out differences in models ----
# 
# for(s in unique(plot_data$scenario)) {
#   if(scenario == "historic") next
#   for(m in unique(plot_data$measure)){
#     d <- plot_data %>% filter(scenario == s & measure == m, 
#                               year %in% c(2015:2025))
#     aov(mean ~ year * gcm, data = d)
#     
#   }
# }
# 
# 
# #### Loop through and make the same figures for each measure
# next_10 <- plot_data %>% group_by(SitRecID, year_group, measure) %>%
#   summarize(mean_sd = sd(mean), mean_index = mean(mean)) %>% filter(year_group < 3) %>% 
#   pivot_wider(id_cols = c(SitRecID, measure), names_from = year_group, values_from = c("mean_sd", "mean_index")) %>% 
#   mutate(diff_mean = mean_index_2 - mean_index_1, 
#          climate_threat = ifelse(SitRecID %in% ids, "Y", "N"))
# 
# for(index in unique(next_10$measure)) {
#   
#   ### Average of index for every 10 years (about)
#   next_10_i <- next_10 %>% filter(measure == index)
#   next_10_i <- left_join(next_10_i, kba_geometry) %>% st_set_geometry("geometry")
#   
#   print(ggplot(data = next_10_i) +
#     ggtitle(paste("Change in Index (avg 2036-2026 minus 2025-2015) ", index)) +
#     geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
#     geom_sf(data = next_10_i, size = 0.0002, aes(fill = diff_mean)) +
#     coord_sf(ylim = c(-22, -35)) +
#     scale_fill_continuous(na.value = "grey") +
#     labs(fill = "Change in Index") +
#     theme_bw())
#   
#   print(ggplot(data = next_10_i) +
#     ggtitle(paste("SD 2025-2015 ", index)) +
#     geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
#     geom_sf(data = next_10_i, size = 0.0002, aes(fill = mean_sd_1)) +
#     coord_sf(ylim = c(-22, -35)) +
#     scale_fill_continuous(na.value = "grey") +
#     labs(fill = "Standard Deviation") +
#     theme_bw())
#   
#   print(ggplot(data = next_10_i) +
#     ggtitle(paste("SD 2026-2036 ", index)) +
#     geom_sf(data = world, size = 0.002, fill = "#e1e7ce") +
#     geom_sf(data = next_10_i, size = 0.0002, aes(fill = mean_sd_2)) +
#     coord_sf(ylim = c(-22, -35)) +
#     scale_fill_continuous(na.value = "grey") +
#     labs(fill = "Standard Deviation") +
#     theme_bw())
#   
# }
# 
# dev.off()
# 
# 
# 
