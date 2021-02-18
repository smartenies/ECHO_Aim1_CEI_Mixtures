#' =============================================================================
#' Project: ECHO Aim 1 
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' Description:
#'  
#' This project examines the relationships between spatially-distributed
#' economic, environmental, and social variables and health outcomes meausred
#' in the Healthy Start cohort (UC Denver)
#' 
#' This script generates maps and plots for the manuscript
#' =============================================================================

library(sf)
library(raster)
library(ggplot2)
library(ggmap)
library(ggsn)
library(ggthemes)
library(ggcorrplot)
library(stringr)
library(tidyverse)
library(lubridate)
library(readxl)
library(viridis)

#' For ggplots
simple_theme <- theme(
  #aspect.ratio = 1,
  text  = element_text(family="Calibri",size = 12, color = 'black'),
  panel.spacing.y = unit(0,"cm"),
  panel.spacing.x = unit(0.25, "lines"),
  panel.grid.minor = element_line(color = "transparent"),
  panel.grid.major = element_line(color = "transparent"),
  panel.border=element_rect(fill = NA),
  panel.background=element_blank(),
  axis.ticks = element_line(colour = "black"),
  axis.text = element_text(color = "black", size=10),
  # legend.position = c(0.1,0.1),
  plot.margin=grid::unit(c(0,0,0,0), "mm"),
  legend.key = element_blank()
)
spaghetti_theme <- theme(
  #aspect.ratio = 1,
  text  = element_text(family="Calibri",size = 12, color = 'black'),
  panel.spacing.y = unit(0,"cm"),
  panel.spacing.x = unit(0.25, "lines"),
  panel.grid.minor = element_line(color = "transparent"),
  panel.grid.major = element_line(color = "transparent"),
  panel.border=element_rect(fill = NA),
  panel.background=element_blank(),
  axis.ticks.y = element_line(colour = "black"),
  axis.text.y = element_text(color = "black", size=10),
  axis.ticks.x = element_blank(),
  axis.text.x = element_blank(),
  # legend.position = c(0.1,0.1),
  plot.margin=grid::unit(c(0,0,0,0), "mm"),
  legend.key = element_blank()
)
map_theme <- theme(
  #aspect.ratio = 1,
  text  = element_text(family="Calibri",size = 12, color = 'black'),
  panel.spacing.y = unit(0,"cm"),
  panel.spacing.x = unit(0.25, "lines"),
  panel.grid.minor = element_line(color = "transparent"),
  panel.grid.major = element_line(color = "transparent"),
  panel.border=element_rect(fill = NA),
  panel.background=element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  # legend.position = c(0.1,0.1),
  plot.margin=grid::unit(c(0,0,0,0), "mm"),
  legend.key = element_blank()
)
windowsFonts(Calibri=windowsFont("TT Calibri"))
options(scipen = 9999) #avoid scientific notation

geo_data <- "T:/Rsch-MRS/ECHO/SEM Large Data/Spatial Data/"
utm_13 <- "+init=epsg:26913"
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
ll_nad83 <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
ll_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#' -----------------------------------------------------------------------------
#' Correlation plot for variables
#' -----------------------------------------------------------------------------

hs_data3 <- read_csv(here::here("Data", "HS_Clean_Data.csv")) %>% 
  filter(!is.na(conception_date))
nrow(hs_data3)

#' Add some additional variables
hs_data3$lat_lon_int <- hs_data3$lat * hs_data3$lon
hs_data3$concep_year <- year(hs_data3$conception_date)

unique(hs_data3$concep_year)

hs_data3$concep_2009 <- ifelse(hs_data3$concep_year == "2009", 1, 0)
hs_data3$concep_2010 <- ifelse(hs_data3$concep_year == "2010", 1, 0)
hs_data3$concep_2011 <- ifelse(hs_data3$concep_year == "2011", 1, 0)
hs_data3$concep_2012 <- ifelse(hs_data3$concep_year == "2012", 1, 0)
hs_data3$concep_2013 <- ifelse(hs_data3$concep_year == "2013", 1, 0)

hs_data3 <- filter(hs_data3, include == 1) %>% 
  select(
    #' Exposures
    mean_pm, mean_o3, mean_temp, pct_tree_cover, pct_impervious,
    mean_aadt_intensity, dist_m_tri:dist_m_mine_well, 
    cvd_rate_adj, res_rate_adj, violent_crime_rate, property_crime_rate,
    pct_less_hs, pct_unemp, pct_limited_eng, pct_hh_pov, pct_poc,
    
    #' Covariates     
    latina_re, black_re, other_re, #' ref is white_re
    ed_no_hs, ed_hs, ed_aa, ed_4yr, #' ref is ed_grad
    low_bmi, ovwt_bmi, obese_bmi, #' ref is norm_bmi
    concep_spring, concep_summer, concep_fall, #' ref is concep_winter
    concep_2010, concep_2011, concep_2012, concep_2013, #ref is concep_2009
    maternal_age, any_smoker, smokeSH, mean_cpss, mean_epsd,
    male, gest_age_w, 
    
    #' Outcomes 
    birth_weight) %>% 
  na.omit()
nrow(hs_data3)

hs_exp <- hs_data3 %>%
  select(mean_pm, mean_o3, mean_temp, pct_tree_cover, pct_impervious,
         mean_aadt_intensity, dist_m_tri:dist_m_mine_well, 
         cvd_rate_adj, res_rate_adj, violent_crime_rate, property_crime_rate,
         pct_less_hs, pct_unemp, pct_limited_eng, pct_hh_pov, pct_poc)

cor_mat <- cor(hs_exp)
colnames(cor_mat)
colnames(cor_mat) <- c("Mean PM2.5", "Mean O3", "Mean Temperature", 
                       "% Tree Cover", "% Impervious", "AADT Intensity", 
                       "Dist. to TRI Sites", "Dist. to NPL Sites", "Dist. to Waste Sites",
                       "Dist. to Major Emitters", "Dist. to CAFOs", "Dist. to Gas Wells",
                       "CVD Hosp. Rate", "Resp. Hosp. Rate", "Violent Crime Rate", "Property Crime Rate",
                       "% Less HS Edu.", "% Unemployment", "% HH Limited English", "% HH Poverty",
                       "% Persons of Color")
rownames(cor_mat)
rownames(cor_mat) <- c("Mean PM2.5", "Mean O3", "Mean Temperature", 
                       "% Tree Cover", "% Impervious", "AADT Intensity", 
                       "Dist. to TRI Sites", "Dist. to NPL Sites", "Dist. to Waste Sites",
                       "Dist. to Major Emitters", "Dist. to CAFOs", "Dist. to Gas Wells",
                       "CVD Hosp. Rate", "Resp. Hosp. Rate", "Violent Crime Rate", "Property Crime Rate",
                       "% Less HS Edu.", "% Unemployment", "% HH Limited English", "% HH Poverty",
                       "% Persons of Color")
ggcorrplot(cor_mat, outline.col = "white", 
           colors = c("#FDE725FF", "white", "#482677FF"), lab = T,
           legend.title = "Pearson's\ncoefficient",
           lab_size = 3, digits = 1) +
  geom_vline(xintercept = 12.5, size = 0.75) +
  geom_hline(yintercept = 12.5, size = 0.75) 
ggsave(filename = here::here("figs", "correlation_matrix.jpeg"), device = "jpeg",
       width = 9, height = 7, units = "in", dpi = 500)




