#' =============================================================================
#' Project: ECHO Aim 1 CEI BKMR Anayysis
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' This project is a follow up the the CEI paper (Env Epi, 2018). Here we're 
#' using BKMR to identify which exposures are driving the associations reported
#' in that original paper
#' 
#' Date Created: October 16, 2019
#' 
#' Description: This code assigns exposures to each HS participant
#' =============================================================================

library(sf)
library(raster)
library(ggplot2)
library(ggmap)
library(ggsn)
library(ggthemes)
library(stringr)
library(tidyverse)
library(lubridate)
library(readxl)

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
windowsFonts(Calibri=windowsFont("TT Calibri"))
options(scipen = 9999) #avoid scientific notation

geo_data <- "R:/RSTOR-Magzamen/Research/Secondary_Data/ACS_Geodatabases/"
utm_13 <- "+init=epsg:26913"
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
ll_nad83 <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
ll_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#' -----------------------------------------------------------------------------
#' Assign AP exposures based on census tract, conception, and delivery
#' -----------------------------------------------------------------------------

load(here::here("Data", "hs_spatial.RData"))
names(hs)

load(here::here("Data", "CT_Air Pollution.RData"))
ct_air_pollution <- st_set_geometry(ct_air_pollution, NULL) %>%
  mutate(week_ending = as.Date(week_ending))

hs_df <- st_set_geometry(hs, NULL) %>%
  select(pid, GEOID, exp_start, exp_end) %>%
  mutate(exp_start = as.Date(exp_start),
         exp_end = as.Date(exp_end)) %>%
  as.data.frame()

#' Add exposure concentrations to hs dataset based on the two-week periods
#' when conception and delivery occurred
for (i in 1:nrow(hs_df)) {
  
  #' subset air pollution data based on which census tract the mom lives in
  ap <- filter(ct_air_pollution, GEOID == as.character(hs_df[i, "GEOID"]))
  
  #' get concentrations for the week-ending dates of her pregnancy
  date_seq <- seq.Date(from = as.Date(hs_df[i, "exp_start"]),
                       to = as.Date(hs_df[i, "exp_end"]),
                       by  = "week")
  pm <- select(ap, week_ending, biweekly_average_pm_pred) %>%
    filter(week_ending %in% date_seq)
  
  o3 <- select(ap, week_ending, biweekly_average_o3_pred) %>%
    filter(week_ending %in% date_seq)
  
  #' Average exposures over the duration of the pregnancy
  hs_df[i,"mean_pm"] <- mean(pm$biweekly_average_pm_pred, na.rm=T)
  hs_df[i,"mean_o3"] <- mean(o3$biweekly_average_o3_pred, na.rm=T)
  hs_df[i,"max_o3"] <- max(o3$biweekly_average_o3_pred, na.rm=T)
  
  #' Average_exposures by trimester
  hs_df[i,"mean_pm_first_tri"] <- mean(pm$biweekly_average_pm_pred[1:7], na.rm=T)
  hs_df[i,"mean_o3_first_tri"] <- mean(o3$biweekly_average_o3_pred[1:7], na.rm=T)
  hs_df[i,"max_o3_first_tri"] <- max(o3$biweekly_average_o3_pred[1:7], na.rm=T)
  
  hs_df[i,"mean_pm_second_tri"] <- mean(pm$biweekly_average_pm_pred[8:14], na.rm=T)
  hs_df[i,"mean_o3_second_tri"] <- mean(o3$biweekly_average_o3_pred[8:14], na.rm=T)
  hs_df[i,"max_o3_second_tri"] <- max(o3$biweekly_average_o3_pred[8:14], na.rm=T)
  
  hs_df[i,"mean_pm_third_tri"] <- mean(pm$biweekly_average_pm_pred[15:nrow(pm)], na.rm=T)
  hs_df[i,"mean_o3_third_tri"] <- mean(o3$biweekly_average_o3_pred[15:nrow(o3)], na.rm=T)
  hs_df[i,"max_o3_third_tri"] <- max(o3$biweekly_average_o3_pred[15:nrow(o3)], na.rm=T)
}

hs_df <- select(hs_df, pid, mean_pm:max_o3_third_tri)
summary(hs_df)

#' Drop participants (N=11) with missing air pollution data 
#' (living the the farthest east census tracts-- too far for kriging)
#' N = 1149 (81.4%) participants with full air pollution data
hs <- left_join(hs, hs_df, by="pid") %>%
  filter(!(is.na(mean_pm)))
summary(hs)

#' -----------------------------------------------------------------------------
#' Assign ENV and SOC exposures based on census tract
#' -----------------------------------------------------------------------------

load(here::here("Data", "CT_Environmental_Updated.RData"))
ct_env <- st_set_geometry(ct_env, NULL)

hs <- left_join(hs, ct_env, by="GEOID") %>% 
  #' reverse code tree cover so that less tree cover is higher ranked
  mutate(pct_no_tree_cover = 100 - pct_tree_cover)

load(here::here("Data", "CT_Social.RData"))
ct_soc <- st_set_geometry(ct_soc, NULL) %>%
  #' Inverse median income so that smaller incomes are ranked higher
  #' (more disadvantaged)
  mutate(inv_med_income = 1/med_income)

hs <- left_join(hs, ct_soc, by="GEOID")

summary(hs)

#' -----------------------------------------------------------------------------
#' Save the exposure data
#' -----------------------------------------------------------------------------

save(hs, file = here::here("Data", "hs_exposures.RData"))