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
#' Description: This code reads in the HS participant data set and makes sure 
#' its good to go (data cleaned, with the correct CT assignment)
#' 
#' Update 7/7/2020: We decided to add census-tract level %POC to the exposures
#' and use lon and lat data for each residential location as covariates to try
#' to account for spatial confounding
#' 
#' We also need a season of conception indicator variable
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
library(haven)
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
#' Read in HS dataset from previous study
#' -----------------------------------------------------------------------------

load(here::here("Data", "hs_data_for_cei.RData"))
names(hs_data)

#' Get conception and delivery dates, as well as "week_ending" dates
hs_data <- select(hs_data, pid, conception_date, di1) %>%
  mutate(conception_date = format(as.Date(conception_date, origin="1960-01-01"),"%Y-%m-%d")) %>%
  mutate(conception_date = as.Date(conception_date)) %>%
  mutate(exp_start = ceiling_date(conception_date, unit = "week"),
         exp_end = ceiling_date(di1, unit = "week")) %>%
  mutate(pregnancy_weeks = difftime(di1, conception_date, unit="week")) %>%
  distinct()

summary(as.numeric(hs_data$pregnancy_weeks))
hist(as.numeric(hs_data$pregnancy_weeks))

#' Drop negative pregnancy weeks (conception date is wrong)
#' and NAs for pregnancy duration (no delivery date)
#' N = 1359 (96.4%) with plausible pregnancy durations
hs_data <- filter(hs_data, !(is.na(pregnancy_weeks))) %>%
  filter(pregnancy_weeks > 0) %>%
  filter(pregnancy_weeks < 50)

summary(as.numeric(hs_data$pregnancy_weeks))
hist(as.numeric(hs_data$pregnancy_weeks))
nrow(hs_data)

#' Get rid of the SAS formats and labels
hs_data <- zap_formats(hs_data) %>%
  zap_labels()

#' Read in geocoded participant addresses
#' N = 1334 geocoded addresses (95%)
load(here::here("Data", "hs_geocode.RData"))

#' Merge HS data with HS locations
hs_all <- left_join(hs_geocode, hs_data, by="pid") %>%
  filter(!(is.na(pregnancy_weeks)))

#' Assign a CT value
load(here::here("Data", "dm_tracts.RData"))
dm_tracts <- select(dm_tracts, GEOID)
head(dm_tracts)

plot(dm_tracts)
nrow(dm_tracts)

#' N = 1202 participants (85%) have geocoded addresses AND live in
#' census tracts within the area
hs <- hs_all %>%
  st_intersection(dm_tracts)
head(hs)
names(hs)

plot(st_geometry(dm_tracts))
plot(st_geometry(hs), col="red", add=T)

summary(as.numeric(hs$pregnancy_weeks))
hist(as.numeric(hs$pregnancy_weeks))

save(hs, file = here::here("Data", "hs_spatial.RData"))








