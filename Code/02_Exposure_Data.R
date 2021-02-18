#' =============================================================================
#' Project: ECHO Aim 1 CEI BKMR Analysis
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' This project is a follow up the the CEI paper (Env Epi, 2018). Here we're 
#' using BKMR to identify which exposures are driving the associations reported
#' in that original paper
#' 
#' Date Created: October 16, 2019
#' 
#' Description: This code reads in the original census-tract level exposure data 
#' from the CEI paper, generates some new variables, and summarizes the dataset 
#' 
#' Update 7/7/2020: We decided to add census-tract level %POC to the exposures
#' and use lon and lat data for each residential location as covariates to try
#' to account for spatial confounding
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
#' Read in the original data sets (copied over from the other project folder)
#' -----------------------------------------------------------------------------

load(here::here("Data", "CT_Air Pollution.RData"))
load(here::here("Data", "CT_Temperature.RData"))
load(here::here("Data", "CT_Environmental.RData"))
load(here::here("Data", "CT_Social.RData"))

names(ct_air_pollution)
names(ct_temperature)
names(ct_env)
names(ct_soc)

# plot(st_geometry(ct_air_pollution))

summary(ct_air_pollution)
summary(ct_temperature)
summary(ct_env)
summary(ct_soc)

#' -----------------------------------------------------------------------------
#' Add some additional measures (because there are a lot of 0's for some of the
#' built environment variables)
#' 
#' New metrics include:
#'     - Min distance to TRI site
#'     - Min distance to NPL site
#'     - Min distnace to waste site
#'     - Min distance to major emitters
#'     - Min distance to CAFO
#'     - Min distance to mine/well
#' -----------------------------------------------------------------------------

#' Census tracts and centroids
load(here::here("Data", "dm_tracts.RData"))
load(here::here("Data", "dm_centroids.RData"))

names(dm_cent)

dm_cent <- select(dm_cent, GEOID)
head(dm_cent)

plot(st_geometry(dm_tracts), border = "blue")
plot(st_geometry(dm_cent), fill = "black", add = T)

#' TRI Inventory
load(here::here("Data", "tri_inventory.RData"))
names(tri)
plot(st_geometry(tri_by_facility))

#' NPL Sites
load(here::here("Data", "npl.RData"))
plot(st_geometry(npl))

#' Remove "NPL Deleted" sites
#' N = 3 that are filtered out
npl <- filter(npl, SiteStatus != "NPL DELETED")
plot(st_geometry(npl))

#' Waste sites
load(here::here("Data", "wwtf.RData"))
load(here::here("Data", "lf.RData"))
load(here::here("Data", "compost.RData"))

wwtf <- select(wwtf, Permit_Nam) %>%
  rename(facility_id = Permit_Nam)

lf <- select(lf, PGM_SITE_N) %>%
  rename(facility_id = PGM_SITE_N)

compost <- select(compost, FACILITY_N) %>%
  rename(facility_id = FACILITY_N)

waste_sites <- rbind(wwtf, lf, compost)
plot(st_geometry(waste_sites))

#' Major emittors
load(here::here("Data", "emissions_inventory.RData"))

inventory_major <- filter(inventory_criteria, major == 1) %>%
  select(site_id) %>%
  distinct(site_id, .keep_all = T)

#' CAFOs
load(here::here("Data", "cafo.RData"))

#' Mines and wells
load(here::here("Data", "mines.RData"))
load(here::here("Data", "wells.RData"))

mines <- select(mines, rec_no) %>%
  rename(id = rec_no) %>%
  distinct

wells <- select(wells, API) %>%
  rename(id = API) %>%
  distinct

mines_wells <- rbind(mines, wells)

#' Calculate minimum distance from the census tract centroid

# Empty data frame for new metrics
new_metrics <- data.frame()

for (i in 1:nrow(dm_cent)) {
  print(paste("Location", i, "of", nrow(dm_cent)))
  point <- dplyr::slice(dm_cent, i)
  
  #'     - Min distance to TRI site
  #'     - Min distance to NPL site
  #'     - Min distnace to waste site
  #'     - Min distance to major emittors
  #'     - Min distance to CAFO
  #'     - Min distance to mine/well
  
  met_temp <- point %>%  
    mutate(dist_m_tri = min(unclass(st_distance(point, tri))),
           dist_m_npl = min(unclass(st_distance(point, npl))),
           dist_m_waste_site = min(unclass(st_distance(point, waste_sites))),
           dist_m_major_emit = min(unclass(st_distance(point, inventory_major))),
           dist_m_cafo = min(unclass(st_distance(point, cafo))),
           dist_m_mine_well = min(unclass(st_distance(point, mines_wells)))) %>%
    st_set_geometry(NULL)
  
  #' Add this location's new metrics to the data frame
  new_metrics <- rbind(new_metrics, met_temp)
  rm(met_temp, point)
}

#' drop the geometry from the new metrics
names(new_metrics)

#' Add these new metrics to ct_env
ct_env <- left_join(ct_env, new_metrics, by = "GEOID") 

names(ct_env)

save(ct_env, file = here::here("Data", "CT_Environmental_Updated.RData"))




