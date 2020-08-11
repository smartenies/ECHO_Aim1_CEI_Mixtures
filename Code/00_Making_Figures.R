#' =============================================================================
#' Project: ECHO Aim 1 CEI BKMR Anayysis
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' Date Created: October 1, 2019
#' 
#' Description: This script summarizes the ACS variables used in the analysis
#'  
#' This project is a follow up the the CEI paper (Env Epi, 2018). Here we're 
#' using BKMR to identify which exposures are driving the associations reported
#' in that original paper
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

utm_13 <- "+init=epsg:26913"
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
ll_nad83 <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
ll_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#' -----------------------------------------------------------------------------
#' Mapping some key demographic and SES variables for Colorado
#' -----------------------------------------------------------------------------
library(viridis)

load(file="./Data/Spatial Data/acs.RData")
acs <- mutate_if(acs, is.integer, as.numeric)

ses_vars <- c("pop_dens", "pct_under5", "pct_over64", "pct_poc", "pct_fb",
              "pct_less_hs", "pct_unemp", "pct_limited_eng", "pct_hh_pov",
              "med_income", "med_year_blt_housing")
titles <- c("Population density (persons per km2)",
            "Population under 5 years of age",
            "Population 65 years of age or older",
            "Persons of color", "Foreign born population",
            "Adults 25 or older with less than a high school diploma",
            "Percentage of the civilan workforce (ages 16+) that are unemployed",
            "Households that speak limited English",
            "Households with last year income below the poverty level",
            "Median household income (2014$)",
            "Housing median year built")

#' First, across all of CO, then, in the three county area
ses_plot_fn <- function(df, var_string, name_string) {
  ggplot() +
    ggtitle(titles[i]) +
    geom_sf(data = df, aes_string(fill = var_string), col=NA) +
    scale_fill_viridis(name = name_string) +
    xlab("") + ylab("") +
    theme(legend.position = "right") +
    simple_theme
}

for (i in 1:length(ses_vars)) {
  ses_var <- ses_vars[i]
  ses_name <- ifelse(ses_var=="med_income", "2014$", 
                     ifelse(ses_var=="pop_dens", "Persons per km\u00B2",
                            ifelse(ses_var=="med_year_blt_housing", "Year",
                                   "Percentage\nof census tract\npopulation:")))
  
  ses_plot_fn(df = acs, var_string = ses_var, name_string = ses_name)
  ggsave(filename = paste("./Figures/CEI Figures/ACS Variables/", "co ", ses_var, " ", years, ".jpeg", sep=""), 
         device = "jpeg", dpi=600)
}

load("./Data/Spatial Data/dm_tracts.RData")
dm_ct_list <- unique(dm_tracts$GEOID)
acs_dm <- acs %>%
  filter(GEOID %in% dm_ct_list)

for (i in 1:length(ses_vars)) {
  ses_var <- ses_vars[i]
  ses_name <- ifelse(ses_var=="med_income", "2014$", 
                     ifelse(ses_var=="pop_dens", "Persons per km\u00B2",
                            ifelse(ses_var=="med_year_blt_housing", "Year",
                                   "Percentage\nof census tract\npopulation:")))
  
  ses_plot_fn(df = acs_dm, var_string = ses_var, name_string = ses_name)
  ggsave(filename = paste("./Figures/CEI Figures/ACS Variables/", "dm ", ses_var, " ", years, ".jpeg", sep=""), 
         device = "jpeg", dpi=600)
}

#' -----------------------------------------------------------------------------
#' Map the crime rates
#' -----------------------------------------------------------------------------

crime_rates <- read_csv(here::here("Data", "Crime_AEA.csv")) %>% 
  st_as_sf(wkt = "WKT", crs = albers)

ggplot(data = filter(crime_rates, GEOID != "LAKESIDE")) +
  ggtitle("Crime rate per 1,000 (all reported incidents)") +
  geom_sf(aes(fill = crime_rate), col=NA) +
  scale_fill_viridis(name = "All crimes\nRate per 1,000") +
  xlab("") + ylab("") +
  theme(legend.position = "right") +
  simple_theme
ggsave(filename = here::here("Figs/CEI_Inputs", "All_Crime_Rates_2010_2014.jpeg"), 
       device = "jpeg", dpi=600)

ggplot(data = filter(crime_rates, GEOID != "LAKESIDE")) +
  ggtitle("Crime rate per 1,000 (violent crimes)") +
  geom_sf(data = crime_rates, aes(fill = violent_rate), col=NA) +
  scale_fill_viridis(name = "Violent crimes\nRate per 1,000") +
  xlab("") + ylab("") +
  theme(legend.position = "right") +
  simple_theme
ggsave(filename = here::here("Figs/CEI_Inputs", "Violent_Crime_Rates_2010_2014.jpeg"), 
       device = "jpeg", dpi=600)

ggplot(data = filter(crime_rates, GEOID != "LAKESIDE")) +
  ggtitle("Crime rate per 1,000 (property and non-violent crimes)") +
  geom_sf(data = crime_rates, aes(fill = property_rate), col=NA) +
  scale_fill_viridis(name = "Property and\nnon-violent crimes\nRate per 1,000") +
  xlab("") + ylab("") +
  theme(legend.position = "right") +
  simple_theme
ggsave(filename = here::here("Figs/CEI_Inputs", "Property_Crime_Rates_2010_2014.jpeg"), 
       device = "jpeg", dpi=600)