#' =============================================================================
#' Project: ECHO Aim 1 
#' Date created: February 10, 2021
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' Description:
#' 
#' This project examines the relationships between spatially-distributed
#' economic, environmental, and social variables and health outcomes meausred
#' in the Healthy Start cohort (UC Denver)
#' 
#' This script generates kriged estiamtes of biweekly AQS temperature at the 
#' census tract level to be used to assign exposures during pregnancy
#' 
#' NOTE: don't forget the ./ before the directory when reading in files!
#' =============================================================================

library(sf)
library(raster)
library(ggplot2)
library(ggthemes)
library(stringr)
library(tidyverse)
library(lubridate)
library(readxl)
library(writexl)
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
windowsFonts(Calibri=windowsFont("TT Calibri"))
options(scipen = 9999) #avoid scientific notation

utm_13 <- "+init=epsg:26913"
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
ll_nad83 <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
ll_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#' -----------------------------------------------------------------------------
#' Set up objects for kriging
#' #' Note: gstat (for kriging) requires sp objects (not sf)
#' -----------------------------------------------------------------------------

#' Census tract centroids for kriging
load(here::here("data", "dm_tracts.RData"))
ct_temp <- dplyr::select(dm_tracts, GEOID) %>%
  arrange(GEOID) %>%
  st_centroid()

plot(st_geometry(ct_temp))

ct_temperature <- dplyr::select(dm_tracts, GEOID) %>%
  arrange(GEOID)

#' Converting the census tract centroids from sf to sp 
ct_sp <- as(ct_temp, "Spatial")
plot(ct_sp)

#' Monitoring data
#' Read in files
#' Calculate 2-week averages
#' Use daily means for temperature
temp_files <- list.files(here::here("Data/AQS_Data"), pattern = "daily_TEMP",
                       full.names = T)
temp_monitors_raw <- do.call(bind_rows, lapply(temp_files, read_csv))
colnames(temp_monitors_raw) <- str_replace(colnames(temp_monitors_raw), " ", "_")

temp_monitors <- filter(temp_monitors_raw, State_Code == "08") %>% 
#test_monitors <- filter(temp_monitors_raw, State_Code == "08") %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = ll_wgs84) %>% 
  st_transform(crs = albers) %>% 
  mutate(County_Code = str_pad(County_Code, width = 3, side = "left", pad = "0"),
         Site_Num = str_pad(Site_Num, width = 4, side = "left", pad = "0")) %>% 
  mutate(monitor_id = paste0(State_Code, County_Code, Site_Num)) %>% 
  #filter(POC == 1) %>% 
  mutate(year = strftime(Date_Local, format = "%Y"),
         week_starting = as.Date(floor_date(Date_Local, unit = "week")),
         week_ending = as.Date(ceiling_date(Date_Local, unit = "week"))) %>%
  mutate(Date_Local = as.Date(Date_Local)) %>%
  arrange(Date_Local, monitor_id) %>%
  #distinct() %>% 
  dplyr::group_by(monitor_id, week_ending) %>%
  summarize(weekly_average = mean(Arithmetic_Mean)) %>%
  mutate(pollutant = "temp")

head(temp_monitors)
tail(temp_monitors)

summary(temp_monitors)
hist(temp_monitors$weekly_average)

plot(st_geometry(dm_tracts), border="grey50")
plot(st_geometry(ct_ap), add=T, col="black")
plot(st_geometry(temp_monitors), add=T, col="green")

save(temp_monitors, file = here::here("data", "temp_monitor_metrics.RData"))

#' Choosing a kriging cutoff distance:
#'     - At 30 km we lose some of the eastern census tracts
#'     - At 40 km we lose the two eastern-most census tracts, but they only 
#'       have a couple participants
#'     - Going to go with 40 km for now and see how it looks

#' Check date lists

start <- ceiling_date(as.Date("2009-01-01"), unit="week")  + 7
ap_dates <- data.frame(week_ending = seq.Date(from = start, 
                                              to = as.Date("2017-12-31"),
                                              by="2 weeks"))
week_list <- sort(unique(ap_dates$week_ending))

#' -----------------------------------------------------------------------------
#' Biweekly mean temperature (F)
#' ------------------------------------------------------------------------------

library(gstat)
library(automap)

show.vgms()
all_models <- c("Exp", "Sph", "Gau", "Mat", "Ste", "Lin")

#' cutoff distance is 40 km
c_dist = 40000

temp_ct_data <- data.frame()
temp_cv_results <- data.frame()
temp_diagnostics <- data.frame()

for (i in 1:length(week_list)) {
  print(paste("Week", i, "of", length(week_list)))
  
  #' use second week ending date to identify concentrations
  week_start <- week_list[i] - 7
  week_end <- week_list[i]
  
  #' Weekly concentration at monitors
  #' Drop rows with NA values
  temp_week <- filter(temp_monitors, week_ending %in% c(week_start, week_end)) %>%
    filter(!is.na(weekly_average))
  
  if(nrow(temp_week) == 0) {next}
  
  #' Biweekly average for each monitor
  temp_week <- temp_week %>%
    dplyr::select(-week_ending) %>%
    group_by(monitor_id) %>%
    summarize(pollutant = "temp",
              biweekly_average = mean(weekly_average),
              week_ending = week_end)

  #' Converting the monitor points from sf to sp 
  temp_week <- as(temp_week, "Spatial")

  #' Kriging using gstat
  #' First, fit the empirical variogram
  vgm <- variogram(biweekly_average ~ 1, temp_week, cutoff = c_dist)
  #plot(vgm)

  #' Second, fit the model
  vgm_fit <- fit.variogram(vgm, model=vgm(all_models), fit.kappa = seq(.3,5,.01))
  model <- as.character(vgm_fit$model)[nrow(vgm_fit)]
  #plot(vgm, vgm_fit)
  
  #' Third, krige
  ok_result <- krige(biweekly_average ~ 1, temp_week, ct_sp, vgm_fit,
                     maxdist = c_dist)

  #' Fourth, leave-one out cross validation
  cv_result <- krige.cv(biweekly_average ~ 1, temp_week, vgm_fit)
  summary(cv_result)
  
  #' me_mean = mean error divided by mean observed
  #' MSNE = mean square normalized error (mean of squared z scores)
  #' RMSE = root mean squared error (same units as pollutant)
  #' cor_obs_pred = correlation between observed and predicted, should be 1
  #' cor_pred_res = correlation between predicted and residual, should be 0
  cv_compare <- compare.cv(list(krige.cv_output = cv_result))

  #' Data frame of results
  temp <- data.frame(GEOID = ct_sp@data$GEOID,
                     week_ending = week_end,
                     biweekly_average_pred = ok_result$var1.pred,
                     biweekly_average_var = ok_result$var1.var)
  temp_ct_data <- rbind(temp_ct_data, temp)

  #' Data frame of cross-validation results
  cv_result <- as.data.frame(cv_result)
  cv_result$week_ending <- week_end
  temp_cv_results <- rbind(temp_cv_results, cv_result)
  
  temp2 <- data.frame(week_ending = week_end,
                      monitor_n = nrow(temp_week),
                      monitor_min = min(temp_week$biweekly_average, na.rm=T),
                      monitor_max = max(temp_week$biweekly_average, na.rm=T),
                      monitor_mean = mean(temp_week$biweekly_average, na.rm=T),
                      model = model,
                      modeled_min = min(ok_result$var1.pred, na.rm=T),
                      modeled_max = max(ok_result$var1.pred, na.rm=T),
                      modeled_mean = mean(ok_result$var1.pred, na.rm=T),
                      mean_error = unname(unlist(cv_compare[1,1])),
                      me_mean = unname(unlist(cv_compare[2,1])),
                      msne = unname(unlist(cv_compare[5,1])),
                      rmse = unname(unlist(cv_compare[8,1])),
                      cor_obs_pred = unname(unlist(cv_compare[6,1])),
                      cor_pred_res = unname(unlist(cv_compare[7,1])))
  temp_diagnostics <- rbind(temp_diagnostics, temp2)
  
  rm(temp_week, vgm, vgm_fit, model, ok_result, cv_result, cv_compare,
     temp, temp2, week_end)
}

temp_ct_data <- temp_ct_data %>%
  rename(biweekly_average_temp_pred = biweekly_average_pred,
         biweekly_average_temp_var = biweekly_average_var) %>%
  mutate_if(is.factor, as.character)

temp_ct_data_full <- full_join(ap_dates, temp_ct_data, by="week_ending")

save(temp_ct_data, temp_ct_data_full, temp_cv_results, temp_diagnostics,
     file = here::here("data/Air Quality", "temp kriging results.RData"))

ct_temperature <- temp_ct_data_full
save(ct_temperature, file = here::here("Data", "CT_Temperature.RData"))

load("./Data/Air Quality/temp kriging results.RData")
write_csv(temp_ct_data_full, here::here("data/Air Quality", "Census_Tract_Temp.csv"))

