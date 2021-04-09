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
ggsave(filename = here::here("Manuscripts", "Figure_S1.jpeg"), device = "jpeg",
       width = 9, height = 7, units = "in", dpi = 500)

#' -----------------------------------------------------------------------------
#' GAM plots
#' -----------------------------------------------------------------------------

library(mgcv)
library(mgcViz)
library(directlabels)

#'------------
#' Figure 1: Plot showing results from full model w/ axes restricted to middle 95% 
#'------------
m_ozone <- 48.0
sd_ozone <- 3.1

m_temp <- 52.6
sd_temp <- 4.5

load(file = here::here("Results", "BW_GAM_v4.rdata"))
bw_gam

#' In mgcv
plot(bw_gam)

gam_b <- getViz(bw_gam)
b <- sm(gam_b, 1)

#' 2D plot
gam_b_2d <- plot(b) + 
  l_fitRaster() + l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") + 
  l_points() + l_rug() + 
  labs(title = "A) 2D Smoothed Term", x = "Ozone (ppb)", y = "Temperature (F)") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_ozone - 2*sd_ozone, m_ozone - sd_ozone, m_ozone, 
                                m_ozone + sd_ozone, m_ozone + 2*sd_ozone),
                     limits = c(-2, 2)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_temp - 2*sd_temp, m_temp - sd_temp, m_temp, 
                                m_temp + sd_temp, m_temp + 2*sd_temp),
                     limits = c(-2, 2)) +
  guides(fill=guide_legend(title="Change in\nbirth weight (g)"))
gam_b_2d <- direct.label(gam_b_2d$ggObj, list("bottom.pieces", colour='black'))
gam_b_2d

#' Accumulated Local Effects
temp_ale <- plot(ALE(bw_gam, "mean_temp")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_temp - 2*sd_temp, m_temp - 1*sd_temp, m_temp, 
                                m_temp + 1*sd_temp, m_temp + 2*sd_temp),
                     limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-500, 300)) +
  labs(title = "B) ALE for temperature", x = "Temperature (F)", y = "Change in birth weight (g)")
temp_ale

o3_ale <- plot(ALE(bw_gam, "mean_o3")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_ozone - 2*sd_ozone, m_ozone - 1*sd_ozone, m_ozone, 
                                m_ozone + 1*sd_ozone, m_ozone + 2*sd_ozone),
                     limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-500, 300)) +
  labs(title = "C) ALE for ozone", x = "Ozone (ppb)", y = "Change in birth weight (g)")
o3_ale

jpeg(file = here::here("Manuscripts", "Figure_1.jpeg"),
     width = 12, height = 4, units = "in", res = 500)
gridPrint(gam_b_2d, temp_ale, o3_ale, ncol = 3)
dev.off()

#'------------
#' Figure S2: Plot showing results from full model 
#'------------
m_ozone <- 48.0
sd_ozone <- 3.1

m_temp <- 52.6
sd_temp <- 4.5

load(file = here::here("Results", "BW_GAM_v4.rdata"))
bw_gam

#' In mgcv
plot(bw_gam)

gam_b <- getViz(bw_gam)
b <- sm(gam_b, 1)

#' 2D plot
gam_b_2d <- plot(b) + 
  l_fitRaster() + l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") + 
  l_points() + l_rug() + 
  labs(title = "A) 2D Smoothed Term", x = "Ozone (ppb)", y = "Temperature (F)") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_ozone - 2*sd_ozone, m_ozone - sd_ozone, m_ozone, 
                                m_ozone + sd_ozone, m_ozone + 2*sd_ozone)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_temp - 2*sd_temp, m_temp - sd_temp, m_temp, 
                                m_temp + sd_temp, m_temp + 2*sd_temp)) +
  guides(fill=guide_legend(title="Change in\nbirth weight (g)"))
gam_b_2d <- direct.label(gam_b_2d$ggObj, list("bottom.pieces", colour='black', cex = 0.75))
gam_b_2d

#' Accumulated Local Effects
temp_ale <- plot(ALE(bw_gam, "mean_temp")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_temp - 2*sd_temp, m_temp - sd_temp, m_temp, 
                                m_temp + sd_temp, m_temp + 2*sd_temp)) +
  labs(title = "B) ALE for temperature", x = "Temperature (F)", y = "Change in birth weight (g)")
temp_ale

o3_ale <- plot(ALE(bw_gam, "mean_o3")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_ozone - 2*sd_ozone, m_ozone - sd_ozone, m_ozone, 
                                m_ozone + sd_ozone, m_ozone + 2*sd_ozone)) +
  labs(title = "C) ALE for ozone", x = "Ozone (ppb)", y = "Change in birth weight (g)")
o3_ale

jpeg(file = here::here("Manuscripts", "Figure_S2.jpeg"),
     width = 12, height = 4, units = "in", res = 500)
gridPrint(gam_b_2d, temp_ale, o3_ale, ncol = 3)
dev.off()

#'------------
#' Figure S3: Plot showing results from model where data were restricted
#'------------

m_ozone <- 48.0
sd_ozone <- 3.1

m_temp <- 52.6
sd_temp <- 4.5

load(file = here::here("Results", "BW_GAM_Sensitivity_v4.rdata"))
bw_gam2

#' In mgcv
plot(bw_gam2)

gam_b2 <- getViz(bw_gam2)
b2 <- sm(gam_b2, 1)

#' 2D plot
gam_b2_2d <- plot(b2) + 
  l_fitRaster() + l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") + 
  l_points() + l_rug() + 
  labs(title = "A) 2D Smoothed Term", x = "Ozone (ppb)", y = "Temperature (F)") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_ozone - 2*sd_ozone, m_ozone - sd_ozone, m_ozone, 
                                m_ozone + sd_ozone, m_ozone + 2*sd_ozone)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_temp - 2*sd_temp, m_temp - sd_temp, m_temp, 
                                m_temp + sd_temp, m_temp + 2*sd_temp)) +
  guides(fill=guide_legend(title="Change in\nbirth weight (g)"))
gam_b2_2d <- direct.label(gam_b2_2d$ggObj, list("bottom.pieces", colour='black', cex = 0.75))
gam_b2_2d

#' Accumulated Local Effects
temp_ale2 <- plot(ALE(bw_gam2, "mean_temp")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_temp - 2*sd_temp, m_temp - sd_temp, m_temp, 
                                m_temp + sd_temp, m_temp + 2*sd_temp)) +
  labs(title = "B) ALE for temperature", x = "Temperature (F)", y = "Change in birth weight (g)")
temp_ale2

o3_ale2 <- plot(ALE(bw_gam2, "mean_o3")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(m_ozone - 2*sd_ozone, m_ozone - sd_ozone, m_ozone, 
                                m_ozone + sd_ozone, m_ozone + 2*sd_ozone)) +
  labs(title = "C) ALE for ozone", x = "Ozone (ppb)", y = "Change in birth weight (g)")
o3_ale2

jpeg(file = here::here("Manuscripts", "Figure_S3.jpeg"),
     width = 12, height = 4, units = "in", res = 500)
gridPrint(gam_b2_2d, temp_ale2, o3_ale2, ncol = 3)
dev.off()

#'------------
#' Figure S4: Plot showing results from model where data were stratified by race/ethnicity
#' Note that the plots are showing results for all of the observations
#'------------

nhw_m_ozone <- 47.8
nhw_sd_ozone <- 3.0

nhw_m_temp <- 52.5
nhw_sd_temp <- 4.4

other_m_ozone <- 48.1
other_sd_ozone <- 3.1

other_m_temp <- 52.7
other_sd_temp <- 4.7

load(file = here::here("Results", "BW_GAM_v4b.rdata"))
bw_gam_nhw <- bw_gam
summary(bw_gam_nhw)

load(file = here::here("Results", "BW_GAM_v4c.rdata"))
bw_gam_other <- bw_gam
summary(bw_gam_other)

#' NHW
gam_b_nhw <- getViz(bw_gam_nhw)
b_nhw <- sm(gam_b_nhw, 1)

#' 2D plot
gam_b_2d_nhw <- plot(b_nhw) + 
  l_fitRaster() + l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") + 
  l_points() + l_rug() + 
  labs(title = "A) 2D Smoothed Term (NHW Mothers)", x = "Ozone (ppb)", y = "Temperature (F)") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(nhw_m_ozone - 2*nhw_sd_ozone, nhw_m_ozone - nhw_sd_ozone, nhw_m_ozone, 
                                nhw_m_ozone + nhw_sd_ozone, nhw_m_ozone + 2*nhw_sd_ozone)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(nhw_m_temp - 2*nhw_sd_temp, nhw_m_temp - nhw_sd_temp, nhw_m_temp, 
                                nhw_m_temp + nhw_sd_temp, nhw_m_temp + 2*nhw_sd_temp)) +
  guides(fill=guide_legend(title="Change in\nbirth weight (g)"))
gam_b_2d_nhw <- direct.label(gam_b_2d_nhw$ggObj, list("bottom.pieces", colour='black', cex = 0.75))
gam_b_2d_nhw

#' Accumulated Local Effects
temp_ale_nhw <- plot(ALE(bw_gam_nhw, "mean_temp")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(nhw_m_temp - 2*nhw_sd_temp, nhw_m_temp - nhw_sd_temp, nhw_m_temp, 
                                nhw_m_temp + nhw_sd_temp, nhw_m_temp + 2*nhw_sd_temp)) +
  labs(title = "B) ALE for temperature (NHW Mothers)", x = "Temperature (F)", y = "Change in birth weight (g)")
temp_ale_nhw

o3_ale_nhw <- plot(ALE(bw_gam_nhw, "mean_o3")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(nhw_m_ozone - 2*nhw_sd_ozone, nhw_m_ozone - nhw_sd_ozone, nhw_m_ozone, 
                                nhw_m_ozone + nhw_sd_ozone, nhw_m_ozone + 2*nhw_sd_ozone)) +
  labs(title = "C) ALE for ozone (NHW Mothers)", x = "Ozone (ppb)", y = "Change in birth weight (g)")
o3_ale_nhw

#' All other race/ethnicity
gam_b_other <- getViz(bw_gam_other)
b_other <- sm(gam_b_other, 1)

#' 2D plot
gam_b_2d_other <- plot(b_other) + 
  l_fitRaster() + l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") + 
  l_points() + l_rug() + 
  labs(title = "D) 2D Smoothed Term (non-NHW Mothers)", x = "Ozone (ppb)", y = "Temperature (F)") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(other_m_ozone - 2*other_sd_ozone, other_m_ozone - other_sd_ozone, other_m_ozone, 
                                other_m_ozone + other_sd_ozone, other_m_ozone + 2*other_sd_ozone)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(other_m_temp - 2*other_sd_temp, other_m_temp - other_sd_temp, other_m_temp, 
                                other_m_temp + other_sd_temp, other_m_temp + 2*other_sd_temp)) +
  guides(fill=guide_legend(title="Change in\nbirth weight (g)"))
gam_b_2d_other <- direct.label(gam_b_2d_other$ggObj, list("bottom.pieces", colour='black', cex = 0.65))
gam_b_2d_other

#' Accumulated Local Effects
temp_ale_other <- plot(ALE(bw_gam_other, "mean_temp")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(other_m_temp - 2*other_sd_temp, other_m_temp - other_sd_temp, other_m_temp, 
                                other_m_temp + other_sd_temp, other_m_temp + 2*other_sd_temp)) +
  labs(title = "E) ALE for temperature (non-NHW Mothers)", x = "Temperature (F)", y = "Change in birth weight (g)")
temp_ale_other

o3_ale_other <- plot(ALE(bw_gam_other, "mean_o3")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(other_m_ozone - 2*other_sd_ozone, other_m_ozone - other_sd_ozone, other_m_ozone, 
                                other_m_ozone + other_sd_ozone, other_m_ozone + 2*other_sd_ozone)) +
  labs(title = "F) ALE for ozone (non-NHW Mothers)", x = "Ozone (ppb)", y = "Change in birth weight (g)")
o3_ale_other

jpeg(file = here::here("Manuscripts", "Figure_S4.jpeg"),
     width = 13, height = 9, units = "in", res = 500)
gridPrint(gam_b_2d_nhw, temp_ale_nhw, o3_ale_nhw,
          gam_b_2d_other, temp_ale_other, o3_ale_other, ncol = 3)
dev.off()

#'------------
#' Figure S5: Plot showing results from model where data were stratified by race/ethnicity
#' Note that the plots are showing results for models restricted to the middle 95% of observations
#'------------

nhw_m_ozone <- 47.8
nhw_sd_ozone <- 3.0

nhw_m_temp <- 52.5
nhw_sd_temp <- 4.4

other_m_ozone <- 48.1
other_sd_ozone <- 3.1

other_m_temp <- 52.7
other_sd_temp <- 4.7

load(file = here::here("Results", "BW_GAM_Sensitivity_v4b.rdata"))
bw_gam2_nhw <- bw_gam2
summary(bw_gam2_nhw)

load(file = here::here("Results", "BW_GAM_Sensitivity_v4c.rdata"))
bw_gam2_other <- bw_gam2
summary(bw_gam2_other)

#' NHW
gam_b2_nhw <- getViz(bw_gam2_nhw)
b2_nhw <- sm(gam_b2_nhw, 1)

#' 2D plot
gam_b2_2d_nhw <- plot(b2_nhw) + 
  l_fitRaster() + l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") + 
  l_points() + l_rug() + 
  labs(title = "A) 2D Smoothed Term (NHW Mothers)", x = "Ozone (ppb)", y = "Temperature (F)") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(nhw_m_ozone - 2*nhw_sd_ozone, nhw_m_ozone - nhw_sd_ozone, nhw_m_ozone, 
                                nhw_m_ozone + nhw_sd_ozone, nhw_m_ozone + 2*nhw_sd_ozone)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(nhw_m_temp - 2*nhw_sd_temp, nhw_m_temp - nhw_sd_temp, nhw_m_temp, 
                                nhw_m_temp + nhw_sd_temp, nhw_m_temp + 2*nhw_sd_temp)) +
  guides(fill=guide_legend(title="Change in\nbirth weight (g)"))
gam_b2_2d_nhw <- direct.label(gam_b2_2d_nhw$ggObj, list("bottom.pieces", colour='black', cex = 0.75))
gam_b2_2d_nhw

#' Accumulated Local Effects
temp_ale2_nhw <- plot(ALE(bw_gam2_nhw, "mean_temp")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(nhw_m_temp - 2*nhw_sd_temp, nhw_m_temp - nhw_sd_temp, nhw_m_temp, 
                                nhw_m_temp + nhw_sd_temp, nhw_m_temp + 2*nhw_sd_temp)) +
  labs(title = "B) ALE for temperature (NHW Mothers)", x = "Temperature (F)", y = "Change in birth weight (g)")
temp_ale2_nhw

o3_ale2_nhw <- plot(ALE(bw_gam2_nhw, "mean_o3")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(nhw_m_ozone - 2*nhw_sd_ozone, nhw_m_ozone - nhw_sd_ozone, nhw_m_ozone, 
                                nhw_m_ozone + nhw_sd_ozone, nhw_m_ozone + 2*nhw_sd_ozone)) +
  labs(title = "C) ALE for ozone (NHW Mothers)", x = "Ozone (ppb)", y = "Change in birth weight (g)")
o3_ale2_nhw

#' All other race/ethnicity
gam_b2_other <- getViz(bw_gam2_other)
b2_other <- sm(gam_b2_other, 1)

#' 2D plot
gam_b2_2d_other <- plot(b2_other) + 
  l_fitRaster() + l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") + 
  l_points() + l_rug() + 
  labs(title = "D) 2D Smoothed Term (non-NHW Mothers)", x = "Ozone (ppb)", y = "Temperature (F)") +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(other_m_ozone - 2*other_sd_ozone, other_m_ozone - other_sd_ozone, other_m_ozone, 
                                other_m_ozone + other_sd_ozone, other_m_ozone + 2*other_sd_ozone)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(other_m_temp - 2*other_sd_temp, other_m_temp - other_sd_temp, other_m_temp, 
                                other_m_temp + other_sd_temp, other_m_temp + 2*other_sd_temp)) +
  guides(fill=guide_legend(title="Change in\nbirth weight (g)"))
gam_b2_2d_other <- direct.label(gam_b2_2d_other$ggObj, list("bottom.pieces", colour='black', cex = 0.65))
gam_b2_2d_other

#' Accumulated Local Effects
temp_ale2_other <- plot(ALE(bw_gam2_other, "mean_temp")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(other_m_temp - 2*other_sd_temp, other_m_temp - other_sd_temp, other_m_temp, 
                                other_m_temp + other_sd_temp, other_m_temp + 2*other_sd_temp)) +
  labs(title = "E) ALE for temperature (non-NHW Mothers)", x = "Temperature (F)", y = "Change in birth weight (g)")
temp_ale2_other

o3_ale2_other <- plot(ALE(bw_gam2_other, "mean_o3")) + l_fitLine() + l_ciLine() + l_rug() +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2),
                     labels = c(other_m_ozone - 2*other_sd_ozone, other_m_ozone - other_sd_ozone, other_m_ozone, 
                                other_m_ozone + other_sd_ozone, other_m_ozone + 2*other_sd_ozone)) +
  labs(title = "F) ALE for ozone (non-NHW Mothers)", x = "Ozone (ppb)", y = "Change in birth weight (g)")
o3_ale2_other

jpeg(file = here::here("Manuscripts", "Figure_S5.jpeg"),
     width = 13, height = 9, units = "in", res = 500)
gridPrint(gam_b2_2d_nhw, temp_ale2_nhw, o3_ale2_nhw,
          gam_b2_2d_other, temp_ale2_other, o3_ale2_other, ncol = 3)
dev.off()