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
#' Description: Data cleaning 
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
#' Data cleaning and summary statistics
#' -----------------------------------------------------------------------------

#' Read in HS I data set
load(here::here("Data", "hs_data.RData"))

#' -----------------------------------------------------------------------------
#' Subsetting to relevant variables
#' Formatting some covariates
#' 
#' Data Dictionary of Derived Variables provided by Brandy Ringham 
#' Information on SHS exposure from Jenny Allborg (emailed on 6/25/18)
#' -----------------------------------------------------------------------------

hs_data2 <- hs_data %>% 
  select(pid, maternal_age, epi5, epi7, 
         pre_preg_bmi, pre_preg_BMI_cat,
         mrdata17, pp1,
         mra_databwt, birth_weight, di_srbirthweight,
         infant_sex, gravidity, HH_income, 
         ipv1smokeyn, ipv2smokeyn, ipv3smokeyn,
         gest_age_d_old, gest_age_d, gest_age_w_old, gest_age_w, 
         bm, ffm, fm, fmp,
         ipv3_pp_bm, ipv3_pp_bv, ipv3_pp_length, 
         ipv3_pp_fm, ipv3_pp_ffm, 
         pct_fm, ipv3_pp_fm_pct, ipv3_pp_ffm_pct,
         ZWFL, excflg,
         ipv1_epsdsum, ipv2_epsdsum, pni_epsdsum,
         ipv1_cpsssum, ipv2_cpsssum, pni_cpsssum,
         ipv1_cohenperstress_complete, ipv2_cohenperstress_complete,
         pni_cohenperstress_complete,
         ipv1_depstressstaff, ipv2_depstressstaff, pni_depstressstaff,
         smokeSH_IPV1, smokeSH_IPV2, smokeSH_IPV3, smokeSH) %>%
  rename(maternal_ed = epi5, 
         maternal_race_eth = epi7,
         wfl_sex_z_score = ZWFL) %>% 
  mutate(infant_race_eth = maternal_race_eth)

#' -----------------------------------------------------------------------------
#' Add in the exposure variables
#' Indicator variable if the participant can be included
#' -----------------------------------------------------------------------------

load(here::here("Data", "hs_exposures.RData"))
hs_data2 <- left_join(hs_data2, hs, by = "pid")

exp_pids <- unique(hs$pid)

hs_data2$include <- ifelse(hs_data2$pid %in% exp_pids, 1, 0)
table(hs_data2$include)

#' -----------------------------------------------------------------------------
#' Additional data cleaning
#' -----------------------------------------------------------------------------

#' Check out variables 
#' N = 1410
names(hs_data2)
class(hs_data2$conception_date)

table(hs_data2$maternal_race_eth, hs_data2$infant_race_eth)

#' Clean up data and format some additional covariates
hs_data2 <- hs_data2 %>%
  
  #' set negative gestational ages to NA
  mutate(gest_age_d_old = ifelse(gest_age_d_old < 0, NA, gest_age_d_old),
         gest_age_w_old = ifelse(gest_age_w_old < 0, NA, gest_age_w_old)) %>%
  
  #' Term birth defined by ACOG
  #' https://www.acog.org/Clinical-Guidance-and-Publications/Committee-Opinions/Committee-on-Obstetric-Practice/Definition-of-Term-Pregnancy?IsMobileSet=false
  mutate(term_cat = ifelse(gest_age_w < 37, 1,
                           ifelse(gest_age_w >= 37 & gest_age_w < 39, 2,
                                  ifelse(gest_age_w >= 39 & gest_age_w < 41, 3,
                                         ifelse(gest_age_w >= 41 & gest_age_w < 42, 4,
                                                ifelse(gest_age_w >= 42, 5, NA)))))) %>%
  
  mutate(preterm_birth = ifelse(term_cat == 1, 1, 0),
         earlyterm_birth = ifelse(term_cat == 2, 1, 0),
         fullterm_birth = ifelse(term_cat == 3, 1, 0),
         lateterm_birth = ifelse(term_cat == 4, 1, 0),
         postterm_birth = ifelse(term_cat == 5, 1, 0)) %>%
  
  #' Low birth weight (< 2500 g)
  mutate(low_birth_weight = ifelse(birth_weight < 2500, 1, 0)) %>% 
  
  #' participant reported any smoking during their pregnancy
  mutate(any_smoker = ifelse(ipv1smokeyn == 1 | 
                             ipv2smokeyn == 1 | 
                             ipv3smokeyn == 1, 1, 0)) %>%

  #' Date, month, and season of delivery
  mutate(del_date = as.Date(mrdata17, origin="1960-01-01")) %>%
  mutate(del_month = as.numeric(substr(as.character(del_date),
                                        start = 6, stop = 7))) %>%
  mutate(del_season = ifelse(del_month %in% c(12, 1, 2), 1, 
                             ifelse(del_month %in% c(3, 4, 5), 2,
                                    ifelse(del_month %in% c(6, 7, 8), 3, 4)))) %>%
  
  #' Date, month, and season of conception
  mutate(concep_month = month(conception_date)) %>%
  mutate(concep_season = ifelse(concep_month %in% c(12, 1, 2), 1, 
                             ifelse(concep_month %in% c(3, 4, 5), 2,
                                    ifelse(concep_month %in% c(6, 7, 8), 3, 4)))) %>%
  
  #' peapod date
  rename(pea_pod_date = pp1) %>% 
  mutate(days_to_peapod = difftime(pea_pod_date, del_date, units = "days")) %>% 
  mutate(days_to_peapod = as.numeric(days_to_peapod)) %>% 
  
  #' Indicator variables for season
  mutate(del_winter = ifelse(del_season == 1, 1, 0),
         del_spring = ifelse(del_season == 2, 1, 0),
         del_summer = ifelse(del_season == 3, 1, 0),
         del_fall = ifelse(del_season == 4, 1, 0)) %>%
  
  mutate(concep_winter = ifelse(concep_season == 1, 1, 0),
         concep_spring = ifelse(concep_season == 2, 1, 0),
         concep_summer = ifelse(concep_season == 3, 1, 0),
         concep_fall = ifelse(concep_season == 4, 1, 0)) %>%

  #' Indicator variables for race/ethnicity, income, education, 
  #' pre-pregnancy BMI, and infant sex
  mutate(latina_re = ifelse(maternal_race_eth == 1, 1, 0),
         white_re = ifelse(maternal_race_eth == 2, 1, 0),
         black_re = ifelse(maternal_race_eth == 3, 1, 0),
         other_re = ifelse(maternal_race_eth == 4, 1, 0),
         low_income = ifelse(HH_income == 1, 1, 0),
         mid_income = ifelse(HH_income == 2, 1, 0),
         high_income = ifelse(HH_income == 3, 1, 0),
         other_income = ifelse(HH_income == 4, 1, 0),
         ed_no_hs = ifelse(maternal_ed == 1, 1, 0),
         ed_hs = ifelse(maternal_ed == 2, 1, 0),
         ed_aa = ifelse(maternal_ed == 3, 1, 0),
         ed_4yr = ifelse(maternal_ed == 4, 1, 0),
         ed_grad = ifelse(maternal_ed == 5, 1, 0),
         low_bmi = ifelse(pre_preg_BMI_cat == 1, 1, 0),
         norm_bmi = ifelse(pre_preg_BMI_cat == 2, 1, 0),
         ovwt_bmi = ifelse(pre_preg_BMI_cat == 3, 1, 0),
         obese_bmi = ifelse(pre_preg_BMI_cat == 4, 1, 0),
         male = ifelse(infant_sex == 2, 1, 0)) %>%
         
  #' convert body mass and fat mass from kg to g
  mutate(body_mass_g = (ipv3_pp_bm * 1000),
         fat_mass_g = ipv3_pp_fm * 1000,
         fat_free_mass_g = ipv3_pp_ffm,
         bm_g = bm * 1000,
         fm_g = fm * 1000,
         ffm_g = ffm * 1000) %>%
  
  #' adiposity = infant fat mass / infant body mass *100
  mutate(adiposity_old = (fat_mass_g / body_mass_g) * 100,
         adiposity = (fm_g / bm_g) * 100)
  
#' Check out new variables
#' Lots of missing data, but can work with that later
summary(hs_data2)
glimpse(hs_data2)

table(hs_data2$concep_month)
table(hs_data2$del_month)

#' Appropriate to average the depression and stress scores across pregnancy?
#' Distributions of Cohen Perceived Stress Scale scores are similar
ggplot(hs_data2) +
  geom_density(aes(x = ipv1_cpsssum), col="red") +
  geom_density(aes(x = ipv2_cpsssum), col="blue") +
  geom_density(aes(x = pni_cpsssum), col="green") +
  simple_theme

#' Depresson scores skew lower for the last prenatal visit
#' Should be OK to average, though, for now
ggplot(hs_data2) +
  geom_density(aes(x = ipv1_epsdsum), col="red") +
  geom_density(aes(x = ipv2_epsdsum), col="blue") +
  geom_density(aes(x = pni_epsdsum), col="green") +
  simple_theme

hs_data2 <- hs_data2 %>% 
  mutate(mean_cpss = rowMeans(subset(hs_data2,
                                     select = c(ipv1_cpsssum, ipv2_cpsssum, pni_cpsssum)), 
                              na.rm=T),
         mean_epsd = rowMeans(subset(hs_data2,
                                     select = c(ipv1_epsdsum, ipv2_epsdsum, pni_epsdsum)),
                              na.rm=T))

check <- select(hs_data2, pid, ipv1_cpsssum, ipv2_cpsssum, pni_cpsssum, mean_cpss)

hs_data2$mean_cpss <- ifelse(is.nan(hs_data2$mean_cpss), NA, hs_data2$mean_cpss)
hs_data2$mean_epsd <- ifelse(is.nan(hs_data2$mean_epsd), NA, hs_data2$mean_epsd)

#' -----------------------------------------------------------------------------
#' Add the lat-lon variables
#' -----------------------------------------------------------------------------

load(here::here("Data", "hs_spatial.RData"))
head(hs)

hs_ll <- st_transform(hs, crs = ll_wgs84)

hs_coords <- do.call(rbind, st_geometry(hs_ll)) %>% as_tibble()
names(hs_coords) <- c("lon", "lat")

hs_ll <- bind_cols(hs_ll, hs_coords) %>% 
  st_set_geometry(NULL) %>% 
  select(pid, lon, lat)
head(hs_ll)

hs_data2 <- left_join(hs_data2, hs_ll, by = "pid")

#' -----------------------------------------------------------------------------
#' Save clean data
#' -----------------------------------------------------------------------------

write_csv(hs_data2, here::here("Data", "HS_Clean_Data.csv"))

#' -----------------------------------------------------------------------------
#' Compare included and excluded participants
#' -----------------------------------------------------------------------------

#' Do geocoded HS participants with exposure data differ from those excluded?
#' Chi square test for categorical variables, Student's t for continuous
hs_inc <- filter(hs_data2, include == 1)
hs_exc <- filter(hs_data2, include == 0)

#' Maternal race/ethnicity
table(hs_data2$maternal_race_eth, hs_data2$include)
chisq.test(hs_data2$maternal_race_eth, hs_data2$include)
pairwise.prop.test(table(hs_data2$maternal_race_eth, hs_data2$include))

#' Maternal age at delivery
t.test(hs_inc$maternal_age, hs_exc$maternal_age)

#' Maternal perceived stress (mean across pregnancy)
t.test(hs_inc$mean_cpss, hs_exc$mean_cpss)

#' Maternal depression score (mean across pregnancy)
t.test(hs_inc$mean_epsd, hs_exc$mean_epsd)

#' Maternal pre-pregnancy BMI
table(hs_data2$pre_preg_BMI_cat, hs_data2$include)
chisq.test(hs_data2$pre_preg_BMI_cat, hs_data2$include)

#' Maternal education
table(hs_data2$maternal_ed, hs_data2$include)
chisq.test(hs_data2$maternal_ed, hs_data2$include)
pairwise.prop.test(table(hs_data2$maternal_ed, hs_data2$include))

#' Past year household income
table(hs_data2$HH_income, hs_data2$include)
chisq.test(hs_data2$HH_income, hs_data2$include)
pairwise.prop.test(table(hs_data2$HH_income, hs_data2$include))

#' Smoking during preganacy
table(hs_data2$any_smoker, hs_data2$include)
chisq.test(hs_data2$any_smoker, hs_data2$include)

#' SHS exposure during preganacy
table(hs_data2$smokeSH, hs_data2$include)
chisq.test(hs_data2$smokeSH, hs_data2$include)

#' Infant sex
table(hs_data2$infant_sex, hs_data2$include)
chisq.test(hs_data2$infant_sex, hs_data2$include)

#' Season of birth
table(hs_data2$del_season, hs_data2$include)
chisq.test(hs_data2$del_season, hs_data2$include)
pairwise.prop.test(table(hs_data2$del_season, hs_data2$include))

#' Gestational age (weeks)
t.test(hs_inc$gest_age_w, hs_exc$gest_age_w)

#' Term births
table(hs_data2$term_cat, hs_data2$include)
chisq.test(hs_data2$term_cat, hs_data2$include)
pairwise.prop.test(table(hs_data2$del_season, hs_data2$include))

#' birth_weight
t.test(hs_inc$birth_weight, hs_exc$birth_weight)

#' Low birth weight
table(hs_data2$low_birth_weight, hs_data2$include)
chisq.test(hs_data2$low_birth_weight, hs_data2$include)

#' Body mass
t.test(hs_inc$bm_g, hs_exc$bm_g)

#' Fat mass
t.test(hs_inc$fm_g, hs_exc$fm_g)

#' Fat free mass
t.test(hs_inc$ffm_g, hs_exc$ffm_g)

#' Adiposity
t.test(hs_inc$adiposity, hs_exc$adiposity)

#' weight-for-length
t.test(hs_inc$wfl_sex_z_score, hs_exc$wfl_sex_z_score)

