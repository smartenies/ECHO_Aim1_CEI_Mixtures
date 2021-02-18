#' =============================================================================
#' Project: ECHO Aim 1 
#' Date created: June 4, 2018
#' Date updated: June 12, 2018
#' Date updated: June 25, 2018
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' Description:
#'  
#' This project examines the relationships between spatially-distributed
#' economic, environmental, and social variables and health outcomes meausred
#' in the Healthy Start cohort (UC Denver)
#' 
#' This script includes the data cleaning and summary statistics for exposures,
#' outcomes, and covariates
#' 
#' Looks at correlations between potential covariates and exposures/outcomes
#' 
#' June 25, 2018: read in a new dataset with corrected income variables
#' Need to get in touch with Brandy about how they've addressed missing
#' income variables 
#' 
#' NOTE: don't forget the ./ before the directory when reading in files!
#' NOTE: Must be connected to the LEAD server (P:/ on my computer)
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
library(IC2)
library(haven)

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

geo_data <- "T:/Rsch-MRS/ECHO/SEM Large Data/Spatial Data/"
utm_13 <- "+init=epsg:26913"
albers <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
ll_nad83 <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
ll_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#' -----------------------------------------------------------------------------
#' Data cleaning and summary statistics
#' Table 1
#' -----------------------------------------------------------------------------

#' Read in HS I SAS dataset
#' Note! This takes a while
#' n = 1410 records

#' p_path <- "P:/LEAD/Outside-Users/Martenies/DataPull180318/"
#' hs_data <- read_sas(data_file = paste(p_path, "p0137.sas7bdat", sep=""))
#' length(unique(hs_data$pid))
#' nrow(hs_data)
#' 
#' #' rename old gestational age variable (using one from new data set!)
#' hs_data <- rename(hs_data, gest_age_d_old = Gest_age_birth_days,
#'                            gest_age_w_old = Gest_age_birth_weeks)
#' 
#' #' Add the secondhand smoke dataset from Jenny Aalborg
#' shs_data <- read_sas(data_file = paste(p_path, "p001.sas7bdat", sep=""))
#' length(unique(shs_data$pid))
#' nrow(shs_data)
#' glimpse(shs_data)
#' 
#' hs_data <- left_join(hs_data, shs_data, by="pid")
#' nrow(hs_data)
#' 
#' #' Add the anthro1 dataset from Robyn Harte
#' anthro_data <- read_sas(data_file = paste(p_path, "anthro1.sas7bdat", sep=""))
#' length(unique(anthro_data$pid))
#' nrow(anthro_data)
#' glimpse(anthro_data)
#' 
#' unique(anthro_data$source)
#' 
#' ipv3_anthro_data <- filter(anthro_data, source == "IPV3Visit") %>% 
#'   filter(excflg == 0)
#' length(unique(ipv3_anthro_data$pid))
#' nrow(ipv3_anthro_data)
#' summary(ipv3_anthro_data)
#' 
#' colnames(ipv3_anthro_data) <- gsub("_", "", colnames(ipv3_anthro_data))
#' hs_data <- left_join(hs_data, ipv3_anthro_data, by="pid")
#' 
#' nrow(hs_data)
#' 
#' #' Add the gestational age dataset
#' gest_data <- read_sas(data_file = paste(p_path, "gestationalage.sas7bdat",
#'                                           sep="")) %>%
#'   rename(gest_age_d = GA_at_birth_in_days) %>%
#'   mutate(gest_age_w = gest_age_d / 7)
#' summary(gest_data)
#' 
#' hs_data <- left_join(hs_data, gest_data, by="pid")
#' nrow(hs_data)
#' 
#' save(hs_data, file="./Data/CEI Data/hs_data.RData")
#' names(hs_data)

#' Read in the saved HS I dataset
load(file="./Data/CEI Data/hs_data.RData")

p1 <- filter(hs_data, pid %in% c(10704, 10533)) %>% 
  select(pid, full_name)
View(p1)


#' Load CEI dataset
#' N = 1151 records (81.5% of total)
load(file="./Data/CEI Data/hs_cei.RData")
length(unique(hs_cei$pid))
length(unique(hs_data$pid))

#' Fraction included
(length(unique(hs_cei$pid)) / length(unique(hs_data$pid))) * 100

#' Subset of participants with exposure data
pids <- unique(hs_cei$pid)

#' Indicator variable for participants who could be geocoded
#' and assigned exposure variables (1 = included)
#' N = 259 participants are excluded
hs_data$include <- ifelse(hs_data$pid %in% pids, 1, 0)
table(hs_data$include)

#' -----------------------------------------------------------------------------
#' Subsetting to relevant variables
#' Formatting some covariates
#' 
#' Data Dictionary of Derived Variables provided by Brandy Ringham 
#' Information on SHS exposure from Jenny Allborg (emailed on 6/25/18)
#' -----------------------------------------------------------------------------

hs_data2 <- hs_data %>% 
  select(pid, include, maternal_age, epi5, epi7, 
         pre_preg_bmi, pre_preg_BMI_cat,
         conception_date, mrdata17, pp1,
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

table(hs_data2$maternal_race_eth, hs_data2$infant_race_eth)

test <- select(hs_data2, mra_databwt, birth_weight, di_srbirthweight)

#' Check out variables 
#' N = 1410
summary(hs_data2)

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

  #' Conception date from SAS to R; SAS dates use Jan 1, 1960 as origin
  mutate(conception_date = format(as.Date(conception_date, 
                                          origin="1960-01-01"),"%Y-%m-%d")) %>%
      
  #' Date, month, and season of delivery
  mutate(del_date = as.Date(mrdata17, origin="1960-01-01")) %>%
  mutate(del_month = as.numeric(substr(as.character(del_date),
                                        start = 6, stop = 7))) %>%
  mutate(del_season = ifelse(del_month %in% c(12, 1, 2), 1, 
                             ifelse(del_month %in% c(3, 4, 5), 2,
                                    ifelse(del_month %in% c(6, 7, 8), 3, 4)))) %>%
  
  #' peapod date
  rename(pea_pod_date = pp1) %>% 
  mutate(days_to_peapod = difftime(pea_pod_date, del_date, units = "days")) %>% 
  mutate(days_to_peapod = as.numeric(days_to_peapod)) %>% 
  
  #' Indicator variables for season
  mutate(winter = ifelse(del_season == 1, 1, 0),
         spring = ifelse(del_season == 2, 1, 0),
         summer = ifelse(del_season == 3, 1, 0),
         fall = ifelse(del_season == 4, 1, 0)) %>%

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

#' Summary of infant sex and race/ethnicity
enrollment_table <- select(hs_data2, male, infant_race_eth) %>% 
  #filter(!is.na(male)) %>% 
  count(male, infant_race_eth) %>% 
  mutate(proportion = prop.table(n))
enrollment_table

#' Create "Table 1"
#' Summary of outcomes and covariates
library(qwraps2)
options(qwraps2_markup = "markdown")

tab1_summary <- 
  list("Maternal race/enthicity" =
         list("Hispanic/Latina" = ~ qwraps2::n_perc0(maternal_race_eth == 1, na_rm=T),
              "White non-Hispanic" = ~ qwraps2::n_perc0(maternal_race_eth == 2, na_rm=T),
              "Black" = ~ qwraps2::n_perc0(maternal_race_eth == 3, na_rm=T),
              "Other Race/Ethnicity" = ~ qwraps2::n_perc0(maternal_race_eth == 4, na_rm=T)),
       "Maternal age at delivery (years)" =
         list("Maternal age Mean (SD)" = ~ qwraps2::mean_sd(maternal_age, na_rm=T, 
                                                            denote_sd = "paren",
                                                            show_n = "always")),
       "Maternal perceived stress score" =
         list("CPSS Mean (SD)" = ~ qwraps2::mean_sd(mean_cpss, na_rm=T, 
                                                            denote_sd = "paren",
                                                            show_n = "always")),
       "Maternal depression score" =
         list("EPSD Mean (SD)" = ~ qwraps2::mean_sd(mean_epsd, na_rm=T, 
                                                            denote_sd = "paren",
                                                            show_n = "always")),
       "Maternal pre-pregnancy BMI category" =
         list("BMI Underweight" = ~ qwraps2::n_perc0(pre_preg_BMI_cat == 1, na_rm=T),
              "BMI Normal" = ~ qwraps2::n_perc0(pre_preg_BMI_cat == 2, na_rm=T),
              "BMI Overweight" = ~ qwraps2::n_perc0(pre_preg_BMI_cat == 3, na_rm=T),
              "BMI Obese" = ~ qwraps2::n_perc0(pre_preg_BMI_cat == 4, na_rm=T)),
       "Maternal education level" =
         list("Less than high school" = ~ qwraps2::n_perc0(maternal_ed == 1, na_rm=T),
              "High school or GED" = ~ qwraps2::n_perc0(maternal_ed == 2, na_rm=T),
              "Some college/Associate's" = ~ qwraps2::n_perc0(maternal_ed == 3, na_rm=T),
              "Bachelor's Degree" = ~ qwraps2::n_perc0(maternal_ed == 4, na_rm=T),
              "Graduate Degree" = ~ qwraps2::n_perc0(maternal_ed == 5, na_rm=T)),
       "Household income in the past year" =
         list("Income Less than $40,000" = ~ qwraps2::n_perc0(HH_income == 1, na_rm=T),
              "Income $40,000 - $70,000" = ~ qwraps2::n_perc0(HH_income == 2, na_rm=T),
              "Income Greater than $70,000" = ~ qwraps2::n_perc0(HH_income == 3, na_rm=T),
              "Income Missing or Don't Know" = ~ qwraps2::n_perc0(HH_income == 4, na_rm=T)),
       "Maternal smoking during pregnancy" =
         list("Smoking Yes" = ~ qwraps2::n_perc0(any_smoker == 1, na_rm=T),
              "Smoking No" = ~ qwraps2::n_perc0(any_smoker == 0, na_rm=T)),
       "Secondhand smoking exposure during pregnancy" =
         list("SHS Exposure Yes" = ~ qwraps2::n_perc0(smokeSH == 1, na_rm=T),
              "SHS Exposure No" = ~ qwraps2::n_perc0(smokeSH == 0, na_rm=T)),
       "Infant sex" =
         list("Infant Sex: Male" = ~ qwraps2::n_perc0(infant_sex == 2, na_rm=T),
              "Infant Sex: Female" = ~ qwraps2::n_perc0(infant_sex == 1, na_rm=T)),
       "Birth season" =
         list("Birth Season: Winter" = ~ qwraps2::n_perc0(del_season == 1, na_rm=T),
              "Birth Season: Spring" = ~ qwraps2::n_perc0(del_season == 2, na_rm=T),
              "Birth Season: Summer" = ~ qwraps2::n_perc0(del_season == 3, na_rm=T),
              "Birth Season: Fall" = ~ qwraps2::n_perc0(del_season == 4, na_rm=T)),
       "Gestational age at birth (weeks)" =
         list("Gest age Mean (SD)" = ~ qwraps2::mean_sd(gest_age_w, na_rm=T, 
                                                        denote_sd = "paren",
                                                        show_n = "always")),
       "Term births" =
         list("Preterm" = ~ qwraps2::n_perc0(term_cat == 1, na_rm=T),
              "Early Term" = ~ qwraps2::n_perc0(term_cat == 2, na_rm=T),
              "Full Term" = ~ qwraps2::n_perc0(term_cat == 3, na_rm=T),
              "Late Term" = ~ qwraps2::n_perc0(term_cat == 4, na_rm=T),
              "Postterm" = ~ qwraps2::n_perc0(term_cat == 5, na_rm=T)),
       "Birth weight (g)" =
         list("Birth weight Mean (SD)" = ~ qwraps2::mean_sd(birth_weight, na_rm=T, 
                                                            denote_sd = "paren",
                                                            show_n = "always")),
       "Low birth weight" =
         list("LBW" = ~ qwraps2::n_perc0(low_birth_weight == 1, na_rm=T),
              "Not LBW" = ~ qwraps2::n_perc0(low_birth_weight == 0, na_rm=T)),
       "Days between delivery and PeaPod (n)" =
         list("Days to PeaPod mean (SD)" = ~ qwraps2::mean_sd(days_to_peapod, na_rm=T, 
                                                            denote_sd = "paren",
                                                            show_n = "always")),
       "Infant body mass (g)" =
         list("Body mass Mean (SD)" = ~ qwraps2::mean_sd(bm_g, na_rm=T, 
                                                         denote_sd = "paren",
                                                         show_n = "always")),
       "Infant fat mass (g)" =
         list("Fat mass Mean (SD)" = ~ qwraps2::mean_sd(fm_g, na_rm=T, 
                                                        denote_sd = "paren",
                                                        show_n = "always")),
       "Infant fat free mass (g)" =
         list("Far free mass Mean (SD)" = ~ qwraps2::mean_sd(ffm_g, na_rm=T, 
                                                             denote_sd = "paren",
                                                             show_n = "always")),
       "Infant adiposity (%)" =
         list("Adiposity Mean (SD)" = ~ qwraps2::mean_sd(adiposity, na_rm=T, 
                                                         denote_sd = "paren",
                                                         show_n = "always")),
       "Infant WFL z-score (%)" =
         list("WFL z-score Mean (SD)" = ~ qwraps2::mean_sd(wfl_sex_z_score, na_rm=T, 
                                                           denote_sd = "paren",
                                                           show_n = "always"))
       )

table1a <- summary_table(hs_data2, tab1_summary)
table1a
table1a_df <- as.data.frame(table1a) %>%
  mutate(variable = rownames(.)) %>%
  select(2,1)

table1b <- summary_table(dplyr::group_by(hs_data2, include), tab1_summary)
table1b
table1b_df <- as.data.frame(table1b) %>%
  mutate(variable = rownames(.)) %>%
  select(3,1:2)

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

#' Create "Table 2"
#' Summary of exposures

hs_data3 <- left_join(hs_data2, hs_cei, by="pid") %>%
  filter(include == 1)

exp_vars <- c("mean_pm", "mean_o3", "tri_tpy", "pct_tree_cover", "pct_impervious", 
              "sum_aadt_intensity", "npl_count", "waste_site_count", "major_emit_count", 
              "cafo_count", "mine_well_count", 
              "cvd_rate_adj", "res_rate_adj", 
              "violent_crime_rate", "property_crime_rate",
              "pct_poc", "pct_less_hs", "pct_unemp", "pct_limited_eng", 
              "pct_hh_pov", "med_income",
              "env", "soc", "cei", "soc_2", "cei_2")

table2 <- data.frame() 
for (i in 1:length(exp_vars)) {
  var <- as.data.frame(select(hs_data3, exp_vars[i]))
  temp <- data.frame(variable = exp_vars[i],
                     mean_sd = paste(round(mean(var[,1], na.rm=T), 2),
                                     " (",
                                     round(sd(var[,1], na.rm=T), 2),
                                     ")", sep=""),
                     min = round(min(var[,1], na.rm=T), 2),
                     q25 = round(quantile(var[,1], probs = 0.25, na.rm=T), 2),
                     median = round(median(var[,1], na.rm=T), 2),
                     q75 = round(quantile(var[,1], probs = 0.75, na.rm=T), 2),
                     q95 = round(quantile(var[,1], probs = 0.95, na.rm=T), 2),
                     max = round(max(var[,1], na.rm=T), 2),
                     cv = round(sd(var[,1], na.rm=T) / mean(var[,1], na.rm=T), 2)*100)
  table2 <- rbind(table2, temp)
  rm(temp)
}


summary(hs_data3$env)
summary(hs_data3$soc)
summary(hs_data3$cei)


#' Save cleaned dataset
save(hs_data, hs_data2, hs_data3, file="./Data/CEI Data/cleaned data.RData")

names(hs_data3)

#' Histograms of continuous exposures, outcomes, and covariates
hist(hs_data3$maternal_age)
hist(hs_data3$ipv1_cpsssum)
hist(hs_data3$ipv2_cpsssum)
hist(hs_data3$pni_cpsssum)
hist(hs_data3$mean_cpss)
hist(hs_data3$ipv1_epsdsum)
hist(hs_data3$ipv2_epsdsum)
hist(hs_data3$pni_epsdsum)
hist(hs_data3$mean_epsd)
hist(hs_data3$birth_weight)
hist(hs_data3$gest_age_w)
hist(hs_data3$pre_preg_bmi)
hist(hs_data3$wfl_sex_z_score)
hist(hs_data3$days_to_peapod)
hist(hs_data3$bm_g)
hist(hs_data3$fm_g)
hist(hs_data3$adiposity)
hist(hs_data3$mean_pm)
hist(hs_data3$mean_o3)
hist(hs_data3$env)
hist(hs_data3$soc)
hist(hs_data3$soc_2)
hist(hs_data3$cei)
hist(hs_data3$cei_2)

#' Correlations between variables
library(corrplot)
hs_data_num <- select(hs_data3, maternal_age, maternal_ed, any_smoker, smokeSH,
                      del_season, days_to_peapod,  term_cat,
                      ipv1_cpsssum, ipv2_cpsssum, pni_cpsssum, mean_cpss, 
                      ipv1_epsdsum, ipv2_epsdsum, pni_epsdsum, mean_epsd,
                      birth_weight, gest_age_w, pre_preg_bmi, 
                      wfl_sex_z_score, birth_weight, bm_g, fm_g, adiposity, 
                      mean_pm, mean_o3, env, soc, soc_2, cei, cei_2)

hs_cor <- cor(hs_data_num, use = "complete")
hs_cor
write_xlsx(as.data.frame(hs_cor),
           path="./Results/CEI Results./Variable Correlations.xlsx")

jpeg("./Figures/CEI Figures/Indices/correlation matrix.jpeg")
corrplot(hs_cor, type="upper")
dev.off()
