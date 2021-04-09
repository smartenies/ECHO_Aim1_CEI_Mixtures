#' =============================================================================
#' Project: ECHO Aim 1 CEI BKMR Anayysis
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' Date Created: October 1, 2019
#' 
#' Description:
#'  
#' This project is a follow up the the CEI paper (Env Epi, 2018). Here we're 
#' using BKMR to identify which exposures are driving the associations reported
#' in that original paper
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


#' Read in the cleaned
load(file=here::here("Data/CEI Data", "cleaned data.RData"))

#' Indicator variable for participants who are included in the final 
#' linear regression models
#' N = 897 participants included in the final dataset

load("./Data/CEI Data/env_ids.RData")
hs_data2$regression <- ifelse(hs_data2$pid %in% env_ids$pid, 1, 0)
table(hs_data2$regression)

#' complete and missing data
hs_data_complete <- filter(hs_data2, regression == 1)
hs_data_missing <- filter(hs_data2, regression == 0)
summary(hs_data_missing)

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
write_xlsx(as.data.frame(table1a_df), 
           path="./Results/CEI Results./Table 1 Full Cohort.xlsx")

table1b <- summary_table(dplyr::group_by(hs_data2, regression), tab1_summary)
table1b
table1b_df <- as.data.frame(table1b) %>%
  mutate(variable = rownames(.)) %>%
  select(3,1:2)
write_xlsx(as.data.frame(table1b_df),
           path="./Results/CEI Results./Table 1 Cohort by Regression Status.xlsx")

table1c <- summary_table(dplyr::group_by(hs_data2, include), tab1_summary)
table1c
table1c_df <- as.data.frame(table1c) %>%
  mutate(variable = rownames(.)) %>%
  select(3,1:2)
write_xlsx(as.data.frame(table1c_df),
           path="./Results/CEI Results./Table 1 Cohort by Exposure Status.xlsx")

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
pairwise.prop.test(table(hs_data2$term_cat, hs_data2$include))

#' birth_weight
t.test(hs_inc$birth_weight, hs_exc$birth_weight)

#' LBW
table(hs_data2$low_birth_weight, hs_data2$include)
chisq.test(hs_data2$low_birth_weight, hs_data2$include)
pairwise.prop.test(table(hs_data2$low_birth_weight, hs_data2$include))

#' days to PeaPod
t.test(hs_inc$days_to_peapod, hs_exc$days_to_peapod)

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

#' Do HS participants in the final analytical cohort differ from those
#' excluded based on missing exposure or covariate data?
#' Chi square test for categorical variables, Student's t for continuous
hs_inc <- filter(hs_data2, regression == 1)
hs_exc <- filter(hs_data2, regression == 0)

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

#' days to peapod assessment
t.test(hs_inc$days_to_peapod, hs_exc$days_to_peapod)

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
load("./Data/CEI Data/hs_cei.RData")
hs_data3 <- left_join(hs_data2, hs_cei, by="pid") %>%
  filter(regression == 1)

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
                     mean_sd = paste(round(mean(var[,1], na.rm=T), 1),
                                     " (",
                                     round(sd(var[,1], na.rm=T), 1),
                                     ")", sep=""),
                     min = round(min(var[,1], na.rm=T), 1),
                     q25 = round(quantile(var[,1], probs = 0.25, na.rm=T), 1),
                     median = round(median(var[,1], na.rm=T), 1),
                     q75 = round(quantile(var[,1], probs = 0.75, na.rm=T), 1),
                     q95 = round(quantile(var[,1], probs = 0.95, na.rm=T), 1),
                     max = round(max(var[,1], na.rm=T), 1),
                     cv = round(sd(var[,1], na.rm=T) / mean(var[,1], na.rm=T), 1)*100)
  table2 <- rbind(table2, temp)
  rm(temp)
}

table2
write_xlsx(as.data.frame(table2),
           path="./Results/CEI Results./Table 2 Exposures.xlsx")

#' Create "eTable 2"
#' Summary of indices by maternal race/ethnicity
#' subset of mothers

exp_vars <- c("env", "soc", "cei")

e_table2_re <- hs_data3 %>% 
  group_by(maternal_race_eth) %>% 
  summarize(counts = n(),
            mean_sd_env = paste(round(mean(env, na.rm=T), 1),
                                " (",
                                round(sd(env, na.rm=T), 1),
                                ")", sep=""),
            mean_sd_soc = paste(round(mean(soc, na.rm=T), 1),
                                " (",
                                round(sd(soc, na.rm=T), 1),
                                ")", sep=""),
            mean_sd_cei = paste(round(mean(cei, na.rm=T), 1),
                                " (",
                                round(sd(cei, na.rm=T), 1),
                                ")", sep="")) %>% 
  mutate(metric = "race-ethnicity")

e_table2_ed <- hs_data3 %>% 
  group_by(maternal_ed) %>% 
  summarize(counts = n(),
            mean_sd_env = paste(round(mean(env, na.rm=T), 1),
                                " (",
                                round(sd(env, na.rm=T), 1),
                                ")", sep=""),
            mean_sd_soc = paste(round(mean(soc, na.rm=T), 1),
                                " (",
                                round(sd(soc, na.rm=T), 1),
                                ")", sep=""),
            mean_sd_cei = paste(round(mean(cei, na.rm=T), 1),
                                " (",
                                round(sd(cei, na.rm=T), 1),
                                ")", sep="")) %>% 
  mutate(metric = "education")

e_table2 <- bind_rows(e_table2_re, e_table2_ed)

e_table2
write_xlsx(as.data.frame(e_table2),
           path="./Results/CEI Results./eTable 2 Exposures by RE and Ed.xlsx")

#' P values

