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
#' Description: Preliminary Analyses
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
#' Preliminary Analyses
#' Using birth weight to define the included population
#' -----------------------------------------------------------------------------

#' Read in the clean dataset

hs_data <- read_csv(here::here("Data", "HS_Clean_Data.csv"))
nrow(hs_data)

hs_data2 <- hs_data %>% 
  select(
    #' identifier
    pid, 
    
    #' Exposures
    mean_pm, mean_o3, pct_tree_cover, pct_impervious,
    mean_aadt_intensity, dist_m_tri:dist_m_mine_well, 
    cvd_rate_adj, res_rate_adj, violent_crime_rate, property_crime_rate,
    pct_less_hs, pct_unemp, pct_limited_eng, pct_hh_pov,
    
    #' Covariates     
    latina_re, black_re, other_re, 
    ed_no_hs, ed_hs, ed_aa, ed_4yr, 
    low_bmi, ovwt_bmi, obese_bmi,
    maternal_age, any_smoker, smokeSH, mean_cpss, mean_epsd,
    male, gest_age_w,
    
    #' Categoricals for comparisons
    maternal_race_eth, pre_preg_BMI_cat, maternal_ed, infant_sex, 
    term_cat, low_birth_weight, del_season, days_to_peapod,
    
    #' Outcomes 
    birth_weight, adiposity) %>% 
  filter(!is.na(birth_weight)) %>% 
  filter(!is.na(mean_pm)) %>% 
  filter(!is.na(smokeSH)) %>% 
  filter(!is.na(mean_cpss)) %>% 
  filter(!is.na(infant_sex))
nrow(hs_data2)

#' What percentage of participants do we have in this analysis?
nrow(hs_data2)/ nrow(hs_data)

#' How many participants were excluded?
sum(is.na(hs_data$birth_weight))
sum(is.na(hs_data$mean_pm))
sum(is.na(hs_data$smokeSH))
sum(is.na(hs_data$infant_sex))

#' excluded_data
hs_dataex <- filter(hs_data, !(pid %in% hs_data2$pid))

#' Indicator variable for hs_data
hs_data <- mutate(hs_data, include = ifelse(pid %in% hs_data2$pid, 1, 0))

#' -----------------------------------------------------------------------------
#' check out distributions of exposures and outcomes
#' Based on distributions of the built environment variables, going to go with
#' distance metrics rather than count metrics
#' -----------------------------------------------------------------------------

names(hs_data2)

#' birth weight
hist(hs_data2$birth_weight)

#' Air pollution exposures
hist(hs_data2$mean_pm)
hist(hs_data2$mean_o3)

#' Built environment exposures
hist(hs_data2$pct_tree_cover)
hist(hs_data2$pct_impervious)

hist(hs_data2$mean_aadt_intensity)

hist(hs_data2$dist_m_npl)
hist(hs_data2$dist_m_tri)
hist(hs_data2$dist_m_waste_site)
hist(hs_data2$dist_m_major_emit)

hist(hs_data$dist_m_mine_well)

#' Social determinants of health
hist(hs_data2$cvd_rate_adj)
hist(hs_data2$res_rate_adj)

hist(hs_data2$violent_crime_rate)
hist(hs_data2$property_crime_rate)

hist(hs_data2$pct_less_hs)
hist(hs_data2$pct_unemp)
hist(hs_data2$pct_limited_eng)
hist(hs_data2$pct_hh_pov)

hist(hs_data2$med_income)
hist(hs_data2$pct_poc)

#' -----------------------------------------------------------------------------
#' TABLE 1
#' check out distributions of outcomes and covariates
#' Starting with those used in the last paper
#' -----------------------------------------------------------------------------

#' Birth weight
paste0(round(mean(hs_data$birth_weight, na.rm = T), 1), " (",
       round(sd(hs_data$birth_weight, na.rm = T), 1), ")")

paste0(round(mean(hs_data2$birth_weight, na.rm = T), 1), " (",
       round(sd(hs_data2$birth_weight, na.rm = T), 1), ")")

paste0(round(mean(hs_dataex$birth_weight, na.rm = T), 1), " (",
       round(sd(hs_dataex$birth_weight, na.rm = T), 1), ")")

t.test(hs_data2$birth_weight, hs_dataex$birth_weight)

#' adiposity
nrow(hs_data) - sum(is.na(hs_data$adiposity))
paste0(round(mean(hs_data$adiposity, na.rm = T), 1), " (",
       round(sd(hs_data$adiposity, na.rm = T), 1), ")")

nrow(hs_data2) - sum(is.na(hs_data2$adiposity))
paste0(round(mean(hs_data2$adiposity, na.rm = T), 1), " (",
       round(sd(hs_data2$adiposity, na.rm = T), 1), ")")

nrow(hs_dataex) - sum(is.na(hs_dataex$adiposity))
paste0(round(mean(hs_dataex$adiposity, na.rm = T), 1), " (",
       round(sd(hs_dataex$adiposity, na.rm = T), 1), ")")

t.test(hs_data2$adiposity, hs_dataex$adiposity)

#' Race/ethnicity (1 = Latina, 2 = White, 3 = Black, 4 = Other)
table(hs_data$maternal_race_eth)
round(prop.table(table(hs_data$maternal_race_eth)), 2)

table(hs_data2$maternal_race_eth)
round(prop.table(table(hs_data2$maternal_race_eth)), 2)

table(hs_dataex$maternal_race_eth)
round(prop.table(table(hs_dataex$maternal_race_eth)), 2)

chisq.test(hs_data$maternal_race_eth, hs_data$include)
pairwise.prop.test(table(hs_data$maternal_race_eth, hs_data$include))

#' Maternal age at delivery (years)
paste0(round(mean(hs_data$maternal_age, na.rm = T), 1), " (",
       round(sd(hs_data$maternal_age, na.rm = T), 1), ")")

paste0(round(mean(hs_data2$maternal_age, na.rm = T), 1), " (",
       round(sd(hs_data2$maternal_age, na.rm = T), 1), ")")

paste0(round(mean(hs_dataex$maternal_age, na.rm = T), 1), " (",
       round(sd(hs_dataex$maternal_age, na.rm = T), 1), ")")

t.test(hs_data2$maternal_age, hs_dataex$maternal_age)

#' CPSS Score
paste0(round(mean(hs_data$mean_cpss, na.rm = T), 1), " (",
       round(sd(hs_data$mean_cpss, na.rm = T), 1), ")")

paste0(round(mean(hs_data2$mean_cpss, na.rm = T), 1), " (",
       round(sd(hs_data2$mean_cpss, na.rm = T), 1), ")")

paste0(round(mean(hs_dataex$mean_cpss, na.rm = T), 1), " (",
       round(sd(hs_dataex$mean_cpss, na.rm = T), 1), ")")

t.test(hs_data2$mean_cpss, hs_dataex$mean_cpss)

#' EPDS Score
paste0(round(mean(hs_data$mean_epsd, na.rm = T), 1), " (",
       round(sd(hs_data$mean_epsd, na.rm = T), 1), ")")

paste0(round(mean(hs_data2$mean_epsd, na.rm = T), 1), " (",
       round(sd(hs_data2$mean_epsd, na.rm = T), 1), ")")

paste0(round(mean(hs_dataex$mean_epsd, na.rm = T), 1), " (",
       round(sd(hs_dataex$mean_epsd, na.rm = T), 1), ")")

t.test(hs_data2$mean_epsd, hs_dataex$mean_epsd)

#' Pre-pregnancy BMI (1 = Under, 2 = Normal, 3 = Over, 4 = Obese)
table(hs_data$pre_preg_BMI_cat)
round(prop.table(table(hs_data$pre_preg_BMI_cat)), 2)

table(hs_data2$pre_preg_BMI_cat)
round(prop.table(table(hs_data2$pre_preg_BMI_cat)), 2)

table(hs_dataex$pre_preg_BMI_cat)
round(prop.table(table(hs_dataex$pre_preg_BMI_cat)), 2)

chisq.test(hs_data$pre_preg_BMI_cat, hs_data$include)
pairwise.prop.test(table(hs_data$pre_preg_BMI_cat, hs_data$include))

#' Maternal education (1 = <HS, 2 = HS/GED, 3 = Some College/AA, 4 = College, 5 = Grad Deg)
table(hs_data$maternal_ed)
round(prop.table(table(hs_data$maternal_ed)), 2)

table(hs_data2$maternal_ed)
round(prop.table(table(hs_data2$maternal_ed)), 2)

table(hs_dataex$maternal_ed)
round(prop.table(table(hs_dataex$maternal_ed)), 2)

chisq.test(hs_data$maternal_ed, hs_data$include)
pairwise.prop.test(table(hs_data$maternal_ed, hs_data$include))

#' Smoking during pregnancy
table(hs_data$any_smoker)
round(prop.table(table(hs_data$any_smoker)), 2)

table(hs_data2$any_smoker)
round(prop.table(table(hs_data2$any_smoker)), 2)

table(hs_dataex$any_smoker)
round(prop.table(table(hs_dataex$any_smoker)), 2)

chisq.test(hs_data$any_smoker, hs_data$include)
pairwise.prop.test(table(hs_data$any_smoker, hs_data$include))

#' SHS exposure during prenancy
table(hs_data$smokeSH)
round(prop.table(table(hs_data$smokeSH)), 2)

table(hs_data2$smokeSH)
round(prop.table(table(hs_data2$smokeSH)), 2)

table(hs_dataex$smokeSH)
round(prop.table(table(hs_dataex$smokeSH)), 2)

chisq.test(hs_data$smokeSH, hs_data$include)
pairwise.prop.test(table(hs_data$smokeSH, hs_data$include))

#' Infant sex (1 = girl, 2 = boy)
table(hs_data$male)
round(prop.table(table(hs_data$male)), 2)

table(hs_data2$male)
round(prop.table(table(hs_data2$male)), 2)

table(hs_dataex$male)
round(prop.table(table(hs_dataex$male)), 2)

chisq.test(hs_data$male, hs_data$include)
pairwise.prop.test(table(hs_data$male, hs_data$include))

#' Gestational age (weeks)
paste0(round(mean(hs_data$gest_age_w, na.rm = T), 1), " (",
       round(sd(hs_data$gest_age_w, na.rm = T), 1), ")")

paste0(round(mean(hs_data2$gest_age_w, na.rm = T), 1), " (",
       round(sd(hs_data2$gest_age_w, na.rm = T), 1), ")")

paste0(round(mean(hs_dataex$gest_age_w, na.rm = T), 1), " (",
       round(sd(hs_dataex$gest_age_w, na.rm = T), 1), ")")

t.test(hs_data2$gest_age_w, hs_dataex$gest_age_w)

#' Term status (1 = <37, 2 = 37 - <39, 3 = 39 - <41, 4 = 41 - <42, 5 = 42+)
table(hs_data$term_cat)
round(prop.table(table(hs_data$term_cat)), 2)

table(hs_data2$term_cat)
round(prop.table(table(hs_data2$term_cat)), 2)

table(hs_dataex$term_cat)
round(prop.table(table(hs_dataex$term_cat)), 2)

chisq.test(hs_data$term_cat, hs_data$include)
pairwise.prop.test(table(hs_data$term_cat, hs_data$include))

#' Low birth weight (< 2500 g)
table(hs_data$low_birth_weight)
round(prop.table(table(hs_data$low_birth_weight)), 2)

table(hs_data2$low_birth_weight)
round(prop.table(table(hs_data2$low_birth_weight)), 2)

table(hs_dataex$low_birth_weight)
round(prop.table(table(hs_dataex$low_birth_weight)), 2)

chisq.test(hs_data$low_birth_weight, hs_data$include)
pairwise.prop.test(table(hs_data$low_birth_weight, hs_data$include))

#' Days between delivery and PEA POD measurement
paste0(round(mean(hs_data$days_to_peapod, na.rm = T), 1), " (",
       round(sd(hs_data$days_to_peapod, na.rm = T), 1), ")")

paste0(round(mean(hs_data2$days_to_peapod, na.rm = T), 1), " (",
       round(sd(hs_data2$days_to_peapod, na.rm = T), 1), ")")

paste0(round(mean(hs_dataex$days_to_peapod, na.rm = T), 1), " (",
       round(sd(hs_dataex$days_to_peapod, na.rm = T), 1), ")")

t.test(hs_data2$days_to_peapod, hs_dataex$days_to_peapod)

#' -----------------------------------------------------------------------------
#' TABLE 2
#' check out distributions of exposures
#' JUST INCLUDED PARTICIPANTS
#' -----------------------------------------------------------------------------

#' Environmental exposures
paste0(round(mean(hs_data2$mean_pm, na.rm = T), 1), " (",
       round(sd(hs_data2$mean_pm, na.rm = T), 1), ")")
round(summary(hs_data2$mean_pm), 1)

paste0(round(mean(hs_data2$mean_o3, na.rm = T), 1), " (",
       round(sd(hs_data2$mean_o3, na.rm = T), 1), ")")
round(summary(hs_data2$mean_o3), 1)

paste0(round(mean(hs_data2$pct_tree_cover, na.rm = T), 1), " (",
       round(sd(hs_data2$pct_tree_cover, na.rm = T), 1), ")")
round(summary(hs_data2$pct_tree_cover), 1)

paste0(round(mean(hs_data2$pct_impervious, na.rm = T), 1), " (",
       round(sd(hs_data2$pct_impervious, na.rm = T), 1), ")")
round(summary(hs_data2$pct_impervious), 1)

paste0(round(mean(hs_data2$mean_aadt_intensity, na.rm = T), 1), " (",
       round(sd(hs_data2$mean_aadt_intensity, na.rm = T), 1), ")")
round(summary(hs_data2$mean_aadt_intensity), 1)

paste0(round(mean(hs_data2$dist_m_tri/1000, na.rm = T), 1), " (",
       round(sd(hs_data2$dist_m_tri/1000, na.rm = T), 1), ")")
round(summary(hs_data2$dist_m_tri/1000), 1)

paste0(round(mean(hs_data2$dist_m_npl/1000, na.rm = T), 1), " (",
       round(sd(hs_data2$dist_m_npl/1000, na.rm = T), 1), ")")
round(summary(hs_data2$dist_m_npl/1000), 1)

paste0(round(mean(hs_data2$dist_m_waste_site/1000, na.rm = T), 1), " (",
       round(sd(hs_data2$dist_m_waste_site/1000, na.rm = T), 1), ")")
round(summary(hs_data2$dist_m_waste_site/1000), 1)

paste0(round(mean(hs_data2$dist_m_major_emit/1000, na.rm = T), 1), " (",
       round(sd(hs_data2$dist_m_major_emit/1000, na.rm = T), 1), ")")
round(summary(hs_data2$dist_m_major_emit/1000), 1)

paste0(round(mean(hs_data2$dist_m_cafo/1000, na.rm = T), 1), " (",
       round(sd(hs_data2$dist_m_cafo/1000, na.rm = T), 1), ")")
round(summary(hs_data2$dist_m_cafo/1000), 1)

paste0(round(mean(hs_data2$dist_m_mine_well/1000, na.rm = T), 1), " (",
       round(sd(hs_data2$dist_m_mine_well/1000, na.rm = T), 1), ")")
round(summary(hs_data2$dist_m_mine_well/1000), 1)

#' Social exposures
#"cvd_rate_adj"        "res_rate_adj"        "violent_crime_rate" 
#"property_crime_rate" "pct_less_hs"         "pct_unemp"           "pct_limited_eng"     
#"pct_hh_pov"  

paste0(round(mean(hs_data2$cvd_rate_adj, na.rm = T), 1), " (",
       round(sd(hs_data2$cvd_rate_adj, na.rm = T), 1), ")")
round(summary(hs_data2$cvd_rate_adj), 1)

paste0(round(mean(hs_data2$res_rate_adj, na.rm = T), 1), " (",
       round(sd(hs_data2$res_rate_adj, na.rm = T), 1), ")")
round(summary(hs_data2$res_rate_adj), 1)

paste0(round(mean(hs_data2$violent_crime_rate, na.rm = T), 1), " (",
       round(sd(hs_data2$violent_crime_rate, na.rm = T), 1), ")")
round(summary(hs_data2$violent_crime_rate), 1)

paste0(round(mean(hs_data2$property_crime_rate, na.rm = T), 1), " (",
       round(sd(hs_data2$property_crime_rate, na.rm = T), 1), ")")
round(summary(hs_data2$property_crime_rate), 1)

paste0(round(mean(hs_data2$pct_less_hs, na.rm = T), 1), " (",
       round(sd(hs_data2$pct_less_hs, na.rm = T), 1), ")")
round(summary(hs_data2$pct_less_hs), 1)

paste0(round(mean(hs_data2$pct_unemp, na.rm = T), 1), " (",
       round(sd(hs_data2$pct_unemp, na.rm = T), 1), ")")
round(summary(hs_data2$pct_unemp), 1)

paste0(round(mean(hs_data2$pct_hh_pov, na.rm = T), 1), " (",
       round(sd(hs_data2$pct_hh_pov, na.rm = T), 1), ")")
round(summary(hs_data2$pct_hh_pov), 1)

paste0(round(mean(hs_data2$pct_limited_eng, na.rm = T), 1), " (",
       round(sd(hs_data2$pct_limited_eng, na.rm = T), 1), ")")
round(summary(hs_data2$pct_limited_eng), 1)

#' Coefficients of variation
exp_df <- select(hs_data2, pid, mean_pm:pct_hh_pov) %>% 
  pivot_longer(-pid, names_to = "exp", values_to = "value") %>% 
  group_by(exp) %>% 
  summarize(mean = mean(value),
            sd = sd(value)) %>% 
  mutate(cv = round((sd / mean) * 100,0))

#' -----------------------------------------------------------------------------
#' Correlation plot for exposures
#' -----------------------------------------------------------------------------

exp_cor <- select(hs_data2, mean_pm:pct_hh_pov)

ggcorrplot(cor(exp_cor), type = "upper", show.diag = T,
           ggtheme = simple_theme, lab = T, digits = 1,
           show.legend = T, colors = c("#FDE725FF", "white", "#440154FF"),
           legend.title = "Correlation\nCoefficient") +
  ggplot2::theme(legend.position = c(0.8, 0.3)) +
  ggplot2::scale_x_discrete(labels = c("Mean PM\u2082.\u2085", "Mean O\u2083",
                                       "%Tree Cover", "%Impervious", "AADT",
                                       "Dist. to TRI site", "Dist. to NPL site", 
                                       "Dist to waste site",
                                       "Dist. to major emitter", "Dist. to CAFO",
                                       "Dist to O&G well", "CDV Hosp. Rate", 
                                       "Res Hosp. Rate", "Violent Crime Rate", 
                                       "Property Crime Rate", "%Less than HS", 
                                       "% Unemployed", "%HH Ltd. English",
                                       "%HH Poverty")) +
  ggplot2::scale_y_discrete(labels = c("Mean PM\u2082.\u2085", "Mean O\u2083",
                                       "%Tree Cover", "%Impervious", "AADT",
                                       "Dist. to TRI site", "Dist. to NPL site", 
                                       "Dist to waste site",
                                       "Dist. to major emitter", "Dist. to CAFO",
                                       "Dist to O&G well", "CDV Hosp. Rate", 
                                       "Res Hosp. Rate", "Violent Crime Rate", 
                                       "Property Crime Rate", "%Less than HS", 
                                       "% Unemployed", "%HH Ltd. English",
                                       "%HH Poverty")) 
ggsave(here::here("Figs", "Exp_Correlations.jpeg"), device = "jpeg",
       width = 7, height = 7, units = "in", dpi = 500)
  
  
#' -----------------------------------------------------------------------------
#' Simple linear regressions between the exposures and birth weight
#' -----------------------------------------------------------------------------

#' Air pollution exposures
ggplot(hs_data2) +
  geom_point(aes(x = mean_pm, y = birth_weight)) +
  simple_theme
pm_bw_lm <- lm(birth_weight ~ mean_pm, data = hs_data2)
summary(pm_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = mean_o3, y = birth_weight)) +
  simple_theme
o3_bw_lm <- lm(birth_weight ~ mean_o3, data = hs_data2)
summary(o3_bw_lm)

#' Built environment exposures
ggplot(hs_data2) +
geom_point(aes(x = pct_tree_cover, y = birth_weight)) +
  simple_theme
tc_bw_lm <- lm(birth_weight ~ pct_tree_cover, data = hs_data2)
summary(tc_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_impervious, y = birth_weight)) +
  simple_theme
imp_bw_lm <- lm(birth_weight ~ pct_impervious, data = hs_data2)
summary(imp_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = mean_aadt_intensity, y = birth_weight)) +
  simple_theme
aadt_bw_lm <- lm(birth_weight ~ mean_aadt_intensity, data = hs_data2)
summary(aadt_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_npl, y = birth_weight)) +
  simple_theme
npl_bw_lm <- lm(birth_weight ~ dist_m_npl, data = hs_data2)
summary(npl_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_tri, y = birth_weight)) +
  simple_theme
tri_bw_lm <- lm(birth_weight ~ dist_m_tri, data = hs_data2)
summary(tri_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_waste_site, y = birth_weight)) +
  simple_theme
waste_bw_lm <- lm(birth_weight ~ dist_m_waste_site, data = hs_data2)
summary(waste_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_major_emit, y = birth_weight)) +
  simple_theme
major_emit_bw_lm <- lm(birth_weight ~ dist_m_major_emit, data = hs_data2)
summary(major_emit_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_mine_well, y = birth_weight)) +
  simple_theme
mine_well_bw_lm <- lm(birth_weight ~ dist_m_mine_well, data = hs_data2)
summary(mine_well_bw_lm)

#' Social determinants of health
ggplot(hs_data2) +
  geom_point(aes(x = cvd_rate_adj, y = birth_weight)) +
  simple_theme
cvd_bw_lm <- lm(birth_weight ~ cvd_rate_adj, data = hs_data2)
summary(cvd_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = res_rate_adj, y = birth_weight)) +
  simple_theme
res_bw_lm <- lm(birth_weight ~ res_rate_adj, data = hs_data2)
summary(res_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = violent_crime_rate, y = birth_weight)) +
  simple_theme
vcrime_bw_lm <- lm(birth_weight ~ violent_crime_rate, data = hs_data2)
summary(vcrime_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = property_crime_rate, y = birth_weight)) +
  simple_theme
pcrime_bw_lm <- lm(birth_weight ~ property_crime_rate, data = hs_data2)
summary(pcrime_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_less_hs, y = birth_weight)) +
  simple_theme
less_hs_bw_lm <- lm(birth_weight ~ pct_less_hs, data = hs_data2)
summary(less_hs_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_unemp, y = birth_weight)) +
  simple_theme
unemp_bw_lm <- lm(birth_weight ~ pct_unemp, data = hs_data2)
summary(unemp_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_limited_eng, y = birth_weight)) +
  simple_theme
limited_eng_bw_lm <- lm(birth_weight ~ pct_limited_eng, data = hs_data2)
summary(limited_eng_bw_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_hh_pov, y = birth_weight)) +
  simple_theme
hh_pov_bw_lm <- lm(birth_weight ~ pct_hh_pov, data = hs_data2)
summary(hh_pov_bw_lm)

#' -----------------------------------------------------------------------------
#' Simple linear regressions between the exposures and adiposity
#' -----------------------------------------------------------------------------

#' Air pollution exposures
ggplot(hs_data2) +
  geom_point(aes(x = mean_pm, y = adiposity)) +
  simple_theme
pm_ad_lm <- lm(adiposity ~ mean_pm, data = hs_data2)
summary(pm_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = mean_o3, y = adiposity)) +
  simple_theme
o3_ad_lm <- lm(adiposity ~ mean_o3, data = hs_data2)
summary(o3_ad_lm)

#' Built environment exposures
ggplot(hs_data2) +
  geom_point(aes(x = pct_tree_cover, y = adiposity)) +
  simple_theme
tc_ad_lm <- lm(adiposity ~ pct_tree_cover, data = hs_data2)
summary(tc_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_impervious, y = adiposity)) +
  simple_theme
imp_ad_lm <- lm(adiposity ~ pct_impervious, data = hs_data2)
summary(imp_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = mean_aadt_intensity, y = adiposity)) +
  simple_theme
aadt_ad_lm <- lm(adiposity ~ mean_aadt_intensity, data = hs_data2)
summary(aadt_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_npl, y = adiposity)) +
  simple_theme
npl_ad_lm <- lm(adiposity ~ dist_m_npl, data = hs_data2)
summary(npl_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_tri, y = adiposity)) +
  simple_theme
tri_ad_lm <- lm(adiposity ~ dist_m_tri, data = hs_data2)
summary(tri_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_waste_site, y = adiposity)) +
  simple_theme
waste_ad_lm <- lm(adiposity ~ dist_m_waste_site, data = hs_data2)
summary(waste_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_major_emit, y = adiposity)) +
  simple_theme
major_emit_ad_lm <- lm(adiposity ~ dist_m_major_emit, data = hs_data2)
summary(major_emit_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = dist_m_mine_well, y = adiposity)) +
  simple_theme
mine_well_ad_lm <- lm(adiposity ~ dist_m_mine_well, data = hs_data2)
summary(mine_well_ad_lm)

#' Social determinants of health
ggplot(hs_data2) +
  geom_point(aes(x = cvd_rate_adj, y = adiposity)) +
  simple_theme
cvd_ad_lm <- lm(adiposity ~ cvd_rate_adj, data = hs_data2)
summary(cvd_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = res_rate_adj, y = adiposity)) +
  simple_theme
res_ad_lm <- lm(adiposity ~ res_rate_adj, data = hs_data2)
summary(res_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = violent_crime_rate, y = adiposity)) +
  simple_theme
vcrime_ad_lm <- lm(adiposity ~ violent_crime_rate, data = hs_data2)
summary(vcrime_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = property_crime_rate, y = adiposity)) +
  simple_theme
pcrime_ad_lm <- lm(adiposity ~ property_crime_rate, data = hs_data2)
summary(pcrime_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_less_hs, y = adiposity)) +
  simple_theme
less_hs_ad_lm <- lm(adiposity ~ pct_less_hs, data = hs_data2)
summary(less_hs_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_unemp, y = adiposity)) +
  simple_theme
unemp_ad_lm <- lm(adiposity ~ pct_unemp, data = hs_data2)
summary(unemp_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_limited_eng, y = adiposity)) +
  simple_theme
limited_eng_ad_lm <- lm(adiposity ~ pct_limited_eng, data = hs_data2)
summary(limited_eng_ad_lm)

ggplot(hs_data2) +
  geom_point(aes(x = pct_hh_pov, y = adiposity)) +
  simple_theme
hh_pov_ad_lm <- lm(adiposity ~ pct_hh_pov, data = hs_data2)
summary(hh_pov_ad_lm)