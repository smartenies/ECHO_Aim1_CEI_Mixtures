#' =============================================================================
#' Project: ECHO Aim 1 CEI BKMR Analysis
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' This project is a follow up the the CEI paper (Env Epi, 2018). Here we're 
#' using BKMR to identify which exposures are driving the associations reported
#' in that original paper
#' 
#' Date Created: October 21, 2019
#' 
#' Description: Using Lauren's mmpack to run NBP on the HS data set using birth
#' weight as the outcome of interest
#' =============================================================================

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
#' Read in the dataset and get the exposures, outcome and covariates
#' -----------------------------------------------------------------------------

#' Read in the clean dataset
hs_data <- read_csv(here::here("Data", "HS_Clean_Data.csv"))
nrow(hs_data)

hs_data2 <- filter(hs_data, include == 1) %>% 
  select(
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
    
    #' Outcomes 
         birth_weight) %>% 
  na.omit()
nrow(hs_data2)

#' What percentage of participants do we have in this analysis?
nrow(hs_data2)/ nrow(hs_data)

#' -----------------------------------------------------------------------------
#' Nonparametric Bayesian Shrinkage (NPB): Birth weight
#' -----------------------------------------------------------------------------

# library(devtools)
# install_github("lvhoskovec/mmpack", build_vignettes = TRUE)
# browseVignettes("mmpack")

library(mmpack)

#' Exposure data
X <- select(hs_data2, mean_pm, mean_o3, pct_tree_cover, pct_impervious,
            mean_aadt_intensity, dist_m_tri:dist_m_mine_well, 
            cvd_rate_adj, res_rate_adj, violent_crime_rate, property_crime_rate,
            pct_less_hs, pct_unemp, pct_limited_eng, pct_hh_pov) %>% 
  as.matrix()
head(X)

#' Take a look at the covariance in the exposure data
var(X) 
X.scaled <- apply(X, 2, scale)
var(X.scaled)

#' Covariate data
#' References: white_re, ed_grad, norm_bmi
W <- select(hs_data2, latina_re, black_re, other_re, 
            ed_no_hs, ed_hs, ed_aa, ed_4yr, 
            low_bmi, ovwt_bmi, obese_bmi,
            maternal_age, any_smoker, smokeSH, mean_cpss, mean_epsd,
            male, gest_age_w) %>% 
  as.matrix()
head(W)

#' Scale covariates
W.s <- apply(W[,c(11, 14, 15, 17)], 2, scale) #' just the continuous ones

W.scaled <- cbind(W[,1:10], W.s[,1],
                  W[,12:13], W.s[,2:3],
                  W[,16], W.s[,4])
colnames(W.scaled)
colnames(W.scaled)[c(11,16,17)] <- c("maternal_age", "male", "gest_age_w")
head(W.scaled)

# response data- birthweight
Y <- select(hs_data2, birth_weight) %>% 
  as.matrix()
head(Y)

#' Distribution of our response variable
hist(Y, breaks = 20)

#' -----------------------------------------------------------------------------
#' Finging the NPB priors-- start with Lauren's from the example in the vignette
#'     Lauren's list of NPB priors (a list)
#' -----------------------------------------------------------------------------

set.seed(123)

#' Trying with default priors
priors.npb.1 <- list(alpha.pi = 1, beta.pi = 1, alpha.pi2 = 9, beta.pi2 = 1)

#' Run with standardized X and standardized W
fit.npb.1 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
               priors = priors.npb.1, interact = TRUE)
npb.sum.1 <- summary(fit.npb.1)
npb.sum.1$main.effects
plot(fit.npb.1$beta[,1], type = "l")

#' Try making a.phi1 bigger
priors.npb.2 <- list(alpha.pi = 1, beta.pi = 1, alpha.pi2 = 9, beta.pi2 = 1,
                    a.phi1 = 10)

#' Run with standardized X and standardized W
fit.npb.2 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.2, interact = TRUE)
npb.sum.2 <- summary(fit.npb.2)
npb.sum.2$main.effects
plot(fit.npb.2$beta[,1], type = "l")

#' Try making alpha.pi and beta.pi bigger
priors.npb.3 <- list(alpha.pi = 3, beta.pi = 3, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 50)

#' Run with standardized X and standardized W
fit.npb.3 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.3, interact = TRUE)
npb.sum.3 <- summary(fit.npb.3)
npb.sum.3$main.effects
plot(fit.npb.3$beta[,1], type = "l")

#' Try making a.phi1 much, much bigger
priors.npb.4 <- list(alpha.pi = 3, beta.pi = 3, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 100)

#' Run with standardized X and standardized W
fit.npb.4 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.4, interact = TRUE)
npb.sum.4 <- summary(fit.npb.4)
npb.sum.4$main.effects
plot(fit.npb.4$beta[,1], type = "l")

#' Try increasing alpha.pi and beta.pi one more time
priors.npb.5 <- list(alpha.pi = 10, beta.pi = 10, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 100)

#' Run with standardized X and standardized W
fit.npb.5 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.5, interact = TRUE)
npb.sum.5 <- summary(fit.npb.5)
npb.sum.5$main.effects
plot(fit.npb.5$beta[,1], type = "l")

#' -----------------------------------------------------------------------------
#' Fitting the nonparametric Bayesian shrinkage model:
#' Going to go with the fifth set of priors as described above
#' Because this is a screening model, going to set PIP criterion to 0.40
#' -----------------------------------------------------------------------------

priors.npb <- priors.npb.5

fit.npb <- npb(niter = 10000, nburn = 1000, X = X.scaled, Y = Y, W = W.scaled, 
               scaleY = TRUE, 
               priors = priors.npb, interact = TRUE)
npb.sum <- summary(fit.npb)

#' First, main effect regression coefficients with PIP
#' 8 variables had PIPs above 0.40
#' Also going to include PM2.5 (PIP = 0.38)
npb.sum$main.effects
selected_exp <- c(1,2,4,11,12,13,15,16,17,19)

#' Which variables are these?
exp_names <- colnames(X)[selected_exp]
exp_names

#' Next, all of the interactions
#' No interactions showed up here
npb.sum$interactions

#' Predict fitted values for each individual
pred.npb <- predict(fit.npb)
fittedvals <- pred.npb$fitted.vals

#' Plot predicted outcomes against "measured" outcomes
plot(fittedvals, Y)
abline(a = 0, b = 1, col = "red")

#' Trace plots for the iterations
plot(fit.npb$alpha, type = "l")
plot(fit.npb$mu, type = "l")
plot(fit.npb$sig2inv, type = "l")
plot(fit.npb$beta[,1], type = "l", main = colnames(X)[1])
plot(fit.npb$beta[,2], type = "l", main = colnames(X)[2])
plot(fit.npb$beta[,4], type = "l", main = colnames(X)[4])
plot(fit.npb$beta[,11], type = "l", main = colnames(X)[11])
plot(fit.npb$beta[,12], type = "l", main = colnames(X)[12])
plot(fit.npb$beta[,13], type = "l", main = colnames(X)[13])
plot(fit.npb$beta[,15], type = "l", main = colnames(X)[15])
plot(fit.npb$beta[,16], type = "l", main = colnames(X)[16])
plot(fit.npb$beta[,17], type = "l", main = colnames(X)[17])
plot(fit.npb$beta[,19], type = "l", main = colnames(X)[19])

save(priors.npb, fit.npb, npb.sum, exp_names,
     file = here::here("Results", "NPB_Models_BW.rdata"))
