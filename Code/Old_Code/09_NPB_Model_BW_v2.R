#' =============================================================================
#' Project: ECHO Aim 1 CEI BKMR Analysis
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' This project is a follow up the the CEI paper (Env Epi, 2018). Here we're 
#' using BKMR to identify which exposures are driving the associations reported
#' in that original paper
#' 
#' Date Created: June 3, 2020
#' 
#' Description: After discussions with Ander, Lauren, and Sheryl, we decided to
#' drop the less interesting BKMR analysis and stick with the NPB. One of our
#' concerns was the potential interaction between individual level covariates 
#' and the neighborhood-level environmental and social factors. In response to
#' this, Lauren modified the mmpack package to allow for interactions between
#' exposures and covariates. 
#' 
#' This script follows her guidance on how to fit the NPB model. In this script,
#' we are using adiposity (%fat mass) as the outcome of interest
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
#' Install the new version of the package
#' -----------------------------------------------------------------------------

# library(devtools)
# install_github("lvhoskovec/mmpack", build_vignettes = TRUE, force = TRUE)

library(mmpack)

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
    male, gest_age_w, days_to_peapod,
    
    #' Outcomes 
    birth_weight) %>% 
  na.omit()
nrow(hs_data2)

#' What percentage of participants do we have in this analysis?
nrow(hs_data2)/ nrow(hs_data)

#' -----------------------------------------------------------------------------
#' Nonparametric Bayesian Shrinkage (NPB): Birth weight
#' -----------------------------------------------------------------------------

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
            male, gest_age_w, days_to_peapod) %>% 
  as.matrix()
head(W)

#' Scale covariates
W.s <- apply(W[,c(11, 14, 15, 17, 18)], 2, scale) #' just the continuous ones

W.scaled <- cbind(W[,1:10], W.s[,1],
                  W[,12:13], W.s[,2:3],
                  W[,16], W.s[,4:5])
colnames(W.scaled)
colnames(W.scaled)[c(11,16)] <- c("maternal_age", "male")
head(W.scaled)

# response data- adiposity
Y <- select(hs_data2, birth_weight) %>% 
  as.matrix()
head(Y)

#' Distribution of our response variable
hist(Y, breaks = 20)

#' Scatter plots of the relationship between scaled(Y) and each exposure
#' What do linear relationships between the outcome and each predictor look like?
df <- as.data.frame(cbind(scale(Y), X.scaled))

# par(mfrow=c(5,4))
sapply(2:length(df), function(x){
  lm.x <- lm(birth_weight ~ df[,x], data = df) 
  plot(df[,c(x, 1)],
       xlab = paste0(colnames(df)[x], " beta: ", 
                     round(summary(lm.x)$coef[2,1],4),
                     "; p = ", round(summary(lm.x)$coef[2,4],5)))
  abline(lm.x)
})

#' -----------------------------------------------------------------------------
#' Finding the NPB priors-- start with Lauren's from the example in the vignette
# In an email from April 29, Lauren provided me with some additional guidance
# on finding the NPB priors:
# 1) Keep alpha.pi and beta.pi set to 1, and then let a.phi1 take values 1, 10,
#    and 100 and see how the results change.
# 2) Keep a.phi1 set to 1 (or 10 or 100), and mess with alpha.pi and beta.pi.
#    Run the following code:
#        alpha.pi=1
#        beta.pi=1
#        plot(density(rbeta(10000, alpha.pi, beta.pi)))
#    and then change alpha.pi and beta.pi and see how it changes the prior
#    distribution. This is the distribution of the probability of a main effect
#    regression coefficient being 0 (aka exclusion probability). We don't want
#    this to be too informative (you don't want high mass around just a few
#    values). Also alpha.pi and beta.pi don't have to be the same value. You
#    might try alpha.pi = 1 and beta.pi = 2 to get a slightly lower prior
#    probability of exclusion.
#    Try not to change all three (alpha.pi, beta.pi, and a.phi1) at once.
# 3) When playing with the priors, set "interact=FALSE" and just fit the model
# with the main effects. Most of the interactions were null anyway so it
# shouldn't change the results too much and it will make the code run a lot
# faster. Then when you find a set of priors you like, you can add in
# "interact=TRUE" and "XWinteract=TRUE."
#' -----------------------------------------------------------------------------

set.seed(123)

#' Try with with default priors from the original vignette
priors.npb.1 <- list(alpha.pi = 1, beta.pi = 1, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 1)

fit.npb.1 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.1, interact = F)
npb.sum.1 <- summary(fit.npb.1)
npb.sum.1$main.effects
plot(fit.npb.1$beta[,1], type = "l")

#' Adjust a.phi1 (keep all other priors the same)
#' Try making a.phi1 = 10
#' PIPs go up slightly, the effects are in the expected direction
priors.npb.2 <- list(alpha.pi = 1, beta.pi = 1, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 10)

fit.npb.2 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.2, interact = F)
npb.sum.2 <- summary(fit.npb.2)
npb.sum.2$main.effects
plot(fit.npb.2$beta[,1], type = "l")

#' Try making a.phi1 = 100
#' PIPs again increase slightly
priors.npb.3 <- list(alpha.pi = 1, beta.pi = 1, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 100)

fit.npb.3 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.3, interact = F)
npb.sum.3 <- summary(fit.npb.3)
npb.sum.3$main.effects
plot(fit.npb.3$beta[,1], type = "l")

#' Set a.phi1 = 10, try adjusting alpha.pi and beta.pi
#' Try making alpha.pi = 2
#' PIPs went down
plot(density(rbeta(10000, 2, 1)))

priors.npb.4 <- list(alpha.pi = 2, beta.pi = 1, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 10)

fit.npb.4 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.4, interact = F)
npb.sum.4 <- summary(fit.npb.4)
npb.sum.4$main.effects
plot(fit.npb.4$beta[,1], type = "l")

#' Try making beta.pi = 2
#' PIPs increase
plot(density(rbeta(10000, 1, 2)))
priors.npb.5 <- list(alpha.pi = 1, beta.pi = 2, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 10)

fit.npb.5 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.5, interact = F)
npb.sum.5 <- summary(fit.npb.5)
npb.sum.5$main.effects
plot(fit.npb.5$beta[,1], type = "l")

#' Try making beta.pi = 3
plot(density(rbeta(10000, 1, 3)))
priors.npb.6 <- list(alpha.pi = 1, beta.pi = 3, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 10)

fit.npb.6 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.6, interact = F)
npb.sum.6 <- summary(fit.npb.6)
npb.sum.6$main.effects
plot(fit.npb.6$beta[,1], type = "l")

#' Try making beta.pi = 2.5
plot(density(rbeta(10000, 1, 2.5)))
priors.npb.7 <- list(alpha.pi = 1, beta.pi = 2.5, alpha.pi2 = 9, beta.pi2 = 1,
                     a.phi1 = 10)

fit.npb.7 <- npb(niter = 1000, nburn = 100, X = X.scaled, Y = Y, W = W.scaled, 
                 scaleY = TRUE, 
                 priors = priors.npb.7, interact = F)
npb.sum.7 <- summary(fit.npb.7)
npb.sum.7$main.effects
plot(fit.npb.7$beta[,1], type = "l")

#' -----------------------------------------------------------------------------
#' Fitting the nonparametric Bayesian shrinkage model:
#' Going to go with the 5th set of priors as described above
#' 
#' NOTE: I'll need to check with Lauren to see if the rbeta density plot is too 
#' narrow. 
#' 
#' Also note: based on Lauren's feedback, I think this set of priors would be
#' interpreted as we expect about half the exposures to show up. I'll need to 
#' check to see if this interpretation makes sense
#' -----------------------------------------------------------------------------

priors.npb <- priors.npb.5

fit.npb <- npb(niter = 10000, nburn = 1000, X = X.scaled, Y = Y, W = W.scaled, 
               scaleY = TRUE, 
               priors = priors.npb, interact = TRUE, XWinteract = T)
save(fit.npb, file = here::here("Results", "NPB_Birth_Weight.rdata"))

load(here::here("Results", "NPB_Birth_Weight.rdata"))
npb.sum <- summary(fit.npb)

#' First, main effect regression coefficients with PIPs
#' Look for PIPs above 0.50
npb.sum$main.effects
selected_exp <- which(npb.sum$main.effects[,5] >= 0.5)
selected_exp

#' Which variables are these?
exp_names <- colnames(X)[selected_exp]
exp_names

#' Next, all of the interactions between exposures or between exposures and covariates
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
     file = here::here("Results", "NPB_Models_AD_Full.rdata"))
