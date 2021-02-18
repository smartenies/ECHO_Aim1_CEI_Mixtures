#' =============================================================================
#' Project: ECHO Aim 1 CEI BKMR Analysis
#' Author: Sheena Martenies
#' Contact: Sheena.Martenies@colostate.edu
#' 
#' This project is a follow up the the CEI paper (Env Epi, 2018). Here we're 
#' using BKMR to identify which exposures are driving the associations reported
#' in that original paper
#' 
#' Date Created: April 15, 2020
#' 
#' Description: Using Jennifer Bobb's package to run the BKMR 
#' =============================================================================

library(sf)
library(raster)
library(ggplot2)
library(ggthemes)
library(stringr)
library(tidyverse)
library(lubridate)
library(readxl)
library(viridis)
library(bkmr)

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
         male, gest_age_w, days_to_peapod,
    
    #' Outcome 
         adiposity) %>% 
  na.omit()
nrow(hs_data2)

#' What percentage of participants do we have in this analysis?
nrow(hs_data2)/ nrow(hs_data)

#' -----------------------------------------------------------------------------
#' Bayesian Kernel Machine Regression: Birth weight
#' Use just the variables selected by the NPB model (PIP > .40)
#' -----------------------------------------------------------------------------

#' Which variables should we consider?
#' Remember to include days to peapod measurement in adiposity model!

load(here::here("Results", "NPB_Models_AD.rdata"))
exp_names

hs_data3 <- select(hs_data2, exp_names,
              #' Covariates     
                  latina_re, black_re, other_re, 
                  ed_no_hs, ed_hs, ed_aa, ed_4yr, 
                  low_bmi, ovwt_bmi, obese_bmi,
                  maternal_age, any_smoker, smokeSH, mean_cpss, mean_epsd,
                  male, gest_age_w, days_to_peapod,
                  
              #' Outcome 
                  adiposity)
names(hs_data3)

ncol(hs_data2)
ncol(hs_data3)


set.seed(123)

#' Basic version of the model with component-wise variable selection (no heirarchy)
#' Fit the BKMR model using the kmbayes() function
#' Uses an Markov chain Monte Carlo algorithm, so you need to specify the no.
#' of iterations (10,000 is pretty standard)
#' y is vector of outcomes
#' Z is matrix of exposures
#' X is matrix of covariates (each column = one variable)
#' verbose = T means print each interim output
#' varsel = T means conduct variable selection

y <- hs_data3$adiposity
Z <- as.matrix(select(hs_data3, exp_names))
X <- as.matrix(select(hs_data3, latina_re, black_re, other_re, 
                      ed_no_hs, ed_hs, ed_aa, ed_4yr, 
                      low_bmi, ovwt_bmi, obese_bmi,
                      maternal_age, any_smoker, smokeSH, mean_cpss, mean_epsd,
                      male, gest_age_w, days_to_peapod))

#' Scale these variables
y.scaled <- scale(y)
Z.scaled <- apply(Z, 2, scale)

X.s <- apply(X[,c(11, 14, 15, 17, 18)], 2, scale) #' just the continuous ones

X.scaled <- cbind(X[,1:10], X.s[,1],
                  X[,12:13], X.s[,2:3],
                  X[,16], X.s[,4:5])
colnames(X.scaled)
colnames(X.scaled)[c(11,16)] <- c("maternal_age", "male")

#' Fit the BKMR model using the scaled variables
fitkm.1 <- kmbayes(y = y.scaled, Z = Z.scaled, X = X.scaled, iter = 10000, 
                 verbose = FALSE, varsel = TRUE)

#' Save current results
save(fitkm.1, y, y.scaled, Z, Z.scaled, X, X.scaled,
     file = here::here("Results", "BKMR_Model_AD.rdata"))

#' Investigate model convergence
#' See how different parameters change with each iteration
#' beta = 
#' sigsq.ep = 
#' r = 

load(here::here("Results", "BKMR_Model_AD.rdata"))

fitkm.1
names(fitkm.1)

TracePlot(fit = fitkm.1, par = "beta")
TracePlot(fit = fitkm.1, par = "sigsq.eps")
TracePlot(fit = fitkm.1, par = "r", comp = 1)

#' Show the estimated posterior inclusion probabilities
#' PIPs are the proportion of iterations that include the predictor
#' Higher PIP means higher chance of it being in the model

arrange(ExtractPIPs(fitkm.1), desc(PIP))

#' Estimate h(z) using three different approaches
#' Need to specify at what level of exposure we are estimating h(z)
#' Using median exposures for now

med_vals <- apply(Z.scaled, 2, median) #' medians for all exposures
med_vals

Znew <- matrix(med_vals, nrow = 1) #' matrix of exposure medians
Znew

#' Approach 1: estimates the posterior mean as μh(θ^) and the posterior variance 
#' as Vh(θ^) where θ^ is the posterior mean estimate of θ. 

h_est1 <- ComputePostmeanHnew(fitkm.1, Znew = Znew, method = "approx")
h_est1

#' Approach 2: estimates the posterior mean as E[μh(θ)] by taking the mean of 
#' posterior samples of μh(θ), and the posterior variance as E[Vh(θ)]+Var[μh(θ)] 
#' by taking the mean of the posterior samples of Vh(θ) and the variance of the 
#' posterior samples of μh(θ) and them summing these quantities. 

h_est2 <- ComputePostmeanHnew(fitkm.1, Znew = Znew, method = "exact")
h_est2

#' Approach 3: Approach 3: generates posterior samples of h(z) by sampling 
#' h∼N(μh(θ),Vh(θ)), given particular samples of θ from the fitted BKMR model
#' NOTE: This is the slowest version but it allows you to estimate the full
#' posterior distribution of h(z)

Xnew <- matrix(0, nrow = 1, ncol = ncol(X))
samps3 <- SamplePred(fitkm.1, Znew = Znew, Xnew = Xnew)

#' Compare the three results
h_est_compare <- data.frame(
  method = c(1:3),
  post_mean = c(h_est1$postmean, h_est2$postmean, mean(samps3)),
  post_sd = c(sqrt(h_est1$postvar), sqrt(h_est2$postvar), sd(samps3))
)
h_est_compare

#' Summarize output

#' 1) Plot the predictor-response function
#' This is a cross-section of the surface that is generated by the model
#' Looking at one or two exposures and holding all others at a specific value
#' (Typically hold all others to a particular percentile, 50th is default)

#' First, get the univariate relationships
pred.resp.univar.50 <- PredictorResponseUnivar(fit = fitkm.1)
pred.resp.univar.25 <- PredictorResponseUnivar(fit = fitkm.1, q.fixed = 0.25)
pred.resp.univar.75 <- PredictorResponseUnivar(fit = fitkm.1, q.fixed = 0.75)

pred.resp.univar.50
pred.resp.univar.25
pred.resp.univar.75

#' Then, plot using ggplot
ggplot(pred.resp.univar.50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z) holding other variables at 50th percentile")
ggsave(filename = here::here("Figs", "AD h(x) holding other vars at 50th.jpeg"),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred.resp.univar.25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z) holding other variables at 25th percentile")
ggsave(filename = here::here("Figs", "AD h(x) holding other vars at 25th.jpeg"),
       device = "jpeg", height = 5, width = 5, units = "in")

ggplot(pred.resp.univar.75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z) holding other variables at 75th percentile")
ggsave(filename = here::here("Figs", "AD h(x) holding other vars at 75th.jpeg"),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 2) Look at bivariate relationships with two predictors
#' We can use ggplot to create contour plots for this
#' Note again that in this plot, all of the other predictors besides the two 
#' of interest are held at the 50th percentile (as the default)

pred.resp.bivar.50 <- PredictorResponseBivar(fit = fitkm.1, min.plot.dist = 1)
pred.resp.bivar.50

ggplot(pred.resp.bivar.50, aes(z1, z2, fill = est)) + 
  geom_raster() + 
  facet_grid(variable2 ~ variable1) +
  scale_fill_viridis() +
  # scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  xlab("expos1") +
  ylab("expos2") +
  ggtitle("h(expos1, expos2)")
ggsave(filename = here::here("Figs", "AD h(x1, x2) holding other vars at 50th.jpeg"),
       device = "jpeg", height = 5, width = 5, units = "in")

pred.resp.bivar.25 <- PredictorResponseBivar(fit = fitkm.1, min.plot.dist = 1,
                                             q.fixed = 0.25)
pred.resp.bivar.25

ggplot(pred.resp.bivar.25, aes(z1, z2, fill = est)) + 
  geom_raster() + 
  facet_grid(variable2 ~ variable1) +
  scale_fill_viridis() +
  # scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  xlab("expos1") +
  ylab("expos2") +
  ggtitle("h(expos1, expos2)")
ggsave(filename = here::here("Figs", "AD h(x1, x2) holding other vars at 25th.jpeg"),
       device = "jpeg", height = 5, width = 5, units = "in")

pred.resp.bivar.75 <- PredictorResponseBivar(fit = fitkm.1, min.plot.dist = 1,
                                             q.fixed = 0.75)
pred.resp.bivar.75

ggplot(pred.resp.bivar.75, aes(z1, z2, fill = est)) + 
  geom_raster() + 
  facet_grid(variable2 ~ variable1) +
  scale_fill_viridis() +
  # scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  xlab("expos1") +
  ylab("expos2") +
  ggtitle("h(expos1, expos2)")
ggsave(filename = here::here("Figs", "AD h(x1, x2) holding other vars at 75th.jpeg"),
       device = "jpeg", height = 5, width = 5, units = "in")

#' From the tutorial: Because it can be hard to see what’s going on in these 
#' types of image plots, an alternative approach is to investigate the 
#' predictor-response function of a single predictor in Z for the second 
#' predictor in Z fixed at various quantiles (and for the remaining predictors 
#' fixed to a particular value). These can be obtained using the 
#' PredictorResponseBivarLevels function, which takes as input the bivariate 
#' exposure-response function outputted from the previous command, where the 
#' argument qs specifies a sequence of quantiles at which to fix the second 
#' predictor.

#' Looking at the second predictor at the 10th 50th and 90th percentile and 
#' holding other predictors at the 50th percentile (default from previous 
#' function)
#' NOTE: need to fit PredictorResponseBivar() first!!

pred.resp.bivar.levels.50 <- PredictorResponseBivarLevels(
  pred.resp.df = pred.resp.bivar.50, 
  Z = Z.scaled, qs = c(0.1, 0.5, 0.9))

ggplot(pred.resp.bivar.levels.50, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  scale_color_viridis(discrete = T) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")
ggsave(filename = here::here("Figs", "AD h(x1, quant x2) holding other vars at 50th.jpeg"),
       device = "jpeg", height = 8, width = 8, units = "in")

#' 3) Summary statistics for the predictor-response function
#' What is the overall effect size when the exposures are a given percentile
#' compared to when they are at the 50th percentile?

risks.overall <- OverallRiskSummaries(fit = fitkm.1, y = y.scaled, 
                                      Z = Z.scaled, X = X.scaled, 
                                      qs = seq(0.25, 0.75, by = 0.05), 
                                      q.fixed = 0.5, method = "exact")
risks.overall

#' 4) Plot the risk summaries
ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  xlab("Exposure quantile") +
  ylab("Effect estimate") +
  geom_pointrange()
ggsave(filename = here::here("Figs", "AD Overall risk relative to 50th.jpeg"),
       device = "jpeg", height = 5, width = 5, units = "in")

#' 5) Estimate the risks for the individual variables the exposure is at the 75th
#' percentile relative to the 25th percentile (essentially an IQR increase), 
#' holding all other exposures at the 25th, 50th, and 75th percentiles  
risks.singvar <- SingVarRiskSummaries(fit = fitkm.1, y = y.scaled, 
                                      Z = Z.scaled, X = X.scaled, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "exact")
risks.singvar

ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
  scale_color_viridis(discrete = T) +
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
ggsave(filename = here::here("Figs", "AD Individual risk 25th to 75th.jpeg"),
       device = "jpeg", height = 5, width = 5, units = "in")

#' Interpreting this plot: We see that the predictors z3 and z4 do not 
#' contribute to the risk, and that higher values of z1 and z2 are associated 
#' with higher values of the h function. In addition, the plot suggests that for 
#' z1, as the remaining predictors increase in value from their 25th to their 
#' 75th percentile, the risk of the outcome associated with z1 increases. A 
#' similar pattern occurs for z2. This indicates the potential for interaction 
#' of z1 and z2.

#' 6) Checking for interaction 
#' To make this notion a bit more formal, we may wish to compute specific 
#' ‘interaction’ parameters. For example, we may wish to compare the 
#' single-predictor health risks when all of the other predictors in Z are fixed 
#' to their 75th percentile to when all of the other predictors in Z are fixed 
#' to their 25th percentile. In the previous plot, this corresponds to 
#' substracting the estimate represented by the purple (red) circle from the 
#' estimate represented by the yellow (blue) circle. This can be done using the 
#' function SingVarIntSummaries.

risks.int <- SingVarIntSummaries(fit = fitkm.1, y = y.scaled, 
                                 Z = Z.scaled, X = X.scaled, 
                                 qs.diff = c(0.25, 0.75), 
                                 qs.fixed = c(0.25, 0.75),
                                 method = "exact")
risks.int

ggplot(risks.int, aes(variable, est, ymin = est - 1.96*sd, 
                      ymax = est + 1.96*sd, col = variable)) + 
  scale_color_viridis(discrete = T) +
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()
ggsave(filename = here::here("Figs", "AD Interactions 25th to 75th.jpeg"),
       device = "jpeg", height = 5, width = 5, units = "in")

#' Interpreting the interaction: We see that single pollutant risk associated 
#' with a change in z1 from its 25th to 75th percentile changes very little when 
#' z2 to z4 are fixed at their 25th percentile, as compared to when z2 to z4 are 
#' fixed at their 75th percentile.

#' Saving results
save(h_est_compare, 
     pred.resp.univar.25, pred.resp.univar.50, pred.resp.univar.75,
     pred.resp.bivar.25, pred.resp.bivar.50, pred.resp.bivar.75,
     pred.resp.bivar.levels.50,
     risks.overall, risks.singvar, risks.int,
     file = here::here("Results", "BKMR_Model_AD_Results.rdata"))

