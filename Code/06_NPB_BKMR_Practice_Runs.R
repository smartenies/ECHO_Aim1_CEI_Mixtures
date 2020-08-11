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
#' Description: NPB and BKMR Practice
#' Following Lauren Hoskovec's tutorial from 'mmpack'
#' Following Jennifer Bobb's tutorial from 'bkmr'
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
#' Running through Lauren's tutorial on Nonparametric Bayesian Shrinkage (NPB)
#' -----------------------------------------------------------------------------

# library(devtools)
# install_github("lvhoskovec/mmpack", build_vignettes = TRUE)

library(mmpack)
browseVignettes("mmpack")

head(Xdat) #' simulated data for exposures

#' First, simulate the exposure-response data for 200 participants
set.seed(12345)
dat <- simexpodat(n = 200, Xdat = Xdat)

X <- dat$X # exposure data
W <- dat$W # covariate data
Y <- dat$Y # response data

#' h is the true exposure-response function
#' active = the active components of the mixture
#' active.ints = active interactions between mixture components
h <- dat$h
active <- dat$active
active.ints <- dat$active.ints

#' Distribution of our response variables
hist(Y, breaks = 20)

#' Take a look at the covariance in the exposure data
#' Some have pretty poor variance (e.g., MeBr) due to lots of zeros
var(X) 

X <- X[,-2] #' Drop MeBr due to low variance

#' Fitting the nonparametric Bayesian shrinkage model:
#' Nonparametric Bayes shrinkage (NPB) is a Bayesian linear model that places a 
#' Dirichlet Process (DP) prior on the regression coefficients to set some 
#' coefficients exactly to 0, effectively excluding them from the model, and 
#' clusters correlated exposures to reduce variance of the estimator.

#' NPB using a MCMC, so you have to specify the number of iterations (niter) and
#' the number of burn-in iterations (nburn).
#' X = matrix of exposures/predictors
#' Y = vector of continuous response data
#' W = matrix of covariate data
#' priors = optional list of hyperparameters (not specifying uses default priors)
#' scaleY = scale response prior to fit?
#' 
#' If we only want main effects, interact = F
#' If we want to include all pairwise multiplicative interactions, interact = T
#' 
#' If we want to estiamte an overall intercept, intercept = T

#' Lauren's list of NPB priors (a list)
priors.npb <- list(alpha.pi = 2, beta.pi = 2, alpha.pi2 = 2, beta.pi2 = 2)

#' Fit the NPB model
fit.npb <- npb(niter = 200, nburn = 100, X = X, Y = Y, W = W, scaleY = TRUE, 
               priors = priors.npb, interact = TRUE)

#' Take a look at the output
npb.sum <- summary(fit.npb)

#' First, main effect regression coefficients with PIP
npb.sum$main.effects

#' Next, summary of individual-level risks
head(apply(npb.sum$risk.summary, 2, FUN = function(x) round(x, 2)))

#' Plot these estimated risks against the "true" risk used to simulate the data
#' NPB does a nice job here
plot(npb.sum$risk, h)

#' Predict fitted values for each individual
pred.npb <- predict(fit.npb)
fittedvals <- pred.npb$fitted.vals

#' Plot predicted outcomes against "measured" outcomes
plot(fittedvals, Y)

#' Trace plots for the iterations
plot(fit.npb$alpha, type = "l")
plot(fit.npb$mu, type = "l")
plot(fit.npb$sig2inv, type = "l")



#' -----------------------------------------------------------------------------
#' Running through Jennifer Bobb's tutorial on Bayesian Kernel Machine 
#' Regression (BKMR)
#' Examples: 'https://jenfb.github.io/bkmr/overview.html'
#' -----------------------------------------------------------------------------

library(bkmr)

#' Basic example with component-wise variable selection

#' Create a data set using built-in functions
set.seed(111)
dat <- SimData(n = 50, M = 4)
y <- dat$y
Z <- dat$Z
X <- dat$X

#' Vizualize the underlying E-R used to generate these data
z1 <- seq(min(dat$Z[, 1]), max(dat$Z[, 1]), length = 20)
z2 <- seq(min(dat$Z[, 2]), max(dat$Z[, 2]), length = 20)
hgrid.true <- outer(z1, z2, function(x,y) apply(cbind(x,y), 1, dat$HFun))

res <- persp(z1, z2, hgrid.true, theta = 30, phi = 20, expand = 0.5, 
             col = "lightblue", xlab = "", ylab = "", zlab = "")

#' Fit the BKMR model using the kmbayes() function
#' Uses an Markov chain Monte Carlo algorithm, so you need to specify the no.
#' of iterations (10,000 is pretty standard)
#' y is vector of outcomes
#' Z is matrix of exposures
#' X is matrix of covariates (each column = one variable)
#' verbose = T means print each interim output
#' varsel = T means conduct variable selection

fitkm <- kmbayes(y = y, Z = Z, X = X, iter = 10000, 
                 verbose = FALSE, varsel = TRUE)

#' Investigate movel convergence
#' See how different parameters change with each iteration
#' beta = 
#' sigsq.ep = 
#' r = 

fitkm
names(fitkm)

TracePlot(fit = fitkm, par = "beta")
TracePlot(fit = fitkm, par = "sigsq.eps")
TracePlot(fit = fitkm, par = "r", comp = 1)

#' Show the estimated posterior inclusion probabilities
#' PIPs are the proportion of iterations that include the predictor
#' Higher PIP means higher chance of it being in the model

ExtractPIPs(fitkm)

#' Estimate h(z) using three different approaches
#' Need to specify at what level of exposure we are estimating h(z)
#' Using median exposures for now

med_vals <- apply(Z, 2, median) #' medians for all exposures
med_vals

Znew <- matrix(med_vals, nrow = 1) #' matrix of exposure medians
Znew

h_true <- dat$HFun(Znew) #' True h() based on functions used to simulate data
h_true

#' Approach 1: estimates the posterior mean as μh(θ^) and the posterior variance 
#' as Vh(θ^) where θ^ is the posterior mean estimate of θ. 

h_est1 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "approx")

#' Approach 2: estimates the posterior mean as E[μh(θ)] by taking the mean of 
#' posterior samples of μh(θ), and the posterior variance as E[Vh(θ)]+Var[μh(θ)] 
#' by taking the mean of the posterior samples of Vh(θ) and the variance of the 
#' posterior samples of μh(θ) and them summing these quantities. 

h_est2 <- ComputePostmeanHnew(fitkm, Znew = Znew, method = "exact")

#' Approach 3: Approach 3: generates posterior samples of h(z) by sampling 
#' h∼N(μh(θ),Vh(θ)), given particular samples of θ from the fitted BKMR model
#' NOTE: This is the slowest version but it allows you to estimate the full
#' posterior distribution of h(z)
set.seed(111)
samps3 <- SamplePred(fitkm, Znew = Znew, Xnew = cbind(0))

#' Compare the three results
h_est_compare <- data.frame(
  method = c("truth", 1:3),
  post_mean = c(h_true, h_est1$postmean, h_est2$postmean, mean(samps3)),
  post_sd = c(NA, sqrt(h_est1$postvar), sqrt(h_est2$postvar), sd(samps3))
)
h_est_compare

#' Summarize output

#' 1) Plot the predictor-response function
#' This is a cross-section of the surface that is generated by the model
#' Looking at one or two exposures and holding all others at a specific value
#' (Typically hold all otheres to a particular percentile, 50th is default)

#' First, get the univariate relationships
pred.resp.univar.50 <- PredictorResponseUnivar(fit = fitkm)
pred.resp.univar.25 <- PredictorResponseUnivar(fit = fitkm, q.fixed = 0.25)
pred.resp.univar.75 <- PredictorResponseUnivar(fit = fitkm, q.fixed = 0.75)

pred.resp.univar.50

#' Then, plot using ggplot
library(ggplot2)
ggplot(pred.resp.univar.50, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z) holding other variables at 50th percentile")

#' From the tutorial: Here we see that for z1 and z2, increasing values of the 
#' exposures is associated with increasing values of the outcome Y (for each of 
#' the other predictors in Z fixed to their 50th percentile, and for the 
#' covariates in x held constant). It also looks like z1 may have a nonlinear 
#' relationship with Y, with a similar suggestion of nonlinearity for z2.

ggplot(pred.resp.univar.25, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z) holding other variables at 25th percentile")

ggplot(pred.resp.univar.75, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable) +
  ylab("h(z) holding other variables at 75th percentile")

#' 2) Look at bivariate relationships with two predictors
#' We can use ggplot to create contour plots for this
#' Note again that in this plot, all of the other predictors besides the two 
#' of interest are held at the 50th percentile (as the default)

pred.resp.bivar.50 <- PredictorResponseBivar(fit = fitkm, min.plot.dist = 1)

pred.resp.bivar.50
View(pred.resp.bivar.50)

ggplot(pred.resp.bivar.50, aes(z1, z2, fill = est)) + 
  geom_raster() + 
  facet_grid(variable2 ~ variable1) +
  scale_fill_viridis() +
  # scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
  xlab("expos1") +
  ylab("expos2") +
  ggtitle("h(expos1, expos2)")

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

pred.resp.bivar.levels <- PredictorResponseBivarLevels(
  pred.resp.df = pred.resp.bivar.50, 
  Z = Z, qs = c(0.1, 0.5, 0.9))

ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1) +
  scale_color_viridis(discrete = T) +
  ggtitle("h(expos1 | quantiles of expos2)") +
  xlab("expos1")

#' 3) Summary statistics for the predictor-response function
#' What is the overall effect size when the exposures are a given percentile
#' compared to when they are at the 50th percentile?

risks.overall <- OverallRiskSummaries(fit = fitkm, y = y, Z = Z, X = X, 
                                      qs = seq(0.25, 0.75, by = 0.05), 
                                      q.fixed = 0.5, method = "exact")
risks.overall

#' 4) Plot the risk summaries
ggplot(risks.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  xlab("effect estimate when ")
  geom_pointrange()

#' 5) Estimate the risks for the individual variables the exposure is at the 75th
#' percentile relative to the 25th percentile (essentially an IQR increase), 
#' holding all other exposures at the 25th, 50th, and 75th percentiles  
risks.singvar <- SingVarRiskSummaries(fit = fitkm, y = y, Z = Z, X = X, 
                                      qs.diff = c(0.25, 0.75), 
                                      q.fixed = c(0.25, 0.50, 0.75),
                                      method = "exact")
risks.singvar

ggplot(risks.singvar, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = q.fixed)) + 
  scale_color_viridis(discrete = T) +
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()

#' Interpreting this plot: We see that the predictors z3 and z4 do not 
#' contribute to the risk, and that higher values of z1 and z2 are associated 
#' with higher values of the h function. In addition, the plot suggests that for 
#' z1, as the remaining predictors increase in value from their 25th to their 
#' 75th percentile, the risk of the outcome associated with z1 increases. A 
#' similar pattern occurs for z2. This indicates the potential for interaction 
#' of z1 and z2.

#' 6) Checking for interaction 
#' To make this notion a bit more formal, we may wish to compute specific 
#' ‘interaction’ parameters. For example, we may which to compare the 
#' single-predictor health risks when all of the other predictors in Z are fixed 
#' to their 75th percentile to when all of the other predictors in Z are fixed 
#' to their 25th percentile. In the previous plot, this corresponds to 
#' substracting the estimate represented by the purple (red) circle from the 
#' estimate represented by the yellow (blue) circle. This can be done using the 
#' function SingVarIntSummaries.

risks.int <- SingVarIntSummaries(fit = fitkm, y = y, Z = Z, X = X, 
                                 qs.diff = c(0.25, 0.75), 
                                 qs.fixed = c(0.25, 0.75),
                                 method = "exact")
risks.int

ggplot(risks.int, aes(variable, est, ymin = est - 1.96*sd, 
                          ymax = est + 1.96*sd, col = variable)) + 
  scale_color_viridis(discrete = T) +
  geom_pointrange(position = position_dodge(width = 0.75)) + 
  coord_flip()

#' Interpreting the interaction: We see that single pollutant risk associated 
#' with a change in z1 from its 25th to 75th percentile increases by 0.7 when 
#' z2 to z4 are fixed at their 25th percentile, as compared to when z2 to z4 are 
#' fixed at their 75th percentile.
