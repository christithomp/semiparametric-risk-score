## This script performs simulations using JAGS
## PA and non-PA data are sampled without replacement from real NHANES dataset
## True parameters are specified following the procedure in S.1 of Ma et al. 2017

rm(list = ls())

set.seed(100)

## Load package
library(rlang)
library(ggplot2)
library(corrplot)
library(mgcv)
library(caret)
library(survey)
library(HRW)
library(gridExtra)
library(MASS)
library(rjags)
library(dplyr)
library(mvtnorm)

##########################################################################
# Organize real NHANES data
##########################################################################

## Load NHANES data
data_analysis <- read.csv("./final_code/data/NHANES_data.csv")
nhanes_data <- data_analysis

## Windsorize PA vairables by upper 95th percentile and 5th percentile
nhanes_data[which(data_analysis$MVPA > quantile(data_analysis$MVPA, 0.95)), "MVPA"] <- quantile(data_analysis$MVPA, 0.95)
nhanes_data[which(data_analysis$ASTP > quantile(data_analysis$ASTP, 0.95)), "ASTP"] <- quantile(data_analysis$ASTP, 0.95)
nhanes_data[which(data_analysis$MVPA < quantile(data_analysis$MVPA, 0.05)), "MVPA"] <- quantile(data_analysis$MVPA, 0.05)
nhanes_data[which(data_analysis$ASTP < quantile(data_analysis$ASTP, 0.05)), "ASTP"] <- quantile(data_analysis$ASTP, 0.05)

## Scale PA variables
nhanes_data[,c("MVPA_reg", "ASTP_reg")] <- nhanes_data[,c("MVPA", "ASTP")]
nhanes_data[,c("MVPA", "ASTP")] <- apply(nhanes_data[,c("MVPA", "ASTP")] , 2, scale)

## Surival outcome
dep_vars <- "yr9_mort"
## remove people who are censored before interested survival time -- 0 subject
nhanes_data <- nhanes_data[which(!is.na(nhanes_data[,dep_vars])),]

## Reorganize the covariates 
# Divide age by 100 for numeric stability
nhanes_data$Age <- nhanes_data$Age / 100

# Create indicator for smoker
nhanes_data$Smoker <- ifelse(nhanes_data$SmokeCigs == "Current", 1, 0)

## Separate data by gender
dataM <- nhanes_data[which(nhanes_data$Gender == 'Male'), ]
dataW <- nhanes_data[which(nhanes_data$Gender == 'Female'), ]

## Derive components of single index model
### Response
yM <- as.matrix(dataM[,dep_vars], ncol = 1)
yW <- as.matrix(dataW[,dep_vars], ncol = 1)
colnames(yM) <- colnames(yW) <- "mortstat"

### Non-accelerometry covariates
zM <- as.matrix(cbind(dataM$Age, dataM$Smoker), ncol = 2)
zW <- as.matrix(cbind(dataW$Age, dataW$Smoker), ncol = 2)
colnames(zM) <- colnames(zW) <- varZ <- c("age", "smoker")

### Unscaled accelerometry variables
varX_reg <- c("MVPA_reg", "ASTP_reg")
xM_reg <- as.matrix(dataM[,varX_reg], ncol = length(varX_reg))
xW_reg <- as.matrix(dataW[,varX_reg], ncol = length(varX_reg))
colnames(xM_reg) <- colnames(xW_reg) <- varX_reg

### Accelerometry variables
varX <- c("MVPA", "ASTP")
# varX <- c("MVPA")
xM <- as.matrix(dataM[,varX], ncol = length(var))
xW <- as.matrix(dataW[,varX], ncol = length(var))
colnames(xM) <- colnames(xW) <- varX

## remove unnecessary variables
df_real_m <- data.frame(yM, xM, zM)
df_real_w <- data.frame(yW, xW, zW)
X.all_reg <- rbind(xM_reg, xW_reg)
rm(dep_vars, data_analysis, nhanes_data, dataM, dataW)
rm(xM, yM, zM, xW, yW, zW)

## Specify true values of parameters
sMVPA <- function(x){ -(-0.2*(x-0.8)^3-0.4) }
sASTP <- function(x){ -(0.3*exp(x)-1.25) }
beta1_m <- -1 ## identifiability constraint
set.seed(100)
beta1_w <- -1*runif(1, 0.5, 1.5)
theta_age_m <- runif(1, 7, 10)
theta_age_w <- runif(1, 7, 10)
theta_smoking_m <- runif(1, 0.5, 0.8)
theta_smoking_w <- runif(1, 0.5, 0.8)
beta0_w <- runif(1, -8, -5)
beta0_m <- runif(1, -8, -5)

##########################################################################
# START SIMULATION
##########################################################################

nsim <- 200
est.sim <- matrix(NA, nrow = 7, ncol = nsim) ## point estimates
lower.sim <- matrix(NA, nrow = 7, ncol = nsim) ## estimated lower bounds
upper.sim <- matrix(NA, nrow = 7, ncol = nsim) ## estimated upper bounds
curve.sim <- list() ## estimates of curves
time <- rep(0, nsim) ## running time of each simulation

for(iter in 1:nsim){
  ptm <- proc.time()
  
  ##########################################################################
  # Generate simulated data
  ##########################################################################
  
  ## Sample X and Z from 1000 men and 1000 women without replacement
  nsample <- 1000
  sample_ind_m <- sample(1:nrow(df_real_m), nsample, replace = FALSE)
  sample_ind_w <- sample(1:nrow(df_real_w), nsample, replace = FALSE)
  df_sample_m <- df_real_m[sample_ind_m,2:ncol(df_real_m)]
  df_sample_w <- df_real_w[sample_ind_w,2:ncol(df_real_w)]
  
  ## Generate binary responses y
  ## impose a centering constraint of each true functions
  eta_m <- beta0_m + beta1_m*(sMVPA(df_sample_m[,"MVPA"]) - mean(sMVPA(c(df_real_m$MVPA, df_real_w$MVPA))) +
                              sASTP(df_sample_m[,"ASTP"]) - mean(sASTP(c(df_real_m$ASTP, df_real_w$ASTP)))) + 
    theta_age_m * df_sample_m[,"age"] + theta_smoking_m *  df_sample_m[,"smoker"]
  eta_w <- beta0_w + beta1_w*(sMVPA(df_sample_w[,"MVPA"]) - mean(sMVPA(c(df_real_m$MVPA, df_real_w$MVPA))) +
                              sASTP(df_sample_w[,"ASTP"]) - mean(sASTP(c(df_real_m$ASTP, df_real_w$ASTP)))) +
    theta_age_w * df_sample_w[,"age"] + theta_smoking_w *  df_sample_w[,"smoker"]
  p_m <- plogis(eta_m)
  df_sample_m$mortstat <- rbinom(length(eta_m), size = 1, prob = p_m)
  p_w <- plogis(eta_w)
  df_sample_w$mortstat <- rbinom(length(eta_w), size = 1, prob = p_w)
  rm(eta_m, p_m, eta_w, p_w)
  
  ##########################################################################
  # Fit the model using simulated data
  ##########################################################################
  
  ## Reorganize the data to allow different coefficients between populations
  ## "df_sample" is a dataset with the form suitable for model fitting
  df_sample <- data.frame(mortstat = c(df_sample_m$mortstat, df_sample_w$mortstat),
                          int_m = c(rep(1, nrow(df_sample_m)), rep(0, nrow(df_sample_w))),
                          int_w = c(rep(0, nrow(df_sample_m)), rep(1, nrow(df_sample_w))),
                          age_m = c(df_sample_m$age, rep(0, nrow(df_sample_w))),
                          smoker_m = c(df_sample_m$smoker, rep(0, nrow(df_sample_w))),
                          age_w = c(rep(0, nrow(df_sample_m)), df_sample_w$age),
                          smoker_w = c(rep(0, nrow(df_sample_m)), df_sample_w$smoker))
  df_sample[,varX] <- rbind(df_sample_m, df_sample_w)[,varX]
  
  ## Specify the additive model formula -- to obtain spline basis and penalty matrices from mgcv
  form <- paste0("mortstat ~ -1 + ", paste(colnames(df_sample)[2:(ncol(df_sample)-length(varX))], collapse = " + "))
  for(pa in varX){
    form <- paste0(form, " + s(", pa, ", bs='ps', k=11)")
  }
  
  ## Fit the additive model in mgcv -- to obtain spline basis and corresponding penalty matrices from mgcv
  fit_fake <- gam(as.formula(form), data = df_sample, family = binomial, select = TRUE) 
  X_lp <- predict.gam(fit_fake, type = "lpmatrix") ## derive linear predictor matrix
  S <- list() ## a list storing penalty matrices
  for(i in 1:length(varX)){
    S[[i]] <- cbind(fit_fake$smooth[[i]]$S[[1]], fit_fake$smooth[[i]]$S[[2]])
  }
  rm(form, fit_fake)
  
  ## Build the dataset ready to fit Bayesian model in JAGS
  load.module("glm") 
  df_jags <- list(y = df_sample$mortstat, n = nrow(df_sample), n1 = nrow(df_sample_m),
                    X = X_lp, S1 = S[[1]], S2 = S[[2]], zero = rep(0, ncol(X_lp)))
  ## Fit the Bayesian graphical model using JAGS
  jm <- jags.model("./code_v4/simulation.jags", data = df_jags, n.chains = 1, n.adapt = 500, quiet = TRUE) 
  
  ## Specify number of burn-in
  update(jm, n.burn = 5000)
  
  ## Obtain posterior samples
  samples <- jags.samples(jm, c("b", "c"), n.iter = 2000, thin = 2, quite = TRUE)
  
  ## Derive estimated beta_1_w
  beta1_w_samples <- as.vector(samples$c)
  est.sim[1,iter] <- mean(beta1_w_samples)
  lower.sim[1,iter] <- mean(beta1_w_samples) - 1.96*sd(beta1_w_samples)
  upper.sim[1,iter] <- mean(beta1_w_samples) + 1.96*sd(beta1_w_samples)
  
  ## Derive estimated thetas
  theta_samples <- samples$b[,,1]
  est.sim[2:7,iter] <- apply(theta_samples[1:6,], 1, mean)
  lower.sim[2:7,iter] <- apply(theta_samples[1:6,], 1, mean) - 1.96* apply(theta_samples[1:6,], 1, sd)
  upper.sim[2:7,iter] <- apply(theta_samples[1:6,], 1, mean) + 1.96* apply(theta_samples[1:6,], 1, sd)
  
  ## Derive curve estimates on an equally spaced grid of the range of each variable
  grid_length <- 100 ## number of locations of each curve range
  xind <- list()
  for(i in 1:length(varX)){
    xind[[varX[i]]] <- seq(min(c(df_real_m[,varX[i]], df_real_w[,varX[i]])), 
                           max(c(df_real_m[,varX[i]], df_real_w[,varX[i]])),
                           length.out = grid_length)
  }
  
  curve <- list()
  nbasis <- (ncol(X_lp) - ncol(df_sample) + length(varX) + 1) / length(varX) ## number of basis functions per curve
  for(i in 1:length(varX)){
    ind_basis <- ncol(X_lp) - nbasis*length(varX) + ((i-1)*nbasis+1):(i*nbasis)
    s_samples <- df_jags$X[,ind_basis] %*% samples$b[ind_basis,,1]
    shat <- apply(s_samples, 1, mean) ## posterior means
    shat_lower <- apply(s_samples, 1, function(x) quantile(x, 0.025)) ## pointwise CI
    shat_upper <- apply(s_samples, 1, function(x) quantile(x, 0.975))
    
    ## smooth curves on a specified grid
    x <- df_sample[,varX[i]]
    sm_mean <- gam(shat ~ s(x), method = "REML") ## smooth the mean
    s_x <- unname(predict(sm_mean, newdata = data.frame(x = xind[[varX[i]]]))) ## get smooth estimates on specified grid
    s_x_mean <- mean(s_x) ## average value across all locations of the grid
    s_x <- s_x - s_x_mean ## identifiability constraint
    sm_lower <- gam(shat_lower ~ s(x), method = "REML") ## smooth the lower bound
    s_lower <- unname(predict(sm_lower, newdata = data.frame(x = xind[[varX[i]]]))) - s_x_mean
    sm_upper <- gam(shat_upper ~ s(x), method = "REML") ## smooth the upper bound
    s_upper <- unname(predict(sm_upper, newdata = data.frame(x = xind[[varX[i]]]))) - s_x_mean
  
    curve[[varX[i]]] <- data.frame(x = xind[[varX[i]]], s_x = s_x, lower = s_lower, 
                                   upper = s_upper, s_x_mean = s_x_mean)
  }
  curve.sim[[iter]] <- curve
  
  print(iter)
  
  time[iter] <- (proc.time() - ptm)[3]
}

##########################################################################
# FINISH SIMULATION
##########################################################################

##########################################################################
# Calculate coverages of all parameters
##########################################################################

save(est.sim, lower.sim, upper.sim, curve.sim, time, file = "./results/simulation_jags_results.rda")
load(file = "./final_code/results/simulation_jags_results.rda")

true_param <- c(beta1_w, beta0_m, beta0_w, theta_age_m, theta_smoking_m, theta_age_w, theta_smoking_w)
apply(est.sim, 1, mean)
true_param

## coverage of each parameter
coverage <- rep(0, 7)
for(i in 1:7){
  coverage[i] <- length(which(upper.sim[i,] >= true_param[i] & lower.sim[i,] <= true_param[i]))/nsim
}
coverage

## coverage of curve estimates
grid_length <- 100
coverage.curve <- list()
coverage_curve_result <- matrix(0, nrow = grid_length, ncol = length(varX))
colnames(coverage_curve_result) <- varX
for(var in varX){
  coverage.curve[[var]] <- matrix(0, nrow = grid_length, ncol = nsim)
  coverage.joint <- rep(0, nsim)
  for(i in 1:nsim){
    if(var == "MVPA"){
      ind_cover <- which(curve.sim[[i]][[var]]$lower <= sMVPA(curve.sim[[i]][[var]]$x) - mean(sMVPA(curve.sim[[i]][[var]]$x) ) & 
                           curve.sim[[i]][[var]]$upper >= sMVPA(curve.sim[[i]][[var]]$x) - mean(sMVPA(curve.sim[[i]][[var]]$x)) )
    }else if(var == "ASTP"){
      ind_cover <- which(curve.sim[[i]][[var]]$lower <= sASTP(curve.sim[[i]][[var]]$x) - mean(sASTP(curve.sim[[i]][[var]]$x) ) & 
                           curve.sim[[i]][[var]]$upper >= sASTP(curve.sim[[i]][[var]]$x) - mean(sASTP(curve.sim[[i]][[var]]$x)) )
    }
    coverage.curve[[var]][ind_cover, i] <- 1
  }
  coverage_curve_result[,var] <- apply(coverage.curve[[var]], 1, mean)
}
cov <- apply(coverage_curve_result, 2, mean)
cov

##########################################################################
# Visualize the curve simulation results
##########################################################################
# List to store plot
simplot <- list()
j <- 1
# Loop over pa variables
for(var in varX){
  curves_all <- c()
  # Loop over datasets
  for(i in 1:nsim){
    
    # Combine all curves into dataframe
    curves_all <- rbind(curves_all, data.frame(x = curve.sim[[i]][[var]]$x*sd(X.all_reg[,j]) + mean(X.all_reg[,j]),
                                               s_x = curve.sim[[i]][[var]]$s_x,
                                               sim = i))
  }
  # Get true curve values 
  if(var == "MVPA"){
    xlab <- "Min/day"
    curve_truth <- data.frame(x = curve.sim[[1]][[var]]$x *sd(X.all_reg[,j]) + mean(X.all_reg[,j]),
                              s_x = sMVPA(curve.sim[[1]][[var]]$x) - mean(sMVPA(curve.sim[[1]][[var]]$x)))
  }else if(var == "ASTP"){
    xlab <- "Transition probability"
    curve_truth <- data.frame(x = curve.sim[[1]][[var]]$x *sd(X.all_reg[,j]) + mean(X.all_reg[,j]),
                              s_x = sASTP(curve.sim[[1]][[var]]$x) - mean(sASTP(curve.sim[[1]][[var]]$x)))
  }
  # Obtain mean of estimated curves
  temp_curves <- c()
  for(i in 1:nsim){
    temp_curves <- cbind(temp_curves, curve.sim[[i]][[var]]$s_x)
  }
  curve_mean <- data.frame(x = curve.sim[[1]][[var]]$x*sd(X.all_reg[,j]) + mean(X.all_reg[,j]),
                           s_x = apply(temp_curves, 1, mean))
  
  # Find only the simulations completely within the middle 50th percentile
  curves <- matrix(curves_all$s_x, ncol = nsim, nrow = grid_length, byrow = FALSE)
  # bounds <- apply(curves, 1, function(x) quantile(x, c(0.25, 0.75)))
  # inner.curves <- apply(curves, 2, function(x) which(x >= bounds[1, ] & x <= bounds[2, ]))
  # sim.kept <- which(lapply(inner.curves, function(x) length(x) >= 50) == TRUE)
  # 
  # Find only the simulations with the squared difference between the estimate and the true fit 
  # in the top 50th percentile
  m <- apply(curves, 2, function(x) sum((x - curve_truth$s_x)^2))
  sim.kept <- which(m <= quantile(m, 0.5))
  
  j <- j + 1
  rm(temp_curves)
  
  # Organize results into plots
  simplot[[var]] <- ggplot() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_line(aes(x = x, y = s_x, group = sim, color = "Estimates"), data = curves_all[which(curves_all$sim %in% sim.kept),], alpha = 0.2) +
    geom_line(aes(x = x, y = s_x, color = "Truth"), data = curve_truth, lwd = 1) +
    geom_line(aes(x = x, y = s_x, color = "Average Estimate"), data = curve_mean, lty = 5) +
    scale_colour_manual(name="",
                        values=c("Truth"="orangered2", "Average Estimate"="black", "Estimates"="gray")) +
    # labs(x = paste0("Scaled ",var), y = paste0("s(", var, ")"), title = paste0(var, " Average Coverage = ", round( cov[var], 2))) +
    labs(x = xlab, y = paste0("s(", var, ")"), title = paste0(var, " Average Coverage = ", round( cov[var], 2)),
         subtitle = "Bayesian Method") +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title=element_blank(),
          legend.position = c(0.99, 0.99),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(0, 0, 0, 0),
          legend.background = element_rect(fill=alpha('white', 0))) 
}
simplot[["MVPA"]] <- simplot[["MVPA"]] + ylim(-1.5, 1.5)
simplot[["ASTP"]] <- simplot[["ASTP"]] + ylim(-2.2, 1.2)
# Plot simulation results for curve estimates with true curve
pdf(file = "./final_code/figures/simulation_jags_MVPA_ASTP_spaghetti_rescaled.pdf", width = 10, height = 5)
do.call("grid.arrange", c(simplot, nrow = 1))
dev.off()

simplot_jags <- simplot
save(simplot_jags, file = "./final_code/results/simplot_jags_list.rda")
