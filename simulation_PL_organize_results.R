#########################################################################
# This script calculates the parameter variance and curve estimates of 
# the simulations using the PL method. 
#########################################################################

##########################################################################
# Organize real NHANES data
##########################################################################

library(mgcv)
## Load NHANES data
data_analysis <- read.csv("/users/ecui/riskscore/data/NHANES_data.csv")
nhanes_data <- data_analysis

## Windsorize PA vairables by upper 95th percentile and 5th percentile
nhanes_data[which(data_analysis$MVPA > quantile(data_analysis$MVPA, 0.95)), "MVPA"] <- quantile(data_analysis$MVPA, 0.95)
nhanes_data[which(data_analysis$ASTP > quantile(data_analysis$ASTP, 0.95)), "ASTP"] <- quantile(data_analysis$ASTP, 0.95)
nhanes_data[which(data_analysis$MVPA < quantile(data_analysis$MVPA, 0.05)), "MVPA"] <- quantile(data_analysis$MVPA, 0.05)
nhanes_data[which(data_analysis$ASTP < quantile(data_analysis$ASTP, 0.05)), "ASTP"] <- quantile(data_analysis$ASTP, 0.05)

## Scale PA variables
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

### Accelerometry variables
varX <- c("MVPA", "ASTP")
xM <- as.matrix(dataM[,varX], ncol = length(var))
xW <- as.matrix(dataW[,varX], ncol = length(var))
colnames(xM) <- colnames(xW) <- varX

## remove unnecessary variables
df_real_m <- data.frame(yM, xM, zM)
df_real_w <- data.frame(yW, xW, zW)
rm(dep_vars, data_analysis, nhanes_data, dataM, dataW)
rm(xM, yM, zM, xW, yW, zW)

# Load the simulation results
sim_data <- readRDS("/users/ecui/riskscore/result/simulation_PL_results.rds")
nsim <- length(sim_data) # number of simulated datasets
nbs <- length(sim_data$Dataset1) #number of bs samples per dataset (301)
nsample <- length(sim_data$Dataset1$BS1$MVPA$s.x)

##############################################################
# Organize simulation results
##############################################################

# Matrices to store iterations
est.sim <- matrix(NA, nrow = 7, ncol = nsim) # point estimates of parameters
param.sd.sim <- matrix(NA, nrow = 7, ncol = nsim) # sd of parameters
curve.sim <- list() ## estimates of curves

# Loop over datasets
for(iter in 1:nsim){
  
  # Store each bootstrap sample estimates
  param.est <- matrix(NA, nrow = 7, ncol = (nbs-1))
  curve.est <- list()
  for(pa in varX){
    curve.est[[pa]]$curvemat <- matrix(NA, nrow = nsample, ncol = (nbs-1))
  }
  
  # Loop over bootstrap samples 1-300
  for(bs in 1:(nbs - 1)){
    bs.samp <- sim_data[[paste0("Dataset", iter)]][[paste0("BS", bs)]]
    param.est[, bs] <- bs.samp$param.est
    for(pa in varX){
      curve.est[[pa]]$curvemat[, bs] <- bs.samp[[pa]]$s.x
    }
  }
  # Obtain inference from bootstrap samples
  est.sim[, iter] <- sim_data[[paste0("Dataset", iter)]][[paste0("BS", 301)]]$param.est #point estimate on bs sample 301
  param.sd.sim[, iter] <- apply(param.est, 1, sd)
  
  ## Derive curve estimates on an equally spaced grid of the range of each variable
  grid_length <- 100 ## number of locations of each curve range
  xind <- list()
  for(i in 1:length(varX)){
    xind[[varX[i]]] <- seq(min(c(df_real_m[,varX[i]], df_real_w[,varX[i]])), 
                           max(c(df_real_m[,varX[i]], df_real_w[,varX[i]])),
                           length.out = grid_length)
  }
  # Obtain smoothed curve estimates 
  curve <- list()
  for(i in 1:length(varX)){
    curve.bs <- curve.est[[varX[i]]]$curvemat
    # Obtain curve estimate from bootstrap samples
    shat <- apply(curve.bs, 1, mean) ## 
    shat_lower <- apply(curve.bs, 1, function(x) quantile(x, 0.025)) ## pointwise CI
    shat_upper <- apply(curve.bs, 1, function(x) quantile(x, 0.975))
    
    ## smooth curves on a specified grid
    x <- sim_data[[paste0("Dataset", iter)]][[paste0("BS", bs)]][[varX[i]]]$x
    sm_mean <- gam(shat ~ s(x), method = "REML") ## smooth the mean
    s_x <- unname(predict(sm_mean, newdata = data.frame(x = xind[[varX[i]]]))) ## get smooth estimates on specified grid
    s_x_mean <- mean(s_x) ## average value across all locations of the grid
    s_x <- s_x - s_x_mean ## identifiability constraint
    sm_lower <- gam(shat_lower ~ s(x), method = "REML") ## smooth the lower bound
    s_lower <- unname(predict(sm_lower, newdata = data.frame(x = xind[[varX[i]]]))) - s_x_mean ## identifiability constraint
    sm_upper <- gam(shat_upper ~ s(x), method = "REML") ## smooth the upper bound
    s_upper <- unname(predict(sm_upper, newdata = data.frame(x = xind[[varX[i]]]))) - s_x_mean ## identifiability constraint
    
    # Organize curve results into list
    curve[[varX[i]]] <- data.frame(x = xind[[varX[i]]], s_x = s_x, lower = s_lower, 
                                   upper = s_upper, s_x_mean = s_x_mean)
  }
  # Save estimates for the curve
  curve.sim[[iter]] <- curve
  
  print(iter)
}

# Save organized results of pointwise, variance and curve estimates from simulation
save(est.sim, param.sd.sim, curve.sim, file = "/users/ecui/riskscore/result/simulation_PL_cleaned_results.rda")
