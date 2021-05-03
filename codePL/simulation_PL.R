#####################################################
# This code perform simulations for the profile 
# likelihood estimation procedure, following the idea 
# in Ma et al.(2017)
#####################################################
# Read the command line arguments
command_args <- as.numeric(commandArgs(trailingOnly = TRUE))

# Define command line arguments as a variable
dataset <- command_args[1] #seed for dataset
bootstrap <- command_args[2] #seed for bootstrap

# Source code
source("/users/ecui/riskscore/code/application_PL_update_parameters.R")

## Load packages
library(rlang)
library(mgcv)
library(HRW)

##########################################################################
# Organize real data
##########################################################################
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
colnames(yM) <- colnames(yW) <- "All-Cause Mortality"

### Non-accelerometry covariates
zM <- as.matrix(cbind(dataM$Age, dataM$Smoker), ncol = 2)
zW <- as.matrix(cbind(dataW$Age, dataW$Smoker), ncol = 2)
colnames(zM) <- colnames(zW) <- c("Age", "Smoker")

### Accelerometry variables
var <- c("MVPA", "ASTP")
xM <- as.matrix(dataM[,var], ncol = length(var))
xW <- as.matrix(dataW[,var], ncol = length(var))
colnames(xM) <- colnames(xW) <- var

## remove unnecessary variables
rm(dep_vars, data_analysis, nhanes_data, dataM, dataW)

# Copies of original matrices
yMO <- yM; yWO <- yW
xMO <- xM; xWO <- xW
zMO <- zM; zWO <- zW

# Define spline type and knots
spline <- "os"
nknots <- 8

## Specify true values of parameters
s.MVPA <- function(x){ -0.2*(x-0.8)^3-0.4 }
s.ASTP <- function(x){ 0.3*exp(x)-1.25 }
beta1.m <- -1 ## identifiability constraint
theta.age.m <- 8
theta.age.w <- 9
theta.smoking.m <- 0.6
theta.smoking.w <- 0.7
beta1.w <- 0.8
beta0.m <- -6
beta0.w <- -7

##########################################################################
# Generate simulated dataset
##########################################################################
set.seed(dataset)
n.sample = 1000
## Sample X and Z from n.sample men and n.sample women without replacement
sample.ind.m <- sample(1:nrow(xMO), n.sample, replace = FALSE)
sample.ind.w <- sample(1:nrow(xWO), n.sample, replace = FALSE)
xM.sample <- xMO[sample.ind.m,]
zM.sample <- zMO[sample.ind.m,]
xW.sample <- xWO[sample.ind.w,]
zW.sample <- zWO[sample.ind.w,]
Xall.sample <- rbind(xM.sample, xW.sample)

## Generate binary responses
eta.m <- beta0.m + beta1.m*(s.MVPA(xM.sample[,1]) + s.ASTP(xM.sample[,2])) +
  theta.age.m * zM.sample[,1] + theta.smoking.m * zM.sample[,2]
p.m <- plogis(eta.m)
yM.sample <- matrix(rbinom(length(eta.m), size = 1, prob = p.m), ncol = 1)
eta.w <- beta0.w + beta1.w*(s.MVPA(xW.sample[,1]) + s.ASTP(xW.sample[,2])) +
  theta.age.w * zW.sample[,1] + theta.smoking.w * zW.sample[,2]
p.w <- plogis(eta.w)
yW.sample <- matrix(rbinom(length(eta.w), size = 1, prob = p.w), ncol = 1)

## remove unnecessary variables
rm(eta.m, eta.w, p.m, p.w)

##########################################################################
# Fit the model with bootstrap sample on dataset
##########################################################################
set.seed(bootstrap)

# Get number of men and women for sampling
nrM <- nrow(yM.sample)
nrW <- nrow(yW.sample)

# Obtain indices of bootstrap samples
nM <- sample(1:nrM, nrM, replace = TRUE)
nW <- sample(1:nrW, nrW, replace = TRUE)

# Create new data with indices
yM <- yM.sample[nM, ]; yW <- yW.sample[nW, ]
xM <- xM.sample[nM, ]; xW <- xW.sample[nW, ]
zM <- zM.sample[nM, ]; zW <- zW.sample[nW, ]

## Organize iteration-invariant variables to fit single index model
betaM <- -1
betaW <- 0.5
Y <- matrix(c(yM, yW), ncol = 1) ## response
X.all = rbind(xM, xW) ## PA variables
X01 <- matrix(c(rep(1, length(yM)), rep(0, length(yW))), ncol = 1) ## design matrix of intercept term
X02 <- matrix(c(rep(0, length(yM)), rep(1, length(yW))), ncol = 1)
Z1 <- rbind(zM, matrix(0, nrow = length(yW), ncol = ncol(zM))) ## design matrix of offset
Z2 <- rbind(matrix(0, nrow = length(yM), ncol = ncol(zW)), zW)
dummyID <- factor(rep(1, length(Y)))

# Get spline basis matrix
B <- c()
for(i in 1:ncol(X.all)){
  x <- matrix(X.all[,i], ncol = 1)
  if(spline == "os"){
    # O'Sullivan splines
    # Set range of X values and interior knots
    a <- 1.01*min(x) - 0.01*max(x)
    b <- 1.01*max(x) - 0.01*min(x)
    intKnots <- quantile(unique(x), seq(0, 1, length = (nknots + 2))[-c(1,(nknots + 2))])

    # Get O'Sullivan spline basis functions
    Bg <- ZOSull(x, range.x = c(a,b), intKnots = intKnots)
    B <- cbind(B, Bg)
  }else{
    # Spline from mgcv
    fit <- gam(Y ~ -1 + s(x, k = nknots, fx = TRUE, bs = spline), family = "gaussian") ## fit a fake model
    B <- cbind(B, predict(fit, type = "lpmatrix"))
    rm(fit)
  }
}
rm(x)

# Do iteration
threshold <- 1e-6 ## change in log-likelihood
loglike_old <- 1 ## initial value of log likelihood
nc <- 0
epsilon <- 1
while(nc < 1000 & epsilon > threshold){
  # Reweight B and X by beta
  B.new <- rbind(betaM * B[1:nrow(xM),], betaW * B[(nrow(xM)+1):nrow(X.all),])
  X.all.new <- rbind(betaM * matrix(X.all[1:nrow(xM),], ncol = ncol(xM)),
                     betaW * matrix(X.all[(nrow(xM)+1):nrow(X.all),], ncol = ncol(xW)) )
  # Obtain updated estimates
  updated_estimates = update_parameters()
  alpha <- updated_estimates$alphahat
  alpha_int <- updated_estimates$alphahat_int
  beta01 <- updated_estimates$beta01hat
  beta02 <- updated_estimates$beta02hat
  thetaM <- updated_estimates$thetaMhat
  thetaW <- updated_estimates$thetaWhat
  betaM <- updated_estimates$betaMhat
  betaW <- updated_estimates$betaWhat
  eta <- updated_estimates$eta
  score <- updated_estimates$score
  scoreM <- score[1:nrow(xM)]
  scoreW <- score[(nrow(xM)+1):(nrow(xM) + nrow(xW))]
  curve = updated_estimates$curve
  sigma_a = updated_estimates$sigma_a

  # Update log likelihood
  loglike <- sum(Y*log(plogis(eta)) + (1-Y)*log(1-plogis(eta))) ## log likelihood
  epsilon <- abs((loglike-loglike_old)/loglike_old)
  loglike_old <- loglike
  nc <- nc + 1
}
# Get results from simulation
result <- list()
for(pa in var){
  X.samp <- X.all[, pa]
  X <- Xall.sample[, pa]
  result[[pa]] <- data.frame(x = X[order(X)], s.x = curve[[pa]]$s.x[order(X.samp)])
}
result$param.est <- c(beta01, beta02, betaW, thetaM, thetaW)

# Create job-specific output so we don't overwrite things!
filename <- sprintf("/users/ecui/riskscore/sim_output/dataset_%s_simulation_%s.RData", dataset, bootstrap)
save("result", file = filename)