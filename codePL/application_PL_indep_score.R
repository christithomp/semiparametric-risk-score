#####################################################################
# This program fits the additive model using PL under NHANES 
# application for when the score is independent of and dependent on k. 
#####################################################################
rm(list = ls())
## Derive variables from NHANES to fit single index model
source("./final_code/nhanes_preprocessing.R")
source("./final_code/application_PL_update_parameters.R")

#####################################################################
### Score independent of k
#####################################################################
## Organize iteration-invariant variables to fit additive model
betaM <- -1 ## identifiability constraint
betaW <- 0.5 ## starting value of betaW
Y <- rbind(yM, yW) ## response
X.all = rbind(xM, xW) ## PA variables
X.all_reg <- rbind(xM_reg, xW_reg)
X01 <- matrix(c(rep(1, length(yM)), rep(0, length(yW))), ncol = 1) ## design matrix of intercept term
X02 <- matrix(c(rep(0, length(yM)), rep(1, length(yW))), ncol = 1)
Z1 <- rbind(zM, matrix(0, nrow = length(yW), ncol = ncol(zM))) ## design matrix of offset
Z2 <- rbind(matrix(0, nrow = length(yM), ncol = ncol(zW)), zW)
dummyID <- factor(rep(1, length(Y)))

## Fit a penalized spline SIM
spline <- "os"
nknots <- 8

# Get spline basis matrix
B <- c()

for(i in 1:ncol(X.all)){
  x <- matrix(X.all[,i], ncol = 1)
  x.mean <- matrix(mean(X.all[,i]), nrow = nrow(X.all), ncol = 1)
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
  # Obtain new B and X matrices weighted by beta
  B.new <- rbind(betaM * B[1:nrow(xM),], betaW * B[(nrow(xM)+1):nrow(X.all),]) 
  X.all.new <- rbind(betaM * matrix(X.all[1:nrow(xM),], ncol = ncol(xM)), 
                     betaW * matrix(X.all[(nrow(xM)+1):nrow(X.all),], ncol = ncol(xW)) )
  
  # Get updated parameters
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
  # sigma_a = updated_estimates$sigma_a
  
  # Update log likelihood
  loglike <- sum(Y*log(plogis(eta)) + (1-Y)*log(1-plogis(eta))) ## log likelihood
  epsilon <- abs((loglike-loglike_old)/loglike_old)
  loglike_old <- loglike
  nc <- nc + 1
  print(paste0("change in log likelihood: ", epsilon))
}
print("Convergence complete!")
loglike_indep = loglike_old
# number of parameters = 1 slope, 2 intercept, ncol(theta_1), ncol(theta_2), spline
k_indep <- 3 + length(thetaW) + length(thetaM) + length(alpha) + length(alpha_int) 



#############################################################################
### Score depends on k
#############################################################################

## Fit model for men
BM <- B[1:nrow(xM),]
dummyID.m <- rep(1, length = nrow(xM)) # for random component
suppressMessages(fitM <- gamm(yM ~ zM + xM, random = list(dummyID.m = pdIdent(~-1+BM)), family = "binomial"))

# Obtain coefficients
alphaM <- matrix(unname(unlist(fitM$lme$coefficients$random)), ncol = 1)
beta01 <- unname(fitM$gam$coefficients)[1]
thetaM <- unname(fitM$gam$coefficients)[2:(ncol(zM) + 1)]
alpha.intM <- unname(fitM$gam$coefficients)[(ncol(zM) + 2):length(fitM$gam$coefficients)]

# Calculate score for men
scoreM <- BM %*% alphaM + xM %*% alpha.intM


## Fit model for women
BW <- B[(1+nrow(xM)):nrow(B),]
dummyID.w <- rep(1, length = nrow(xW)) # for random component
suppressMessages(fitW <- gamm(yW ~ zW + xW, random = list(dummyID.w = pdIdent(~-1+BW)), family = "binomial"))

# Obtain coefficients
alphaW <- matrix(unname(unlist(fitW$lme$coefficients$random)), ncol = 1)
beta02 <- unname(fitW$gam$coefficients)[1]
thetaW <- unname(fitW$gam$coefficients)[2:(ncol(zM) + 1)]
alpha.intW <- unname(fitW$gam$coefficients)[(ncol(zM) + 2):length(fitW$gam$coefficients)]

# Calculate score for women
scoreW <- BW %*% alphaW + xW %*% alpha.intW


## Calculate log likelihood for dependent score
eta_dep <- X01 %*% beta01 + X02 %*% beta02 + Z1 %*% thetaM + Z2 %*% thetaW + 
  rbind(matrix(scoreM, ncol = 1), matrix(scoreW, ncol = 1))
loglike_dep <- sum(Y*log(plogis(eta_dep)) + (1-Y)*log(1-plogis(eta_dep)))
k_dep <- length(fitM$gam$coefficients) + length(alphaM) + length(fitW$gam$coefficients) + length(alphaW)


#########################################################################
# Calculate AIC and BIC for each model
##########################################################################
## Calculate AIC for each model
# AIC = 2k - 2loglike
AIC_indep <- 2 * k_indep - 2* loglike_indep
AIC_dep <- 2 * k_dep - 2* loglike_dep

## Calculate BIC for each model
# BIC = klog(n) - 2loglike
BIC_indep <- log(nrow(X.all)) * k_indep - 2* loglike_indep
BIC_dep <- log(nrow(X.all)) * k_dep - 2* loglike_dep

# Print results
results <- cbind(c(AIC_indep, AIC_dep), c(BIC_indep, BIC_dep))
rownames(results) <- c("Score is independent of population", 
                       "Score is dependent of population")
colnames(results) <- c("AIC", "BIC")
print(results)
