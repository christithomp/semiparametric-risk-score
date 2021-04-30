#####################################################################
# This program fits the additive model using JAGS under NHANES 
# application for when the score is independent of and dependent on k. 
#####################################################################
rm(list = ls())
set.seed(100)

## Derive variables from NHANES to fit single index model
source("./final_code/nhanes_preprocessing.R")

#####################################################################
### Score independent of k
#####################################################################
# Organize variables to fit model
Y <- rbind(yM, yW) ## response
X = rbind(xM, xW) ## PA measures
X01 <- matrix(c(rep(1, length(yM)), rep(0, length(yW))), ncol = 1) ## design matrix of intercept term
X02 <- matrix(c(rep(0, length(yM)), rep(1, length(yW))), ncol = 1)
Z1 <- rbind(zM, matrix(0, nrow = length(yW), ncol = ncol(zM))) ## design matrix of other predictors
Z2 <- rbind(matrix(0, nrow = length(yM), ncol = ncol(zW)), zW)

## dataset containing variables of interest
data.nhanes <- data.frame(Y, X01, X02, Z1, Z2, X)
colnames(data.nhanes) <- c("mort", "X0.m", "X0.w", paste0(colnames(zM), ".m"),
                           paste0(colnames(zW), ".w"), var)

## specify the additive model formula -- to obtain spline basis and corresponding penalty matrices from mgcv
form <- paste0("mort ~ -1 + ", paste(colnames(data.nhanes)[2:(ncol(data.nhanes)-ncol(X))], collapse = " + "))
for(pa in var){
  form <- paste0(form, " + s(", pa, ", bs = 'tp', k = 11)")
}
## fit the additive model in mgcv
fit <- gam(as.formula(form), data = data.nhanes, family = binomial, select = TRUE) 
X.mat <- predict.gam(fit, type = "lpmatrix") ## derive linear predictor matrix
S1 <- cbind(fit$smooth[[1]]$S[[1]], fit$smooth[[1]]$S[[2]]) ## penalty matrices for MVPA
S2 <- cbind(fit$smooth[[2]]$S[[1]], fit$smooth[[2]]$S[[2]]) ## penalty matrices for ASTP

## build the dataset to fit Bayesian model
data.jags <- list(y = as.vector(Y), n = length(Y), n1 = length(yM),
                  X = X.mat, S1 = S1, S2 = S2, zero = rep(0, ncol(X.mat)))

## fit the Bayesian graphical model using JAGS
load.module("glm") 
jm <- jags.model("./final_code/application.jags", data = data.jags, n.chains = 1, n.adapt = 1000) 

## specify number of burn-in
update(jm, n.burn = 10000)

## obtain posterior samples
samples <- jags.samples(jm, c("b", "c"), n.iter = 5000, thin = 5)

## estimate beta_{12}
betaW.sample <- as.vector(samples$c)
betaW.hat <- mean(betaW.sample)
# betaW.hat + 2 * sd(betaW.sample)

## estimate of predictors
theta.sample <- samples$b[,,1]
theta.hat <- apply(theta.sample, 1, mean)

# Calculate log likelihood
n1 <- length(yM)
n <- length(Y)
ind.spl <- (3 + ncol(zM) + ncol(zW)):length(theta.hat)
eta <- matrix(0, nrow = n, ncol = 1)
eta[1:n1] <- X.mat[1:n1, 1:(ncol(zM) + ncol(zW) + 2)] %*% theta.hat[1:(ncol(zM) + ncol(zW) + 2)]+
  X.mat[1:n1, ind.spl] %*% theta.hat[ind.spl] * -1
eta[(n1+1):n] <- X.mat[(n1+1):n, 1:(ncol(zM) + ncol(zW) + 2)] %*% theta.hat[1:(ncol(zM) + ncol(zW) + 2)] + 
  X.mat[(n1+1):n, ind.spl] %*% theta.hat[ind.spl] * betaW.hat[1]
loglike_indep <- sum(Y*log(plogis(eta)) + (1-Y)*log(1-plogis(eta))) ## log likelihood
k_indep <- length(betaW.hat) + length(theta.hat)


#####################################################################
### Score dependent of k
#####################################################################
## Create dataset for males 
dataM.nhanes <- data.frame(yM, zM, xM)
colnames(dataM.nhanes) <- c("mort", paste0(colnames(zM), ".m"), var)

## specify the additive model formula -- to obtain spline basis and corresponding penalty matrices from mgcv
formM <- paste0("mort ~ ", paste(colnames(dataM.nhanes)[2:(ncol(dataM.nhanes)-ncol(X))], collapse = " + "))
for(pa in var){
  formM <- paste0(formM, " + s(", pa, ", bs = 'tp', k = 11)")
}
## fit the additive model in mgcv
fitM <- gam(as.formula(formM), data = dataM.nhanes, family = binomial, select = TRUE) 
xM.mat <- predict.gam(fitM, type = "lpmatrix") ## derive linear predictor matrix
S1M <- cbind(fitM$smooth[[1]]$S[[1]], fitM$smooth[[1]]$S[[2]]) ## penalty matrices for MVPA
S2M <- cbind(fitM$smooth[[2]]$S[[1]], fitM$smooth[[2]]$S[[2]]) ## penalty matrices for ASTP

## build the dataset to fit Bayesian model
dataM.jags <- list(y = as.vector(yM), n = length(yM),
                   X = xM.mat, S1 = S1M, S2 = S2M, zero = rep(0, ncol(xM.mat)))

## fit the Bayesian graphical model using JAGS
jmM <- jags.model("./final_code/application_dependent_score.jags", data = dataM.jags, n.chains = 1, n.adapt = 1000) 

## specify number of burn-in
update(jmM, n.burn = 10000)

## obtain posterior samples
samples <- jags.samples(jmM, c("b"), n.iter = 5000, thin = 5)

## estimate of predictors
thetaM.sample <- samples$b[,,1]
thetaM.hat <- apply(thetaM.sample, 1, mean)



## Create dataset for women 
dataW.nhanes <- data.frame(yW, zW, xW)
colnames(dataW.nhanes) <- c("mort", paste0(colnames(zM), ".w"), var)

## specify the additive model formula -- to obtain spline basis and corresponding penalty matrices from mgcv
formW <- paste0("mort ~ ", paste(colnames(dataW.nhanes)[2:(ncol(dataW.nhanes)-ncol(X))], collapse = " + "))
for(pa in var){
  formW <- paste0(formW, " + s(", pa, ", bs = 'tp', k = 11)")
}
## fit the additive model in mgcv
fitW <- gam(as.formula(formW), data = dataW.nhanes, family = binomial, select = TRUE) 
xW.mat <- predict.gam(fitW, type = "lpmatrix") ## derive linear predictor matrix
S1W <- cbind(fitW$smooth[[1]]$S[[1]], fitW$smooth[[1]]$S[[2]]) ## penalty matrices for MVPA
S2W <- cbind(fitW$smooth[[2]]$S[[1]], fitW$smooth[[2]]$S[[2]]) ## penalty matrices for ASTP

## build the dataset to fit Bayesian model
dataW.jags <- list(y = as.vector(yW), n = length(yW),
                   X = xW.mat, S1 = S1W, S2 = S2W, zero = rep(0, ncol(xW.mat)))

## fit the Bayesian graphical model using JAGS
jmW <- jags.model("./final_code/application_dependent_score.jags", data = dataW.jags, n.chains = 1, n.adapt = 1000) 

## specify number of burn-in
update(jmW, n.burn = 10000)

## obtain posterior samples
samplesW <- jags.samples(jmW, c("b"), n.iter = 5000, thin = 5)

## estimate of predictors
thetaW.sample <- samplesW$b[,,1]
thetaW.hat <- apply(thetaW.sample, 1, mean)


# Obtain log likelihood for dependent model
etaM <- xM.mat %*% thetaM.hat
etaW <- xW.mat %*% thetaW.hat
eta_all <- rbind(etaM, etaW)
loglike_dep <- sum(Y*log(plogis(eta_all)) + (1-Y)*log(1-plogis(eta_all))) ## log likelihood
k_dep <- length(thetaM.hat) + length(thetaW.hat)


#########################################################################
# Calculate AIC and BIC for each model
##########################################################################
## Calculate AIC for each model
# AIC = 2k - 2loglike
AIC_indep <- 2 * k_indep - 2* loglike_indep
AIC_dep <- 2 * k_dep - 2* loglike_dep

## Calculate BIC for each model
# BIC = klog(n) - 2loglike
BIC_indep <- log(n) * k_indep - 2* loglike_indep
BIC_dep <- log(n) * k_dep - 2* loglike_dep

# Print results
results <- cbind(c(AIC_indep, AIC_dep), c(BIC_indep, BIC_dep))
rownames(results) <- c("Score is independent of population", 
                       "Score is dependent of population")
colnames(results) <- c("AIC", "BIC")
print(results)
