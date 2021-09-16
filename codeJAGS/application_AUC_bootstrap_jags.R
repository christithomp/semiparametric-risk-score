## This script obtain bootstrap samples of the change in AUC and perform a t-test on these changes

rm(list = ls())

## Derive variables from NHANES to fit single index model
source("./final_code/nhanes_preprocessing.R")

####################################################################
# Bootstrap
####################################################################

xM0 <- xM
xW0 <- xW
yM0 <- yM
yW0 <- yW
zM0 <- zM
zW0 <- zW
nhanes0 <- nhanes_data

nboot <- 100 ## number of bootstrap iterations
auc_delta <- rep(0, nboot)
auc_boot <- matrix(0, nrow = nboot, ncol = 5)

for(b in 1:nboot){
  set.seed(b)
  indM <- sample(1:nrow(xM), nrow(xM), replace = TRUE)
  indW <- sample(1:nrow(xW), nrow(xW), replace = TRUE)
  
  xM <- xM0[indM,]
  yM <- yM0[indM,,drop = F]
  zM <- zM0[indM,]
  xW <- xW0[indW,]
  yW <- yW0[indW,,drop = F]
  zW <- zW0[indW,]
  
  nhanes_data[which(nhanes0$Gender == "Male"),] <- nhanes0[which(nhanes0$Gender == "Male")[indM],]
  nhanes_data[which(nhanes0$Gender == "Female"),] <- nhanes0[which(nhanes0$Gender == "Female")[indW],]
  
  Y <- rbind(yM, yW) ## response
  X <- rbind(xM, xW) ## PA measures
  X01 <- matrix(c(rep(1, length(yM)), rep(0, length(yW))), ncol = 1) ## design matrix of intercept term
  X02 <- matrix(c(rep(0, length(yM)), rep(1, length(yW))), ncol = 1)
  Z1 <- rbind(zM, matrix(0, nrow = length(yW), ncol = ncol(zM))) ## design matrix of other predictors
  Z2 <- rbind(matrix(0, nrow = length(yM), ncol = ncol(zW)), zW)
  
  ## Build a dataset containing variables of interest
  data.nhanes <- data.frame(Y, X01, X02, Z1, Z2, X)
  colnames(data.nhanes) <- c("mort", "X0.m", "X0.w", paste0(colnames(zM), ".m"),
                             paste0(colnames(zW), ".w"), var)
  
  ## Specify the additive model formula to obtain spline basis and corresponding penalty matrices from mgcv
  form <- paste0("mort ~ -1 + ", paste(colnames(data.nhanes)[2:(ncol(data.nhanes)-ncol(X))], collapse = " + "))
  k <- 11 # number of interior knots for each smooth function
  for(pa in var){
    form <- paste0(form, " + s(", pa, ", bs = 'tp', k = ", k, ")")
  }
  
  ### fit the additive model in mgcv (just for spline basis and penalty matrices purpose)
  fit <- gam(as.formula(form), data = data.nhanes, family = binomial, select = TRUE) 
  X.mat <- predict.gam(fit, type = "lpmatrix") ## derive linear predictor matrix
  S1 <- cbind(fit$smooth[[1]]$S[[1]], fit$smooth[[1]]$S[[2]]) ## penalty matrices for MVPA
  S2 <- cbind(fit$smooth[[2]]$S[[1]], fit$smooth[[2]]$S[[2]]) ## penalty matrices for ASTP
  
  ## Build the final dataset to fit Bayesian model
  data.jags <- list(y = as.vector(Y), n = length(Y), n1 = length(yM), 
                    k = k, nc = (ncol(Z1) + ncol(Z2) + 2),
                    X = X.mat, S1 = S1, S2 = S2, zero = rep(0, ncol(X.mat)))
  
  ####################################################################
  # Model Fitting 
  ####################################################################
  
  ## Fit the Bayesian graphical model using JAGS
  load.module("glm") 
  jm <- jags.model("./final_code/application.jags", data = data.jags, 
                   n.chains = 1, n.adapt = 1000) 
  ### specify number of burn-in
  update(jm, n.burn = 10000)
  ### obtain posterior samples
  samples <- jags.samples(jm, c("b", "c"), n.iter = 2000, thin = 2)
  
  ####################################################################
  # Organize the Results
  ####################################################################

  ## Get total score for rescaling
  scores.raw <- c()
  i <- 1
  nb <- (ncol(data.jags$X)-2-ncol(zM)-ncol(zW))/length(var) ## number of basis per term
  for(i in 1:length(var)){
    ### Get index of pa variable spline bases
    ind.spl <- 2+ncol(zM)+ncol(zW) + ((i-1)*nb+1):(i*nb) 
    i <- i + 1
    ### Obtain estimate for curve
    fHat <- data.jags$X[,ind.spl] %*% samples$b[ind.spl,,1]
    
    ### Store score for pa variable
    scores.raw <- cbind(scores.raw, apply(fHat, 1, mean))
  }
  
  ## Obtain rescaled scores
  scores.adj <- apply(scores.raw, 2, function(x) x - min(x))
  total <- sum(apply(scores.adj, 2, max))
  total.scores.adj <- (100/total)*rowSums(scores.adj)
  range(total.scores.adj)
  
  ####################################################################
  # Calculate Cross-validated AUC
  ####################################################################
  
  score_reorder <- rep(NA, length(total.scores.adj))
  score_reorder[which(nhanes_data$Gender == "Male")] <- total.scores.adj[1:nrow(dataM)]
  score_reorder[which(nhanes_data$Gender == "Female")] <- total.scores.adj[(nrow(dataM)+1):length(total.scores.adj)]
  nhanes_data$score <- score_reorder
  source("./final_code/application_AUC.R")
  # print(auc.result)
  auc_delta[b] <- auc.result$Cross.Validated.AUC[1] - auc.result$Cross.Validated.AUC[5]
  auc_boot[b,] <- auc.result$Cross.Validated.AUC
  
  print(b)
}

print(auc_boot)
pvals <- 2 * (1 - pnorm(abs(mean(auc_delta) / sd(auc_delta) ) ) )
print(pvals)

save(auc_boot, auc_delta, pvals, file = "./final_code/bootstrap_AUC.rda")




