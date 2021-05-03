############################################################
# This code performs nonparametric boostrap to fit the model
# using the profile likelihood (PL) method.
############################################################
rm(list = ls())
set.seed(4428967)

## Source code to pre-process the NHANES data
source("./final_code/nhanes_preprocessing.R")
## Source code to update parameter estimates for PL method
source("./final_code/application_PL_update_parameters.R")

############################################################
# Setting up the bootstrap sampler
############################################################
## Specify knots and type of spline for score components
spline <- "os" # O'Sullivan splines
nknots <- 8 # number of interior knots

iter = 300 #number of bootstrap samples

boot <- list() # list to store bootstrap samples 
boot$paramMat <- matrix(0, nrow = iter, ncol = (3 + ncol(zM) + ncol(zW))) #store betas and thetas
boot$alphaMat <- matrix(0, nrow = iter, ncol = ncol(xM)*(nknots + 2) + 2) #store alpha
for(pa in var){
  boot[[pa]]$fdmat <- matrix(0, nrow = iter, ncol = (nrow(xM) + nrow(xW))) #store first derivative
  boot[[pa]]$curvemat <- matrix(0, nrow = iter, ncol = (nrow(xM) + nrow(xW))) #store curve
}
# store bootstrap samples of each fd spline basis
for(i in 1:(2*(nknots+3))){
  boot[[paste0("basis",i)]]$basismat <- matrix(0, nrow = iter, ncol = (nrow(xM) + nrow(xW)))
}

# Get number of men and women for sampling
nrM <- nrow(yM)
nrW <- nrow(yW)

# Copies of original matrices 
yMO <- yM; yWO <- yW
xMO <- xM; xWO <- xW
zMO <- zM; zWO <- zW


############################################################
# Nonparametric bootstrap sampler
############################################################

for(t in 1:iter){
  
  # Obtain indices of bootstrap samples
  nM <- sample(1:nrM, nrM, replace = TRUE)
  nW <- sample(1:nrW, nrW, replace = TRUE)
  
  # Create new data with indices
  yM <- yMO[nM, ]; yW <- yWO[nW, ]
  xM <- xMO[nM, ]; xW <- xWO[nW, ]
  zM <- zMO[nM, ]; zW <- zWO[nW, ]
  
  ## Organize iteration-invariant variables to fit model
  betaM <- -1 # identifiability constraint
  betaW <- 0.5 # starting value of betaW
  Y <- matrix(c(yM, yW), ncol = 1) # response
  X.all <- rbind(xM, xW) ##PA variables
  X01 <- matrix(c(rep(1, length(yM)), rep(0, length(yW))), ncol = 1) # design matrix of intercept term for men
  X02 <- matrix(c(rep(0, length(yM)), rep(1, length(yW))), ncol = 1) # design matrix of intercept term for women
  Z1 <- rbind(zM, matrix(0, nrow = length(yW), ncol = ncol(zM))) # design matrix of covariates for men
  Z2 <- rbind(matrix(0, nrow = length(yM), ncol = ncol(zW)), zW) # design matrix of covariates for women
  dummyID <- factor(rep(1, length(Y))) # full index set for GLMM
  
  # Get epsilon shift in X.all for first derivative calculation
  eps <- 1e-7
  X.alleps <- X.all + eps
  
  ######################################################################
  # Get spline basis matrix for PA variables
  ######################################################################
  B <- c() # spline bases
  Beps <- c() # epsilon shifted spline bases
  for(i in 1:ncol(X.all)){
    # Obtain ith PA variable 
    x <- matrix(X.all[,i], ncol = 1)
    Xeps <- matrix(X.alleps[,i], ncol = 1) #epsilon shift
    
    if(spline == "os"){
      # O'Sullivan splines
      # Set range of X values and interior knots
      a <- 1.01*min(x) - 0.01*max(x)
      b <- 1.01*max(x) - 0.01*min(x)
      intKnots <- quantile(unique(x), seq(0, 1, length = (nknots + 2))[-c(1,(nknots + 2))])
      
      # Get O'Sullivan spline basis functions
      Bg <- ZOSull(x, range.x = c(a,b), intKnots = intKnots)
      B <- cbind(B, Bg)
      
      # Get O'Sullivan spline basis functions for epsilon shift
      Bepsg <- ZOSull(Xeps, range.x = c(a,b), intKnots = intKnots)
      Beps <- cbind(Beps, Bepsg)
    }else{
      # Spline from mgcv
      fit <- gam(Y ~ -1 + s(x, k = nknots, fx = TRUE, bs = spline), family = "gaussian") ## fit a fake model
      B <- cbind(B, predict(fit, type = "lpmatrix"))
      Beps <- cbind(Beps, predict(fit, data.frame(x = Xeps), type = "lpmatrix"))
      rm(fit)
    }
  }
  rm(x)
  
  #########################################################################
  # Fit model using iterative PL algorithm Steps 1-3
  #########################################################################
  threshold <- 1e-6 # change in log-likelihood
  loglike_old <- 1 # initial value of log likelihood
  nc <- 0 # number of max iterations
  epsilon <- 1 #initial change in log-likelihood
  
  # Start iterative algorithm
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
    sigma_a = updated_estimates$sigma_a
    
    # Update log likelihood
    loglike <- sum(Y*log(plogis(eta)) + (1-Y)*log(1-plogis(eta))) # log likelihood
    epsilon <- abs((loglike-loglike_old)/loglike_old) # change in log-likelihoog
    loglike_old <- loglike
    nc <- nc + 1
    print(paste0("change in log likelihood: ", epsilon))
  }
  print("Convergence complete!")
  
  ########################################################################
  # Collect results from PL algorithm 
  ########################################################################
  # Save parametric coefficients for bootstrap sample t
  boot$paramMat[t, ] <- c(beta01, beta02, betaW, thetaM, thetaW)
  boot$alphaMat[t, ] <- c(alpha_int[1], alpha[1:(nknots + 2)], alpha_int[2], alpha[(nknots + 3):length(alpha)])
  
  
  ########################################################################
  # Calculate curve and first derivative estimated
  #########################################################################
  # Finite difference approximation of first derivative design matrix
  Bp <- (Beps - B) / eps
  Xp <- (X.alleps - X.all) / eps
  
  # Save first derivative estimate bases
  bases <- cbind(Xp[,1], Bp[,1:(nknots + 2)], Xp[,2], Bp[,(nknots + 3):length(alpha)])
  for(j in 1:(2*(nknots+3))){
    boot[[paste0("basis",j)]]$basismat[t, ] <- bases[,j]
  }
  
  # Get curve and first derivative estimate
  nb <- ncol(Bp) / ncol(X.all) # number of basis functions
  for(i in 1:ncol(X.all)){
    # Indices of corresponding basis functions for variable i
    ind <- ((i-1)*nb + 1):(i*nb)
    
    # Isolate the parts of X and B corresponding to the component
    Bi <- Bp * 0
    Bi[, ind] <- Bp[, ind]
    Xi <- Xp * 0
    Xi[, i] <- Xp[, i]
    
    # Compute first derivative approximation
    fd <- Bi %*% alpha + Xi %*% alpha_int
    X <- X.all[, i]
    
    # Save results
    boot[[var[i]]]$fdmat[t, ] <- fd[order(X)] # first derivative estimate
    boot[[var[i]]]$curvemat[t, ] <- curve[[var[i]]]$s.x[order(X)] # curve estimate
  }
  
  # Print number of iteration ever 10 iterations
  cat(paste0("Iteration: ", t, "\n"))
  
}

# Save bootstrap results to rda file
save(boot, file = "./results/application_PL_bootstrap_results_updated.rda")
