## This script fit the model using Bayesian method (JAGS)

rm(list = ls())
set.seed(100)

####################################################################
# Pre-process the NHANES Data
####################################################################

## Derive variables from NHANES to fit single index model
source("./final_code/nhanes_preprocessing.R")
Y <- rbind(yM, yW) ## response
X = rbind(xM, xW) ## PA measures
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
for(pa in var){
  form <- paste0(form, " + s(", pa, ", bs = 'tp', k = 11)")
}
### fit the additive model in mgcv (just for spline basis and penalty matrices purpose)
fit <- gam(as.formula(form), data = data.nhanes, family = binomial, select = TRUE) 
X.mat <- predict.gam(fit, type = "lpmatrix") ## derive linear predictor matrix
S1 <- cbind(fit$smooth[[1]]$S[[1]], fit$smooth[[1]]$S[[2]]) ## penalty matrices for MVPA
S2 <- cbind(fit$smooth[[2]]$S[[1]], fit$smooth[[2]]$S[[2]]) ## penalty matrices for ASTP


## Add epsilon difference to X.mat to get spline basis matrix for first derivative (fd) calculation
eps <- 1e-7
Xnew.mat <- X + eps
newdat <- cbind(Xnew.mat, data.nhanes[, c("X0.m", "X0.w", paste0(colnames(zM), ".m"),
                                          paste0(colnames(zW), ".w"))])
Xeps <- predict(fit, data.frame(newdat), type = 'lpmatrix')


## Build the final dataset to fit Bayesian model
data.jags <- list(y = as.vector(Y), n = length(Y), n1 = length(yM),
                  X = X.mat, S1 = S1, S2 = S2, zero = rep(0, ncol(X.mat)))


####################################################################
# Model Fitting
####################################################################

## Fit the Bayesian graphical model using JAGS
load.module("glm") 
jm <- jags.model("./final_code/application.jags", data = data.jags, n.chains = 1, n.adapt = 1000) 
### specify number of burn-in
update(jm, n.burn = 10000)
### obtain posterior samples
samples <- jags.samples(jm, c("b", "c"), n.iter = 5000, thin = 5)


####################################################################
# Organize the Results
####################################################################

## Estimate of beta_{12}
betaW.sample <- as.vector(samples$c)
betaW.hat <- mean(betaW.sample)
betaW.sd <- sd(betaW.sample)

## Estimate of predictors
theta.sample <- samples$b[,,1]
theta.hat <- apply(theta.sample, 1, mean)
theta.sd <- apply(theta.sample, 1, sd)
theta.lower <- theta.hat - 2*apply(theta.sample, 1, sd)
theta.upper <- theta.hat + 2*apply(theta.sample, 1, sd)
theta.est <- rbind(theta.hat, theta.lower, theta.upper)
colnames(theta.est) <- colnames(X.mat)

## Organize and print results
est <- c(betaW.hat, theta.hat[3:(2+ncol(zM) + ncol(zW))])
se <- c(betaW.sd, theta.sd[3:(2+ncol(zM) + ncol(zW))])
pvals <-  2 * (1 - pnorm(abs(est / se)))
results <- cbind(est, se, pvals)
rownames(results) <- c("beta1.w", colnames(X.mat)[3:(2+ncol(zM) + ncol(zW))])
print(results, digits = 2)


## Get finite difference approx for each sample
Xp <- (Xeps - X.mat) 
FD <- list()
nb <- (ncol(data.jags$X)-2-ncol(zM)-ncol(zW))/length(var)
for(i in 1:length(var)){
  ### Index of spline bases for term
  ind <- 2+ncol(zM)+ncol(zW) + ((i-1)*nb+1):(i*nb)
  
  ### Get finite difference matrix
  Bpg <- (Xp[,ind]) / eps
  FDg <- c()
  for(j in 1:ncol(theta.sample)){
    ### Get fd estimate
    FDg <- cbind(FDg, Bpg %*% matrix(samples$b[ind,j,1], ncol = 1))
  }
  
  FD[[var[i]]] <- FDg
}

## Get fd basis matrix
Cg <- Xp / eps
Cg <- Cg[, (3 + ncol(zM) + ncol(zW)):ncol(X.mat)]


## Derive covariance matrix of spline coefficients from bootstrap
ind <- (3+ncol(zM)+ncol(zW)):ncol(X.mat)
alphaest <- as.matrix(samples$b[ind,,1])
Vb <- cov(t(alphaest))

## Generate 1,000 simulations from MVN
alphahat <- rowMeans(alphaest) # estimate of alpha
betahat <- mvrnorm(n = 1000, mu = alphahat, Sigma = Vb)

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


## Rescale beta coefficients
C <- total / 100
beta.rs <- c(1, results["beta1.w", "est"]) * C
se.rs <- c(NA, results["beta1.w", "se"] * C)
pval.rs <- c(NA, results["beta1.w", "pvals"] )
beta.df <- data.frame(est = beta.rs, se = se.rs, pval = pval.rs)
rownames(beta.df) <- c("betaM", "betaW")
print(beta.df)


## Estimate s_j(X) and first derivative, as well as its confidence band
X_reg <- rbind(as.matrix(dataM[,var_reg], ncol = length(var_reg)),
                   as.matrix(dataW[,var_reg], ncol = length(var_reg)))
curve <- list()
nb <- (ncol(data.jags$X)-2-ncol(zM)-ncol(zW))/length(var) ## number of basis per term
for(i in 1:length(var)){
  ind.spl <- 2+ncol(zM)+ncol(zW) + ((i-1)*nb+1):(i*nb) ## indices of corresponding spline basis in X.spl
  
  ### Obtain estimate and CI for curve
  fHat.all <- data.jags$X[,ind.spl] %*% samples$b[ind.spl,,1]
  fHat.mean <- apply(fHat.all, 1, mean)
  lower <- apply(fHat.all, 1, function(x) quantile(x, 0.025))
  upper <- apply(fHat.all, 1, function(x) quantile(x, 0.975))
  
  ### Rescale to sum to 100
  fHat.mean.rs <- sapply(fHat.mean + abs(min(fHat.mean)), function(x) x * (100 / total))
  lower.rs <- sapply(lower + abs(min(fHat.mean)), function(x) x * (100 / total))
  upper.rs <- sapply(upper + abs(min(fHat.mean)), function(x) x * (100 / total))
  
  ### Obtain estimate and PWCI for fd
  fd.mean <- apply(FD[[var[i]]], 1, mean)
  fd.se <- apply(FD[[var[i]]], 1, sd)
  fd.lower <- fd.mean + qnorm(0.025) * fd.se
  fd.upper <- fd.mean + qnorm(0.975) * fd.se
  
  ### Calculate simultaneous CI
  ## following blog post by Gavin Simpson
  ## Isolate basis components for current PA variable
  ind <- ((i-1)*nb + 1):(i*nb)
  Ci <- Cg * 0
  Ci[, ind] <- Cg[, ind] # only nonzero for bases of current PA variable
  
  ### Calculate fhat(x) - f(x)
  simDev <- fd.mean - Ci %*% t(betahat) 
  ### Find the absolute values of the standardized deviations from the true mode
  absDev <- abs(sweep(simDev, 1, fd.se, FUN = "/"))
  ### Find maximum of the absolute standardized deviations at the grid of x values for each simulation 
  maxsd <- apply(absDev, 2L, max)
  ### Find the critical value used to scale the standard errors to yield the simultaneous interval
  crit <- quantile(maxsd, prob = 0.95, type = 8)
  print(crit)
  
  ### Get simultaneous CI
  fd.lower.s <- fd.mean - (crit * fd.se)
  fd.upper.s <- fd.mean + (crit * fd.se)
  
  ### Smooth estimates and save in list
  x <- X[,i]
  x.reg <- X_reg[,i]
  curve[[var[i]]] <- data.frame(x = x, x.reg = x.reg, s.x = fHat.mean, s.x.rs = fHat.mean.rs)
  curve[[var[i]]]$lower <- gam(lower ~ s(x), method = "REML")$fitted.values
  curve[[var[i]]]$upper <- gam(upper ~ s(x), method = "REML")$fitted.values
  curve[[var[i]]]$lower.rs <- gam(lower.rs ~ s(x.reg), method = "REML")$fitted.values
  curve[[var[i]]]$upper.rs <- gam(upper.rs ~ s(x.reg), method = "REML")$fitted.values
  curve[[var[i]]]$fd <- gam(fd.mean ~ s(x.reg), method = "REML")$fitted.values # smooth fd estimate
  curve[[var[i]]]$fd.lw <- gam(fd.lower ~ s(x.reg), method = "REML")$fitted.values # smooth fd PWCI lower bound
  curve[[var[i]]]$fd.up <- gam(fd.upper ~ s(x.reg), method = "REML")$fitted.values # smooth fd PWCI upper bound
  curve[[var[i]]]$fd.lw.s <- gam(fd.lower.s ~ s(x.reg), method = "REML")$fitted.values # smooth fd SCI lower bound
  curve[[var[i]]]$fd.up.s <- gam(fd.upper.s ~ s(x.reg), method = "REML")$fitted.values # smooth fd SCI upper bound
}

####################################################################
# Calculate Cross-validated AUC
####################################################################

score_reorder <- rep(NA, length(total.scores.adj))
score_reorder[which(nhanes_data$Gender == "Male")] <- total.scores.adj[1:nrow(dataM)]
score_reorder[which(nhanes_data$Gender == "Female")] <- total.scores.adj[(nrow(dataM)+1):length(total.scores.adj)]
nhanes_data$score <- score_reorder
source("./final_code/application_AUC.R")
print(auc.result)

####################################################################
# Plot the Additive Terms and Their First Derivatives
####################################################################

## Plot the estimate rescaled to sum to 100
glist <- list()
for(pa in var){
  if(pa == "MVPA"){
    xlab <- "Min/day"
  }
  else{
    xlab <- "Transition probability"
  }
  curve.order <- curve[[pa]][order(curve[[pa]][,1]),]
  curve.order <- distinct(curve.order) ## keep only rows with distinct values
  glist[[pa]] <- ggplot(curve.order, aes(x = x.reg, y = s.x.rs)) +
    geom_line(col = "blue") +
    geom_line(aes(x = x.reg, y = lower.rs), col = "blue", linetype = "dashed") +
    geom_line(aes(x = x.reg, y = upper.rs), col = "blue", linetype = "dashed") +
    geom_rug(sides = "b") +
    theme_bw() +
    labs(x = xlab, y = "Score", title = pa, subtitle = "Bayesian Method") +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
          text = element_text(size = 15), plot.subtitle = element_text(hjust = 0.5))
}
glist[["MVPA"]] <- glist[["MVPA"]] + ylim(-10, 62)
glist[["ASTP"]] <- glist[["ASTP"]] + ylim(-30, 75)
pdf(file = "./final_code/figures/Estimates_JAGS.pdf", width = 10, height = 4.5)
do.call("grid.arrange", c(glist, nrow = floor(sqrt(length(glist)))))
dev.off()


## Plot the first derivative
fdlist <- list()
for(pa in var){
  if(pa == "MVPA"){
    xlab <- "Min/day"
  }
  else{
    xlab <- "Transition probability"
  }
  curve.order <- curve[[pa]][order(curve[[pa]][,1]),]
  curve.order <- distinct(curve.order) ## keep only rows with distinct values
  fdlist[[pa]] <- ggplot(curve.order, aes(x = x.reg, y = fd)) +
    geom_line(aes(x = x.reg, y = fd.lw), col = "blue", linetype = "dashed") +
    geom_line(aes(x = x.reg, y = fd.up), col = "blue", linetype = "dashed") +
    geom_line(col = "blue") + geom_hline(yintercept = 0) + 
    geom_rug(sides = "b") + theme_bw() + 
    labs(x = xlab, y = "First derivative", title = paste0("First Derivative of ",pa), 
         subtitle = "Bayesian Method") +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
          text = element_text(size = 15), plot.subtitle = element_text(hjust = 0.5))
}
fdlist[["MVPA"]] <- fdlist[["MVPA"]] + ylim(-1, 1.9)
fdlist[["ASTP"]] <- fdlist[["ASTP"]] + ylim(-1.6, 1.3)
pdf(file = "./final_code/figures/First_Deriv_JAGS.pdf", width = 10, height = 4.5)
do.call("grid.arrange", c(fdlist, nrow = floor(sqrt(length(fdlist)))))
dev.off()

glist_jags <- glist
fdlist_jags <- fdlist
save(glist_jags, fdlist_jags, file = "./final_code/results/est_deriv_plotlist_jags.rda")


