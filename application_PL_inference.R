###############################################################################
# This script fits the estimated curve and first derivative of the model using 
# the PL estimation procedure
###############################################################################
rm(list = ls())
set.seed(4428967)

# Source pre-processed NHANES data
source("./final_code/nhanes_preprocessing.R")

# Load bootstrap data from rda files with 300 samples
load(file = "./final_code/results/application_PL_results.rda")

##############################################################################
# Rescale total score to sum between 0 and 100
##############################################################################
X.all = rbind(xM, xW) # scaled PA variables
X.all_reg <- rbind(as.matrix(dataM[,var_reg], ncol = length(var_reg)),
                   as.matrix(dataW[,var_reg], ncol = length(var_reg))) # original PA variables

# Get total score 
scores.raw <- c()
for(pa in var){
  # Obtain s.x and combine in dataframe
  score.pa <- colMeans(boot[[pa]]$curvemat)
  data.curve <- data.frame(cbind(X.all[,pa][order(X.all[,pa])], score.pa))
  # Current order of s.x
  X.order <- as.integer(rownames(data.curve))
  # Original order of x
  X.orig <- as.integer(rownames(X.all))   
  # Reorganize s.x in original order of x
  score.pa.ordered <- data.curve[match(X.orig, X.order), 2]   
  
  scores.raw <- cbind(scores.raw, score.pa.ordered)
}
colnames(scores.raw) <- var
total.score.raw <- rowSums(scores.raw)

# Obtain rescaled total score
scores.adj <- apply(scores.raw, 2, function(x) x - min(x))
total <- sum(apply(scores.adj, 2, max))
total.scores.adj <- (100/total) * rowSums(scores.adj)
range(total.scores.adj)


############################################################################
# Organize nonparametric bootstrap results for parametric coefficients
############################################################################
# Get estimates and 95% CI of linear parameters
est = colMeans(boot$paramMat)
bands = apply(boot$paramMat, 2, function(x) quantile(x, c(0.025, 0.975)))
estmat <- data.frame(est = est, lower = bands[1,], upper = bands[2,])
rownames(estmat) <- c("int.m", "int.w","beta_12", paste0(colnames(zM), ".m"),  paste0(colnames(zW), ".w"))
se = (estmat$upper - estmat$lower)/ (2 * qnorm(0.975)) # compute SE
pvals = 2 * (1 - pnorm(abs(est / se))) # compute p-values for coefficients != 0

# Print the estimate results
estmat = data.frame(list(est = est, se = se, lower = bands[1, ], upper = bands[2, ], pval = pvals))
rownames(estmat) <- c("int.m", "int.w","beta_12", paste0(colnames(zM), ".m"),  paste0(colnames(zW), ".w"))
print(round(estmat, 2))

# Rescale beta coefficients
C <- (max(total.score.raw) - min(total.score.raw)) / 100
beta.rs <- c(1, estmat["beta_12", "est"]) * C
se.rs <- c(NA, estmat["beta_12", "se"] * C)
pval.rs <- c(NA, estmat["beta_12", "pval"] )
beta.df <- data.frame(est = beta.rs, se = se.rs, pval = pval.rs)
rownames(beta.df) <- c("betaM", "betaW")
print(beta.df)


#############################################################################
# Calculate estimate of first derivative bases and simulation for
# simultaneous CI (SCI)
#############################################################################
# Get estimated spline bases matrix
nknots <- 8
Cg <- c()
for(i in 1:(2*(nknots+3))){
  Cg <- cbind(Cg, apply(boot[[paste0("basis",i)]]$basismat, 2, mean))
}

# Derive covariance matrix of spline coefficients from bootstrap
Vb <- cov(boot$alphaMat)

# Get 1000 MVN samples from N(alpha, Var(alpha))
alphahat <- colMeans(boot$alphaMat) # estimate of alpha
betahat <- mvrnorm(n = 1000, mu = alphahat, Sigma = Vb)


##############################################################################
# Get mean and CI of bootstrap estimate and first derivative samples for each 
# PA component
##############################################################################
nb <- ncol(boot$alphaMat) / ncol(X.all) # number of basis functions
i = 1
est.bs = list() # list to hold estimates for curve

for(pa in var){
  # Isolate the PA component
  est.var <- boot[[pa]]
  
  # Identify column of X values
  X <- X.all[, pa]
  x <- X[order(X)]
  X_reg <- X.all_reg[, var_reg[i]]
  x.reg <- X_reg[order(X_reg)]
  
  #############################################################
  # Bootstrap estimate and CI
  #############################################################
  # Get estimates for smooth function of estimate and point-wise confidence bands
  curve.est <- colMeans(est.var$curvemat) # mean of curve
  lower.est <- apply(est.var$curvemat, 2, function(x) quantile(x, 0.025)) # lower CI estimate of curve
  upper.est <- apply(est.var$curvemat, 2, function(x) quantile(x, 0.975)) # upper CI estimate of curve
  
  # Rescale to sumfrom 0 to 100
  curve.rs <- sapply(curve.est + abs(min(curve.est)), function(x) x * (100 / total))
  lower.rs <- sapply(lower.est + abs(min(curve.est)), function(x) x * (100 / total))
  upper.rs <- sapply(upper.est + abs(min(curve.est)), function(x) x * (100 / total))
  
  #############################################################
  # Bootstrap first derivative and CI
  #############################################################
  # Get estimates for smooth function of derivative with PWCI and SCI
  # Isolate basis components for current PA variable
  ind <- ((i-1)*nb + 1):(i*nb)
  Ci <- Cg * 0
  Ci[, ind] <- Cg[, ind] # only nonzero for bases of current PA variable
  i <- i + 1
  
  # Derive first derivative estimate and se
  fd.est <- colMeans(est.var$fdmat) #estimate
  se.fit <- apply(est.var$fdmat, 2, sd) # sd along curve from bootstrap estimates
  
  # Point-wise confidence intervals 
  fd.lower.pw <- fd.est + qnorm(0.025) * se.fit
  fd.upper.pw <- fd.est + qnorm(0.975) * se.fit
  
  # Calculate simultaneous CI
  ## following blog post by Gavin Simpson
  # Calculate fhat(x) - f(x)
  simDev <- fd.est - Ci %*% t(betahat) 
  # Find the absolute values of the standardized deviations from the true model
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  # Find maximum of the absolute standardized deviations at the grid of x values for each simulation 
  maxsd <- apply(absDev, 2L, max)
  # Find the critical value used to scale the standard errors to yield the simultaneous interval
  crit <- quantile(maxsd, prob = 0.95, type = 8)
  
  # Get simultaneous CI
  fd.lower.s <- fd.est - (crit * se.fit)
  fd.upper.s <- fd.est + (crit * se.fit)
  
  
  # Smooth estimates and save in list
  est.bs[[pa]] <- data.frame(x = x, x.reg = x.reg, s.x = curve.est)
  est.bs[[pa]]$lower <- gam(lower.est ~ s(x), method = "REML")$fitted.values
  est.bs[[pa]]$upper <- gam(upper.est ~ s(x), method = "REML")$fitted.values
  est.bs[[pa]]$s.x.rs <- gam(curve.rs ~ s(x), method = "REML")$fitted.values
  est.bs[[pa]]$lower.rs <- gam(lower.rs ~ s(x), method = "REML")$fitted.values
  est.bs[[pa]]$upper.rs <- gam(upper.rs ~ s(x), method = "REML")$fitted.values
  est.bs[[pa]]$fd <- gam(fd.est ~ s(x), method = "REML")$fitted.values 
  est.bs[[pa]]$fd.lw <- gam(fd.lower.pw ~ s(x), method = "REML")$fitted.values # smooth fd PWCI lower bound
  est.bs[[pa]]$fd.up <- gam(fd.upper.pw ~ s(x), method = "REML")$fitted.values # smooth fd PWCI upper bound
  est.bs[[pa]]$fd.lw.s <- gam(fd.lower.s ~ s(x), method = "REML")$fitted.values # smooth fd SCI lower bound
  est.bs[[pa]]$fd.up.s <- gam(fd.upper.s ~ s(x), method = "REML")$fitted.values # smooth fd SCI upper bound
  
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
# Plot the bootstrap estimates for the curves and first derivatives
####################################################################

# Plot the bootstrap estimates with CI rescaled to sum to 100
glist <- list()
for(pa in var){
  if(pa == "MVPA"){
    xlab <- "Min/day"
  }
  else{
    xlab <- "Transition probability"
  }
  curve.order <- est.bs[[pa]]
  curve.order <- distinct(curve.order) ## keep only rows with distinct values
  glist[[pa]] <- ggplot(curve.order, aes(x = x.reg, y = s.x.rs)) +
    geom_line(col = "blue") +
    geom_line(aes(x = x.reg, y = lower.rs), col = "blue", linetype = "dashed") +
    geom_line(aes(x = x.reg, y = upper.rs), col = "blue", linetype = "dashed") +
    geom_rug(sides = "b") +
    theme_bw() +
    labs(x = xlab, y = "Score", title = pa, subtitle = "Profile Likelihood") +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
          text = element_text(size = 15), plot.subtitle = element_text(hjust = 0.5))
}
glist[["MVPA"]] <- glist[["MVPA"]] + ylim(-10, 62)
glist[["ASTP"]] <- glist[["ASTP"]] + ylim(-30, 75)
pdf(file = "./final_code/figures/Estimates_PL.pdf", width = 10, height = 4.5)
do.call("grid.arrange", c(glist, nrow = floor(sqrt(length(glist)))))
dev.off()


# Plot the bootstrap derivatives with CI
fdlist <- list()
for(pa in var){
  if(pa == "MVPA"){
    xlab <- "Min/day"
  }
  else{
    xlab <- "Transition probability"
  }
  curve.order <- est.bs[[pa]]
  curve.order <- distinct(curve.order) ## keep only rows with distinct values
  fdlist[[pa]] <- ggplot(curve.order, aes(x = x.reg, y = fd)) +
    #geom_ribbon(aes(ymin = fd.lw.s, ymax = fd.up.s), fill = "darkgray") +
    #geom_ribbon(aes(ymin = fd.lw, ymax = fd.up), fill = "lightgray") +
    geom_line(aes(x = x.reg, y = fd.lw), col = "blue", linetype = "dashed") +
    geom_line(aes(x = x.reg, y = fd.up), col = "blue", linetype = "dashed") +
    geom_line(col = "blue") + geom_hline(yintercept = 0) + 
    geom_rug(sides = "b") + theme_bw() + 
    labs(x = xlab, y = "First derivative", title = paste0("First Derivative of ",pa), 
         subtitle = "Profile Likelihood") +
    theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
          text = element_text(size = 15), plot.subtitle = element_text(hjust = 0.5))
}
fdlist[["MVPA"]] <- fdlist[["MVPA"]] + ylim(-1, 1.9)
fdlist[["ASTP"]] <- fdlist[["ASTP"]] + ylim(-1.6, 1.3)
pdf(file = "./final_code/figures/First_Deriv_PL.pdf", width = 10, height = 4.5)
do.call("grid.arrange", c(fdlist, nrow = floor(sqrt(length(fdlist)))))
dev.off()

glist_PL <- glist
fdlist_PL <- fdlist
save(glist_PL, fdlist_PL, file = "./final_code/results/est_deriv_plotlist_PL.rda")
