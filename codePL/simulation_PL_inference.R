##########################################################################
# This script organizes simulations using the PL method.
##########################################################################
rm(list = ls())
set.seed(100)

library(ggplot2)
library(gridExtra)

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

### Scaled accelerometry variables
varX <- c("MVPA", "ASTP")
xM <- as.matrix(dataM[,varX], ncol = length(var))
xW <- as.matrix(dataW[,varX], ncol = length(var))
colnames(xM) <- colnames(xW) <- varX

## remove unnecessary variables
df_real_m <- data.frame(yM, xM, zM, xM_reg)
df_real_w <- data.frame(yW, xW, zW, xW_reg)
X.all_reg <- rbind(xM_reg, xW_reg)
rm(dep_vars, data_analysis, nhanes_data, dataM, dataW)
rm(xM, yM, zM, xW, yW, zW)

#####################################################################
# Specify true values of parameters
#####################################################################
sMVPA <- function(x){ -(-0.2*(x-0.8)^3-0.4) }
sASTP <- function(x){ -(0.3*exp(x)-1.25) }
beta1.m <- -1 ## identifiability constraint
set.seed(100)
beta1.w <- -1 * runif(1, 0.5, 1.5)
theta.age.m <- runif(1, 7, 10)
theta.age.w <- runif(1, 7, 10)
theta.smoking.m <- runif(1, 0.5, 0.8)
theta.smoking.w <- runif(1, 0.5, 0.8)
beta0.w <- runif(1, -8, -5)
beta0.m <- runif(1, -8, -5)
true.param = c(beta0.m, beta0.w, beta1.w, theta.age.m, theta.smoking.m, 
               theta.age.w, theta.smoking.w)

# Load simulation results
nsim <- 200
load(file = "./final_code/results/simulation_PL_results_cleaned.rda")


##########################################################################
# Calculate coverages of linear parameters
##########################################################################
# Obtain estimates of each parameter
lower.sim <- est.sim - qnorm(0.975) * param.sd.sim
upper.sim <- est.sim + qnorm(0.975) * param.sd.sim

est.param <- apply(est.sim, 1, mean)
est.se <- apply(est.sim, 1, sd) #/ sqrt(200)
cbind(est.param - 1.96*est.se, est.param + 1.96*est.se)

## coverage of each parameter 
coverage <- rep(0, 7)
for(i in 1:7){
  coverage[i] <- length(which(upper.sim[i,] >= true.param[i] & lower.sim[i,] <= true.param[i]))/nsim
}
param.res <- cbind(true.param, est.param, coverage)
rownames(param.res) <- c("beta0.m", "beta0.w", "beta1.w", "age.m", "smoking.m", "age.w", "smoking.w")
print(param.res, digits = 2)


############################################################################
# Calculate coverage of curve estimates
############################################################################
grid_length <- 100 # evaluate coverage over 100 grid points

# Store coverage estimates
coverage.curve <- list()
coverage_curve_result <- matrix(0, nrow = grid_length, ncol = length(varX))
colnames(coverage_curve_result) <- varX

# Calculate coverage for each PA variable
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

# Obtain average coverage for each curve
cov <- apply(coverage_curve_result, 2, mean)
cov

##########################################################################
# Visualize the simulation results for each curve
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
         subtitle = "Profile Likelihood") +
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

# Plot simulation results for curve estimates with true curves
pdf(file = "./final_code/figures/simulation_PL_MVPA_ASTP_spaghetti_rescaled.pdf", width = 10, height = 5)
do.call("grid.arrange", c(simplot, nrow = 1))
dev.off()

simplot_PL <- simplot
save(simplot_PL, file = "./final_code/results/simplot_PL_list.rda")
