
#' Estimation of parameters for additive spline model
#'
#' This function will update all parameters needed for the profile likelihood (PL) algorithm.
#'
#' @return A list of the score (score), the weights (alphahat), eta (likelihood), coefficients on the covariates (thetaMhat and thetaWhat), 
#' coefficients of the intercept (beta01hat and beta02hat), and the coefficients on the PA variables (betaMhat and betaWhat).
#' @export
#'
#' @examples

update_parameters <- function(){
  
  #############################################################
  # First step: Estimate all other parameters except beta_12
  #############################################################
  ## Fit a joint GLMM 
  suppressMessages( fit.sp <- gamm(Y ~ -1 + X01 + X02 + Z1 + Z2 + X.all.new, random = list(dummyID = pdIdent(~-1+B.new)), family = "binomial") )
  
  # Obtain coefficient estimates
  alpha.spline <- matrix(unname(unlist(fit.sp$lme$coefficients$random)), ncol = 1)
  beta01 <- unname(fit.sp$gam$coefficients)[1]
  beta02 <- unname(fit.sp$gam$coefficients)[2]
  thetaM <- matrix(unname(fit.sp$gam$coefficients)[3:(2+ncol(zM))], ncol = 1)
  thetaW <- matrix(unname(fit.sp$gam$coefficients)[(3+ncol(zM)):(2+ncol(zM)+ncol(zW))], ncol = 1)
  betaHat <- matrix(unname(fit.sp$gam$coefficients)[(3+ncol(zM)+ncol(zW)):length(fit.sp$gam$coefficients)], ncol = 1)
  sigma_alpha = VarCorr(fit.sp$lme)[2]
  
  ## Derive the score estimates
  score <- B %*% alpha.spline + X.all %*% betaHat
  scoreM <- score[1:nrow(xM),] # score for men
  scoreW <- score[(nrow(xM)+1):nrow(X.all),] # score for women
  
  # Obtain curve of each PA component
  curve <- list() 
  for(j in 1:ncol(xM)){
    # Inded for each PA variable
    ind <- (ncol(B)/ncol(xM) * (j-1) + 1):(ncol(B)/ncol(xM) * j)
    # Dataframe of x and s(x)
    curve[[j]] <- data.frame(x = X.all[ ,j], 
                             s.x = B[ ,ind]%*%alpha.spline[ind, 1] + X.all[ ,j] * betaHat[j])
  }
  names(curve) <- var
  
  ########################################################################
  # Second step: Estimating beta_12 with all other parameters fixed.
  ########################################################################
  ## Fit a glm for each population minus K = 1
  fit.beta <- glm(yW ~ -1 + scoreW, offset = zW %*% thetaW + beta02, family = "binomial")
  # Obtain coefficient estimate
  betaW <- unname(fit.beta$coefficients)
  
  
  ########################################################################
  # Third step: Calculate likelihood
  ########################################################################
  eta <- X01 %*% beta01 + X02 %*% beta02 + Z1 %*% thetaM + Z2 %*% thetaW + 
    rbind(betaM*matrix(scoreM, ncol = 1), betaW*matrix(scoreW, ncol = 1))
  
  # Return estimates
  theresult <- list(score = c(score), beta01hat = beta01, thetaMhat = c(thetaM), 
                    beta02hat = beta02, thetaWhat = c(thetaW), betaMhat = betaM, betaWhat = betaW,
                    eta = eta, alphahat = c(alpha.spline), curve = curve, alphahat_int = betaHat, sigma_a = sigma_alpha)
  return(theresult)
}