## The code calculates cross-validated survey-weighted AUC for selected predictors in the NHANES dataset

## function which will get cut-points to use in AUC calculation
## returns all unique predicted values of the model
get_ctpts <- function(model, type='link',...) c(Inf,sort(unique(predict(model,type=type,...)),decreasing=TRUE))
## calculate weighted AUC using the appropriate weights
calc_weighted_AUC <- function(response, labels, cutpts, weights=NULL){
  stopifnot(length(response) == length(labels))
  stopifnot(all(labels %in% c(0,1)))
  N      <- length(response)
  u_resp <- unique(sort(response))
  n_cuts <- length(cutpts)
  x1 <- x2 <- rep(NA, n_cuts)
  if(is.null(weights)) weights = rep(1/N,N)
  f_0    <- weights[labels==0]/mean(weights[labels==0])
  resp_0 <- response[labels==0]
  f_1    <- weights[labels==1]/mean(weights[labels==1])
  resp_1 <- response[labels==1]
  for(j in 1:n_cuts){
    x1[j] <- mean((resp_1 >= cutpts[j])*f_1)
    x2[j] <- mean((resp_0 >= cutpts[j])*f_0)
  }
  sum((x1[1:(n_cuts-1)] + x1[2:(n_cuts)])/2 * (x2[2:(n_cuts)] - x2[1:(n_cuts-1)]))
}

## build the survey weighted dataset
nhanes_data_svy <- svydesign(id= ~SDMVPSU, strata = ~SDMVSTRA, 
                             weights = ~wtmec4yr_adj_norm, data = nhanes_data, nest = TRUE)


## set the seed so cross-validation results are reproducible
set.seed(1244)
## get the training and testing datasets for 10-fold cross validation
n_folds <- 10
## split the data to have an (approximately) equal number of alive/died in each training/test dataset
inx_id_alive <- which(nhanes_data[,dep_vars]==0)
inx_id_died  <- which(nhanes_data[,dep_vars]==1)
nid_alive    <- length(inx_id_alive)
nid_died     <- length(inx_id_died)
inx_ls_alive <- split(sample(inx_id_alive, size=nid_alive, replace=FALSE),
                      rep(1:n_folds,each=ceiling(nid_alive/n_folds))[1:nid_alive])
inx_ls_died <- split(sample(inx_id_died, size=nid_died, replace=FALSE),
                     rep(1:n_folds,each=ceiling(nid_died/n_folds))[1:nid_died])
inx_ls <- lapply(1:n_folds, function(x) c(inx_ls_alive[[x]], inx_ls_died[[x]]))
rm(list=c("inx_id_alive","inx_id_died","nid_alive","nid_died","inx_ls_alive","inx_ls_died"))

ind_vars <- c("Age", "Age + MVPA", "Age + ASTP", "Age + MVPA + ASTP", "Age + score")
auc_ijk_adj <-  matrix(NA, nrow=length(ind_vars), ncol=n_folds)
auc_mat_adj  <-  data.frame("Variable" = rep(NA_character_,length(ind_vars)),
                            "Cross-Validated AUC" = rep(NA_real_,length(ind_vars)),
                            "SE" = rep(NA_real_, length(ind_vars)),
                            stringsAsFactors = FALSE)
for(i in 1:length(ind_vars)){
  form_adj <- ind_vars[i] ## name of covariates

  ## get cross-validated AUC
  for(k in 1:n_folds){
    ## subset test and training data sets
    SEQN_train <- nhanes_data$SEQN[-inx_ls[[k]]]
    SEQN_test  <- nhanes_data$SEQN[inx_ls[[k]]]
    data_test  <- subset(nhanes_data, SEQN %in% SEQN_test)
    data_train <- subset(nhanes_data, SEQN %in% SEQN_train)
    ## Fit the appropriate models by subsetting the data.
    ## By subsetting the existing svydesign objects instead of creating new svydesign objects,
    ## we retain information on the number of PSU/strata in the original study.
    fit_adj_cv <- svyglm(as.formula(paste(dep_vars, " ~", form_adj)),
                         design=subset(nhanes_data_svy, SEQN %in% SEQN_train),
                         family=quasibinomial())
    ## Calculate weighted AUC using the appropriate weights
    auc_ijk_adj[i,k] <- calc_weighted_AUC(response=predict(fit_adj_cv, newdata=data_test, type='link'),
                                          cutpts=get_ctpts(fit_adj_cv,type='link',newdata=data_test),
                                          labels=data_test[,dep_vars],
                                          weights=data_test$wtmec4yr_adj_norm)
    rm(list=c("data_train","data_test","SEQN_train","SEQN_test",
              paste0("fit_adj_cv")))
  }
  # mean(auc_ijk_adj[i,])
  # sd(auc_ijk_adj[i,])/sqrt(10)

  auc_mat_adj[i,1] <- form_adj
  auc_mat_adj[i,2] <- mean(auc_ijk_adj[i,])
  auc_mat_adj[i,3] <- sd(auc_ijk_adj[i,]) / sqrt(n_folds)
}
auc.result <- auc_mat_adj[order(-auc_mat_adj[,2]),]

rm(nhanes_data_svy, auc_mat_adj, form_adj, auc_ijk_adj, inx_ls)


