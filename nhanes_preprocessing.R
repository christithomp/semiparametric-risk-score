## This code pre-processes the NHANES data to derive variables for the analysis

## Load package
library(rlang)
library(ggplot2)
library(corrplot)
library(mgcv)
library(caret)
library(survey)
library(HRW)
library(gridExtra)
library(MASS)
library(scam)
library(rjags)
library(dplyr)

## Load the data
data_analysis <- read.csv("./final_code/data/NHANES_data.csv")
nhanes_data <- data_analysis

## Winsorize PA variables by upper 95th percentile and 5th percentile
nhanes_data[which(data_analysis$MVPA > quantile(data_analysis$MVPA, 0.95)), "MVPA"] <- quantile(data_analysis$MVPA, 0.95)
nhanes_data[which(data_analysis$ASTP > quantile(data_analysis$ASTP, 0.95)), "ASTP"] <- quantile(data_analysis$ASTP, 0.95)
nhanes_data[which(data_analysis$MVPA < quantile(data_analysis$MVPA, 0.05)), "MVPA"] <- quantile(data_analysis$MVPA, 0.05)
nhanes_data[which(data_analysis$ASTP < quantile(data_analysis$ASTP, 0.05)), "ASTP"] <- quantile(data_analysis$ASTP, 0.05)

## Scale PA variables
nhanes_data[,c("ASTP_reg","MVPA_reg")] <- nhanes_data[,c("ASTP","MVPA")] ## original scale
nhanes_data[,c("ASTP","MVPA")] <- apply(nhanes_data[,c("ASTP", "MVPA")] , 2, scale)

## Survival outcome
dep_vars <- "yr9_mort"
## remove people who are censored before interested survival time
nhanes_data <- nhanes_data[which(!is.na(nhanes_data[,dep_vars])),]

## Reorganize the covariates 
# Divide age by 100 for numeric stability
nhanes_data$Age <- nhanes_data$Age / 100

# Create indicator for obese BMI
nhanes_data$BMI_Obese <- ifelse(nhanes_data$BMI_cat == "Obese", 1, 0)

# Create indicator for nondrinker
nhanes_data$Nondrinker <- ifelse(nhanes_data$DrinkStatus == "Non-Drinker", 1, 0)

# Create indicator for white
nhanes_data$White <- ifelse(nhanes_data$Race == "White", 1, 0)

# Create indicator for education more than high school
nhanes_data$EduMoreThanHS <- ifelse(nhanes_data$EducationAdult == "More than high school", 1, 0)

# Create indicator for smoker
nhanes_data$Smoker <- ifelse(nhanes_data$SmokeCigs == "Current", 1, 0)

# Create indicators for Overall health categories
nhanes_data$Health_verygood <- ifelse(nhanes_data$Overall_health == "Very good", 1, 0)
nhanes_data$Health_excellent <- ifelse(nhanes_data$Overall_health == "Excellent", 1, 0)
nhanes_data$Health_good <- ifelse(nhanes_data$Overall_health == "Good", 1, 0)
nhanes_data$Health_fair <- ifelse(nhanes_data$Overall_health == "Fair", 1, 0)
nhanes_data$Health_poor <- ifelse(nhanes_data$Overall_health == "Poor", 1, 0)
nhanes_data$Health_MoreThanGood <- ifelse(nhanes_data$Health_excellent + nhanes_data$Health_verygood >= 1, 1, 0)
nhanes_data$Health_LessThanGood <- ifelse(nhanes_data$Health_fair + nhanes_data$Health_poor >= 1, 1, 0)

## Separate data by gender
indM <- which(nhanes_data$Gender == 'Male')
indW <- which(nhanes_data$Gender == 'Female')
dataM <- nhanes_data[indM, ]
dataW <- nhanes_data[indW, ]

## Derive components of single index model
### Response
yM <- as.matrix(dataM[,dep_vars], ncol = 1)
yW <- as.matrix(dataW[,dep_vars], ncol = 1)
colnames(yM) <- colnames(yW) <- "All-Cause Mortality"

### Other covariates
zM <- as.matrix(cbind(dataM$Age, dataM$Smoker, dataM$White, 
                      dataM$EduMoreThanHS, dataM$BMI_Obese,
                      dataM$Health_MoreThanGood, dataM$Health_LessThanGood), ncol = 7)
zW <- as.matrix(cbind(dataW$Age, dataW$Smoker, dataW$White, 
                      dataW$EduMoreThanHS, dataW$BMI_Obese,
                      dataW$Health_MoreThanGood, dataW$Health_LessThanGood), ncol = 7)
colnames(zM) <- colnames(zW) <- c("Age", "Smoker", "White", "EduMoreThanHS", "BMI_Obese", 
                                  "HealthMoreThanGood", "HealthLessThanGood")

### PA variables
var <- c("MVPA", "ASTP")
var_reg <- c("MVPA_reg", "ASTP_reg")
xM <- as.matrix(dataM[,var], ncol = length(var))
xW <- as.matrix(dataW[,var], ncol = length(var))
colnames(xM) <- colnames(xW) <- var

