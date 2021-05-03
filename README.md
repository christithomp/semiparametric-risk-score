These R programs accompany the paper "A Semiparametric Risk Score for Physical Activity." They let a user reproduce the results from Section 4 and 5 of the paper.

The NHANES data is available in data/NHANES_preprocessed.csv. To preprocess the data for analysis, use nhanes_preprocessing.R. The code to calculate the 10-fold survey-weighted cross-validated AUC for both methods is given in calculation_AUC.R. After estimation and inference, the figures for both methods are combined into a single file using combine_images.R.


PROFILE LIKELIHOOD:
--------------------------------------------------------------------------------------------------------
The code in the codePL folder gives the results for the profile likelihood (PL) method from Section 3.1. 

APPLICATION:

Implementation of Steps 1-3 of Section 3.1 are programmed in application_PL_update_parameters.R. The program application_PL_bootstrap.R loops through Steps 1-3 for each nonparametric bootstrap sample. Inference on these samples from Sections 4.2.1 and 4.2.2 is implemented in application_PL_inference.R. For PL results in 4.2.3, use application_PL_indep_score.R. 

SIMULATION:

For simulation results from Section 5.2, use simulation_PL.R to run in parallel on a server. These results are combined and saved into an rds file on the server. They are then cleaned using simulation_PL_organize_results.R. Inference is performed in simulation_PL_inference.R.


BAYESIAN METHOD:
----------------------------------------------------------------------------------------------
The code in the codeJAGS folder provides the results for the Bayesian method from Section 3.2. 

APPLICATION:

The JAGS code for the model and prior is given in application.jags. Implementation and results for Sections 4.2.1 and 4.2.2 in R is provided in application_jags.R. For Bayesian method results in 4.2.3, use application_dependent_score.jags for the dependent JAGS model and application_jags_indep_score.R in R for results.

SIMULATION:

For simulation results from Section 5.2, use simulation_jags.R for estimation and inference. 


