##==============================================================================
## Script for extracting information for results
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

######################
## S-TS description
######################


######################
## LB-TS description
######################


#############################################
## Persistence and critical flow thresholds
#############################################
## Import stan model fit
stan_model_output_STS <- readRDS("./rds files/stan_6riv_output_AR_2022_02_22.rds")








