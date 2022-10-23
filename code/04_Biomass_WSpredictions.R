##==============================================================================
## Script for within-sample predictions
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"
PM_Gompertz.col <- "#1C474D"

## Import stan fits - simulate one at a time
stan_model_output_AR <- readRDS("./rds files/stan_6riv_output_AR_2022_02_22.rds")
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker_2022_02_27.rds")
#stan_model_output_Gompertz <- readRDS("./rds files/stan_6riv_output_Gompertz_2022_01_23.rds")

##########################
## Model 1 Output - AR
#########################
names(df); names(stan_model_output_AR)
AR_list <- Map(c, stan_model_output_AR, df)

AR_sim_fxn <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars1 <- extract(output, c("phi","alpha","beta","sig_p","sig_o"))
  simmat1<-matrix(NA,length(df$GPP),length(unlist(pars1$phi)))
  rmsemat1<-matrix(NA,length(df$GPP),1)
  
  # Simulate
  for (i in 1:length(pars1$phi)){
    simmat1[,i]<-PM_AR(pars1$phi[i],pars1$alpha[i],pars1$beta[i],pars1$sig_p[i],pars1$sig_o[i],df)
    rmsemat1[i]<-sqrt(sum((simmat1[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat1, rmsemat1)
  return(l)
  
}

#test <- AR_sim_fxn(AR_list$nwis_01608500)
AR_sim <- lapply(AR_list, function(x) AR_sim_fxn(x))

## Save simulation
saveRDS(AR_sim, "./rds files/Sim_6riv_AR_ws_2022_02_27.rds")


###############################
## Model 2 Output - Ricker
###############################
names(df); names(stan_model_output_Ricker)
Ricker_list <- Map(c, stan_model_output_Ricker, df[(names(stan_model_output_Ricker))])

Ricker_sim_fxn <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars3<-extract(output, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o"))
  simmat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  biomat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  rmsemat3<-matrix(NA,length(df$GPP),1)
  #Simulated
  for (i in 1:length(pars3$r)){
    simmat3[,i]<-PM_Ricker(pars3$r[i],pars3$lambda[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    biomat3[,i]<-PM_Ricker_B(pars3$r[i],pars3$lambda[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
  }

  l <- list(simmat3, rmsemat3, biomat3)
  return(l)
  
}

Ricker_sim <- lapply(Ricker_list, function(x) Ricker_sim_fxn(x))

## Save simulation
saveRDS(Ricker_sim, "./rds files/Sim_6riv_Ricker_ws_2022_02_27.rds")


###############################
## Model 3 Output - Gompertz
###############################
names(df); names(stan_model_output_Gompertz)
Gompertz_list <- Map(c, stan_model_output_Gompertz, df)

Gompertz_sim_fxn <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars5<-extract(output, c("beta_0","beta_1","s","c","B","P","pred_GPP","sig_p"))
  simmat5<-matrix(NA,length(df$GPP),length(unlist(pars5$sig_p)))
  biomat5<-matrix(NA,length(df$GPP),length(unlist(pars5$sig_p)))
  rmsemat5<-matrix(NA,length(df$GPP),1)
  #Simulate
  for (i in 1:length(pars5$sig_p)){
    simmat5[,i]<-PM_Gompertz(pars5$beta_0[i],pars5$beta_1[i],pars5$s[i],pars5$c[i],pars5$sig_p[i],pars5$sig_o[i],df)
    biomat5[,i]<-PM_Gompertz_B(pars5$beta_0[i],pars5$beta_1[i],pars5$s[i],pars5$c[i],pars5$sig_p[i],pars5$sig_o[i],df)
    rmsemat5[i]<-sqrt(sum((simmat5[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat5, rmsemat5, biomat5)
  return(l)
  
}

Gompertz_sim <- lapply(Gompertz_list, function(x) Gompertz_sim_fxn(x))

## Save simulation
saveRDS(Gompertz_sim, "./rds files/Sim_6riv_Gompertz_ws_2022_01_23.rds")
