## Fitting models to data
## JR Blaszczak

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","Metrics","MCMCglmm","tictoc"), require, character.only=T)

## Source data
source("DataSource_6rivers.R")
df <- df[c("nwis_07191222","nwis_01608500")]

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates
#########################################

## PM 1 - Standard time series
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x,chains=3,iter=5000, control=list(max_treedepth=12))) #5000 seconds

## PM 2 - Latent Biomass (Ricker)
init_Ricker <- function(...) {
  list(c = 0.5, s = 200)
}

PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12)))

## PM 3 - Latent Biomass (Ricker) - light modification
PM_outputlist_Ricker_Lmod <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr_Lmod.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12)))

## PM 3 - Latent Biomass (Ricker) - GPP modification
PM_outputlist_Ricker_Gmod <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr_Gmod.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12)))

###################################
## Visualize differences among models
##################################

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p, sig_o

##########################
## Model 1 Output - AR
#########################
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

AR_sim <- lapply(PM_outputlist_AR, function(x) AR_sim_fxn(x))


###############################
## Model 2 Output - Ricker
###############################
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

Ricker_sim <- lapply(PM_outputlist_Ricker, function(x) Ricker_sim_fxn(x))



###############################
## Model 3 Output - Ricker, Lmod
###############################
Ricker_sim_fxn_Lmod <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars3<-extract(output, c("r","lambda","alpha","s","c","B","P","pred_GPP","sig_p","sig_o"))
  simmat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  biomat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  rmsemat3<-matrix(NA,length(df$GPP),1)
  #Simulated
  for (i in 1:length(pars3$r)){
    simmat3[,i]<-PM_Ricker_Lmod(pars3$r[i],pars3$lambda[i],pars3$alpha[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    biomat3[,i]<-PM_Ricker_Lmod_B(pars3$r[i],pars3$lambda[i],pars3$alpha[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat3, rmsemat3, biomat3)
  return(l)
  
}

Ricker_sim_Lmod <- lapply(PM_outputlist_Ricker_Lmod, function(x) Ricker_sim_fxn_Lmod(x))


###############################
## Model 4 Output - Ricker, Gmod
###############################
Ricker_sim_fxn_Gmod <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars3<-extract(output, c("r","lambda","alpha","s","c","B","P","pred_GPP","sig_p","sig_o"))
  simmat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  biomat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  rmsemat3<-matrix(NA,length(df$GPP),1)
  #Simulated
  for (i in 1:length(pars3$r)){
    simmat3[,i]<-PM_Ricker_Gmod(pars3$r[i],pars3$lambda[i],pars3$alpha[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    biomat3[,i]<-PM_Ricker_Gmod_B(pars3$r[i],pars3$lambda[i],pars3$alpha[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat3, rmsemat3, biomat3)
  return(l)
  
}

Ricker_sim_Gmod <- lapply(PM_outputlist_Ricker_Gmod, function(x) Ricker_sim_fxn_Gmod(x))












