## Fitting models to data
## JR Blaszczak

# load packages
devtools::install_github("collectivemedia/tictoc")
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","Metrics","MCMCglmm","tictoc"), require, character.only=T)

## Source data
source("DataSource_9rivers.R")

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
test <- stan("Stan_ProductivityModel3_Ricker_fixedinit.stan",
             data=stan_data_l$nwis_14211010,
             chains=3,iter=1000, control=list(max_treedepth=12))

## PM 1 - Phenomenological
PM_outputlist_AR <- lapply(stan_data_l, function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive.stan", data=x,chains=3,iter=5000, control=list(max_treedepth=12))) #5000 seconds
saveRDS(PM_outputlist_AR, "./rds files/stan_9riv_output_AR.rds") 

## PM 3 - Ricker
PM_outputlist_Ricker <- lapply(stan_data_l, function(x) stan("Stan_ProductivityModel3_Ricker.stan", data=x,chains=3,iter=5000, control=list(max_treedepth=12))) #78 seconds
saveRDS(PM_outputlist_Ricker, "./rds files/stan_9riv_output_Ricker.rds")

## PM 4 - Gompertz
PM_outputlist_Gompertz <- lapply(stan_data_l, function(x) stan("Stan_ProductivityModel5_Gompertz.stan", data=x,chains=3,iter=5000, control=list(max_treedepth=12)))
saveRDS(PM_outputlist_Gompertz, "./rds files/stan_3newriv_output_Gompertz.rds")

launch_shinystan(PM_outputlist_Ricker$nwis_01649190)




