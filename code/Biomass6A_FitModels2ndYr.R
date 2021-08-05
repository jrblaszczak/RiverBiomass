## Fitting models to data
## JR Blaszczak

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools","shinystan"), require, character.only=T)

## Source data
source("DataSource_6rivers_oos_StreamLight.R")
df <- dat_oos

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
## Run Stan to get parameter estimates - all sites
#########################################

## PM 1 - Standard time series
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x,chains=3,iter=5000, control=list(max_treedepth=12)))
saveRDS(PM_outputlist_AR, "./rds files/stan_6riv_2ndYr_output_AR_2021_08_04.rds")

## PM 2 - Latent Biomass (Ricker)
init_Ricker <- function(...) {
  list(c = 0.5, s = 100)
}

PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12)))
PM_Ricker_elapsedtime <- lapply(PM_outputlist_Ricker, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_Ricker, "./rds files/stan_6riv_2ndYr_output_Ricker_2021_06_15.rds")
saveRDS(PM_Ricker_elapsedtime, "./rds files/stan_6riv_2ndYr_Ricker_time_2021_06_15.rds")


## View

