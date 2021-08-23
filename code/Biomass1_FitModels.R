##==============================================================================
## Script for fitting stan models to data
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools","shinystan",
         "MCMCglmm"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")

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

##############################################################
## Run Stan to get parameter estimates - initial tests
##############################################################

## Initial tests
#AR
test_ar <- stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
             data=stan_data_l$nwis_01608500,
             chains=3,iter=5000, control=list(max_treedepth=12))
launch_shinystan(test_ar)

#Ricker - P reparameterized
init_Ricker <- function(...) {
  list(c = 0.5, s = 0.5)
}
test_ricker <- stan("Stan_ProductivityModel2_Ricker_s_modification.stan",
                    data=stan_data_l$nwis_01608500,
                    init = init_Ricker,
                    chains=3,iter=5000, control=list(max_treedepth=12,
                                                     adapt_delta=0.95))
launch_shinystan(test_ricker)

#Gompertz
init_Gompertz <- function(...) {
  list(c = 0.5, s = 100)
}
test_Gompertz <- stan("Stan_ProductivityModel3_Gompertz_fixedinit_obserr.stan",
                    data=stan_data_l$nwis_01608500,
                    init = init_Gompertz,
                    chains=3,iter=5000, control=list(max_treedepth=12))
launch_shinystan(test_Gompertz)


###################################################
## Run Stan to get parameter estimates - all sites
###################################################

## PM 1 - Standard time series
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x,chains=3,iter=5000, control=list(max_treedepth=12)))
PM_AR_elapsedtime <- lapply(PM_outputlist_AR, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_AR, "./rds files/stan_6riv_output_AR_2021_08_12.rds")
saveRDS(PM_AR_elapsedtime, "./rds files/stan_6riv_AR_time_2021_08_12.rds")

## PM 2 - Latent Biomass (Ricker)
init_Ricker <- function(...) {
  list(c = 0.5, s = 0.5)
}

PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_s_modification.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12, adapt_delta=0.95)))
PM_Ricker_elapsedtime <- lapply(PM_outputlist_Ricker, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_Ricker, "./rds files/stan_6riv_output_Ricker_2021_08_23.rds")
saveRDS(PM_Ricker_elapsedtime, "./rds files/stan_6riv_Ricker_time_2021_08_23.rds")

#####################################
## View & check acceptance criteria
#####################################
stan_model_output_AR <- readRDS("./rds files/stan_6riv_output_AR_2021_08_12.rds")
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker_2021_08_12.rds")


traceplot(stan_model_output_AR$nwis_01608500, pars = c("phi","alpha","beta","sig_p","sig_o"))
stan_dens(stan_model_output_AR$nwis_01608500, pars = c("phi","alpha","beta","sig_p","sig_o"))

launch_shinystan(stan_model_output_Ricker$nwis_01649190)


