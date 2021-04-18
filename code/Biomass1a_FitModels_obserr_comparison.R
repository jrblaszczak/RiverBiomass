## Fitting models to data
## JR Blaszczak

# load packages
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
## Run Stan to get parameter estimates - initial tests
#########################################

## Initial tests
#AR
test_ar <- stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
             data=stan_data_l$nwis_08180700,
             chains=3,iter=5000, control=list(max_treedepth=12))
launch_shinystan(test_ar)

#Ricker
init_Ricker <- function(...) {
  list(c = 0.5, s = 200)
}
test_ricker <- stan("Stan_ProductivityModel3_Ricker_fixedinit_obserr.stan",
             data=stan_data_l$nwis_08180700,
             init = init_Ricker,
             chains=3,iter=2000, control=list(max_treedepth=12))
launch_shinystan(test_ricker)


#########################################
## Run Stan to get parameter estimates - all sites
#########################################

## PM 1 - Phenomenological
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x,chains=3,iter=5000, control=list(max_treedepth=12))) #5000 seconds
PM_AR_elapsedtime <- lapply(PM_outputlist_AR, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_AR, "./rds files/stan_9riv_output_AR_2021_04_18_obserr.rds")
saveRDS(PM_AR_elapsedtime, "./rds files/stan_9riv_AR_time_2021_04_18_obserr.rds")

## PM 3 - Ricker
init_Ricker <- function(...) {
  list(c = 0.5, s = 200)
}

PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel3_Ricker_fixedinit_obserr.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12)))
PM_Ricker_elapsedtime <- lapply(PM_outputlist_Ricker, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_Ricker, "./rds files/stan_9riv_output_Ricker_2021_04_18_obserr.rds")
saveRDS(PM_Ricker_elapsedtime, "./rds files/stan_9riv_Ricker_time_2021_04_18_obserr.rds")

## PM 4 - Gompertz
init_Gompertz <- function(...) {
  list(c = 0.5, s = 200)
}

PM_outputlist_Gompertz <- lapply(stan_data_l,
                                 function(x) stan("Stan_ProductivityModel5_Gompertz_fixedinit.stan",
                                                  data=x,chains=3,iter=5000, 
                                                  control=list(max_treedepth=12)))
PM_Gompertz_elapsedtime <- lapply(PM_outputlist_Gompertz, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_Gompertz, "./rds files/stan_9riv_output_Gompertz_2021_03_05.rds")
saveRDS(PM_Gompertz_elapsedtime, "./rds files/stan_9riv_Gompertz_time_2021_03_05.rds")




