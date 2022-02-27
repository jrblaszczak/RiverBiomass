##==============================================================================
## Script for fitting stan models to first year of 2 year data
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
options(mc.cores=16)

stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel_PAR, GPP = x$GPP,
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
             chains=4,iter=5000, control=list(max_treedepth=12))
launch_shinystan(test_ar)

#Ricker - P reparameterized
init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}

test_ricker <- stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                    data=stan_data_l$nwis_01608500,
                    init = init_Ricker,chains=4,iter=5000,
                    control=list(max_treedepth=12, adapt_delta=0.95))
launch_shinystan(test_ricker)

#Gompertz
init_Gompertz <- function(...) {
  list(c = 0.5, s = 1.5)
}

test_Gompertz <- stan("Stan_ProductivityModel3_Gompertz.stan",
                    data=stan_data_l$nwis_01608500,
                    init = init_Gompertz,chains=3,iter=5000,
                    control=list(max_treedepth=12, adapt_delta=0.95))
launch_shinystan(test_Gompertz)


###################################################
## Run Stan to get parameter estimates - all sites
###################################################

## PM 1 - Standard time series
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x, chains=4, iter=5000,
                                                   control=list(max_treedepth=12, adapt_delta=0.95)))
saveRDS(PM_outputlist_AR, "./rds files/stan_6riv_output_AR_2022_02_22.rds")


## PM 2 - Latent Biomass (Ricker)
init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}
PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                                                data=x, init = init_Ricker, chains=4, iter=5000,
                                                control=list(max_treedepth=12, adapt_delta=0.95)))
saveRDS(PM_outputlist_Ricker, "./rds files/stan_6riv_output_Ricker_2022_02_27.rds")


## PM 3 - Latent Biomass (Gompertz)
init_Gompertz <- function(...) {
  list(c = 0.5, s = 1.5)
}
PM_outputlist_Gompertz <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel3_Gompertz.stan",
                                                data=x, init = init_Gompertz, chains=3, iter=5000,
                                                control=list(max_treedepth=12, adapt_delta=0.95)))
saveRDS(PM_outputlist_Gompertz, "./rds files/stan_6riv_output_Gompertz_2022_01_23.rds")



PM_outputlist_AR <- readRDS("./rds files/stan_6riv_output_AR_2022_01_23.rds")
PM_outputlist_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker_2022_01_23.rds")
PM_outputlist_Gompertz <- readRDS("./rds files/stan_6riv_output_Gompertz_2022_01_23.rds")

#####################################
## Summary of divergent transitions
#####################################

launch_shinystan(PM_outputlist_AR$nwis_08447300)
## 01608500 (South Branch Potomac) - 0 divergent transitions
## 01649190 (Paint Branch) - 0 divergent transitions
## 02336526 (Proctor Creek) - 0 divergent transitions
## 07191222 (Beatty Creek) - 0 divergent transitions
## 08447300 (Pecos River) - 0 divergent transitions
## 11044000 (Santa Margarita) - 0 divergent transitions

launch_shinystan(PM_outputlist_Ricker$nwis_08447300)
## 01608500 (South Branch Potomac) - 0 divergent transitions
## 01649190 (Paint Branch) - 70 divergent transitions
## 02336526 (Proctor Creek) - 14 divergent transitions
## 07191222 (Beatty Creek) - 0 divergent transitions
## 08447300 (Pecos River) - 298 divergent transitions
## 11044000 (Santa Margarita) - 4 divergent transitions

launch_shinystan(PM_outputlist_Gompertz$nwis_11044000) # 3 chains, 5000 iterations
## 01608500 (South Branch Potomac) - 0 divergent transitions
## 01649190 (Paint Branch) - 0 divergent transitions
## 02336526 (Proctor Creek) - 3 divergent transitions
## 07191222 (Beatty Creek) - 0 divergent transitions
## 08447300 (Pecos River) - 1 divergent transitions - c not converging
## 11044000 (Santa Margarita) - 10 divergent transitions - c not converging





