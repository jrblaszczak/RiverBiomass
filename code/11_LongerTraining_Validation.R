##==============================================================================
## Script for fitting stan models to more years of training data for predictions
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools","shinystan",
         "MCMCglmm"), require, character.only=T)

## Read in data files
Pot_TS <- readRDS("./rds files/SBPotomac_longTS.rds")

## Read in streamlight
Pot_SL <- readRDS("./rds files/SBPotomac_SL.rds")

#######################################################################
## Longer time series data prep - South Branch Potomac River, WV
#######################################################################

## Join TS data and StreamLight
colnames(Pot_SL)[colnames(Pot_SL) == "Date"] <- "date"
data <- left_join(Pot_TS, Pot_SL, by=c("site_name", "date"))

## How many days of data per site per year
data$year <- year(data$date)
data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()
data_siteyears$year
# SB Potomac = 2008 2010 2012 2013 2014 2015 2016 (previous out-of-sample year was 2013 predicted from 2012)

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
data[which(data$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

## Relativize light and discharge
rel_LQT <- function(x){
  x$light_rel_PPFD <- x$light/max(x$light)
  x$light_rel_PAR <- x$PAR_surface/max(x$PAR_surface)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)
  
  x<-x[order(x$date),]
  return(x)
}

Pot_df <- rel_LQT(data)
# remove 2013
Pot_longtrain <- Pot_df[-which(Pot_df$year == "2013"),]
Pot_shorttrain <- Pot_df[which(Pot_df$year == "2012"),]

ggplot(Pot_longtrain, aes(date, GPP))+
  geom_point()+
  geom_point(data = Pot_shorttrain, aes(date, GPP), color = "blue")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
options(mc.cores=4)

stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel_PAR, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

stan_data_Pot_long <- stan_data_compile(Pot_longtrain)
stan_data_Pot_short <- stan_data_compile(Pot_shorttrain)

##############################################################
## Run Stan to get parameter estimates - initial tests
##############################################################

## Initial tests
#AR
test_ar <- stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
             data=stan_data_Pot_short,
             chains=4,iter=5000, 
             control=list(max_treedepth=12, adapt_delta=0.95))
launch_shinystan(test_ar)

#Ricker - P reparameterized
init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}

test_ricker <- stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                    data=stan_data_Pot_short,
                    init = init_Ricker,chains=4,iter=5000,
                    control=list(max_treedepth=12, adapt_delta=0.95))
launch_shinystan(test_ricker)


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





