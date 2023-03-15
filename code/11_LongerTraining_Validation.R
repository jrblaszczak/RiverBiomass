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

# remove 2013
Pot_longtrain <- data[which(data$year %in% c("2012","2014","2015","2016")),]
Pot_shorttrain <- data[which(data$year == "2012"),]

## Relativize light and discharge
rel_LQT <- function(x){
  x$light_rel_PPFD <- x$light/max(x$light)
  x$light_rel_PAR <- x$PAR_surface/max(x$PAR_surface)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)
  
  x<-x[order(x$date),]
  return(x)
}

Pot_longtrain <- rel_LQT(Pot_longtrain)
Pot_shorttrain <- rel_LQT(Pot_shorttrain)

ggplot(Pot_longtrain, aes(date, GPP))+
  geom_point()+
  geom_point(data = Pot_shorttrain, aes(date, GPP), color = "blue")+
  geom_point(data = data[which(data$year == "2013"),], aes(date, GPP), color = "purple")

plot_grid(
  ggplot(Pot_longtrain, aes(date, light_rel_PAR))+
  geom_point(),
  ggplot(Pot_longtrain, aes(date, tQ))+
    geom_point(),
  ncol= 1)

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

rm(test_ar, test_ricker)


###################################################
## Run Stan to get parameter estimates - all sites
###################################################
stan_data_l <- list(stan_data_Pot_long, stan_data_Pot_short)

## PM 1 - Standard time series
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x, chains=4, iter=5000,
                                                   control=list(max_treedepth=12, adapt_delta=0.95)))
names(PM_outputlist_AR) <- c("Pot_long","Pot_short")
saveRDS(PM_outputlist_AR, "./rds files/stan_Pot_output_AR_2023_03_12.rds")


## PM 2 - Latent Biomass (Ricker)
init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}
PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                                                data=x, init = init_Ricker, chains=4, iter=5000,
                                                control=list(max_treedepth=12, adapt_delta=0.99)))
names(PM_outputlist_Ricker) <- c("Pot_long","Pot_short")
launch_shinystan(PM_outputlist_Ricker$Pot_long)
launch_shinystan(PM_outputlist_Ricker$Pot_short)

saveRDS(PM_outputlist_Ricker, "./rds files/stan_Pot_output_Ricker_2023_03_12.rds")


########################################################################
## Compare out of sample predictions from long versus short TS
########################################################################

## if need to import
#PM_outputlist_AR <- readRDS("./rds files/stan_Pot_output_AR_2023_03_12.rds")
#PM_outputlist_Ricker <- readRDS("./rds files/stan_Pot_output_Ricker_2023_03_12.rds")

## Prep out-of-sample data for short versus long TS








# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p














