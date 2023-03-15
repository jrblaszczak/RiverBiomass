##==============================================================================
## Script for fitting stan models to more years of training data for predictions
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork","grid","gridExtra"), require, character.only=T)

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

## visualize years included
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

ggplot(Pot_longtrain, aes(tQ))+
  geom_histogram()+
  geom_histogram(data = Pot_shorttrain, aes(tQ), fill = "blue")

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
## Compare out-of-sample predictions from long versus short TS
########################################################################

## Import model fit if needed
PM_outputlist_AR <- readRDS("./rds files/stan_Pot_output_AR_2023_03_12.rds")
PM_outputlist_Ricker <- readRDS("./rds files/stan_Pot_output_Ricker_2023_03_12.rds")

## Prep out-of-sample data for short versus long TS
# define oos year
Pot_oos <- data[which(data$year == "2013"),]
# extract longtrain and shorttrain max light and discharge by site
long_max <- Pot_longtrain %>% 
  summarise_at(.vars = c("Q","PAR_surface"), .funs = max)
short_max <- Pot_shorttrain %>% 
  summarise_at(.vars = c("Q","PAR_surface"), .funs = max)
# relativize Q and light by the max specific to the long and short training data sets
oos_relativize <- function(prev_max, post_dat){
  
  max.vals <- prev_max
  dat <- post_dat
  
  dat$light_rel_PAR <- dat$PAR_surface/max.vals$PAR_surface
  dat$tQ <- dat$Q/max.vals$Q
  
  dat <- dat[order(dat$date),]
  return(dat)
  
}

Pot_oos_long <- oos_relativize(long_max, Pot_oos)
Pot_oos_short <- oos_relativize(short_max, Pot_oos)

# visualize difference
ggplot(Pot_oos_long, aes(date, tQ))+
  geom_line()+
  geom_line(data = Pot_oos_short, aes(date, tQ), color="purple")


## source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p

## Ricker simulations using long and short posteriors
Ricker_sim_fxn <- function(x, oos_dat){
  #separate data
  output <- x
  df <- oos_dat
  
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

Ricker_sim_Pot_long <- Ricker_sim_fxn(PM_outputlist_Ricker$Pot_long, Pot_oos_long)
Ricker_sim_Pot_short <- Ricker_sim_fxn(PM_outputlist_Ricker$Pot_short, Pot_oos_short)

## Save simulation
saveRDS(Ricker_sim_Pot_long, "./rds files/sim_Pot_long_Ricker_2023_03_15.rds")
saveRDS(Ricker_sim_Pot_short, "./rds files/sim_Pot_short_Ricker_2023_03_15.rds")










