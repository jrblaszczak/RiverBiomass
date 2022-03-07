##==============================================================================
## Script for fitting stan models to second year of 2 year data
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools","shinystan",
         "ggExtra","patchwork","grid","gridExtra"), require, character.only=T)

## Source data
source("DataSource_6rivers_2ndYr_StreamLight.R")

##################
## Examine data ##
##################

plotting_covar <- function(x) {
  
  df <- x
  
  gpp_plot <- ggplot(df, aes(date, GPP))+
    geom_point(color="chartreuse4", size=2)+
    geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$short_name[1])+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12))
  
  Q_plot <- ggplot(df, aes(date, tQ))+
    geom_line(size=1.5, color="deepskyblue4")+
    labs(y="Rel. Q")+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12))
  
  Temp_plot <- ggplot(df, aes(date, temp))+
    geom_line(size=1.5, color="#A11F22")+
    labs(y="Water Temp (C)")+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12))
  
  Light_plot <- ggplot(df, aes(date, light_rel_PAR))+
    geom_point(size=2, color="darkgoldenrod3")+
    labs(y="Rel. Light", x="Date")+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=12))
  
  
  
  p <- gpp_plot/Q_plot/Temp_plot/Light_plot
  
  return(p)
  
}

plotting_covar(df$nwis_11044000)



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

#########################################
## Run Stan to get parameter estimates - all sites
#########################################

## PM 1 - Standard time series
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x, chains=4, iter=5000,
                                                   control=list(max_treedepth=12, adapt_delta=0.95)))
saveRDS(PM_outputlist_AR, "./rds files/stan_6riv_2ndYr_output_AR_2022_03_06.rds")
#launch_shinystan(PM_outputlist_AR$nwis_01608500)


## PM 2 - Latent Biomass (Ricker)
init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}
PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                                                data=x, init = init_Ricker, chains=4, iter=5000,
                                                control=list(max_treedepth=12, adapt_delta=0.95)))
saveRDS(PM_outputlist_Ricker, "./rds files/stan_6riv_2ndYr_output_Ricker_2022_03_06.rds")

PM_outputlist_Ricker <- readRDS("./rds files/stan_6riv_2ndYr_output_Ricker_2022_03_06.rds")
launch_shinystan(PM_outputlist_Ricker$nwis_11044000)








