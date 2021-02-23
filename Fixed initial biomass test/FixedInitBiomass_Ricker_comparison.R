## Comparing Ricker with and without fixed initial biomass

## See simulation_matrix_output.R for code to create figures
# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","rstan","shinystan","MCMCglmm"), require, character.only=T)

#################
## Source data
#################
source("DataSource_9rivers.R")
df <- dat

####################
## Stan data prep 
####################
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ, B_int=log(x$GPP[1]/x$light_rel[1]))
  return(data)
}

stan_data_l <- lapply(dat, function(x) stan_data_compile(x))

########################
## Run Stan comparison
########################

## Limit to just 2 sites for comparison
stan_data_l <- stan_data_l[3:4]

## Ricker - free initial biomass
output_Ricker <- lapply(stan_data_l, function(x) stan("Stan_ProductivityModel3_Ricker.stan",
                                                             data=x,chains=3,iter=5000,
                                                             control=list(max_treedepth=12)))

## Ricker - fixed initial biomass
output_Ricker_fixedinit <- lapply(stan_data_l, function(x) stan("Stan_ProductivityModel3_Ricker_fixedinit.stan",
                                                      data=x,chains=3,iter=5000,
                                                      control=list(max_treedepth=12)))
















