# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools","shinystan",
         "MCMCglmm"), require, character.only=T)

theme_set(theme(legend.position = "none",
                panel.background = element_rect(color = "black", fill=NA, size=1),
                axis.text = element_text(size=12),
                axis.title = element_text(size=12)))


#############################
## Source & visualize data
#############################
source("DataSource_6rivers_StreamLight.R")

## Choose example data
ex <- df$nwis_01608500

## Stan data prep
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

stan_data_compile_PAR <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel_PAR, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

ex_stan <- stan_data_compile_PAR(ex)

test <- stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
             data = ex_stan, chains = 3, iter = 3000)
             #control=list(max_treedepth=12))
launch_shinystan(test)














