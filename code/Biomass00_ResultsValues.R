##==============================================================================
## Script for extracting information for results
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

########################################
## Parameter summaries
########################################
## Import stan model fit
stan_model_output_STS <- readRDS("./rds files/stan_6riv_output_AR_2022_02_22.rds")
stan_model_output_LBTS <- readRDS("./rds files/stan_6riv_output_Ricker_2022_02_27.rds")

stan_psum <- function(x){
  
  s_init <- as.data.frame(summary(x)$summary)
  s_init$pars <- row.names(s_init)
  s_init$n_eff_pct <- s_init$n_eff/20000  ## effective samples are the number of independent samples with the same estimation power as the N autocorrelated samples
  s_init$n_eff_less10pct <- ifelse(s_init$n_eff_pct < 0.10, yes = "true", no = "false") # 10% is often used as a threshold, below which the chains for a parameter did not properly converge
  s <- s_init[,c("pars","50%", "2.5%", "97.5%","Rhat","n_eff","n_eff_less10pct")]

  return(s)
  
}

pSTS <- ldply(lapply(stan_model_output_STS, function(z) stan_psum(z)), data.frame)
#STS params of interest: c("phi","alpha","beta","sig_p","sig_o")
write.csv(pSTS, "./tables/STS_ws_posterior_sum.csv")
pSTS_sub <- pSTS[which(pSTS$pars %in% c("phi","alpha","beta","sig_p","sig_o")),]
write.csv(pSTS_sub, "./tables/STS_ws_posteriorsubset_sum.csv")

pLBTS <- ldply(lapply(stan_model_output_LBTS, function(z) stan_psum(z)), data.frame)
#LB-TS params of interest: c("r","lambda","s","c","sig_p","sig_o")
write.csv(pLBTS, "./tables/LBTS_ws_posterior_sum.csv")
pLBTS_sub <- pLBTS[which(pLBTS$pars %in% c("r","lambda","s","c","sig_p","sig_o")),]
write.csv(pLBTS_sub, "./tables/LBTS_ws_posteriorsubset_sum.csv")

######################
## S-TS description
######################


######################
## LB-TS description
######################


#############################################
## Critical flow thresholds and sensitivity
#############################################
## Import and merge bankfull discharge with site info
RI_2 <- read.csv("../data/RI_2yr_flood_6riv.csv", header=T)
sapply(RI_2, class)
site_info <- merge(site_info, RI_2, by="site_name")













