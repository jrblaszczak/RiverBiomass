## Supplementary Table - Parameter estimate summary
## columns: Site, Parameter, Median, 2.5%, 97.5%

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","reshape2"), require, character.only=T)

## Import stan fits
stan_model_output_AR_yr1 <- readRDS("./rds files/stan_6riv_output_AR_2021_06_01.rds")
stan_model_output_AR_yr2 <- readRDS("./rds files/stan_6riv_2ndYr_output_AR_2021_08_04.rds")

stan_model_output_Ricker_yr1 <- readRDS("./rds files/stan_6riv_output_Ricker_2021_06_01.rds")
stan_model_output_Ricker_yr2 <- readRDS("./rds files/stan_6riv_2ndYr_output_Ricker_2021_06_15.rds")

## List together
AR_outputs <- list(stan_model_output_AR_yr1, stan_model_output_AR_yr2)
Ricker_outputs <- list(stan_model_output_Ricker_yr1, stan_model_output_Ricker_yr2)
names(AR_outputs) <- c("AR_yr1", "AR_yr2")
names(Ricker_outputs) <- c("Ricker_yr1", "Ricker_yr2")

## Extract parameters
AR_pars <- c("phi","alpha","beta","sig_p","sig_o")
Ricker_pars <- c("r","lambda","s","c","sig_p","sig_o")

AR_parout <- lapply(AR_outputs, function(x) lapply(x, function(y) rstan::extract(y, AR_pars)))
Ricker_parout <- lapply(Ricker_outputs, function(x) lapply(x, function(y) rstan::extract(y, Ricker_pars)))

#######################
## Parameter summary
######################
par_summary <- function(par) {
  
  ## Find the median
  median_par <- ldply(lapply(par, function(x) median(x)), data.frame)
  lower_par <- ldply(lapply(par, function(x) quantile(x, probs = 0.025)), data.frame)
  upper_par <- ldply(lapply(par, function(x) quantile(x, probs = 0.975)), data.frame)
  
  ## Compile in list and return
  par_vals <- cbind(median_par, lower_par$X..i.., upper_par$X..i..)
  names(par_vals) <- c("parameter","median","lowerCI","upperCI")
  return(par_vals)
}

sum_AR <- lapply(AR_parout, function(x) ldply(lapply(x, function(y) par_summary(y)), data.frame))
sum_Ricker <- lapply(Ricker_parout, function(x) ldply(lapply(x, function(y) par_summary(y)), data.frame))

## Save as separate tables
write.csv(sum_AR$AR_yr1, "./rds files/AR_medianpar_Yr1.csv")
write.csv(sum_AR$AR_yr2, "./rds files/AR_medianpar_Yr2.csv")

write.csv(sum_Ricker$Ricker_yr1, "./rds files/Ricker_medianpar_Yr1.csv")
write.csv(sum_Ricker$Ricker_yr2, "./rds files/Ricker_medianpar_Yr2.csv")















