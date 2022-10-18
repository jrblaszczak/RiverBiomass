##==============================================================================
## Script for within-sample RRMSE & MRE
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"
PM_Gompertz.col <- "#1C474D"

##########################
## Model 1 Output - AR
#########################
simmat1_list <- readRDS("./rds files/Sim_6riv_AR_ws_2022_02_27.rds")

# For every day extract median and CI
median_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
lower_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
upper_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)

## Plot simulated GPP
dat1 <- ldply(df, data.frame)
df_sim1 <- as.data.frame(cbind(dat1$site_name, as.character(dat1$date), dat1$GPP, median_simmat1$X..i.., lower_simmat1$X..i.., upper_simmat1$X..i..))
colnames(df_sim1) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
df_sim1$Date <- as.POSIXct(as.character(df_sim1$Date), format="%Y-%m-%d")
df_sim1[,3:6] <- apply(df_sim1[,3:6],2,function(x) as.numeric(as.character(x)))

## Arrange rivers by river order
df_sim1 <- left_join(df_sim1, site_info[,c("site_name","short_name")])
df_sim1$short_name <- factor(df_sim1$short_name, levels=site_order_list)