##==============================================================================
## Script to compare RMSE within and out-of-sample and among S-TS and LB-TS
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","patchwork"), require, character.only=T)

PM_AR.col <- "#d95f02" # AR
PM_Ricker.col <- "#1C474D" # Ricker
PM_Gompertz.col <- "#743731" # Gompertz


################################
## Import and combine by site
################################
## Import within sample data
source("DataSource_6rivers_StreamLight.R")
dat_ws <- df

# Within-sample predictions
ws_AR <- readRDS("./rds files/Sim_6riv_AR_ws_2022_02_27.rds")
ws_Ricker <- readRDS("./rds files/Sim_6riv_Ricker_ws_2022_02_27.rds")

## Import within sample data
source("DataSource_6rivers_oos_StreamLight.R")

# Out-of-sample predictions
oos_AR <- readRDS("./rds files/Sim_6riv_AR_oos_2022_02_27.rds")
oos_Ricker <- readRDS("./rds files/Sim_6riv_Ricker_oos_2022_02_27.rds")

########################################
## Extract median daily GPP predictions
########################################
preds <- list(ws_AR, ws_Ricker, oos_AR, oos_Ricker)
names(preds) <- c("ws_AR", "ws_Ricker", "oos_AR", "oos_Ricker")
med_preds <- lapply(preds, function(x) lapply(x, function(z) apply(z[[1]], 1, function(x) median(x))))

## pair predictions to data sources and calculate rmse

rmse_calcs <- function(dat_source, p, pred_type){
  rmse_list <- rep(NA, 6)
  for(i in 1:length(names(dat_source))){
    
    data <- dat_source[[i]]
    GPP.pred <- p[[i]]
    
    rmse <- sqrt(sum((GPP.pred - data$GPP)^2)/length(data$GPP))
    rmse_list[[i]] <- rmse
  }
  
  rmse_df <- as.data.frame(cbind(names(dat_ws), rmse_list))
  colnames(rmse_df) <- c("site_name", "rmse")
  rmse_df$pred_type <- pred_type
  rmse_df$rmse <- as.numeric(rmse_df$rmse)
  return(rmse_df)
}

rmse_values <- rbind(rmse_calcs(dat_ws, med_preds$ws_AR, "Within-Sample S-TS"),
                     rmse_calcs(dat_ws, med_preds$ws_Ricker, "Within-Sample LB-TS"),
                     rmse_calcs(dat_oos, med_preds$oos_AR, "Out-of-Sample S-TS"),
                     rmse_calcs(dat_oos, med_preds$oos_Ricker, "Out-of-Sample LB-TS"))


rmse_spread <- spread(data = rmse_values, key = pred_type, value = rmse)

ggplot(rmse_spread, aes(`Within-Sample S-TS`, `Within-Sample LB-TS`, color = factor(site_name)))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size=4)+
  theme_bw()

ggplot(rmse_spread, aes(`Out-of-Sample S-TS`, `Out-of-Sample LB-TS`, color = factor(site_name)))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size=4)+
  theme_bw()


