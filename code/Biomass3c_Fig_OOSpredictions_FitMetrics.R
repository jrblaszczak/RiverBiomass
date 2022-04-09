##==============================================================================
## Script for out-of-sample prediction metric comparison
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork","grid","gridExtra"), require, character.only=T)

## Source data
source("DataSource_6rivers_oos_StreamLight.R")
df <- dat_oos

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p
#source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"
PM_Gompertz.col <- "#1C474D"

## Check rivers
site_info[,c("site_name","long_name","NHD_STREAMORDE")]

################################
## Predicted time series - compile original GPP data with simulated GPP based on median parameter estimates
################################
GPP_oos_preds_ts <- function(preds, df){
  
  simmat_list <- readRDS(preds)
  simmat_list <- simmat_list[names(df)]
  
  # For every day extract median and CI
  mean_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) mean(x))), data.frame)
  lower_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
  upper_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)
  
  ## Plot simulated GPP
  dat <- ldply(df, data.frame)
  df_sim <- as.data.frame(cbind(dat$site_name, as.character(dat$date), dat$GPP, mean_simmat$X..i.., lower_simmat$X..i.., upper_simmat$X..i..))
  colnames(df_sim) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
  df_sim$Date <- as.POSIXct(as.character(df_sim$Date), format="%Y-%m-%d")
  df_sim[,3:6] <- apply(df_sim[,3:6],2,function(x) as.numeric(as.character(x)))
  
  ## Arrange rivers by river order
  df_sim <- left_join(df_sim, site_info[,c("site_name","short_name")])
  df_sim$short_name <- factor(df_sim$short_name, levels=site_order_list)
  
  return(df_sim)
  
}

STS_simdat <- GPP_oos_preds_ts("./rds files/Sim_6riv_AR_oos_2022_02_27.rds",df)
LB_simdat <- GPP_oos_preds_ts("./rds files/Sim_6riv_Ricker_oos_2022_02_27.rds",df)

###################################################################
## Calculate metrics - coverage, RRMSE (accuracy), MRE (bias)
###################################################################

## split dataframes by site
STS_simdat_l <- split(STS_simdat, STS_simdat$site_name)
LB_simdat_l <- split(LB_simdat, LB_simdat$site_name)

## Goodness of fit metrics function
calc_gof_metrics <- function(x, mod){
  
  ts_sub <- x
  
  ## RRMSE
  rrmse_pre <- sum(((ts_sub$sim_GPP-ts_sub$GPP)/ts_sub$GPP)^2)
  rrmse <- sqrt(rrmse_pre*(100/length(ts_sub$GPP)))
  
  ## MRE
  mre_pre <- sum((ts_sub$sim_GPP-ts_sub$GPP)/ts_sub$GPP)
  mre <- mre_pre*(100/length(ts_sub$GPP))
  
  ## Coverage
  ts_sub$c_yn <- ifelse(ts_sub$GPP >= ts_sub$sim_GPP_lower & ts_sub$GPP <= ts_sub$sim_GPP_upper,
                        yes=1, no=0)
  cov_pct <- (sum(ts_sub$c_yn)/length(ts_sub$c_yn))*100
  
  ## Compile
  metrics <- as.data.frame(cbind(rrmse, mre, cov_pct,mod))
  return(metrics)
  
}

STS_gof <- ldply(lapply(STS_simdat_l, function(x) calc_gof_metrics(x,"S-TS")), data.frame)
LB_gof <- ldply(lapply(LB_simdat_l, function(x) calc_gof_metrics(x,"LB-TS")), data.frame)

gof <- rbind(STS_gof, LB_gof)
## export table
write.csv(gof, "../figures and tables/2022 Figures and Tables/tables/TableS_OOS_gof.csv")














