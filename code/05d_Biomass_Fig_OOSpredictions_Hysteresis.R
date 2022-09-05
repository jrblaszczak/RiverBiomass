##==============================================================================
## Script for out-of-sample prediction hysteresis
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

################################
## Latent Biomass Predictions
################################

LB_oos_preds <- function(preds, df){
  
  simmat_list <- readRDS(preds)
  simmat_list <- simmat_list[names(df)]
  
  ## Plot latent biomass
  # For every day extract median and CI
  median_biomat <- ldply(lapply(simmat_list, function(z) apply(z[[3]], 1, function(x) median(x))), data.frame)
  lower_biomat <- ldply(lapply(simmat_list, function(z) apply(z[[3]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
  upper_biomat <- ldply(lapply(simmat_list, function(z) apply(z[[3]], 1, function(x) quantile(x, probs = 0.975))), data.frame)
  
  ## Plot simulated GPP
  dat <- ldply(df, data.frame)
  df_bio <- as.data.frame(cbind(dat$site_name, as.character(dat$date),
                                median_biomat[,2], lower_biomat[,2], upper_biomat[,2]))
  colnames(df_bio) <- c("site_name","Date","B","B_lower","B_upper")
  df_bio$Date <- as.POSIXct(as.character(df_bio$Date), format="%Y-%m-%d")
  df_bio[,3:5] <- apply(df_bio[,3:5],2,function(x) as.numeric(as.character(x)))
  df_bio <- left_join(df_bio, site_info[,c("site_name","short_name")])
  
  return(df_bio)
  
}


LB_LBpreds <- LB_oos_preds("./rds files/Sim_6riv_Ricker_oos_2022_02_27.rds",df)



###############################
## Hysteresis Figures
###############################
oosdat_df <- ldply(dat_oos, data.frame)
colnames(oosdat_df)[which(colnames(oosdat_df) == "date")] <- "Date"

## combine with light
LBpred_light <- merge(LB_simdat, oosdat_df[,c("site_name","Date","light","PAR_surface","light_rel_PAR","Q")], by = c("site_name","Date"))
colnames(LBpred_light)

## combine with latent biomass preds
LBpreds_oos_df <- merge(LBpred_light, LB_LBpreds, by=c("site_name","Date","short_name"))

## Find storm events and recovery curves
l <- split(LBpreds_oos_df, LBpreds_oos_df$short_name)



ggplot(l$`S. Br. Potomac River, WV`, aes(PAR_surface, GPP))+
  geom_point()
ggplot(l$`S. Br. Potomac River, WV`, aes(PAR_surface, B))+
  geom_point()
ggplot(l$`S. Br. Potomac River, WV`, aes(PAR_surface, sim_GPP))+
  geom_point()
ggplot(l$`S. Br. Potomac River, WV`, aes(Date, sim_GPP))+
  geom_point()+
  geom_line(aes(Date, Q/50), color = "blue")+
  geom_line(aes(Date, exp(B)), color="chartreuse4")

View(l$`S. Br. Potomac River, WV`)


hysteresis_fig <- function(site, start, end){
  
  storm <- subset(site, Date >= start & Date <= end)
  
  plot_grid(
    ggplot(storm, aes(Date, Q))+
      geom_line(color="midnightblue", size=1)+
      geom_ribbon(aes(ymin=-Inf,ymax=Q), color="midnightblue", fill="midnightblue")+
      geom_line(aes(Date, PAR_surface), color = "purple")+
      scale_x_datetime(expand = c(0,0))+
      theme_bw(),
    
  ggplot(storm, aes(Date, sim_GPP))+
    geom_point()+
    geom_point(aes(Date, exp(B)), color="chartreuse4")+
    geom_line(aes(Date, exp(B)), color="chartreuse4")+
    scale_x_datetime(expand = c(0,0))+
    theme_bw(),
  
  plot_grid(
  
  ggplot(storm, aes(Q, sim_GPP))+
    geom_point(data = site, aes(Q, sim_GPP), alpha=0.2, color="gray75")+
    geom_point(data = storm, aes(Q, sim_GPP, color=Date), size=2)+
    geom_path()+
    theme_bw()+
    theme(legend.position = "none"),
  
  ggplot(storm, aes(Q, exp(B)))+
    geom_point(data = site, aes(Q, exp(B)), alpha=0.2, color="gray75")+
    geom_point(size=2, aes(color=Date))+
    geom_path()+
    theme_bw()+
    theme(legend.position = "none"),
  
  ggplot(storm, aes(PAR_surface, sim_GPP))+
    geom_point(data = site, aes(PAR_surface, sim_GPP), alpha=0.2, color="gray75")+
    geom_point(data = storm, aes(PAR_surface, sim_GPP, color=Date), size=2)+
    geom_path()+
    theme_bw()+
    theme(legend.position = "none"),
  
  ncol = 3),
  
  ncol = 1)
  
  
}

hysteresis_fig(l$`S. Br. Potomac River, WV`, "2013-05-01","2013-05-31")
#hysteresis_fig(l$`S. Br. Potomac River, WV`, "2013-01-28","2013-03-10")
#hysteresis_fig(l$`S. Br. Potomac River, WV`, "2013-04-09","2013-04-27")











