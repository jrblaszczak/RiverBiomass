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

#original Potomac for short vs long training
saveRDS(LB_simdat[which(LB_simdat$site_name == "nwis_01608500"),], "./rds files/Potomac_orig_oos_simdat.rds")

###################################################################
## Calculate metrics - coverage, RRMSE (accuracy), MRE (bias)
###################################################################

## split dataframes by site
STS_simdat_l <- split(STS_simdat, STS_simdat$site_name)
LB_simdat_l <- split(LB_simdat, LB_simdat$site_name)

## Goodness of fit metrics function
calc_gof_metrics <- function(x, mod){
  
  ts_sub <- x
  
  ##RMSE
  rmse <- sqrt(sum((ts_sub$sim_GPP-ts_sub$GPP)^2)/length(ts_sub$GPP))
  
  ##NRMSE
  nrmse <- rmse/(max(ts_sub$GPP)-min(ts_sub$GPP))
  
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
  metrics <- as.data.frame(cbind(rmse, nrmse, rrmse, mre, cov_pct,mod))
  return(metrics)
  
}

STS_gof <- ldply(lapply(STS_simdat_l, function(x) calc_gof_metrics(x,"S-TS")), data.frame)
LB_gof <- ldply(lapply(LB_simdat_l, function(x) calc_gof_metrics(x,"LB-TS")), data.frame)

gof <- rbind(STS_gof, LB_gof)
## export table
#write.csv(gof, "../figures and tables/2022 Figures and Tables/tables/TableS_OOS_gof.csv")


########################
## GOF figure for JASM
########################
names(gof)
gof_sub <- gof[,c(".id","nrmse","cov_pct","mod")]
colnames(gof_sub)[1] <- "site_name"
gof_sub <- left_join(gof_sub, site_info, by="site_name")
sapply(gof_sub, class)
gof_sub$nrmse <- as.numeric(gof_sub$nrmse)
gof_sub$cov_pct <- as.numeric(gof_sub$cov_pct)


gofplot_nrmse <- ggplot(gof_sub, aes(mod, nrmse, color = mod))+
  geom_boxplot()+
  geom_jitter(width = 0.2, aes(size=ORD_STRA))+
  labs(x="Model",y="NRMSE")+
  scale_color_manual(values = c("S-TS" = PM_AR.col,
                                "LB-TS" = PM_Ricker.col))+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(),axis.title = element_text(size=12),
        axis.text = element_text(size=12))

gofplot_cov_pct <- ggplot(gof_sub, aes(mod, cov_pct, color = mod))+
  geom_boxplot()+
  geom_jitter(width = 0.2, aes(size=ORD_STRA))+
  labs(x="Model",y="Coverage (%)")+
  scale_y_continuous(labels = function(x) paste0(x, '%'))+
  scale_color_manual(values = c("S-TS" = PM_AR.col,
                                "LB-TS" = PM_Ricker.col))+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(),axis.title = element_text(size=12),
        axis.text = element_text(size=12))

plot_grid(gofplot_nrmse+theme(legend.position = "none"),
          gofplot_cov_pct+theme(legend.position = "none"),
          nrow=2, align = "hv")
plot_grid(get_legend(gofplot_nrmse))



