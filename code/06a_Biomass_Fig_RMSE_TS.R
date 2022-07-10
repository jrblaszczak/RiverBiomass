## OOS Prediction and direct comparison figure

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork"), require, character.only=T)

## Source data
source("DataSource_6rivers_oos.R")
df <- dat_oos

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"
PM_Gompertz.col <- "#1C474D"

#################################################
## GPP predicted TS
##################################################

GPP_oos_preds <- function(preds, df, PM.col, PM.title){
  
  simmat_list <- readRDS(preds)
  simmat_list <- simmat_list[names(df)]
  
  # For every day extract median and CI
  median_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
  lower_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
  upper_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)
  
  ## Plot simulated GPP
  dat <- ldply(df, data.frame)
  df_sim <- as.data.frame(cbind(dat$site_name, as.character(dat$date), dat$GPP, median_simmat$X..i.., lower_simmat$X..i.., upper_simmat$X..i..))
  colnames(df_sim) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
  df_sim$Date <- as.POSIXct(as.character(df_sim$Date), format="%Y-%m-%d")
  df_sim[,3:6] <- apply(df_sim[,3:6],2,function(x) as.numeric(as.character(x)))
  
  ## Arrange rivers by river order
  df_sim <- left_join(df_sim, site_info[,c("site_name","short_name")])
  df_sim$short_name <- factor(df_sim$short_name, levels=c("Black Earth Creek, WI",
                                                          "Fatlick Branch, VA",
                                                          "Beaty Creek, OK",
                                                          "Fanno Creek, OR",
                                                          "South Branch Potomac River, WV",
                                                          "San Joaquin River, CA"))
  
  ## plot only a subset of sites
  df_sim <- df_sim[which(df_sim$short_name %in% c("Fanno Creek, OR","South Branch Potomac River, WV")),]
  
  ## Plot
  df_sim_plot <- ggplot(df_sim, aes(Date, GPP))+
    geom_point(size=2, color="black")+
    geom_line(aes(Date, sim_GPP), color=PM.col, size=1.2)+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title=PM.title)+
    geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                fill=PM.col, alpha=0.3, show.legend = FALSE)+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))+
    #coord_cartesian(ylim = c(0,15))+
    facet_wrap(~short_name, scales = "free", ncol = 1)
  
  return(df_sim_plot)
  
}

GPP_TS_Phenom <- GPP_oos_preds("./rds files/Sim_6riv_AR_oos.rds",df, PM_AR.col, "STS")
GPP_TS_Ricker <- GPP_oos_preds("./rds files/Sim_6riv_Ricker_oos.rds",df, PM_Ricker.col, "LB")

GPP_TS_Phenom
GPP_TS_Ricker

###################################
## RMSE Comparison Plot
###################################
RMSE_all <- readRDS("./rds files/RMSE_all_2021_05_17.rds")

RMSE_OOS <- RMSE_all[which(RMSE_all$WS_vs_OOS == "OOS"),]
RMSE_OOS <-left_join(RMSE_OOS, RMSE_prev_post, by="site_name")
RMSE_OOS$quant <- revalue(as.character(RMSE_OOS$quant), replace=c("0.025"="lowerCI",
                                                                  "0.5"="median",
                                                                  "0.975"="upperCI"))

## AR vs. Ricker
oos_biplot_R <- RMSE_OOS[which(RMSE_OOS$model_type == "Ricker"),c("site_name","RMSE","quant","model_type")]
oos_biplot_R <- spread(oos_biplot_R, key = quant, value=RMSE)
colnames(oos_biplot_R) <- c("site_name","model_type","Ricker_lowerCI","Ricker_median","Ricker_upperCI")

oos_biplot_AR <- RMSE_OOS[which(RMSE_OOS$model_type == "AR"),c("site_name","RMSE","quant","model_type")]
oos_biplot_AR <- spread(oos_biplot_AR, key = quant, value=RMSE)
colnames(oos_biplot_AR) <- c("site_name","model_type","AR_lowerCI","AR_median","AR_upperCI")

oos_biplot <- merge(oos_biplot_AR, oos_biplot_R, by="site_name")

oos_biplot$short_name <- revalue(as.character(oos_biplot$site_name), replace = c("nwis_05406457"="Black Earth Creek, WI",
                                                                                 "nwis_01656903"="Fatlick Branch, VA",
                                                                                 "nwis_07191222"="Beaty Creek, OK",
                                                                                 "nwis_14206950"="Fanno Creek, OR",
                                                                                 "nwis_01608500"="South Branch Potomac River, WV",
                                                                                 "nwis_11273400"="San Joaquin River, CA"))


six_RMSE_plot <- ggplot(oos_biplot, aes(AR_median, Ricker_median, color=short_name))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Ricker_lowerCI, ymax=Ricker_upperCI, color=short_name), size=0.5)+
  geom_errorbarh(aes(xmin=AR_lowerCI, xmax=AR_upperCI, color=short_name), size=0.5)+
  coord_cartesian(xlim = c(0,15),y=c(0,5))+
  scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x="Phenomenological model OOS RMSE",
       y="Ricker model OOS RMSE",
       color="Site")+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16),
        legend.position = "bottom")



(GPP_TS_Phenom | GPP_TS_Ricker)/six_RMSE_plot























