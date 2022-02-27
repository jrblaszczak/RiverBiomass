##==============================================================================
## Script for out-of-sample prediction comparison
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
source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"
PM_Gompertz.col <- "#1C474D"

## Check rivers
site_info[,c("site_name","long_name","NHD_STREAMORDE")]

################################
## Model output plot function
################################
GPP_oos_preds <- function(preds, df){

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
  df_sim$short_name <- factor(df_sim$short_name, levels=site_order_list)
  
  return(df_sim)

}

STS_simdat <- GPP_oos_preds("./rds files/Sim_6riv_AR_oos_2022_02_27.rds",df)
LB_simdat <- GPP_oos_preds("./rds files/Sim_6riv_Ricker_oos_2022_02_27.rds",df)

### Overlain
scaleFUN <- function(x) sprintf("%.1f", x)

overlainGPP_plot <- function(y){
  site <- y
  
  STS_simdat_site <- STS_simdat[which(STS_simdat$short_name == site),]
  LB_simdat_site <- LB_simdat[which(LB_simdat$short_name == site),]
  
  plot <- ggplot(STS_simdat_site, aes(Date, GPP))+
    geom_point(size=1.5, color="black")+
    
    geom_line(aes(Date, sim_GPP), color=PM_AR.col, size=1.2)+
    geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                fill=PM_AR.col, alpha=0.2, show.legend = FALSE)+
    
    geom_line(data=LB_simdat_site, aes(Date, sim_GPP), color=PM_Ricker.col, size=1.2)+
    geom_ribbon(data=LB_simdat_site, aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                fill=PM_Ricker.col, alpha=0.2, show.legend = FALSE)+
    
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_blank(), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))+
    #coord_cartesian(ylim=c(0,max(LB_simdat_site$sim_GPP_upper)*1.2))+
    scale_y_continuous(labels=scaleFUN)+
    #labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'))+
    facet_wrap(~short_name, scales = "free", ncol = 2)
  
  return(plot)
  
}

plots <- lapply(site_order_list, function(x) overlainGPP_plot(x))
names(plots) <- site_order_list

grobplots <- lapply(plots, function(x) ggplotGrob(x))

## plot aligned
all_plot <- plot_grid(grobplots$`Proctor Creek, GA`, grobplots$`Paint Branch, MD`,
          grobplots$`Beaty Creek, OK`, grobplots$`S. Br. Potomac River, WV`,
          grobplots$`Santa Margarita River, CA`, grobplots$`Pecos River, TX`,
          ncol = 2, align="hv")

#create common x and y labels
y.grob <- textGrob(expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), 
                   gp=gpar(fontsize=16), rot=90)
x.grob <- textGrob("Date", 
                   gp=gpar(fontsize=16))

#add to plot
grid.arrange(arrangeGrob(all_plot, left = y.grob, bottom = x.grob))


###############################
## RMSE comparison
###############################

## Recalc RMSE for median prediction for each model
STS_simdat_l <- split(STS_simdat, STS_simdat$short_name)
LB_simdat_l <- split(LB_simdat, LB_simdat$short_name)

RMSE_calcs <- function(x){
  rmse_med <- sqrt(sum((x$sim_GPP - x$GPP)^2)/length(x$GPP))
  rmse_lower <- sqrt(sum((x$sim_GPP_lower - x$GPP)^2)/length(x$GPP))
  rmse_upper <- sqrt(sum((x$sim_GPP_upper - x$GPP)^2)/length(x$GPP))
  
  rmse_vals <- as.data.frame(cbind(rmse_med, rmse_lower, rmse_upper))
  return(rmse_vals)
  
}

STS_rmse <- ldply(lapply(STS_simdat_l, function(x) RMSE_calcs(x)), data.frame)
LB_rmse <- ldply(lapply(LB_simdat_l, function(x) RMSE_calcs(x)), data.frame)
STS_rmse$model <- "STS"; LB_rmse$model <- "LB"

rmses <- rbind(STS_rmse, LB_rmse)
rmse_melt <- melt(rmses, .idvars=c("model",".id"))
rmse_melt$.id <- factor(rmse_melt$.id, levels=site_order_list)
rmse_melt$model <- revalue(rmse_melt$model, replace=c("STS"="S-TS","LB"="LB-TS"))
rmse_melt$model <- factor(rmse_melt$model, levels=c("S-TS","LB-TS"))
rmse_melt$variable <- revalue(rmse_melt$variable, replace=c("rmse_med"="Lower CI Pred. RMSE",
                                                            "rmse_lower"="Median Pred. RMSE",
                                                            "rmse_upper"="Upper CI Pred. RMSE"))
rmse_melt$variable <- factor(rmse_melt$variable, levels=c("Lower CI Pred. RMSE",
                                                          "Median Pred. RMSE",
                                                          "Upper CI Pred. RMSE"))
rmse_melt$Site <- revalue(rmse_melt$.id, replace=c("Proctor Creek, GA"="GA (2)",
                                                    "Paint Branch, MD"="MD (2)",
                                                    "Beaty Creek, OK"="OK (3)",
                                                    "S. Br. Potomac River, WV"="WV (5)",
                                                    "Santa Margarita River, CA"="CA (6)",
                                                    "Pecos River, TX"="TX (7)"))


ggplot(rmse_melt, aes(fill=model, y=value, x=Site))+
  geom_bar(stat = "identity",position="dodge",
           color="black",alpha=0.9)+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=10),
        axis.title.y = element_text(size=12), axis.text.x = element_text(angle=35, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=12))+
  labs(y="RMSE",x="Site")+
  facet_wrap(~variable)+
  scale_fill_manual("Model",values = c("S-TS"=PM_AR.col,"LB-TS"=PM_Ricker.col))




###############################
## Latent Biomass predictions
###############################

LatBio_oos_preds <- function(preds, df, PM.col, PM.title){
  
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
  
  ## Arrange rivers by river order
  df_bio <- left_join(df_bio, site_info[,c("site_name","short_name")])
  df_bio$short_name <- factor(df_bio$short_name, levels= site_order_list)
  
  ## plot
  df_bio_plot <- ggplot(df_bio, aes(Date, exp(B)))+
    geom_line(size=1.2, color="chartreuse4")+
    labs(y="Latent Biomass",title=PM.title)+
    geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
                fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))+
    coord_cartesian(ylim=c(0,30))+
    facet_wrap(~short_name, scales = "free_x", ncol = 2)
  
  
  return(df_bio_plot)
  
  
}


LatBio_oos_preds("./rds files/Sim_6riv_Ricker_oos.rds",df, PM_Ricker.col, "PM: LB")
#LatBio_oos_preds("./rds files/Sim_6riv_Gompertz_oos.rds",df, PM_Gompertz.col, "PM: Gompertz")


########################
## 1 to 1 comparison
#######################

ggplot(STS_simdat, aes(GPP, sim_GPP))+
  facet_wrap(~short_name, ncol = 2, scales="free")+
  geom_point(size=2, color=PM_AR.col)+
  geom_abline(slope=1, intercept = 0)+
  geom_point(data=LB_simdat, aes(GPP, sim_GPP), color=PM_Ricker.col)+
  
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_blank(), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  scale_y_continuous(labels=scaleFUN)+
  #labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'))+
  facet_wrap(~short_name, scales = "free", ncol = 2)



