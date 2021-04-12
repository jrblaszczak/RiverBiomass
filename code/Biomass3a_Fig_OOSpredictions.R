## Out-of-sample predictions

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

## Source data
source("DataSource_9rivers_oos.R")
df <- dat_oos
# Subset source data
df <- df[c("nwis_01649500","nwis_02234000","nwis_03058000",
           "nwis_08180700","nwis_10129900","nwis_14211010")]

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"
PM_Gompertz.col <- "#1C474D"

## Change river names to short names
site_info[,c("site_name","long_name","NHD_STREAMORDE")]
site_info <- site_info[which(site_info$site_name %in% names(df)),]
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_01649500"="Anacostia River, MD",
                                                                               "nwis_02234000"="St. John's River, FL",
                                                                               "nwis_03058000"="West Fork River, WV",
                                                                               "nwis_08180700"="Medina River, TX",
                                                                               "nwis_10129900"="Silver Creek, UT",
                                                                               "nwis_14211010"="Clackamas River, OR"))


################################
## Model output plot function
################################
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
  df_sim$short_name <- factor(df_sim$short_name, levels=c("Silver Creek, UT",
                                                            "Medina River, TX",
                                                            "Anacostia River, MD",
                                                            "West Fork River, WV",
                                                            "St. John's River, FL",
                                                            "Clackamas River, OR"))
  
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
    coord_cartesian(ylim = c(0,15))+
    facet_wrap(~short_name, scales = "free_x", ncol = 2)
  
  return(df_sim_plot)

}

GPP_oos_preds("./rds files/Sim_9riv_AR_oos.rds",df, PM_AR.col, "PM: Phenomenological")
GPP_oos_preds("./rds files/Sim_9riv_Ricker_oos.rds",df, PM_Ricker.col, "PM: Ricker")
#GPP_oos_preds("./rds files/Sim_9riv_Gompertz_oos.rds",df, PM_Gompertz.col, "PM: Gompertz")



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
  df_bio$short_name <- factor(df_bio$short_name, levels=c("Silver Creek, UT",
                                                            "Medina River, TX",
                                                            "Anacostia River, MD",
                                                            "West Fork River, WV",
                                                            "St. John's River, FL",
                                                            "Clackamas River, OR"))
  
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
    coord_cartesian(ylim=c(0,20))+
    facet_wrap(~short_name, scales = "free_x", ncol = 2)
  
  
  return(df_bio_plot)
  
  
}


LatBio_oos_preds("./rds files/Sim_9riv_Ricker_oos.rds",df, PM_Ricker.col, "PM: Ricker")
LatBio_oos_preds("./rds files/Sim_9riv_Gompertz_oos.rds",df, PM_Gompertz.col, "PM: Gompertz")


