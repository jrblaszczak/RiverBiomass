## Figure comparing within versus out-of-sample predictions

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
source("DataSource_9rivers.R")
dat_ws <- dat

# Within-sample predictions
ws_AR <- readRDS("./rds files/Sim_6riv_AR.rds")
ws_Ricker <- readRDS("./rds files/Sim_6riv_Ricker_ws.rds")

## Import within sample data
source("DataSource_9rivers_oos.R")

# Out-of-sample predictions
oos_AR <- readRDS("./rds files/Sim_6riv_AR_oos.rds")
oos_Ricker <- readRDS("./rds files/Sim_6riv_Ricker_oos.rds")

#########################################
## Extract and summarize predictions
#########################################

## GPP predictions
GPP_matrix_extract <- function(matrix_list, dat, site_info, WS_or_OOS){
  
  # For every day extract median and CI
  median_predmat <- ldply(lapply(matrix_list, function(z) apply(z[[1]], 1, function(x) median(x,na.rm = TRUE))), data.frame)
  lower_predmat <- ldply(lapply(matrix_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025,na.rm = TRUE))), data.frame)
  upper_predmat <- ldply(lapply(matrix_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975,na.rm = TRUE))), data.frame)
  
  ## Plot simulated GPP
  df <- ldply(dat[names(matrix_list)], data.frame)
  GPP_df <- as.data.frame(cbind(df$site_name, as.character(df$date),
                                df$GPP, median_predmat[,2], lower_predmat[,2], upper_predmat[,2]))
  colnames(GPP_df) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
  GPP_df$Date <- as.POSIXct(as.character(GPP_df$Date), format="%Y-%m-%d")
  GPP_df[,3:6] <- apply(GPP_df[,3:6],2,function(x) as.numeric(as.character(x)))
  
  ## Arrange rivers by river order
  GPP_df <- left_join(GPP_df, site_info[,c("site_name","short_name")])
  GPP_df$short_name <- factor(GPP_df$short_name, levels=c("Silver Creek, UT",
                                                            "Medina River, TX",
                                                            "Anacostia River, MD",
                                                            "West Fork River, WV",
                                                            "St. John's River, FL",
                                                            "Clackamas River, OR"))
  
  ## Specify whether within or out-of-sample predictions
  GPP_df$WS_or_OOS <- WS_or_OOS
  
  return(GPP_df)
  
}

## Latent biomass predictions
Bio_matrix_extract <- function(matrix_list, dat, site_info, WS_or_OOS){
  
  # For every day extract median and CI
  median_predmat <- ldply(lapply(matrix_list, function(z) apply(z[[3]], 1, function(x) median(x,na.rm = TRUE))), data.frame)
  lower_predmat <- ldply(lapply(matrix_list, function(z) apply(z[[3]], 1, function(x) quantile(x, probs = 0.025,na.rm = TRUE))), data.frame)
  upper_predmat <- ldply(lapply(matrix_list, function(z) apply(z[[3]], 1, function(x) quantile(x, probs = 0.975,na.rm = TRUE))), data.frame)
  
  ## Extract latent biomass
  df <- ldply(dat[names(matrix_list)], data.frame)
  LB_df <- as.data.frame(cbind(df$site_name, as.character(df$date),
                                df$GPP, median_predmat[,2], lower_predmat[,2], upper_predmat[,2]))
  colnames(LB_df) <- c("site_name","Date","GPP","Biomass","Biomass_lower","Biomass_upper")
  LB_df$Date <- as.POSIXct(as.character(LB_df$Date), format="%Y-%m-%d")
  LB_df[,3:6] <- apply(LB_df[,3:6],2,function(x) as.numeric(as.character(x)))
  
  ## Arrange rivers by river order
  LB_df <- left_join(LB_df, site_info[,c("site_name","short_name")])
  LB_df$short_name <- factor(LB_df$short_name, levels=c("Silver Creek, UT",
                                                          "Medina River, TX",
                                                          "Anacostia River, MD",
                                                          "West Fork River, WV",
                                                          "St. John's River, FL",
                                                          "Clackamas River, OR"))
  
  ## Specify whether within or out-of-sample predictions
  LB_df$WS_or_OOS <- WS_or_OOS
  
  return(LB_df)
  
}

## Summarize across models
# GPP - AR
GPP_ws_AR <- GPP_matrix_extract(ws_AR, dat_ws, site_info,"WS")
GPP_oos_AR <- GPP_matrix_extract(oos_AR, dat_oos, site_info,"OOS")
GPP_comb_AR <- rbind(GPP_ws_AR, GPP_oos_AR)

# GPP - Ricker
GPP_ws_Ricker <- GPP_matrix_extract(ws_Ricker, dat_ws, site_info,"WS")
GPP_oos_Ricker <- GPP_matrix_extract(oos_Ricker, dat_oos, site_info,"OOS")
GPP_comb_Ricker <- rbind(GPP_ws_Ricker, GPP_oos_Ricker)
# Latent Biomass - Ricker
LB_ws_Ricker <- Bio_matrix_extract(ws_Ricker, dat_ws, site_info,"WS")
LB_oos_Ricker <- Bio_matrix_extract(oos_Ricker, dat_oos, site_info,"OOS")
LB_comb_Ricker <- rbind(LB_ws_Ricker, LB_oos_Ricker)

######################
## Plot together
######################
AR_list <- split(GPP_comb_AR, GPP_comb_AR$site_name)
Ricker_list <- split(GPP_comb_Ricker, GPP_comb_Ricker$site_name)
Ricker_LB_list <- split(LB_comb_Ricker, LB_comb_Ricker$site_name)


AR_vs_Ricker_plot <- function(site, AR_list, Ricker_list, Ricker_LB_list, dat_ws, dat_oos){
  
  ## predictions
  AR <- AR_list[[site]]
  Ricker <- Ricker_list[[site]]
  Ricker_LB <- Ricker_LB_list[[site]]
  ## data
  data_ws <- dat_ws[[site]]
  data_oos <- dat_oos[[site]]
  data <- rbind(data_ws, data_oos)
  
  ## Plot
  AR_plot <- ggplot(AR, aes(Date, GPP))+
    geom_point(size=1.5, color="black")+
    geom_line(aes(Date, sim_GPP), color=PM_AR.col, size=1)+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM1: AR")+
    geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                fill=PM_AR.col, alpha=0.2, show.legend = FALSE)+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=12), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15), plot.title = element_text(hjust = 1))+
    scale_y_continuous(limits=c(0,max(AR$sim_GPP_upper, na.rm = TRUE)))+
    geom_vline(xintercept = AR[which(AR$WS_or_OOS == "OOS"),]$Date[1],
               size=1, linetype="dashed")
  
  Ricker_plot <- ggplot(Ricker, aes(Date, GPP))+
    geom_point(size=1.5, color="black")+
    geom_line(aes(Date, sim_GPP), color=PM_Ricker.col, size=1)+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM2: Ricker")+
    geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                fill=PM_Ricker.col, alpha=0.2, show.legend = FALSE)+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=12), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15), plot.title = element_text(hjust = 1))+
    scale_y_continuous(limits=c(0,max(AR$sim_GPP_upper, na.rm = TRUE)))+
    geom_vline(xintercept = Ricker[which(Ricker$WS_or_OOS == "OOS"),]$Date[1],
               size=1, linetype="dashed")
  
  Ricker_LB_plot <- ggplot(Ricker_LB, aes(Date, exp(Biomass)))+
    geom_line(size=1.5, color="chartreuse4")+
    labs(y="Latent Biomass", title="PM2: Ricker Latent Biomass")+
    geom_ribbon(aes(ymin=exp(Biomass_lower),ymax=exp(Biomass_upper)),
                fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=12), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15), plot.title = element_text(hjust = 1))+
    geom_vline(xintercept = Ricker_LB[which(Ricker_LB$WS_or_OOS == "OOS"),]$Date[1],
               size=1, linetype="dashed")
  
  ratio_QL <- max(data$light)/max(data$Q)
  
  data_plot <- ggplot(data, aes(date, Q*ratio_QL))+
    geom_point(data=data, aes(date, light), size=1.5, color="darkgoldenrod3")+
    geom_line(size=1, color="deepskyblue4")+
    scale_y_continuous(sec.axis = sec_axis(~./ratio_QL, name=expression("Daily Q (cms)")))+
    labs(y=expression('Daily PPFD'), title="External Drivers")+# ('*~mu~mol~ m^-2~d^-1*')'), x="Date")+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y.left = element_text(size=12, color="darkgoldenrod3"),
          axis.title.y.right = element_text(size=12, color="deepskyblue4"),
          axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15), plot.title = element_text(hjust = 1))+
    geom_vline(xintercept = Ricker_LB[which(Ricker_LB$WS_or_OOS == "OOS"),]$Date[1],
               size=1, linetype="dashed")
  
  (AR_plot+theme(axis.text.x = element_blank())) + 
    (Ricker_plot+theme(axis.text.x = element_blank())) +
    (Ricker_LB_plot+theme(axis.text.x = element_blank())) +
    data_plot + plot_layout(ncol=1) + plot_annotation(title=AR$short_name[1])

}

# test plot
AR_vs_Ricker_plot(site = names(AR_list)[5], AR_list, Ricker_list, Ricker_LB_list, dat_ws, dat_oos)

# save figures
setwd("~/GitHub/RiverBiomass/figures/WS OOS comparison by site")

for(i in 1:length(names(AR_list))){
  ggsave(plot = AR_vs_Ricker_plot(site = names(AR_list)[i], AR_list, Ricker_list, Ricker_LB_list, dat_ws, dat_oos),
         filename = paste(names(AR_list)[i],"_ws_oos_TS.jpg",sep = ""))
}










