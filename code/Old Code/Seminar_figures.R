## Figure comparing within versus out-of-sample predictions

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","patchwork"), require, character.only=T)

PM_AR.col <- "#d95f02" # AR
PM_Ricker.col <- "#7654B4"#"#1C474D" # Ricker
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

## Import out-of-sample data
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
  
  AR <- AR[which(AR$WS_or_OOS == "OOS"),]
  Ricker <- Ricker[which(Ricker$WS_or_OOS == "OOS"),]
  
  ## Plot
  AR_plot <- ggplot(AR, aes(Date, GPP))+
    geom_point(size=1.5, color="black")+
    geom_line(aes(Date, sim_GPP), color=PM_AR.col, size=1.5)+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM1: AR")+
    geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                fill=PM_AR.col, alpha=0.2, show.legend = FALSE)+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=12), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15), plot.title = element_text(hjust = 1))+
    scale_y_continuous(limits=c(0,max(AR$sim_GPP_upper, na.rm = TRUE)))
  
  Ricker_plot <- ggplot(Ricker, aes(Date, GPP))+
    geom_point(size=1.5, color="black")+
    geom_line(aes(Date, sim_GPP), color=PM_Ricker.col, size=1.5)+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM2: Ricker")+
    geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                fill=PM_Ricker.col, alpha=0.2, show.legend = FALSE)+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.text = element_text(size=12),
          axis.title.y = element_text(size=12), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15), plot.title = element_text(hjust = 1))+
    scale_y_continuous(limits=c(0,max(AR$sim_GPP_upper, na.rm = TRUE)))
  
  
  (AR_plot+theme(axis.text.x = element_blank())) + 
    (Ricker_plot) +
    plot_layout(ncol=1)
  
}

# test plot
AR_vs_Ricker_plot(site = names(AR_list)[1], AR_list, Ricker_list, Ricker_LB_list, dat_ws, dat_oos)


##################################
## Extract and calculate RMSE
##################################

## extract RMSE summary stats
RMSE_extract <- function(pred_rmse_list,WS_or_OOS,model.type){
  
  rmsemat <- ldply(lapply(pred_rmse_list, function(x) return(x[[2]])), data.frame)
  colnames(rmsemat) <- c("site_name","RMSE")
  
  rmsemat_sum <- rmsemat %>%
    group_by(site_name) %>%
    summarise(RMSE = quantile(RMSE, c(0.025, 0.5, 0.975)), quant = c(0.025, 0.5, 0.975))
  
  # specify whether within or out-of-sample predictions
  rmsemat_sum$WS_vs_OOS <- WS_or_OOS
  
  # specify model
  rmsemat_sum$model_type <- model.type
  
  return(rmsemat_sum)
  
}

## Summarize for each model
#RMSE_ws_AR <- RMSE_extract(ws_AR,"WS","AR")
RMSE_oos_AR <- RMSE_extract(oos_AR,"OOS","AR")
#RMSE_ws_Ricker <- RMSE_extract(ws_Ricker,"WS","Ricker")
RMSE_oos_Ricker <- RMSE_extract(oos_Ricker,"OOS","Ricker")

#RMSE_all <- rbind(RMSE_ws_AR, RMSE_oos_AR, RMSE_ws_Ricker, RMSE_oos_Ricker)

combined <- rbind(RMSE_oos_AR, RMSE_oos_Ricker)

combined$quant <- revalue(as.character(combined$quant), replace=c("0.025"="lowerCI",
                                                                  "0.5"="median",
                                                                  "0.975"="upperCI"))
RMSE_OOS_wide <- spread(combined, key = c("quant"), value="RMSE")

OR <- RMSE_OOS_wide[which(RMSE_OOS_wide$site_name == "nwis_14211010"),]

ggplot(OR, aes(as.factor(model_type), median, color=as.factor(model_type)))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI), width=0.2, size=1)+
  labs(x="Model",y="Out-of-sample RMSE")+
  scale_color_manual(values=c("AR" = PM_AR.col, "Ricker"= PM_Ricker.col))+
  theme_bw()+theme(legend.position = "none")



