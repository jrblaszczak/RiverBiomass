## RMSE comparison within and out-of-sample

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
ws_AR <- readRDS("./rds files/Sim_6riv_AR_ws_2022_01_23.rds")
ws_Ricker <- readRDS("./rds files/Sim_6riv_Ricker_ws_2022_01_23.rds")

## Import within sample data
source("DataSource_6rivers_oos_StreamLight.R")

# Out-of-sample predictions
oos_AR <- readRDS("./rds files/Sim_6riv_AR_oos_2022_02_01.rds")
oos_Ricker <- readRDS("./rds files/Sim_6riv_Ricker_oos_2022_02_01.rds")

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
RMSE_ws_AR <- RMSE_extract(ws_AR,"WS","AR")
RMSE_oos_AR <- RMSE_extract(oos_AR,"OOS","AR")
RMSE_ws_Ricker <- RMSE_extract(ws_Ricker,"WS","Ricker")
RMSE_oos_Ricker <- RMSE_extract(oos_Ricker,"OOS","Ricker")

RMSE_all <- rbind(RMSE_ws_AR, RMSE_oos_AR, RMSE_ws_Ricker, RMSE_oos_Ricker)

saveRDS(RMSE_all, "./rds files/RMSE_all_2022_02_01.rds")

###########################################################################
## Figure out RMSE if just using last year's data to predict this year
##########################################################################

# Previous year GPP
prev_year_GPP <- ldply(lapply(dat_ws[names(oos_AR)], function(x) return(x[,c("date","doy","GPP")])), data.frame)
colnames(prev_year_GPP) <- c("site_name","Date","DOY","Prev_yr_GPP")

# Following year GPP
post_year_GPP <- ldply(lapply(dat_oos[names(oos_AR)], function(x) return(x[,c("date","doy","GPP")])), data.frame)
colnames(post_year_GPP) <- c("site_name","Date","DOY","Post_yr_GPP")

# Create a DOY template to then join both years onto
doy_vec <- seq(from=1, to=366, by=1)

site_name_list <- list()
for(i in 1:length(names(oos_AR))){
  site_name_list[[i]] <- rep(names(oos_AR)[i], times=366)
}
site_name_df <- ldply(site_name_list, data.frame)
site_name_df <- as.data.frame(cbind(site_name_df[,1], rep(doy_vec, times=length(names(oos_AR)))))
colnames(site_name_df) <- c("site_name","DOY")
site_name_df$DOY <- as.integer(site_name_df$DOY)

# join prev and post years
joined <- left_join(site_name_df, prev_year_GPP, by=c("site_name","DOY"))
joined <- left_join(joined, post_year_GPP, by=c("site_name","DOY"))

joined_list <- split(joined, joined$site_name)

## Calculate RMSE
RMSE_prev_yr_fxn <- function(x){
  
  x <- na.omit(x)

  rmse_val <- sqrt(sum((x$Post_yr_GPP-x$Prev_yr_GPP)^2)/length(x$Prev_yr_GPP))

  return(rmse_val)
}

RMSE_prev_post <- ldply(lapply(joined_list, function(x) RMSE_prev_yr_fxn(x)), data.frame)
colnames(RMSE_prev_post) <- c("site_name","RMSE_among_yrs")

###############################
## Visualize
###############################
RMSE_OOS <- RMSE_all[which(RMSE_all$WS_vs_OOS == "OOS"),]
RMSE_OOS <-left_join(RMSE_OOS, RMSE_prev_post, by="site_name")
RMSE_OOS$quant <- revalue(as.character(RMSE_OOS$quant), replace=c("0.025"="lowerCI",
                                                                  "0.5"="median",
                                                                  "0.975"="upperCI"))
RMSE_OOS_wide <- spread(RMSE_OOS, key = c("quant"), value="RMSE")

## merge RMSE wide with sitename
RMSE_OOS_wide <- merge(RMSE_OOS_wide, site_info[,c("site_name","short_name")])



OOS_prev_yr_plot <-ggplot(RMSE_OOS_wide, aes(RMSE_among_yrs, median, color=short_name))+#, shape=model_type))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI, color=short_name))+#, group=model_type))+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
    axis.text = element_text(size=12),
        axis.title = element_text(size=14), 
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  geom_abline(slope=1, intercept = 0)+
  facet_wrap(~model_type)+
  labs(x="Previous year prediction RMSE",y="Out-of-sample prediction RMSE")+
  coord_cartesian(ylim = c(0,12))
OOS_prev_yr_plot


## Within Sample
RMSE_WS <- RMSE_all[which(RMSE_all$WS_vs_OOS == "WS"),]
RMSE_WS <-left_join(RMSE_WS, RMSE_prev_post, by="site_name")
RMSE_WS$quant <- revalue(as.character(RMSE_WS$quant), replace=c("0.025"="lowerCI",
                                                                  "0.5"="median",
                                                                  "0.975"="upperCI"))
RMSE_WS_wide <- spread(RMSE_WS, key = c("quant"), value="RMSE")
## merge RMSE wide with sitename
RMSE_WS_wide <- merge(RMSE_WS_wide, site_info[,c("site_name","short_name")])

WS_plot <- ggplot(RMSE_WS_wide, aes(short_name, median))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI, color=short_name))+#, group=model_type))+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size=14), 
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15),
        legend.position = "none")+
  #geom_abline(slope=1, intercept = 0)+
  facet_wrap(~model_type)+
  labs(x="Previous year prediction RMSE",y="Within-sample prediction RMSE")+
  coord_cartesian(ylim = c(0,20))

plot_grid(WS_plot, OOS_prev_yr_plot+theme(legend.position = "none"), align = "hv", nrow=2)

################################
## AR vs Ricker OOS comparison
################################

oos_biplot_R <- RMSE_OOS[which(RMSE_OOS$model_type == "Ricker"),c("site_name","RMSE","quant","model_type")]
oos_biplot_R <- spread(oos_biplot_R, key = quant, value=RMSE)
colnames(oos_biplot_R) <- c("site_name","model_type","Ricker_lowerCI","Ricker_median","Ricker_upperCI")

oos_biplot_AR <- RMSE_OOS[which(RMSE_OOS$model_type == "AR"),c("site_name","RMSE","quant","model_type")]
oos_biplot_AR <- spread(oos_biplot_AR, key = quant, value=RMSE)
colnames(oos_biplot_AR) <- c("site_name","model_type","AR_lowerCI","AR_median","AR_upperCI")

oos_biplot <- merge(oos_biplot_AR, oos_biplot_R, by="site_name")
oos_biplot <- merge(oos_biplot, site_info, by="site_name")

six_RMSE_plot <- ggplot(oos_biplot, aes(AR_median, Ricker_median, color=short_name))+
  geom_abline(slope = 1, intercept = 0)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=Ricker_lowerCI, ymax=Ricker_upperCI, color=short_name), size=0.5)+
  geom_errorbarh(aes(xmin=AR_lowerCI, xmax=AR_upperCI, color=short_name), size=0.5)+
  coord_cartesian(xlim = c(0,12),y=c(0,12))+
  scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x="Phenomenological model OOS RMSE",
       y="Ricker model OOS RMSE",
       color="Site")+
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16))
  

six_RMSE_plot

















