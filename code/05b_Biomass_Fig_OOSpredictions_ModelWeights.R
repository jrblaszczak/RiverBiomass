##==============================================================================
## Script for out-of-sample prediction model weight comparison
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
## Predicted time series
################################
GPP_oos_preds_ts <- function(preds, df){

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

STS_simdat <- GPP_oos_preds_ts("./rds files/Sim_6riv_AR_oos_2022_02_27.rds",df)
LB_simdat <- GPP_oos_preds_ts("./rds files/Sim_6riv_Ricker_oos_2022_02_27.rds",df)


#############################################
## Predicted time series for model weights
#############################################
GPP_oos_preds_mw <- function(preds, df, mean_mod, se_mod){
  
  simmat_list <- readRDS(preds)
  simmat_list <- simmat_list[names(df)]
  
  # For every day extract median and CI
  mean_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) mean(x))), data.frame)
  se_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) sd(x))), data.frame)
  
  ## Plot simulated GPP
  dat <- ldply(df, data.frame)
  df_sim <- as.data.frame(cbind(dat$site_name, as.character(dat$date), 
                                dat$GPP, dat$GPP_sd, 
                                mean_simmat$X..i.., se_simmat$X..i..))
  colnames(df_sim) <- c("site_name","Date","GPP","GPP_se", mean_mod, se_mod)
  df_sim$Date <- as.POSIXct(as.character(df_sim$Date), format="%Y-%m-%d")
  df_sim[,3:6] <- apply(df_sim[,3:6],2,function(x) as.numeric(as.character(x)))
  df_sim[,3:6] <- log(df_sim[,3:6])
  
  ## Link to site_info
  df_sim <- left_join(df_sim, site_info[,c("site_name","short_name")])
  
  return(df_sim)
  
}

STS_oos_mw <- GPP_oos_preds_mw("./rds files/Sim_6riv_AR_oos_2022_02_27.rds", df, "mean_p_GPP_STS", "se_p_GPP_STS")
LB_oos_mw <- GPP_oos_preds_mw("./rds files/Sim_6riv_Ricker_oos_2022_02_27.rds", df, "mean_p_GPP_LB", "se_p_GPP_LB")

## Combine
oos_preds_mw <- merge(STS_oos_mw, LB_oos_mw[,c("site_name","Date","mean_p_GPP_LB", "se_p_GPP_LB")])

## Split by short_name
oos_pl <- split(oos_preds_mw, oos_preds_mw$short_name)

#################################
## Calculate model weights
##################################

support_calc <- function(short.name){
  
  test <- oos_pl[[short.name]]
  test$p_obs_STS <- dnorm(test$GPP, test$mean_p_GPP_STS, sqrt(test$GPP_se^2 + test$se_p_GPP_STS^2), log = FALSE)
  test$p_obs_LB <- dnorm(test$GPP, test$mean_p_GPP_LB, sqrt(test$GPP_se^2 + test$se_p_GPP_LB^2), log = FALSE)
  
  test2 <- test[-1,]
  
  #probability vectors
  p_STS <- c(0,test2$p_obs_STS); p_LB <- c(0, test2$p_obs_LB)
  #empty vectors
  weight_STS <- numeric(nrow(test2)+1)
  weight_LB <- numeric(nrow(test2)+1)
  #starting values
  weight_STS[1] <- 0.5; weight_LB[1] <- 0.5
  
  for(i in 2:length(weight_STS)){
    weight_STS[i] <- (weight_STS[(i-1)]*p_STS[i])/((weight_STS[(i-1)]*p_STS[i])+(weight_LB[(i-1)]*p_LB[i]))
    weight_LB[i] <- 1-weight_STS[i]
  }
  support <- cbind(test2,"STS_support" = weight_STS[-1], "LB_support" = weight_LB[-1])
  return(support)
}

#############################
## Visualize overlain 
#############################
scaleFUN <- function(x) sprintf("%.1f", x)

overlainGPP_plot <- function(y){
  site <- y
  
  STS_simdat_site <- STS_simdat[which(STS_simdat$short_name == site),]
  LB_simdat_site <- LB_simdat[which(LB_simdat$short_name == site),]
  
  support <- support_calc(site)
  
  inset_plot <- ggplot(support, aes(Date, STS_support))+
    geom_line(color=PM_AR.col,size=0.9)+
    geom_line(aes(Date, LB_support),color=PM_Ricker.col,size=0.9)+
    labs(y="Model Weight")+
    theme_bw(base_size = 8)+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())+
    scale_y_continuous(labels=scaleFUN)
  
  main_plot <- ggplot(STS_simdat_site, aes(Date, GPP))+
    geom_point(size=1, color="black")+
    
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
    coord_cartesian(ylim=c(0,max(LB_simdat_site$sim_GPP_upper)*2))+
    scale_y_continuous(labels=scaleFUN)+
    #labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'))+
    facet_wrap(~short_name, scales = "free", ncol = 2)
  
  comb.plot <- ggdraw()+
    draw_plot(main_plot)+
    draw_plot(inset_plot, x=0.55, y=0.6, width=0.4, height = 0.25)

  return(comb.plot)
  
}

plots <- lapply(site_order_list, function(x) overlainGPP_plot(x))
names(plots) <- site_order_list
plots$`Proctor Creek, GA`
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




