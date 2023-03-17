##==============================================================================
## Script for out-of-sample LB-TS prediction and 1:1 plots
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

#STS_simdat <- GPP_oos_preds_ts("./rds files/Sim_6riv_AR_oos_2022_02_27.rds",df)
LB_simdat <- GPP_oos_preds_ts("./rds files/Sim_6riv_Ricker_oos_2022_02_27.rds",df)

xy.limits <- range(c(LB_simdat[which(LB_simdat$short_name == "Proctor Creek, GA"),]$GPP,
                     LB_simdat[which(LB_simdat$short_name == "Proctor Creek, GA"),]$sim_GPP))

ggplot(LB_simdat[which(LB_simdat$short_name == "Proctor Creek, GA"),], aes(GPP, sim_GPP))+
  geom_point()+
  scale_x_continuous(limits=c(0,xy.limits[2])) + 
  scale_y_continuous(limits=c(0,xy.limits[2])) +
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()



#############################
## Visualize overlain 
#############################
scaleFUN <- function(x) sprintf("%.1f", x)

overlainGPP_plot <- function(y){
  site <- y
  
  LB_simdat_site <- LB_simdat[which(LB_simdat$short_name == site),]
  
  xy.limits <- range(c(LB_simdat_site$GPP,
                       LB_simdat_site$sim_GPP))
  
  inset_plot <- ggplot(LB_simdat_site, aes(GPP, sim_GPP))+
    geom_point(color = PM_Ricker.col, size=0.5)+
    scale_x_continuous(limits=c(0,xy.limits[2])) + 
    scale_y_continuous(limits=c(0,xy.limits[2])) +
    labs(x="GPP Data",y="Predicted\n GPP")+
    geom_abline(slope = 1, intercept = 0)+
    theme_classic()+
    theme(axis.title = element_text(size = 8),
          axis.text = element_text(size = 7))
    
  
  main_plot <- ggplot(LB_simdat_site, aes(Date, GPP))+
    geom_point(size=1, color="black")+
    
    geom_line(data=LB_simdat_site, aes(Date, sim_GPP), color=PM_Ricker.col, size=1)+
    geom_ribbon(data=LB_simdat_site, aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                fill=PM_Ricker.col, alpha=0.4, show.legend = FALSE)+
    
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
    draw_plot(inset_plot, x=0.58, y=0.52, width=0.35, height = 0.35)

  return(comb.plot)
  
}


plots <- lapply(site_order_list, function(x) overlainGPP_plot(x))
names(plots) <- site_order_list
#plots$`Proctor Creek, GA`
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




