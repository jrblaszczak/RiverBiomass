##==============================================================================
## Script for persistence curve plots
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork","grid","gridExtra"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p
#source("Simulated_ProductivityModel5_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"
PM_Gompertz.col <- "#1C474D"

##################################################
## Extract 
###################################################
## Import stan fits - simulate one at a time
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker_2022_02_27.rds")
#2nd yr
#stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_2ndYr_output_Ricker_2021_06_15.rds")
#source("DataSource_6rivers_2ndYr_StreamLight.R")

## Extract and summarize parameters
par_Ricker <- lapply(stan_model_output_Ricker, function(x) rstan::extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o")))
#par_Gompertz <- lapply(stan_model_output_Gompertz, function(x) rstan::extract(x, c("beta_0","beta_1","s","c","B","P","pred_GPP","sig_p")))


##################################
## Persistence plots
##################################

## mean parameter
par_mean <- function(par) {
  ## Find the mean
  mean_par <- lapply(par, function(x) mean(x))
  
  mean_pred_GPP_ts <- apply(par$pred_GPP,2,mean)
  mean_P_ts <- apply(par$P,2,mean)
  mean_B_ts <- apply(par$B,2,mean)
  
  ## Compile in list and return
  mean_par_ts <- list(mean_par, mean_pred_GPP_ts, mean_P_ts, mean_B_ts)
  names(mean_par_ts) <- c("par","pred_GPP","P","B")
  return(mean_par_ts)
}

meanpar_R <- lapply(par_Ricker, function(x) par_mean(x))


## Plot persistence
persistence_list <- function(y, data){
  Ppars <- list()
  
  for(i in 1:length(y)){
    df <- y[[i]]
    dat <- data[[i]]
    Ppars[[i]] <- list("tQ"=dat$tQ,"range"=range(dat$tQ),"c"=df$c,"s"=df$s, "site_name"=dat$site_name[1])
    }
  
  names(Ppars) <- names(y)
  
  return(Ppars)

}

P_R <- persistence_list(par_Ricker, df)
#P_G <- persistence_list(par_Gompertz, df)

## plot
plotting_P_dat <- function(x){
  pq <- seq(x$range[1],x$range[2], length=length(x$s))
  p_median <- numeric()
  p_up <- numeric()
  p_down <- numeric()
  name <- substring(deparse(substitute(x)),7)
  for(i in 1:length(pq)){
    t <- exp(-exp(x$s*100*(pq[i]-x$c)))
    p_median[i] <- median(t)
    p_up[i] <- quantile(t, probs = 0.975)
    p_down[i] <- quantile(t, probs = 0.025)
  }
  df <- as.data.frame(as.matrix(cbind(pq,p_median, p_up, p_down)))
  df$site_name <- x$site_name
  
  return(df)
}

P_dat_R <- ldply(lapply(P_R, function(z) plotting_P_dat(z)), data.frame); P_dat_R$PM <- "Ricker"

P_df <- P_dat_R

#####################
## Visualize & compare to bankfull discharge
#####################

## Import and merge bankfull discharge with site info
RI_2 <- read.csv("../data/RI_2yr_flood_6riv.csv", header=T)
sapply(RI_2, class)
site_info <- merge(site_info, RI_2, by="site_name")
## merge with site_info, convert 2 year flood to cms (divide by 35.314666212661), and calculate differences
#crits <- merge(site_info[,c("site_name","RI_2yr_Q","short_name")],Qc_all, by="site_name")
#crits$Q_2yrRI_cms <- crits$RI_2yr_Q/35.314666212661



## join by river name
P_df <- left_join(P_df, site_info[,c("site_name","short_name")], by="site_name")
P_df$short_name <- factor(P_df$short_name, levels= site_order_list)


Persistence_plots <- function(site, df, site_info, P_df){

  Q_sub <- df[[site]]
  #Q_sub$pq <- Q_sub$tQ
  Q_sub$p_for_q <- 0.5
  
  ## critical Q based on velocity
  crit_Q <- site_info[which(site_info$site_name == site),]$RI_2yr_Q
  crit_Q <- crit_Q/35.314666212661
  
  ## convert relativized Q to original values
  P <- P_df[which(P_df$site_name == site),]
  P$Q <- P$pq*max(Q_sub$Q, na.rm = T)
  
  ## critical Q based on GPP - Q correction needed
  c <- meanpar_R[[site]]$par$c*max(Q_sub$Q, na.rm = T)
  
  scaleFUN <- function(x) sprintf("%.1f", x)
  
  ## Plot
  Persist_plot <- ggplot(P, aes(Q, p_median))+
    scale_x_continuous(trans = "log", labels = scaleFUN)+
    geom_point(data=Q_sub, aes(Q, p_for_q), color="white")+
    geom_line(size=1.5, alpha=0.9, color="chartreuse4")+
    geom_ribbon(data=P, aes(ymin=p_down, ymax=p_up), alpha=0.3, fill="chartreuse4", color=NA)+
    theme(panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12, angle=45, hjust=1),
          axis.title = element_blank(), 
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))+
    #labs(x="Range of Standardized Discharge",y="Persistence")+
    scale_y_continuous(limits=c(0,1), labels = function(x) paste0(x*100, '%'))+
    geom_vline(xintercept = crit_Q, size=1, linetype="dotted", color="grey25")+
    geom_vline(xintercept = c, size=1, linetype="dashed")
  
  
  Persist_plot2 <- ggExtra::ggMarginal(Persist_plot, data=Q_sub, type="histogram",
                                       size=4, x = Q, margins = "x", color="black",
                                       fill="deepskyblue4", xparams = list(alpha=0.8))
  
  return(Persist_plot2)
  
}

site_list <- levels(as.factor(site_info$site_name))

plots <- lapply(site_list, function(x) Persistence_plots(x,df,site_info,P_df))
names(plots) <- site_list

## order based on river order
grid.arrange(
  arrangeGrob(grobs=list(plots$nwis_02336526,
                         plots$nwis_01649190,
                         plots$nwis_07191222,
                         plots$nwis_01608500,
                         plots$nwis_11044000,
                         plots$nwis_08447300),
              ncol = 2,
              bottom=textGrob(expression('Discharge ('*m^3~s^-1*')'), gp=gpar(fontsize=14)), 
              left=textGrob("Day-to-Day Biomass Persistence", gp=gpar(fontsize=14), rot=90)),
  widths=c(9,1)
)





################################################
## View distributions of flow across all sites
###################################################
test <- ldply(df, data.frame)
test <- left_join(test, site_info[,c("site_name","short_name")], by="site_name")
test$short_name <- factor(test$short_name, levels= site_order_list)

library(viridis)

ggplot(test, aes(Q, fill = short_name))+
  geom_density(color="black")+
  scale_fill_viridis("Site", discrete = TRUE, alpha=0.6, option="A") +
  scale_x_continuous(trans="log")+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  theme_bw()


#########################
## S vs. C plots
##########################
## Extract s versus c
sc_R <- ldply(lapply(par_Ricker, function(x) cbind(x$s, x$c)), data.frame); colnames(sc_R) <- c("site_name","s","c")
#sc_G <- ldply(lapply(par_Gompertz, function(x) cbind(x$s, x$c)), data.frame); colnames(sc_G) <- c("site_name","s","c")

## Plot s versus c
sc_plotting <- function(z, t){
  
  ## join by river name
  z <- left_join(z, site_info[,c("site_name","short_name")])
  z$short_name <- factor(z$short_name, levels= site_order_list)
  ggplot(z, aes(s, c))+geom_point()+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))+
    facet_wrap(~short_name, ncol = 2, scales = "free_y")+
    labs(title = t)
  
}

sc_plotting(sc_R, "LB-TS Productivity Model s and c Estimates")
#sc_plotting(sc_G, "PM4: Gompertz")








