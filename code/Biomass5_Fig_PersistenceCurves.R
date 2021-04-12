## Figure - Persistence Curves

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork","grid","gridExtra"), require, character.only=T)

## Source data
source("DataSource_9rivers.R")
# Subset source data
df <- df[c("nwis_01649500","nwis_02234000","nwis_03058000",
           "nwis_08180700","nwis_10129900","nwis_14211010")]
## Change river names to short names
site_info[,c("site_name","long_name","NHD_STREAMORDE")]
site_info <- site_info[which(site_info$site_name %in% names(df)),]
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_01649500"="Anacostia River, MD",
                                                                               "nwis_02234000"="St. John's River, FL",
                                                                               "nwis_03058000"="West Fork River, WV",
                                                                               "nwis_08180700"="Medina River, TX",
                                                                               "nwis_10129900"="Silver Creek, UT",
                                                                               "nwis_14211010"="Clackamas River, OR"))


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
stan_model_output_Ricker <- readRDS("./rds files/stan_9riv_output_Ricker_2021_03_05.rds")
stan_model_output_Ricker <- stan_model_output_Ricker[names(df)] ## limit to the six chosen sites

#stan_model_output_Gompertz <- readRDS("./rds files/stan_6riv_output_Gompertz.rds")

## Extract and summarize parameters
par_Ricker <- lapply(stan_model_output_Ricker, function(x) rstan::extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p")))
#par_Gompertz <- lapply(stan_model_output_Gompertz, function(x) rstan::extract(x, c("beta_0","beta_1","s","c","B","P","pred_GPP","sig_p")))


##################################
## Persistence plots
##################################

## Mean parameter
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
#meanpar_G <- lapply(par_Gompertz, function(x) par_mean(x))

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
    temp <- exp(-exp(x$s*(pq[i]-x$c)))
    p_median[i] <- median(temp)
    p_up[i] <- quantile(temp, probs = 0.975)
    p_down[i] <- quantile(temp, probs = 0.025)
  }
  df <- as.data.frame(as.matrix(cbind(pq,p_median, p_up, p_down)))
  df$site_name <- x$site_name
  
  return(df)
}

P_dat_R <- ldply(lapply(P_R, function(z) plotting_P_dat(z)), data.frame); P_dat_R$PM <- "Ricker"
#P_dat_G <- ldply(lapply(P_G, function(z) plotting_P_dat(z)), data.frame); P_dat_G$PM <- "Gompertz"

P_df <- P_dat_R

#####################
## Visualize & compare to bankfull discharge
#####################

## Import and merge bankfull discharge with site info
RI_2 <- read.csv("../data/RI_2yr_flood.csv", header=T)
sapply(RI_2, class)
site_info <- merge(site_info, RI_2, by="site_name")

## join by river name
P_df <- left_join(P_df, site_info[,c("site_name","short_name")], by="site_name")
P_df$short_name <- factor(P_df$short_name, levels=c("Silver Creek, UT",
                                              "Medina River, TX",
                                              "Anacostia River, MD",
                                              "West Fork River, WV",
                                              "St. John's River, FL",
                                              "Clackamas River, OR"))


Persistence_plots <- function(site, df, site_info, P_df){

  Q_sub <- df[[site]]
  #Q_sub$pq <- Q_sub$tQ
  Q_sub$p_for_q <- 0.5
  
  ## critical Q based on velocity
  crit_Q <- site_info[which(site_info$site_name == site),]$RI_2yr_Q
  
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
    annotate("text", label=as.character(P$short_name[1]),
             x = 1.2*c, 
             y= 0.9, size=3.75, hjust=0)+
    labs(x="Range of Standardized Discharge",y="Persistence")+
    scale_y_continuous(limits=c(0,1))+
    geom_vline(xintercept = crit_Q, size=1, linetype="dotted", color="grey25")+
    geom_vline(xintercept = c, size=1, linetype="dashed")
  
  
  Persist_plot2 <- ggExtra::ggMarginal(Persist_plot, data=Q_sub, type="histogram",
                                       size=4, x = Q, margins = "x", color="black",
                                       fill="deepskyblue4", xparams = list(alpha=0.8))
  
  return(Persist_plot2)
  
}

Anacostia <- Persistence_plots("nwis_01649500", df, site_info, P_df) #, 0.05*max(df$nwis_01649500$Q),0.1)
St_Johns <- Persistence_plots("nwis_02234000", df, site_info, P_df) #,0.05*max(df$nwis_02234000$Q),0.1)
West_Fork <- Persistence_plots("nwis_03058000", df, site_info, P_df) #,0.05*max(df$nwis_03058000$Q),0.1)
Medina <- Persistence_plots("nwis_08180700", df, site_info, P_df) #,0.05*max(df$nwis_08180700$Q),0.1)
Silver <- Persistence_plots("nwis_10129900", df, site_info, P_df) #,0.05*max(df$nwis_10129900$Q),0.1)
Clackamas <- Persistence_plots("nwis_14211010", df, site_info, P_df) #,0.05*max(df$nwis_14211010$Q),0.1)

## order based on river order
grid.arrange(
  arrangeGrob(grobs=list(Silver, Medina, Anacostia, West_Fork, St_Johns, Clackamas),
              ncol = 2,
              bottom=textGrob("Discharge (cms)", gp=gpar(fontsize=16)), 
              left=textGrob("Biomass Persistence", gp=gpar(fontsize=16), rot=90)),
  widths=c(9,1)
)


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
  z$short_name <- factor(z$short_name, levels=c("Silver Creek, UT",
                                                "Medina River, TX",
                                                "Anacostia River, MD",
                                                "West Fork River, WV",
                                                "St. John's River, FL",
                                                "Clackamas River, OR"))
  ggplot(z, aes(s, c))+geom_point()+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))+
    facet_wrap(~short_name, ncol = 2)+
    labs(title = t)
  
}

sc_plotting(sc_R, "Productivity Model: Ricker")
#sc_plotting(sc_G, "PM4: Gompertz")








