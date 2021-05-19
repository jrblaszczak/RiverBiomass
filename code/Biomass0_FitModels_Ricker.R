## Fitting models to data
## JR Blaszczak

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","Metrics","MCMCglmm","tictoc"), require, character.only=T)

## Source data
source("DataSource_6rivers.R")
df <- df[c("nwis_07191222","nwis_01608500")]

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates
#########################################

## PM 1 - Standard time series
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x,chains=3,iter=5000, control=list(max_treedepth=12))) #5000 seconds

## PM 2 - Latent Biomass (Ricker)
init_Ricker <- function(...) {
  list(c = 0.5, s = 200)
}

PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12)))

## PM 3 - Latent Biomass (Ricker) - light modification
PM_outputlist_Ricker_Lmod <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr_Lmod.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12)))

## PM 3 - Latent Biomass (Ricker) - GPP modification
PM_outputlist_Ricker_Gmod <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_obserr_Gmod.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12)))


###################################
## Visualize differences among models
##################################
# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p, sig_o

##########################
## Model 1 Output - AR
#########################
names(df); names(PM_outputlist_AR)
AR_list <- Map(c, PM_outputlist_AR, df)


AR_sim_fxn <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars1 <- extract(output, c("phi","alpha","beta","sig_p","sig_o"))
  simmat1<-matrix(NA,length(df$GPP),length(unlist(pars1$phi)))
  rmsemat1<-matrix(NA,length(df$GPP),1)
  
  # Simulate
  for (i in 1:length(pars1$phi)){
    simmat1[,i]<-PM_AR(pars1$phi[i],pars1$alpha[i],pars1$beta[i],pars1$sig_p[i],pars1$sig_o[i],df)
    rmsemat1[i]<-sqrt(sum((simmat1[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat1, rmsemat1)
  return(l)
  
}

AR_sim <- lapply(AR_list, function(x) AR_sim_fxn(x))


###############################
## Model 2 Output - Ricker
###############################
names(df); names(PM_outputlist_Ricker)
Ricker_list <- Map(c, PM_outputlist_Ricker, df)

Ricker_sim_fxn <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars3<-extract(output, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o"))
  simmat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  biomat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  rmsemat3<-matrix(NA,length(df$GPP),1)
  #Simulated
  for (i in 1:length(pars3$r)){
    simmat3[,i]<-PM_Ricker(pars3$r[i],pars3$lambda[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    biomat3[,i]<-PM_Ricker_B(pars3$r[i],pars3$lambda[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat3, rmsemat3, biomat3)
  return(l)
  
}

Ricker_sim <- lapply(Ricker_list, function(x) Ricker_sim_fxn(x))



###############################
## Model 3 Output - Ricker, Lmod
###############################
names(df); names(PM_outputlist_Ricker_Lmod)
Ricker_list_Lmod <- Map(c, PM_outputlist_Ricker_Lmod, df)

Ricker_sim_fxn_Lmod <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars3<-extract(output, c("r","lambda","alpha","s","c","B","P","pred_GPP","sig_p","sig_o"))
  simmat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  biomat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  rmsemat3<-matrix(NA,length(df$GPP),1)
  #Simulated
  for (i in 1:length(pars3$r)){
    simmat3[,i]<-PM_Ricker_Lmod(pars3$r[i],pars3$lambda[i],pars3$alpha[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    biomat3[,i]<-PM_Ricker_Lmod_B(pars3$r[i],pars3$lambda[i],pars3$alpha[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat3, rmsemat3, biomat3)
  return(l)
  
}

Ricker_sim_Lmod <- lapply(Ricker_list_Lmod, function(x) Ricker_sim_fxn_Lmod(x))


###############################
## Model 4 Output - Ricker, Gmod
###############################
names(df); names(PM_outputlist_Ricker_Gmod)
Ricker_list_Gmod <- Map(c, PM_outputlist_Ricker_Gmod, df)

Ricker_sim_fxn_Gmod <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars3<-extract(output, c("r","lambda","alpha","s","c","B","P","pred_GPP","sig_p","sig_o"))
  simmat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  biomat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  rmsemat3<-matrix(NA,length(df$GPP),1)
  #Simulated
  for (i in 1:length(pars3$r)){
    simmat3[,i]<-PM_Ricker_Gmod(pars3$r[i],pars3$lambda[i],pars3$alpha[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    biomat3[,i]<-PM_Ricker_Gmod_B(pars3$r[i],pars3$lambda[i],pars3$alpha[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat3, rmsemat3, biomat3)
  return(l)
  
}

Ricker_sim_Gmod <- lapply(Ricker_list_Gmod, function(x) Ricker_sim_fxn_Gmod(x))


## Save prediction
saveRDS(AR_sim, "./rds files/AR_ws_comparison.rds")
saveRDS(Ricker_sim, "./rds files/Ricker_ws_comparison.rds")
saveRDS(Ricker_sim_Lmod, "./rds files/Ricker_Lmod_ws_comparison.rds")
saveRDS(Ricker_sim_Gmod, "./rds files/Ricker_Gmod_ws_comparison.rds")

#################################
## Plot within sample pred comparison
#################################

vis_preds <- function(simmat1_list, dataTS, mod_title){
  
  # For every day extract median and CI
  median_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
  lower_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
  upper_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)
  
  ## Plot simulated GPP
  dat1 <- ldply(dataTS, data.frame)
  df_sim1 <- as.data.frame(cbind(dat1$site_name, as.character(dat1$date), dat1$GPP, median_simmat1$X..i.., lower_simmat1$X..i.., upper_simmat1$X..i..))
  colnames(df_sim1) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
  df_sim1$Date <- as.POSIXct(as.character(df_sim1$Date), format="%Y-%m-%d")
  df_sim1[,3:6] <- apply(df_sim1[,3:6],2,function(x) as.numeric(as.character(x)))
  
  ## Arrange rivers by river order
  df_sim1 <- left_join(df_sim1, site_info[,c("site_name","short_name")])
  df_sim1$short_name <- factor(df_sim1$short_name, levels=site_order_list)
  
  ## Plot
  df_sim1_plot <- ggplot(df_sim1, aes(Date, GPP))+
    geom_point(size=2, color="black")+
    geom_line(aes(Date, sim_GPP), color="#1C474D", size=1.2)+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=mod_title)+
    geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                fill="#1C474D", alpha=0.3, show.legend = FALSE)+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))+
    coord_cartesian(ylim = c(0,17))+
    facet_wrap(~short_name, scales = "free", ncol = 2)
  return(df_sim1_plot)
  
}

plot_grid(
  vis_preds(Ricker_sim, df, "LB - Ricker"),
  vis_preds(Ricker_sim_Lmod, df, "LB - Ricker Light"),
  vis_preds(Ricker_sim_Gmod, df, "LB - Ricker GPP"),
  ncol = 1)



## Latent biomass

vis_LB_preds <- function(simmat3_list, dataTS, mod_title){
  
  ## Plot latent biomass
  # For every day extract median and CI
  median_biomat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[3]], 1, function(x) median(x))), data.frame)
  lower_biomat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[3]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
  upper_biomat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[3]], 1, function(x) quantile(x, probs = 0.975))), data.frame)
  
  ## Plot simulated GPP
  dat3 <- ldply(dataTS, data.frame)
  df_bio3 <- as.data.frame(cbind(dat3$site_name, as.character(dat3$date),
                                 median_biomat3[,2], lower_biomat3[,2], upper_biomat3[,2]))
  colnames(df_bio3) <- c("site_name","Date","B","B_lower","B_upper")
  df_bio3$Date <- as.POSIXct(as.character(df_bio3$Date), format="%Y-%m-%d")
  df_bio3[,3:5] <- apply(df_bio3[,3:5],2,function(x) as.numeric(as.character(x)))
  
  ## Arrange rivers by river order
  df_bio3 <- left_join(df_bio3, site_info[,c("site_name","short_name")])
  df_bio3$short_name <- factor(df_bio3$short_name, levels=site_order_list)
  
  df_modB3 <- df_bio3
  
  ## plot
  df_modB3_plot <- ggplot(df_modB3, aes(Date, exp(B)))+
    geom_line(size=1.2, color="chartreuse4")+
    labs(y="Latent Biomass",title=mod_title)+
    geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
                fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
    theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))+
    coord_cartesian(ylim=c(0,25))+
    facet_wrap(~short_name, scales = "free_x", ncol = 2)
  return(df_modB3_plot)
  
}

plot_grid(
  vis_LB_preds(Ricker_sim, df, "LB - Ricker"),
  vis_LB_preds(Ricker_sim_Lmod, df, "LB - Ricker Light"),
  vis_LB_preds(Ricker_sim_Gmod, df, "LB - Ricker GPP"),
  ncol = 1)

##############################
## Out of sample predictions
##############################
## Source data
source("DataSource_6rivers_oos.R")

dat_oos <- dat_oos[names(df)]

## predict the next year
# STS
AR_list_oos <- Map(c, PM_outputlist_AR, dat_oos)
AR_sim_oos <- lapply(AR_list_oos, function(x) AR_sim_fxn(x))
# Ricker
Ricker_list_oos <- Map(c, PM_outputlist_Ricker, dat_oos)
Ricker_sim_oos <- lapply(Ricker_list_oos, function(x) Ricker_sim_fxn(x))
# Ricker Lmod
Ricker_list_Lmod_oos <- Map(c, PM_outputlist_Ricker_Lmod, dat_oos)
Ricker_sim_Lmod_oos <- lapply(Ricker_list_Lmod_oos, function(x) Ricker_sim_fxn_Lmod(x))
# Ricker Gmod
Ricker_list_Gmod_oos <- Map(c, PM_outputlist_Ricker_Gmod, dat_oos)
Ricker_sim_Gmod_oos <- lapply(Ricker_list_Gmod_oos, function(x) Ricker_sim_fxn_Gmod(x))

## Save prediction
saveRDS(AR_sim_oos, "./rds files/AR_oos_comparison.rds")
saveRDS(Ricker_sim_oos, "./rds files/Ricker_oos_comparison.rds")
saveRDS(Ricker_sim_Lmod_oos, "./rds files/Ricker_Lmod_oos_comparison.rds")
saveRDS(Ricker_sim_Gmod_oos, "./rds files/Ricker_Gmod_oos_comparison.rds")


## Plot OOS predictions
plot_grid(
  vis_preds(Ricker_sim_oos, dat_oos, "LB - Ricker - OOS Prediction"),
  vis_preds(Ricker_sim_Lmod_oos, dat_oos, "LB - Ricker Light - OOS Prediction"),
  vis_preds(Ricker_sim_Gmod_oos, dat_oos, "LB - Ricker GPP - OOS Prediction"),
  ncol = 1)

plot_grid(
  vis_LB_preds(Ricker_sim_oos, dat_oos, "LB - Ricker - OOS Prediction"),
  vis_LB_preds(Ricker_sim_Lmod_oos, dat_oos, "LB - Ricker Light - OOS Prediction"),
  vis_LB_preds(Ricker_sim_Gmod_oos, dat_oos, "LB - Ricker GPP - OOS Prediction"),
  ncol = 1)











