## Persistence Curves

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2"), require, character.only=T)

## Source data
source("DataSource_6rivers.R")

# source simulation models
source("Simulated_ProductivityModel1_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Simulated_ProductivityModel2_Logistic.R") # parameters: r, K, s, c, sig_p
source("Simulated_ProductivityModel3_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Simulated_ProductivityModel5_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM1.col <- "#d95f02"
PM2.col <- "#7570b3"
PM3.col <- "#1C474D"
PM4.col <- "#743731"

## Change river names to short names
site_info[,c("site_name","long_name","NHD_STREAMORDE")]
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_01649500"="Anacostia River, MD",
                                                                               "nwis_02234000"="St. John's River, FL",
                                                                               "nwis_03058000"="West Fork River, WV",
                                                                               "nwis_08180700"="Medina River, TX",
                                                                               "nwis_10129900"="Silver Creek, UT",
                                                                               "nwis_14211010"="Clackamas River, OR"))


##################################################
## Extract 
###################################################

## Import stan fits - simulate one at a time
stan_model_output_Logistic <- readRDS("./rds files/stan_6riv_output_Logistic.rds")
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker.rds")
stan_model_output_Gompertz <- readRDS("./rds files/stan_6riv_output_Gompertz.rds")

## Extract and summarize parameters
par_Logistic <- lapply(stan_model_output_Logistic, function(x) rstan::extract(x, c("r","K","s","c","B","P","pred_GPP","sig_p")))
par_Ricker <- lapply(stan_model_output_Ricker, function(x) rstan::extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p")))
par_Gompertz <- lapply(stan_model_output_Gompertz, function(x) rstan::extract(x, c("beta_0","beta_1","s","c","B","P","pred_GPP","sig_p")))

## Extract s versus c
sc_L <- ldply(lapply(par_Logistic, function(x) cbind(x$s, x$c)), data.frame); colnames(sc_L) <- c("site_name","s","c")
sc_R <- ldply(lapply(par_Ricker, function(x) cbind(x$s, x$c)), data.frame); colnames(sc_R) <- c("site_name","s","c")
sc_G <- ldply(lapply(par_Gompertz, function(x) cbind(x$s, x$c)), data.frame); colnames(sc_G) <- c("site_name","s","c")

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

sc_plotting(sc_L, "PM2: Logistic")
sc_plotting(sc_R, "PM3: Ricker")
sc_plotting(sc_G, "PM4: Gompertz")


##################################
## Persistence plots
##################################
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
meanpar_L <- lapply(par_Logistic, function(x) par_mean(x))
meanpar_R <- lapply(par_Ricker, function(x) par_mean(x))
meanpar_G <- lapply(par_Gompertz, function(x) par_mean(x))

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

P_L <- persistence_list(par_Logistic, dat)
P_R <- persistence_list(par_Ricker, dat)
P_G <- persistence_list(par_Gompertz, dat)


## plot
plotting_P_dat <- function(x){
  pq <- seq(x$range[1],x$range[2], length=length(x$s))
  p_mean <- numeric()
  p_up <- numeric()
  p_down <- numeric()
  name <- substring(deparse(substitute(x)),7)
  for(i in 1:length(pq)){
    temp <- exp(-exp(x$s*(pq[i]-x$c)))
    p_mean[i] <- mean(temp)
    p_up[i] <- quantile(temp, probs = 0.975)
    p_down[i] <- quantile(temp, probs = 0.025)
  }
  df <- as.data.frame(as.matrix(cbind(pq,p_mean, p_up, p_down)))
  df$site_name <- x$site_name
  
  return(df)
}

P_dat_L <- ldply(lapply(P_L, function(z) plotting_P_dat(z)), data.frame); P_dat_L$PM <- "Logistic"
P_dat_R <- ldply(lapply(P_R, function(z) plotting_P_dat(z)), data.frame); P_dat_R$PM <- "Ricker"
P_dat_G <- ldply(lapply(P_G, function(z) plotting_P_dat(z)), data.frame); P_dat_G$PM <- "Gompertz"

P_df <- rbind(P_dat_L, P_dat_R)

#####################
## Visualize
#####################

## join by river name
P_df <- left_join(P_df, site_info[,c("site_name","short_name")])
P_df$short_name <- factor(P_df$short_name, levels=c("Silver Creek, UT",
                                              "Medina River, TX",
                                              "Anacostia River, MD",
                                              "West Fork River, WV",
                                              "St. John's River, FL",
                                              "Clackamas River, OR"))

ggplot(P_df, aes(pq, p_mean, color=PM))+
  geom_line(size=1.5, alpha=0.9)+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15), 
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  facet_wrap(~short_name, ncol = 2)+
  scale_color_manual(values = c("Logistic" = PM2.col, "Ricker"= PM3.col))+
  labs(x="Range of Standardized Discharge",y="Persistence")

















