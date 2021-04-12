## Figure - Density-dependence

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork"), require, character.only=T)

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
source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02" # AR
PM_Ricker.col <- "#1C474D" # Ricker
PM_Gompertz.col <- "#743731" # Gompertz


## Import stan fits - simulate one at a time
stan_model_output_Ricker <- readRDS("./rds files/stan_9riv_output_Ricker_2021_03_05.rds")
stan_model_output_Ricker <- stan_model_output_Ricker[names(df)]
#stan_model_output_Gompertz <- readRDS("./rds files/stan_6riv_output_Gompertz.rds")












##################################################
## Extract and summarize parameters
###################################################
## Extract and summarize parameters
par_Ricker <- lapply(stan_model_output_Ricker, function(x) rstan::extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p")))
#par_Gompertz <- lapply(stan_model_output_Gompertz, function(x) rstan::extract(x, c("beta_0","beta_1","s","c","B","P","pred_GPP","sig_p")))


## Median parameter
par_median <- function(par) {
  ## Find the median
  median_par <- lapply(par, function(x) median(x))
  
  median_pred_GPP_ts <- apply(par$pred_GPP,2,median)
  median_P_ts <- apply(par$P,2,median)
  median_B_ts <- apply(par$B,2,median)
  
  ## Compile in list and return
  median_par_ts <- list(median_par, median_pred_GPP_ts, median_P_ts, median_B_ts)
  names(median_par_ts) <- c("par","pred_GPP","P","B")
  return(median_par_ts)
}

medianpar_R <- lapply(par_Ricker, function(x) par_median(x))

####################################################################
## Extract latent biomass values below the median Qcrit threshold
####################################################################

DD_plots <- function(medianpar_R, site_num, site_info) {

  dat <- medianpar_R[[site_num]]
  
  short_name <- site_info[which(site_info$site_name == names(medianpar_R)[site_num]),]$short_name
  
  PB <- as.data.frame(cbind(dat$P, dat$B))
  colnames(PB) <- c("P","B")
  PB$B <- ifelse(PB$P < 0.5, yes = NA, no = PB$B)
  PB$B_exp <- exp(PB$B)
  
  PB <- PB %>%
    mutate(B_diff = B_exp - lag(B_exp, n = 1),
           B_lag = lag(B_exp, n=1))

  PB_plot<- ggplot(PB, aes(B_lag, B_diff))+
    geom_point()+
    theme(panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.text = element_text(size=12),
          axis.title = element_blank())+
    annotate("text", label=as.character(short_name), x = 15, y= 7, size=3.5)+
    scale_x_continuous(limits=c(0,20))+scale_y_continuous(limits=c(-10,10))
  
  return(PB_plot)
  
}

anacostia <- DD_plots(medianpar_R, 1, site_info)
st_johns <- DD_plots(medianpar_R, 2, site_info)
west_fork <- DD_plots(medianpar_R, 3, site_info)
medina <- DD_plots(medianpar_R, 4, site_info)
silver <- DD_plots(medianpar_R, 5, site_info)
clackamas <- DD_plots(medianpar_R, 6, site_info)

## order based on river order

(silver | medina)/
  (anacostia | west_fork)/
  (st_johns | clackamas)






