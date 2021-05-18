## Figure - Density-dependence

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork"), require, character.only=T)

## Source data
source("DataSource_6rivers.R")

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02" # AR
PM_Ricker.col <- "#1C474D" # Ricker
PM_Gompertz.col <- "#743731" # Gompertz

## Import stan fits - simulate one at a time
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker_2021_05_16.rds")
#stan_model_output_Gompertz <- readRDS("./rds files/stan_6riv_output_Gompertz.rds")

## Extract parameters
par_Ricker <- lapply(stan_model_output_Ricker, function(x) rstan::extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o")))
#par_Gompertz <- lapply(stan_model_output_Gompertz, function(x) rstan::extract(x, c("beta_0","beta_1","s","c","B","P","pred_GPP","sig_p")))

#################################################
## r and K distributions
################################################

r_K_func <- function(x) {
  rx <- x$r
  rm <- as.data.frame(cbind(rx[1:2500],rx[2501:5000],rx[5001:7500]))
  rm <- rm %>%
    rowwise() %>% mutate(Avg=mean(c(V1, V2, V3))) 
  
  lx <- x$lambda
  lambda <- as.data.frame(cbind(lx[1:2500],lx[2501:5000],lx[5001:7500]))
  lambda <- lambda %>%
    rowwise() %>% mutate(Avg=mean(c(V1, V2, V3))) 
  
  new <- as.data.frame(cbind(rm$Avg, lambda$Avg))
  colnames(new) <- c("r","lambda")
  new$K <- (-1*new$r)/new$lambda
  
  long <- gather(new)
  
  return(long)
}


rK <- ldply(lapply(par_Ricker, function(x) r_K_func(x)), data.frame)
rK$short_name <- revalue(as.character(rK$.id), replace = c("nwis_05406457"="Black Earth Creek, WI",
                                                           "nwis_01656903"="Fatlick Branch, VA",
                                                           "nwis_07191222"="Beaty Creek, OK",
                                                           "nwis_14206950"="Fanno Creek, OR",
                                                           "nwis_01608500"="South Branch Potomac River, WV",
                                                           "nwis_11273400"="San Joaquin River, CA"))
rK$short_name <- factor(rK$short_name, levels= site_order_list)
rK$key <- factor(rK$key, levels = c("r","lambda","K"))

## Visualize
rK <- rK[-which(rK$short_name == "San Joaquin River, CA"),]


ggplot(rK, aes(value, fill = short_name))+
  geom_histogram()+
  facet_grid(~key, scales = "free_x")
 



##################################################
## Summarize parameters
###################################################

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






