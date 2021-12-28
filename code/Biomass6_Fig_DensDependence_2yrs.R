## Figure - Density-dependence

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork","viridis"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p
#source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02" # AR
PM_Ricker.col <- "#1C474D" # Ricker
PM_Gompertz.col <- "#743731" # Gompertz

## Import stan fits - simulate one at a time
stan_model_output_Ricker_yr1 <- readRDS("./rds files/stan_6riv_output_Ricker_2021_06_01.rds")
stan_model_output_Ricker_yr2 <- readRDS("./rds files/stan_6riv_2ndYr_output_Ricker_2021_06_15.rds")

## Extract parameters
par_Ricker_yr1 <- lapply(stan_model_output_Ricker_yr1, function(x) rstan::extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o")))
par_Ricker_yr2 <- lapply(stan_model_output_Ricker_yr2, function(x) rstan::extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o")))


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


rK_yr1 <- ldply(lapply(par_Ricker_yr1, function(x) r_K_func(x)), data.frame)
rK_yr2 <- ldply(lapply(par_Ricker_yr2, function(x) r_K_func(x)), data.frame)


rk_formatting <- function(x){
  x$short_name <- revalue(as.character(x$.id), replace = c("nwis_02336526"="Proctor Creek, GA",
                                                             "nwis_01649190"="Paint Branch, MD",
                                                             "nwis_07191222"="Beaty Creek, OK",
                                                             "nwis_01608500"="S. Br. Potomac River, WV",
                                                             "nwis_11044000"="Santa Margarita River, CA",
                                                             "nwis_08447300"="Pecos River, TX"))
  x$short_name_SO <- revalue(as.character(x$.id), replace = c("nwis_02336526"="Proctor Creek, GA (2)",
                                                                "nwis_01649190"="Paint Branch, MD (2)",
                                                                "nwis_07191222"="Beaty Creek, OK (3)",
                                                                "nwis_01608500"="S. Br. Potomac River, WV (5)",
                                                                "nwis_11044000"="Santa Margarita River, CA (6)",
                                                                "nwis_08447300"="Pecos River, TX (7)"))
  x$short_name_SO <- factor(x$short_name_SO, levels= c("Proctor Creek, GA (2)",
                                                         "Paint Branch, MD (2)",
                                                         "Beaty Creek, OK (3)",
                                                         "S. Br. Potomac River, WV (5)", 
                                                         "Santa Margarita River, CA (6)",
                                                         "Pecos River, TX (7)"))
  #x$key <- factor(x$key, levels = c("r","lambda","K"))
  return(x)
}

rK1 <- rk_formatting(rK_yr1)  
rK2 <- rk_formatting(rK_yr2)
rK1$year <- "Year 1"; rK2$year <- "Year 2"

#Combine
rK_vals <- rbind(rK1, rK2)


###################
## Visualize
####################

## Distributions
vis.rK <- function(par){
  p <- ggplot(rK_vals[which(rK_vals$key == par),], aes(value, fill = year))+
    geom_density(color="black")+
    scale_fill_viridis("Site (Stream Order)", discrete = TRUE, alpha=0.6, option="A") +
    theme(axis.text.x = element_text(angle=45, hjust=1))+
    facet_wrap(~short_name_SO, ncol=1)+
    theme(panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.title.x = element_blank(), axis.text = element_text(size=12),
          axis.title.y = element_blank(), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))
  
  return(p)
  
}

r_plot <- vis.rK("r")
K_plot <- vis.rK("K")


## Biplot


rK_sumY1 <- rK_yr1[which(rK_yr1$key %in% c("r","K")),] %>%
  group_by(.id, key) %>%
  summarise(quant_med_yr1 = quantile(value, 0.5),
            quant_lower_yr1 = quantile(value, 0.025),
            quant_upper_yr1 = quantile(value, 0.975))

rK_sumY2 <- rK_yr2[which(rK_yr2$key %in% c("r","K")),] %>%
  group_by(.id, key) %>%
  summarise(quant_med_yr2 = quantile(value, 0.5),
            quant_lower_yr2 = quantile(value, 0.025),
            quant_upper_yr2 = quantile(value, 0.975))
rk_yrs <- merge(rK_sumY1, rK_sumY2, by=c(".id","key")) 

# formatting
rk_yrs$short_name_SO <- revalue(as.character(rk_yrs$.id), replace = c("nwis_02336526"="Proctor Creek, GA (2)",
                                                            "nwis_01649190"="Paint Branch, MD (2)",
                                                            "nwis_07191222"="Beaty Creek, OK (3)",
                                                            "nwis_01608500"="S. Br. Potomac River, WV (5)",
                                                            "nwis_11044000"="Santa Margarita River, CA (6)",
                                                            "nwis_08447300"="Pecos River, TX (7)"))
rk_yrs$short_name_SO <- factor(rk_yrs$short_name_SO, levels= c("Proctor Creek, GA (2)",
                                                     "Paint Branch, MD (2)",
                                                     "Beaty Creek, OK (3)",
                                                     "S. Br. Potomac River, WV (5)", 
                                                     "Santa Margarita River, CA (6)",
                                                     "Pecos River, TX (7)"))
r_yrs <- rk_yrs[which(rk_yrs$key == "r"),];K_yrs <- rk_yrs[which(rk_yrs$key == "K"),]

r_biplot <- ggplot(r_yrs, aes(quant_med_yr1, quant_med_yr2, fill=short_name_SO))+
  scale_x_continuous(limits=c(0,0.8))+scale_y_continuous(limits = c(0,0.8))+
  scale_fill_viridis("Site (Stream Order)", discrete = TRUE, alpha=0.6, option="A")+
  labs(x="Year 1",y="Year 2",title=expression("r"))+
  geom_point(shape=21, size=4)+ geom_abline(slope=1, intercept = 0)+
  geom_errorbar(aes(ymin=quant_lower_yr2, ymax=quant_upper_yr2))+
  geom_errorbarh(aes(xmin=quant_lower_yr1, xmax=quant_upper_yr1))
  

K_biplot <- ggplot(K_yrs, aes(quant_med_yr1, quant_med_yr2, fill=short_name_SO))+
  scale_x_continuous(limits=c(0,22))+scale_y_continuous(limits = c(0,22))+
  scale_fill_viridis("Site (Stream Order)", discrete = TRUE, alpha=0.6, option="A")+
  labs(x="Year 1",y="Year 2")+
  geom_point(shape=21, size=4)+ geom_abline(slope=1, intercept = 0)+
  geom_errorbar(aes(ymin=quant_lower_yr2, ymax=quant_upper_yr2))+
  geom_errorbarh(aes(xmin=quant_lower_yr1, xmax=quant_upper_yr1))


plot_grid( plot_grid(r_biplot+theme(legend.position = "none"),
          K_biplot+theme(legend.position = "none"),
          ncol=2, align="hv", labels=c('A','B')),
          get_legend(r_biplot),
          ncol=2, rel_widths = c(0.7,0.25))














## Median parameter
par_summary <- function(par) {

  ## Find the median
  median_par <- ldply(lapply(par, function(x) median(x)), data.frame)
  lower_par <- ldply(lapply(par, function(x) quantile(x, probs = 0.025)), data.frame)
  upper_par <- ldply(lapply(par, function(x) quantile(x, probs = 0.975)), data.frame)
  
  ## Compile in list and return
  par_vals <- cbind(median_par, lower_par$X..i.., upper_par$X..i..)
  names(par_vals) <- c("parameter","median","lowerCI","upperCI")
  return(par_vals)
}

medianpar_yr1 <- ldply(lapply(par_Ricker_yr1, function(x) par_summary(x)), data.frame)
medianpar_yr2 <- ldply(lapply(par_Ricker_yr2, function(x) par_summary(x)), data.frame)
medianpar_yr1$year <- "Year 1"; medianpar_yr2$year <- "Year 2"

#Combine
rK_summary <- rbind(medianpar_yr1, medianpar_yr2)
rsum <- rK_summary[which(rK_summary$parameter == "r"),]
Ksum <- rK_summary[which(rK_summary$parameter == "K"),]

rK_wide <- spread(rK_summary, key = )






















