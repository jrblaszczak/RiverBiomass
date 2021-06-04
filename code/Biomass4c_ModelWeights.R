## Out-of-sample predictions - model weights

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse"), require, character.only=T)

## Source data
source("DataSource_6rivers_oos_StreamLight.R")
df <- dat_oos

################################
## Model output plot function
################################
GPP_oos_preds <- function(preds, df, mean_mod, se_mod){

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

STS_oos <- GPP_oos_preds("./rds files/Sim_6riv_AR_oos.rds", df, "mean_p_GPP_STS", "se_p_GPP_STS")
LB_oos <- GPP_oos_preds("./rds files/Sim_6riv_Ricker_oos.rds", df, "mean_p_GPP_LB", "se_p_GPP_LB")
#Gomp_oos <- GPP_oos_preds("./rds files/Sim_9riv_Gompertz_oos.rds", df, "mean_p_GPP_Gomp", "se_p_GPP_Gomp")

## Combine
oos_preds <- merge(STS_oos, LB_oos[,c("site_name","Date","mean_p_GPP_LB", "se_p_GPP_LB")])
#oos_preds <- merge(oos_preds, Gomp_oos[,c("site_name","Date","mean_p_GPP_Gomp", "se_p_GPP_Gomp")],
#                   by=c("site_name","Date"))

## Split by short_name
oos_pl <- split(oos_preds, oos_preds$short_name)

#######################
## Calc weights
########################

visualize_support <- function(short.name){
  
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
  
  ## visualize
  theme_set(theme_bw())
  
  vis.plot <- plot_grid(
    ggplot(support, aes(Date, exp(GPP)))+geom_point(color="grey75")+
      geom_line(aes(Date, exp(mean_p_GPP_STS)),color="purple",size=0.9)+
      geom_ribbon(aes(Date, ymin=exp(mean_p_GPP_STS) - exp(se_p_GPP_STS),
                      ymax=exp(mean_p_GPP_STS) + exp(se_p_GPP_STS)), fill="purple",alpha=0.3)+
      geom_line(aes(Date, exp(mean_p_GPP_LB)),color="chartreuse4",size=0.9)+
      geom_ribbon(aes(Date, ymin=exp(mean_p_GPP_LB) - exp(se_p_GPP_LB),
                      ymax=exp(mean_p_GPP_LB) + exp(se_p_GPP_LB)), fill="chartreuse4",alpha=0.3)+
      labs(y="GPP",title = short.name),
    
    ggplot(support, aes(Date, STS_support))+geom_line(color="purple",size=0.9)+
      geom_line(aes(Date, LB_support),color="chartreuse4",size=0.9),
    ncol=1, align = "hv")
  
  return(vis.plot)

}

## save images
setwd("../figures/Model support")

for(i in 1:length(site_info$short_name)){
  
  ggsave(filename = paste(site_info$short_name[i],"_model_support.jpeg",sep=""),
         plot = visualize_support(site_info$short_name[i]),device = "jpeg")
  

}

