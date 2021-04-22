## Out-of-sample predictions - model weights

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse"), require, character.only=T)

## Source data
source("DataSource_9rivers_oos.R")
df <- dat_oos
# Subset source data
df <- df[c("nwis_01649500","nwis_02234000","nwis_03058000",
           "nwis_08180700","nwis_10129900","nwis_14211010")]

## Change river names to short names
site_info[,c("site_name","long_name","NHD_STREAMORDE")]
site_info <- site_info[which(site_info$site_name %in% names(df)),]

################################
## Model output plot function
################################
GPP_oos_preds <- function(preds, df, mean_mod, se_mod){

  simmat_list <- readRDS(preds)
  simmat_list <- simmat_list[names(df)]
  
  # For every day extract median and CI
  mean_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) mean(x))), data.frame)
  se_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) sd(x)/sqrt(7500))), data.frame)
  
  ## Plot simulated GPP
  dat <- ldply(df, data.frame)
  df_sim <- as.data.frame(cbind(dat$site_name, as.character(dat$date), 
                                dat$GPP, dat$GPP_sd/sqrt(2000), 
                                mean_simmat$X..i.., se_simmat$X..i..))
  colnames(df_sim) <- c("site_name","Date","GPP","GPP_se", mean_mod, se_mod)
  df_sim$Date <- as.POSIXct(as.character(df_sim$Date), format="%Y-%m-%d")
  df_sim[,3:6] <- apply(df_sim[,3:6],2,function(x) as.numeric(as.character(x)))
  df_sim[,3:6] <- log(df_sim[,3:6])
  
  ## Arrange rivers by river order
  df_sim <- left_join(df_sim, site_info[,c("site_name","short_name")])
  df_sim$short_name <- factor(df_sim$short_name, levels=c("Silver Creek, UT",
                                                          "Medina River, TX",
                                                          "Anacostia River, MD",
                                                          "West Fork River, WV",
                                                          "St. John's River, FL",
                                                          "Clackamas River, OR"))
  
  return(df_sim)
  
}

STS_oos <- GPP_oos_preds("./rds files/Sim_9riv_AR_oos.rds", df, "mean_p_GPP_STS", "se_p_GPP_STS")
LB_oos <- GPP_oos_preds("./rds files/Sim_9riv_Ricker_oos.rds", df, "mean_p_GPP_LB", "se_p_GPP_LB")
Gomp_oos <- GPP_oos_preds("./rds files/Sim_9riv_Gompertz_oos.rds", df, "mean_p_GPP_Gomp", "se_p_GPP_Gomp")

## Combine
oos_preds <- merge(STS_oos, LB_oos[,c("site_name","Date","mean_p_GPP_LB", "se_p_GPP_LB")])
oos_preds <- merge(oos_preds, Gomp_oos[,c("site_name","Date","mean_p_GPP_Gomp", "se_p_GPP_Gomp")],
                   by=c("site_name","Date"))

## Split by short_name
oos_pl <- split(oos_preds, oos_preds$short_name)

test <- oos_pl$`Anacostia River, MD`
test$p_obs_STS <- dnorm(test$GPP, test$mean_p_GPP_STS, sqrt(test$GPP_se^2 + test$se_p_GPP_STS^2), log = FALSE)
test$p_obs_LB <- dnorm(test$GPP, test$mean_p_GPP_LB, sqrt(test$GPP_se^2 + test$se_p_GPP_LB^2), log = FALSE)


ggplot(test, aes(Date, p_obs_STS))+geom_line()+
  geom_line(aes(Date, p_obs_LB), color="blue")


#######################
## Calc weights
########################
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


plot(weight_STS)
plot(weight_LB)
