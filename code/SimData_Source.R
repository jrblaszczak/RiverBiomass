## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- read.csv("NWIS_sub_SiteYears.csv", header=T)
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

## subset to only 2010
data <- data[which(data$year == "2010"),]

## Change GPP sd column name
colnames(data)[which(colnames(data) == "GPP_daily_sd_name")] <- "GPP_sd"

## Assign IDs
data$ID <- NA
data[which(data$site_name == "nwis_06711565"),]$ID <- "SouthPlatte_CO"
data[which(data$site_name == "nwis_08181500"),]$ID <- "Medina_TX"
data[which(data$site_name == "nwis_14211010"),]$ID <- "Clackamas_OR"

## split list by ID
l <- split(data, data$ID)

rel_LQT <- function(x){
  x$light_rel <- x$light/max(x$light)
  x$daylength_rel <- x$daylength_hrs/max(x$daylength_hrs)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- (x$Q-mean(x$Q))/sd(x$Q)
  x<-x[order(x$date),]
  return(x)
}

dat <- lapply(l, function(x) rel_LQT(x))



