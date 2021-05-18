## 6 rivers data source

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2",
         "PerformanceAnalytics","jpeg","grid",
         "rstan","bayesplot","shinystan"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- readRDS("./rds files/NWIS_6site_subset_v2.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("./rds files/NWIS_6siteinfo_subset_v2.rds")

## Change river names to short names
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_02336526"="Proctor Creek, GA",
                                                                               "nwis_01656903"="Fatlick Branch, VA",
                                                                               "nwis_07191222"="Beaty Creek, OK",
                                                                               "nwis_14206950"="Fanno Creek, OR",
                                                                               "nwis_01608500"="South Br. Potomac River, WV",
                                                                               "nwis_04101500"="St. Joseph River, MI"))

## Order for figures by stream order
site_order_list <- c("Proctor Creek, GA",
                     "Fatlick Branch, VA",
                     "Beaty Creek, OK",
                     "Fanno Creek, OR",
                     "South Br. Potomac River, WV",
                     "St. Joseph River, MI")


## How many days of data per site per year
data$year <- year(data$date)
data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()
## Select the first of two years
data <- rbind(data[which(data$site_name == "nwis_02336526" & data$year %in% c(2015)),],
              data[which(data$site_name == "nwis_01656903" & data$year %in% c(2013)),],
              data[which(data$site_name == "nwis_07191222" & data$year %in% c(2009)),],
              data[which(data$site_name == "nwis_14206950" & data$year %in% c(2015)),],
              data[which(data$site_name == "nwis_01608500" & data$year %in% c(2012)),],
              data[which(data$site_name == "nwis_04101500" & data$year %in% c(2012)),])

## Set any GPP < 0 to a small value close to 0
data[which(data$GPP < 0),]$GPP <- sample(exp(-6):exp(-4), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

## visualize
ggplot(data, aes(date, GPP))+
  geom_point()+geom_line()+
  facet_wrap(~site_name,scales = "free_x")

## split list by ID
l <- split(data, data$site_name)

rel_LQT <- function(x){
  x$light_rel <- x$light/max(x$light)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)

  x<-x[order(x$date),]
  return(x)
}

dat <- lapply(l, function(x) rel_LQT(x))
df <- dat

rm(data,l, data_siteyears, dat)
