##==============================================================================
## Script for compiling and formatting data with stream light for stan
## Second Year - out of sample prediction
## Code author: J.R. Blaszczak
##==============================================================================

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- readRDS("./rds files/NWIS_6site_subset_SL.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("./rds files/NWIS_6siteinfo_subset_SL.rds")

# Read in StreamLight processed data (Savoy)
SL <- readRDS("./rds files/StreamLight_daily_6riv_all.rds")
colnames(SL)[colnames(SL) == "Date"] <- "date"

## Join data and StreamLight
data <- left_join(data, SL, by=c("site_name", "date"))

## Change river names to short names
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_02336526"="Proctor Creek, GA",
                                                                               "nwis_01649190"="Paint Branch, MD",
                                                                               "nwis_07191222"="Beaty Creek, OK",
                                                                               "nwis_01608500"="S. Br. Potomac River, WV",
                                                                               "nwis_11044000"="Santa Margarita River, CA",
                                                                               "nwis_08447300"="Pecos River, TX"))

## Order for figures by stream order
site_order_list <- c("Proctor Creek, GA",
                     "Paint Branch, MD",
                     "Beaty Creek, OK",
                     "S. Br. Potomac River, WV",
                     "Santa Margarita River, CA",
                     "Pecos River, TX")

## How many days of data per site per year
data$year <- year(data$date)
data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()

## Set any GPP < 0 to a small value close to 0
data[which(data$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2


#########################################################################
## Subset data to relativize by the first year max discharge and light
#########################################################################
prev_post <- data %>%
  group_by(site_name) %>%
  summarise_at(.vars = "year", .funs = c(min, max))
colnames(prev_post) <- c("site_name","year1","year2")
prev_post$site_year1 <- paste(prev_post$site_name, prev_post$year1)
prev_post$site_year2 <- paste(prev_post$site_name, prev_post$year2)

data$site_year <- paste(data$site_name, data$year)

## subset time series data
prev_dat <- subset(data, site_year %in% prev_post$site_year1)
post_dat <- subset(data, site_year %in% prev_post$site_year2)

## extract first year max light and discharge by site
prev_max <- prev_dat %>% group_by(site_name) %>%
  summarise_at(.vars = c("Q","PAR_surface"), .funs = max)

oos_relativize <- function(prev_max, post_dat, id){
  
  max.vals <- prev_max[which(prev_max$site_name == id),]
  dat <- post_dat[which(post_dat$site_name == id),]
  
  dat$light_rel_PAR <- dat$PAR_surface/max.vals$PAR_surface
  dat$tQ <- dat$Q/max.vals$Q
  
  dat <- dat[order(dat$date),]
  return(dat)
  
}


dat_oos <- list()
for(i in 1:nrow(prev_max)){
  dat_oos[[i]] <- oos_relativize(prev_max, post_dat, prev_max[i,]$site_name)
}
names(dat_oos) <- prev_max$site_name

rm(data, SL, data_siteyears, prev_max, prev_dat, post_dat, prev_post, i, oos_relativize)

#lapply(dat_oos, function(x) ggplot(x, aes(date, tQ))+geom_line())

