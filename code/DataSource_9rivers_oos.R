## 9 rivers seccond year data source
## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot",
         "lubridate","tidyverse", "reshape2"), require, character.only=T)

##############################
## Data Import & Processing ##
##############################
data <- readRDS("./rds files/NWIS_9site_subset.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("./rds files/NWIS_9siteinfo_subset.rds")

## Change river names to short names
site_info$short_name <- revalue(as.character(site_info$site_name),
                                replace = c("nwis_01649500"="Anacostia River, MD",
                                            "nwis_02234000"="St. John's River, FL",
                                            "nwis_03058000"="West Fork River, WV",
                                            "nwis_08180700"="Medina River, TX",
                                            "nwis_10129900"="Silver Creek, UT",
                                            "nwis_14211010"="Clackamas River, OR",
                                            "nwis_01645762"="Little Difficult Run, VA",
                                            "nwis_01649190"="Paint Branch Crrek, MD",
                                            "nwis_04137500"="Au Sable River, MI"))



## How many days of data per site per year
data$year <- year(data$date)
data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()

## Set any GPP < 0 to a small value close to 0
data[which(data$GPP < 0),]$GPP <- sample(exp(-6):exp(-4), 1)

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
  summarise_at(.vars = c("Q","light"), .funs = max)

oos_relativize <- function(prev_max, post_dat, id){
  
  max.vals <- prev_max[which(prev_max$site_name == id),]
  dat <- post_dat[which(post_dat$site_name == id),]
  
  dat$light_rel <- dat$light/max.vals$light
  dat$tQ <- dat$Q/max.vals$Q
  
  dat <- dat[order(dat$date),]
  return(dat)
  
}


dat_oos <- list()
for(i in 1:nrow(prev_max)){
  dat_oos[[i]] <- oos_relativize(prev_max, post_dat, prev_max[i,]$site_name)
}
names(dat_oos) <- prev_max$site_name


rm(data, data_siteyears, prev_max, prev_dat, post_dat, prev_post, i, oos_relativize)
