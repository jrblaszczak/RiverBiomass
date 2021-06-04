## Discharge Data Retrieval

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork","dataRetrieval"), require, character.only=T)

# Site numbers
site_numbers <- c("02336526","01649190","07191222",
                  "01608500","11044000","08447300")

# Site information
NWIS_Info <- lapply(site_numbers, function(x) dataRetrieval::readNWISsite(x))

# Mean daily discharge
parameterCd <- "00060"

# Raw daily data
rawDailyData <- lapply(site_numbers, function(y) readNWISdv(y,parameterCd,
                                                            "1970-01-01","2020-12-31"))
names(rawDailyData) <- site_numbers

# Extract relevant information
DailyQ <- lapply(rawDailyData, function(z) return(z[,c("site_no","Date","X_00060_00003","X_00060_00003_cd")]))

# Remove provisional data
DailyQ_clean <- lapply(DailyQ, function(x) return(x[which(x$X_00060_00003_cd %in% c("A","A e")),]))

################################################
## Calculate 2 year flood recurrence interval
##############################################

two_year_flood <- function(data){
  
  #Grouping maximum average daily discharge for each year
  max.by.year<-data %>% group_by(year=floor_date(Date, "year")) %>% summarize(amount=max(X_00060_00003))
  
  #Recording the maximum discharges by year and removing N.A. values
  maximas<-max.by.year$amount
  maximas<-maximas[!is.na(maximas)]
  
  #Sorting maxima by decreasing order
  sorted.maximas<-as.data.frame(sort(maximas,decreasing=T))
  sorted.maximas$rank <- seq(length=nrow(sorted.maximas))
  colnames(sorted.maximas) <- c("Q_max","rank")
  
  #Fit relationship
  sorted.maximas$ln_Q_max <- log(sorted.maximas$Q_max)
  sorted.maximas$exceed_prob <- sorted.maximas$rank/(length(sorted.maximas$rank)+1)
  
  #visualize
  ggplot(sorted.maximas, aes(ln_Q_max, exceed_prob))+
    geom_point()+
    geom_smooth(method = "lm")+
    scale_y_reverse()
  
  #extract coefficients and calc 
  mod <- lm(exceed_prob ~ ln_Q_max, data = sorted.maximas)
  l <- as.data.frame(t(as.matrix(coefficients(mod))))
  colnames(l) <- c("int","slope")
  yr_2 <- exp((0.5 - l$int)/l$slope)
  
  return(yr_2)

}

RI_two <- ldply(lapply(DailyQ_clean, function(y) two_year_flood(y)), data.frame)
colnames(RI_two) <- c(".id","RI_2yr_Q")
RI_two$site_name <- paste("nwis_",RI_two$.id, sep = "")
RI_two <- RI_two[,c("site_name","RI_2yr_Q")]

## Export - save to data folder
write.csv(RI_two, "../data/RI_2yr_flood_6riv.csv")








