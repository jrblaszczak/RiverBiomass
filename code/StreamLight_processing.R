## StreamLight output

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork"), require, character.only=T)

## Sites of interest that have StreamLight data
sites <- c("nwis_02336526",
           "nwis_01649190",
           "nwis_14206950",
           "nwis_07191222",
           "nwis_01608500")
sites_files <- rep(NA, length(sites))
for(i in 1:length(sites)){
  sites_files[i] <- paste(sites[i],"_input_output.txt", sep="")
}

## Import streamlight
setwd("~/GitHub/RiverBiomass/data/StreamLight model inputs and outputs/individual files")

SL <- ldply(sites_files, function(filename) {
  d <- read.table(filename, header = T, sep = "\t")
  d$file <- filename
  return(d)
})

## take the mean daily incoming PAR at the surface
SL_split <- split(SL, SL$file)
View(SL_split$nwis_14206950_input_output.txt)

meandaily_PAR <- function(y){
  df <- y %>%
  group_by(jday) %>%
  summarize_at(.vars = c("DOY","Year","PAR_surface","PAR_turb"), .funs = mean, na.rm = TRUE)
  
  df$origin <- as.Date(paste0(df$Year, "-01-01"),tz = "UTC") - days(1)
  df$Date <- as.Date(df$DOY, origin = df$origin, tz = "UTC") 
  
  return(df)
}

SL_daily <- lapply(SL_split, function(x) meandaily_PAR(x))
tail(SL_daily$nwis_01608500_input_output.txt)

## visualize
ggplot(SL_daily$nwis_01608500_input_output.txt, aes(Date, PAR_surface))+
  geom_point()+
  geom_point(aes(Date, PAR_turb), color="blue")


###########################
## Subset and evaluate NA
###########################
SF_df <- ldply(SL_daily, data.frame)

## add site_name
SF_df$site_name <- substr(SF_df$.id, 1, nchar(SF_df$.id)-17)

## From NWIS_RiverSelection
site_subset <- rbind(SF_df[which(SF_df$site_name == "nwis_02336526" & SF_df$Year %in% c(2015,2016)),],
                     SF_df[which(SF_df$site_name == "nwis_01649190" & SF_df$Year %in% c(2010,2011)),],
                     SF_df[which(SF_df$site_name == "nwis_14206950" & SF_df$Year %in% c(2013,2014)),],
                     SF_df[which(SF_df$site_name == "nwis_07191222" & SF_df$Year %in% c(2009,2010)),],
                     SF_df[which(SF_df$site_name == "nwis_01608500" & SF_df$Year %in% c(2012,2013)),])
                     #SF_df[which(SF_df$site_name == "nwis_11044000" & SF_df$Year %in% c(2015,2016)),])

site_subset_split <- split(site_subset, site_subset$.id)

lapply(site_subset_split, function(x) sum(is.na(x$PAR_surface))) # all 0
lapply(site_subset_split, function(x) sum(is.na(x$PAR_turb))) # 108, 126, 153, 46, 56

## Save
setwd("~/GitHub/RiverBiomass/data")

saveRDS(site_subset, "StreamLight_daily_6riv.rds")






