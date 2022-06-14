##==============================================================================
## Script to check covariate data quality & relationships
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","reshape2","ggExtra","patchwork","grid","gridExtra"), require, character.only=T)

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

## Add short name
data$short_name <- revalue(as.character(data$site_name), replace = c("nwis_02336526"="Proctor Creek, GA",
                                                                     "nwis_01649190"="Paint Branch, MD",
                                                                     "nwis_07191222"="Beaty Creek, OK",
                                                                     "nwis_01608500"="S. Br. Potomac River, WV",
                                                                     "nwis_11044000"="Santa Margarita River, CA",
                                                                     "nwis_08447300"="Pecos River, TX"))

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
data[which(data$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## split
site_sub_list <- split(data, data$site_name)

## Plot
plotting_covar <- function(x) {
  
  df <- x
  
  gpp_plot <- ggplot(df, aes(date, GPP))+
      geom_point(color="chartreuse4", size=2)+
      geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$short_name[1])+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12))
    
    Q_plot <- ggplot(df, aes(date, Q))+
      geom_line(size=1.5, color="deepskyblue4")+
      labs(y="Q (cms)")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12))
    
    Temp_plot <- ggplot(df, aes(date, temp))+
      geom_line(size=1.5, color="#A11F22")+
      labs(y="Water Temp (C)")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12))
    
    Light_plot <- ggplot(df, aes(date, PAR_surface))+
      geom_point(size=2, color="darkgoldenrod3")+
      labs(y="Water surface light", x="Date")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.text = element_text(size=12),
            axis.title = element_text(size=12))
  
  
  
  p <- gpp_plot/Q_plot/Temp_plot/Light_plot
 
  return(p)
  
}

plotting_covar(site_sub_list$nwis_11044000)

setwd("~/GitHub/RiverBiomass/figures and tables/Site Covariate Plots")
lapply(site_sub_list, function(x) ggsave(plot = plotting_covar(x),filename = paste(x$site_name[1],"covar.jpg",sep = "")))

#######################################
## Light and temperature correlations
#######################################
lapply(site_sub_list, function(x) plot(x$PAR_surface, x$temp))






