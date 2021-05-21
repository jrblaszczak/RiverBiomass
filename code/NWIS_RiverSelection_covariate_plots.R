## Check covariate data quality

##############################
## Data Import & Processing ##
##############################
data <- readRDS("./rds files/NWIS_6site_subset_SL.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("./rds files/NWIS_6siteinfo_subset_SL.rds")

# Read in StreamLight processed data (Savoy)
SL <- readRDS("./rds files/StreamLight_daily_6riv.rds")
colnames(SL)[colnames(SL) == "Date"] <- "date"

## Join data and StreamLight
data <- left_join(data, SL, by=c("site_name", "date"))

## split
site_sub_list <- split(data, data$site_name)

## For only Santa Margarita River, set PAR_surface to light
# because no StreamLight available
site_sub_list$nwis_11044000$PAR_surface <- site_sub_list$nwis_11044000$light

plotting_covar <- function(x) {
  
  df <- x
  p <- plot_grid(
    
    ggplot(df, aes(date, GPP))+
      geom_point(color="chartreuse4", size=2)+
      geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$site_name[1])+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12)),
    
    ggplot(df, aes(date, Q))+
      geom_line(size=1.5, color="deepskyblue4")+
      labs(y="Q (cms)")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12)),
    
    ggplot(df, aes(date, temp))+
      geom_line(size=1.5, color="#A11F22")+
      labs(y="Water Temp (C)")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12)),
    
    ggplot(df, aes(date, PAR_surface))+
      geom_point(size=2, color="darkgoldenrod3")+
      labs(y="Water surface light", x="Date")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.text = element_text(size=12),
            axis.title = element_text(size=12)),
    ncol=1, align="hv")
  
  return(p)
  
}

plotting_covar(site_sub_list$nwis_02336526)

setwd("~/GitHub/RiverBiomass/figures/Site Covariate Plots")
lapply(site_sub_list, function(x) ggsave(plot = plotting_covar(x),filename = paste(x$site_name[1],"covar.jpg",sep = "")))
