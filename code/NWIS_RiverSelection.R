## Subset data from hypoxia database that is already linked to NHD
## JRB

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork"), require, character.only=T)

############################
## To create linked file
############################

## Import site data from Appling
setwd("../data")
site <- fread("site_data.csv")
names(site)

## Import hypoxia dataset and subset to Appling
hyp <- fread("GRDO_GEE_HA_NHD_2021_02_07.csv")
PC <- hyp[which(hyp$DB_Source == "PC"),]
names(PC)

## Merge both based on siteID
head(site$site_name); head(PC$SiteID)
colnames(PC)[which(colnames(PC) == "SiteID")] <- "site_name"

df <- left_join(site, PC, by="site_name")

## Export
write.csv(df, "PC_site_attribs.csv")

##########################################
## Re-import for site selection
##########################################
setwd("../data")
df <- fread("PC_site_attribs.csv")
colnames(df)

## subset
sub <- df[,c("site_name","long_name","ORD_STRA","NHD_STREAMORDE",
              "Start_time","End_time","n_time",
              "US_state")]

## Need multiple years of data
sub$Start_year <- year(sub$Start_time)
sub$End_year <- year(sub$End_time)
sub$n_years <- sub$End_year - sub$Start_year

s <- sub[which(sub$n_years >=2),]

## Instead choose sites based on grouping of NHD_STREAMORDE
## low (1,2,3), mid(4,5,6), high(7,8,9)
s_low <- s[which(s$NHD_STREAMORDE %in% c(1,2,3)),]
s_mid <- s[which(s$NHD_STREAMORDE %in% c(4,5,6)),]
s_high <- s[which(s$NHD_STREAMORDE %in% c(7,8,9)),]


##################################################
## Further select sites based on data quality
##################################################
## Import
NWIS <- read.table("daily_predictions.tsv", sep='\t', header = TRUE)
NWIS$date <- as.POSIXct(as.character(NWIS$date), format="%Y-%m-%d")
head(NWIS)

## Subset columns and sites
NWIS_sub <- NWIS[,c("site_name","date","GPP","GPP.lower","GPP.upper", "GPP.Rhat",
                    "ER","ER.lower","ER.upper","K600","K600.lower","K600.upper",
                    "temp.water","discharge","shortwave","velocity")]
colnames(NWIS_sub) <- c("site_name","date","GPP","GPP.lower","GPP.upper", "GPP.Rhat",
                        "ER","ER.lower","ER.upper","K600","K600.lower","K600.upper",
                        "temp","Q","light","velocity")

## Subset to sites in s (with at least 3 years of data)
NWIS_sub <- NWIS_sub[which(NWIS_sub$site_name %in% s$site_name),]

## Identify which sites have the most continuous data
NWIS_sub$doy <- yday(NWIS_sub$date)
NWIS_sub$year <- year(NWIS_sub$date)

## count days per year
dat_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  count()
## identify the max day gap per year
gap_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  mutate(gap = doy - lag(doy, default=doy[1]))
maxgap <- gap_per_year %>%
  group_by(site_name, year) %>%
  summarize_at(.vars = "gap", .funs = max)
## subset for sites with a max gap of 14 days
sub_by_gap <- maxgap[which(maxgap$gap <= 14),]
length(levels(as.factor(sub_by_gap$site_name)))
## merge with number of days per year
sub_by_gap <- merge(sub_by_gap, dat_per_year, by=c("site_name","year"))
## at least 300 days per year
sub_by_gap <- sub_by_gap[which(sub_by_gap$n > 275),]
sub_by_gap_sum <- sub_by_gap %>% group_by(site_name) %>% count()
high_q <- sub_by_gap_sum[which(sub_by_gap_sum$n >= 2),]

## Subset NWIS_sub
TS <- NWIS_sub[which(NWIS_sub$site_name %in% high_q$site_name),]
TS_site <- s[which(s$site_name %in% high_q$site_name),]

## Assign a stream order classfication
TS_site$order_group <- "NA"
TS_site[which(TS_site$NHD_STREAMORDE %in% c(1,2)),]$order_group <- "small"
TS_site[which(TS_site$NHD_STREAMORDE %in% c(3,4)),]$order_group <- "mid"
TS_site[which(TS_site$NHD_STREAMORDE %in% c(5,6)),]$order_group <- "large"

##################################################
## Choose rivers with higher productivity rates
################################################
TS$GPP_temp <- TS$GPP
TS[which(TS$GPP < 0),]$GPP_temp <- sample(exp(-6):exp(-4), 1)
TS_gpp <- TS %>%
  group_by(site_name) %>%
  summarise_at(.vars = "GPP_temp", .funs = c(mean, max))
colnames(TS_gpp) <- c("site_name","GPP_mean","GPP_max")

TS_site <- left_join(TS_site, TS_gpp, by="site_name")

## View sites by river order
TS_highmean <- TS_site[which(TS_site$GPP_mean > 0.25),c("site_name","long_name","order_group")]
View(sub_by_gap[which(sub_by_gap$site_name %in% TS_highmean$site_name),])

## plot
ggplot(TS[which(TS$site_name == "nwis_01645762" & TS$year %in% c(2012,2013)),], aes(date, GPP))+
  geom_line()

## small: nwis_08180700 2010-2011 (Medina Riv, TX), nwis_10129900 2015-2016 (Silver Creek, UT)
## mid: nwis_03058000 2014-2015 (West Fork, WV), nwis_01649500 2012-2013 (Anacostia Riv, MD)
## large:  nwis_14211010 2009-2010 (Clackamas Riv, OR), nwis_02234000 2013-2014 (St. John Riv, FL)

site_subset <- rbind(TS[which(TS$site_name == "nwis_08180700" & TS$year %in% c(2010,2011)),],
               TS[which(TS$site_name == "nwis_10129900" & TS$year %in% c(2015,2016)),],
               TS[which(TS$site_name == "nwis_03058000" & TS$year %in% c(2014,2015)),],
               TS[which(TS$site_name == "nwis_01649500" & TS$year %in% c(2012,2013)),],
               TS[which(TS$site_name == "nwis_14211010" & TS$year %in% c(2012,2013)),],
               TS[which(TS$site_name == "nwis_02234000" & TS$year %in% c(2013,2014)),],
               TS[which(TS$site_name == "nwis_04137500" & TS$year %in% c(2010,2011)),],
               TS[which(TS$site_name == "nwis_01645762" & TS$year %in% c(2012,2013)),],
               TS[which(TS$site_name == "nwis_01649190" & TS$year %in% c(2010,2011)),])

TS_site_subset <- df[which(df$site_name %in% site_subset$site_name),]

###################################################
## Check other covariate data quality
###################################################
site_sub_list <- split(site_subset, site_subset$site_name)

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
    
    ggplot(df, aes(date, light))+
      geom_point(size=2, color="darkgoldenrod3")+
      labs(y="Incoming Light", x="Date")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.text = element_text(size=12),
            axis.title = element_text(size=12)),
    ncol=1, align="hv")
  
  return(p)
  
}

plotting_covar(site_sub_list$nwis_01645762)

setwd("~/GitHub/RiverBiomass/figures/Site Covariate Plots")
lapply(site_sub_list, function(x) ggsave(plot = plotting_covar(x),filename = paste(x$site_name[1],"covar.jpg",sep = "")))


###########################
## Export
###########################

## NWIS site subset
setwd("~/GitHub/RiverBiomass/code")
saveRDS(site_subset, "./rds files/NWIS_9site_subset.rds")
saveRDS(TS_site_subset, "./rds files/NWIS_9siteinfo_subset.rds")



##############################
## Plot for talk
###########################

df <- site_sub_list$nwis_14211010
ratio_QL <- max(df$light)/max(df$Q)
GPP_plot <- ggplot(df, aes(date, GPP))+
  geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="chartreuse4")+
      geom_point(color="chartreuse4", size=2)+geom_line(color="chartreuse4", size=1)+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'))+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12))+
  geom_vline(xintercept = as.POSIXct("2013-01-01", format="%Y-%m-%d"),size=1, linetype="dashed")
    
data_plot <- ggplot(df, aes(date, Q*ratio_QL))+
      geom_point(data=df, aes(date, light), size=1.5, color="darkgoldenrod3")+
      geom_line(size=1, color="deepskyblue4")+
      scale_y_continuous(sec.axis = sec_axis(~./ratio_QL, name=expression("Daily Q (cms)")))+
      labs(y=expression('Daily PPFD'))+# ('*~mu~mol~ m^-2~d^-1*')'), x="Date")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text = element_text(size=12),
            axis.title.y.left = element_text(size=12, color="darkgoldenrod3"),
            axis.title.y.right = element_text(size=12, color="deepskyblue4"),
            axis.text.x = element_text(angle=25, hjust = 1),
            strip.background = element_rect(fill="white", color="black"),
            strip.text = element_text(size=15))+
  geom_vline(xintercept = as.POSIXct("2013-01-01", format="%Y-%m-%d"),size=1, linetype="dashed")
    
GPP_plot + data_plot + plot_layout(ncol = 1)

ggplot(df, aes(light, GPP))+
  geom_point(size=1.5)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x="Daily PPFD")+
  theme_bw()+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18))


