## StreamLight output

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork","devtools"), require, character.only=T)

## Sites of interest that have StreamLight data
sites <- c("nwis_02336526",
           "nwis_01649190",
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
View(SL_split$nwis_01608500_input_output.txt)

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

## From NWIS_RiverSelection - no light for Santa Margarita or Pecos River
site_subset <- rbind(SF_df[which(SF_df$site_name == "nwis_02336526" & SF_df$Year %in% c(2015,2016)),],
                     SF_df[which(SF_df$site_name == "nwis_01649190" & SF_df$Year %in% c(2010,2011)),],
                     SF_df[which(SF_df$site_name == "nwis_07191222" & SF_df$Year %in% c(2009,2010)),],
                     SF_df[which(SF_df$site_name == "nwis_01608500" & SF_df$Year %in% c(2012,2013)),])

site_subset_split <- split(site_subset, site_subset$.id)

lapply(site_subset_split, function(x) sum(is.na(x$PAR_surface))) # all 0
lapply(site_subset_split, function(x) sum(is.na(x$PAR_turb))) # 108, 126, 153, 46

####################################################
## For remaining sites - Santa Margarita and Pecos
####################################################
#Use the devtools packge to install StreamLightUtils
#devtools::install_github("psavoy/StreamLightUtils")
#devtools::install_github("psavoy/StreamLight")
library("StreamLightUtils")
library("StreamLight")

###########
## NLDAS
###########
#Make a table for the MODIS request
setwd("~/GitHub/RiverBiomass/code")
site_info <- readRDS("./rds files/NWIS_6siteinfo_subset_SL.rds")
requested_sites <- site_info[which(site_info$site_name %in% c("nwis_11044000","nwis_08447300")),]
req_table <- requested_sites[,c("site_name","lat","lon")]
colnames(req_table) <- c("Site_ID","Lat","Lon")

NLDAS_DL_bulk(save_dir = "~/GitHub/RiverBiomass/data/NLDAS",
              site_locs = req_table, startDate = "2011-12-31")

#List of successfully downloaded sites
NLDAS_list <- stringr::str_sub(list.files("~/GitHub/RiverBiomass/data/NLDAS"), 1, -11)

#Processing the downloaded NLDAS data
NLDAS_processed <- StreamLightUtils::NLDAS_proc(read_dir = "~/GitHub/RiverBiomass/data/NLDAS", NLDAS_list)

###########
## MODIS
###########
#First time - Modify table for the MODIS request
#MODIS_req_table <- req_table
#colnames(MODIS_req_table) <- c("ID","Latitude","Longitude")
#write.csv(MODIS_req_table, "NASA_Light_request.csv")

## Submit request following directions here: https://psavoy.github.io/StreamLight/articles/2%20Download%20and%20process%20MODIS%20LAI.html
setwd("~/GitHub/RiverBiomass/code")
MOD_unpack <- AppEEARS_unpack_QC(zip_file = "biomass-rivers.zip",
                   zip_dir = "../data",
                   c("nwis11044000","nwis08447300")
                   )
MOD_processed <- AppEEARS_proc(unpacked_LAI = MOD_unpack,
                               fit_method = "Gu",
                               plot=TRUE)

## Gu methods prevents the interpolation of the final ~30 days of 2016 in nwis11044000
## fill in the missing LAI_proc data with a final simple linear interpolation
MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date == "2017-01-01"),]$LAI_proc <- 0.2
nrow(MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date >= "2016-12-02" & MOD_processed$nwis11044000$Date <= "2017-01-01"),])


interp_LAI <- rev(approx(x = c(MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date == "2016-12-02"),]$LAI_proc,
             MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date == "2017-01-01"),]$LAI_proc),
       y = c(0,0),
       n = 31)$x)
## replace
MOD_processed$nwis11044000[which(MOD_processed$nwis11044000$Date >= "2016-12-02" & MOD_processed$nwis11044000$Date <= "2017-01-01"),]$LAI_proc <- interp_LAI
View(MOD_processed$nwis11044000) ## scroll to bottom to confirm

######################
## Use stream_light
######################
## Prepare driver file
site_locs <- req_table
site_locs$epsg_crs <- 4326
names(MOD_processed) <- c("nwis_08447300","nwis_11044000")

## make driver file
source("../data/NLDAS/make_driver_mod.R")
RivBiomass_driver <- make_driver(req_table, NLDAS_processed, MOD_processed,
                                 write_output = FALSE, save_dir = NULL)

## Extract tree height per site
source("../data/NLDAS/extract_height_mod.R")
Pecos_TH <- extract_height(
  Site_ID = site_locs$Site_ID[1], 
  Lat = site_locs$Lat[1],
  Lon = site_locs$Lon[1],
  site_crs = 4326,
  simard_loc = "../data/NLDAS/simard2011.asc"
)

SantaMar_TH <- extract_height(
  Site_ID = site_locs$Site_ID[2], 
  Lat = site_locs$Lat[2],
  Lon = site_locs$Lon[2],
  site_crs = 4326,
  simard_loc = "../data/NLDAS/simard2011.asc"
)

## Estimate width based on mode of discharge
site_data <- readRDS("./rds files/NWIS_6site_subset_SL.rds")
median_Q <- site_data %>%
  group_by(site_name) %>%
  summarise_at(.vars = "Q", .funs = median, na.rm=TRUE)
colnames(median_Q) <- c("site_name","median_Q")
site_info <- left_join(site_info, median_Q, by="site_name")

site_info$width_calc <- site_info$dvqcoefs.a*(site_info$median_Q^site_info$dvqcoefs.b)
Pecos_width <- site_info[which(site_info$site_name == "nwis_08447300"),]$width_calc
SantaMar_width <- site_info[which(site_info$site_name == "nwis_11044000"),]$width_calc

## Run the model
Pecos_site <- site_locs[which(site_locs$Site_ID == "nwis_08447300"),]
SantaMar_site <- site_locs[which(site_locs$Site_ID == "nwis_11044000"),]

Pecos_modeled <- stream_light(RivBiomass_driver$nwis_08447300,
                              Lat = Pecos_site$Lat,
                              Lon = Pecos_site$Lon,
                              channel_azimuth = 45,
                              bottom_width = Pecos_width,
                              BH = 0.1,
                              BS = 100,
                              WL = 0.1,
                              TH = Pecos_TH$TH,
                              overhang = 0.1*Pecos_TH$TH,
                              overhang_height = NA,
                              x_LAD = 1)
## Pecos missing PAR_surface for 1st 20 days in 2012 because LAI & PAR_bc is missing
## First 20 days = PAR_inc
View(Pecos_modeled[is.na(Pecos_modeled$LAI),])
Pecos_modeled[is.na(Pecos_modeled$LAI),]$PAR_surface <- Pecos_modeled[is.na(Pecos_modeled$LAI),]$PAR_inc

SantaMar_modeled <- stream_light(RivBiomass_driver$nwis_11044000,
                              Lat = SantaMar_site$Lat,
                              Lon = SantaMar_site$Lon,
                              channel_azimuth = 80,
                              bottom_width = SantaMar_width,
                              BH = 0.1,
                              BS = 100,
                              WL = 0.1,
                              TH = SantaMar_TH$TH,
                              overhang = 0.1*SantaMar_TH$TH,
                              overhang_height = NA,
                              x_LAD = 1)

## Convert 0s to NA for PAR_surface
SantaMar_modeled$PAR_surface[which(SantaMar_modeled$PAR_surface == 0)] <- NA
Pecos_modeled$PAR_surface[which(Pecos_modeled$PAR_surface == 0)] <- NA

## Daily mean light using slightly modified function above

meandaily_PARv2 <- function(y){
  df <- y %>%
    group_by(jday) %>%
    summarize_at(.vars = c("DOY","Year","PAR_surface"), .funs = mean, na.rm = TRUE)
  
  df$origin <- as.Date(paste0(df$Year, "-01-01"),tz = "UTC") - days(1)
  df$Date <- as.Date(df$DOY, origin = df$origin, tz = "UTC") 
  
  return(df)
}

SantaMar_daily <- meandaily_PARv2(SantaMar_modeled)
Pecos_daily <- meandaily_PARv2(Pecos_modeled)

## View
ggplot(SantaMar_daily, aes(Date, PAR_surface))+geom_line()
ggplot(Pecos_daily, aes(Date, PAR_surface))+geom_line()

##################
## SAVE
###############
## Add SantaMar and Pecos to the dataframe
colnames(site_subset);colnames(Pecos_daily)
joint_cols <- c("site_name","jday","DOY","Year","Date","PAR_surface")

site_subset <- site_subset[,joint_cols]
Pecos_daily$site_name <- "nwis_08447300"; SantaMar_daily$site_name <- "nwis_11044000"
Pecos_daily <- Pecos_daily[,joint_cols]
SantaMar_daily <- SantaMar_daily[,joint_cols]

all_site_light <- rbind(site_subset, SantaMar_daily, Pecos_daily)

## Save
saveRDS(all_site_light, "./rds files/StreamLight_daily_6riv_all.rds")







