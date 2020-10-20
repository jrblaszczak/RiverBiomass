## Klamath Data Example - No disturbance

## Import
km <- read.table("KlamMetab.txt", header=T)
names(km)
## Fix date
km$date <- as.POSIXct(as.character(km$date), format="%Y-%m-%d")
head(km$date)
## Sites
levels(as.factor(km$site))
l <- split(km, km$site)

## Choose one year and site
d <- subset(km, site == "WE" & year == "2012")
## Start on June 1st
d <- d[which(d$date >= "2012-06-15"),]
ggplot(d, aes(jday, GPP))+geom_point(aes(color=year), alpha=0.5)+
  geom_line(aes(jday, GPP, group = year, color=year), alpha=0.5)

rel_LQT <- function(x){
  x$light_rel <- x$solar_rad/max(x$solar_rad)
  x$temp_rel <- x$mean_wtemp/max(x$mean_wtemp)
  x$tQ <- x$cms/max(x$cms) #(x$cms-mean(x$cms))/sd(x$cms)
  x$Q95 <- 1 #ifelse(x$cms>quantile(x$cms, probs = 0.99), yes=1, no=0)
  x<-x[order(x$date),]
  x <- x[,c("site","date","jday","GPP","light_rel","temp_rel","tQ","Q95")]
  return(x)
}

r <- rel_LQT(d)
df <- na.omit(r)

rm(d,km,l,r)
