## Testing different distributions

# Beta
p <- seq(0,1,length=100)
plot(p,dbeta(p,1,1), ylab="density",type="l",ylim=c(0,10))
lines(p, dbeta(p,10,1), col="blue")

# Gamma
z <- seq(from=0, to=5, by=0.1)
a <- 1
gfunc <- dgamma(z, shape=a)
plot(z,gfunc,
     xlab="z",
     ylab="Probability distribution of [z|alpha,beta]",
     title("Probability Density Function of a Gamma Distribution"))

# Negative binomial
x<-seq(0,25,1)
nb<-dnbinom(x, size=3, prob=0.15)
plot(x, nb)

# Weibull
library(tidyr)
library(dplyr)
library(ggplot2)

x <- seq(0, 3, 0.1)
less_than_one <- dweibull(x, shape = 0.5, scale = 1)
one <- dweibull(x, shape = 1, scale = 1)
more_than_one <- dweibull(x, shape = 1.5, scale = 1)

tibble(less_than_one, one, more_than_one, x) %>% gather("k_parameter", "sample", 1:3) %>% ggplot(aes(x = x, y = sample, col = k_parameter))+
  geom_line(size = 1.4)+theme_classic()+
  scale_color_manual(values = c("#fec44f", "#999999", "#5ab4ac"))+ theme(legend.position = c(.2,.8))+facet_wrap(~k_parameter)

# Rayleigh
library(VGAM)

x <- seq(0, 3, 0.1)
ray_less_than_one <- drayleigh(x, scale = 0.5)
ray_one <- drayleigh(x, scale = 1)
ray_more_than_one <- drayleigh(x, scale = 1.5)

tibble(ray_less_than_one, ray_one, ray_more_than_one, x) %>% 
  gather("k_parameter", "sample", 1:3) %>% ggplot(aes(x = x, y = sample, col = k_parameter))+
  geom_line(size = 1.4)+theme_classic()+
  scale_color_manual(values = c("#fec44f", "#999999", "#5ab4ac"))+
  theme(legend.position = c(.2,.8))+facet_wrap(~k_parameter)+
  scale_x_continuous(limits=c(0,5))






