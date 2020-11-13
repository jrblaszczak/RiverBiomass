## Persistence simulations

# Simulate a oversimplified hydrograph
Q <- rnorm(350, mean=10, sd=0.1)
Q_peaks <- sample(1:350, 10)
Q[Q_peaks] <- rnorm(length(Q_peaks), mean=11.5, sd=0.2)
# Visualize oversimplified hydrograph
plot(Q, ylim = c(9,12))

# Create persistence vector P and set steepness of the curve (s) and critical discharge (c)
P <- numeric(length(Q))
P[1] <- 1
s <- 3
c <- 11
for(i in 2:length(Q)){
  P[i] = exp(-exp(s*(Q[i] - c)))
}
#P #P is never 0 and therefore can serve as refuge biomass but cannot be modeled using a binomial fit

# Repeat with a steeper and shallower curve (higher and lower s values)
P_steeper <- numeric(length(Q)); P_shallower <- numeric(length(Q))
P_steeper[1] <- 1; P_shallower[1] <- 1
s_steep <- 10; s_shallow <- 1

for(i in 2:length(Q)){
  P_steeper[i] = exp(-exp(s_steep*(Q[i] - c)))
  P_shallower[i] = exp(-exp(s_shallow*(Q[i] - c)))
}

## Visualize
library(ggplot2)
library(cowplot)

df <- as.data.frame(cbind(Q, P, P_steeper, P_shallower))

plot_grid(
  ggplot(df, aes(Q, P))+
    geom_line(lwd=2)+
    geom_line(aes(Q, P_steeper), color="purple",lwd=1)+
    geom_line(aes(Q, P_shallower), color="red")+
    theme_classic()+
    scale_y_continuous(limits=c(0,1)),
  
  qplot(seq_along(Q), Q)+
    geom_point(color="blue")+
    theme_classic()+
    geom_hline(yintercept = c), #horizontal line to represent critQ
  
  ncol=1)

# Notes: I tried a range of steeper s values than 10 and nothing budged.
# The curve only seemed to change if the value was lower than the original s
