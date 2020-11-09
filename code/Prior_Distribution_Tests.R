## Testing different distributions

# Beta
p <- seq(0,1,length=100)
plot(p,dbeta(p,1,1), ylab="density",type="l",ylim=c(0,10))
lines(p, dbeta(p,10,1), col="blue")



