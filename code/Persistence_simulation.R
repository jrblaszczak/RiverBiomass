## Persistence simulations

# Simulate a oversimplified hydrograph
Q <- rnorm(350, mean=10, sd=0.1)
Q_peaks <- sample(1:350, 10)
Q[Q_peaks] <- rnorm(length(Q_peaks), mean=11.5, sd=0.3)
# Visualize oversimplified hydrograph
plot(Q, ylim = c(9,12))

# Create persistence vector P and set steepness of the curve (s) and critical discharge (c)
P <- numeric(length(Q))
P[1] <- 1
s <- 10
c <- 11

for(i in 2:length(Q)){
  P[i] = exp(-exp(s*(Q[i] - c)))
}

## Visualize
df <- as.data.frame(cbind(Q,P))
fit <- glm(df$P ~ df$Q, data=df, family = )
newdat <- data.frame(Q=seq(min(df$Q), max(df$Q),len=350))
newdat$P = predict(fit, newdata=newdat, type="response")
plot(P ~ Q, data=df, col="red4")
lines(P ~ Q, newdat, col="green4", lwd=2)



plot(Q, P, ylim = c(0,1))


mtcars <- mtcars
fit = glm(vs ~ hp, data=mtcars, family=binomial)
predicted= predict(fit, newdata=mtcars, type="response")
plot(vs~hp, data=mtcars, col="red4")
lines(mtcars$hp, predicted, col="green4", lwd=2)






