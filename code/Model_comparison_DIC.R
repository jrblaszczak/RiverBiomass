###Calcualte DIC for model fits


CalcDIC<-function(P,Y){
  meanP <- apply(P,2,mean) ## mean model prediction
  logLikelihood1 <- sum(dbinom(Y,1,meanP,log=T)) ## Likelihood of data given mean prediction
  logLikevec2<-numeric(length(P[,1]))
  for (i in 1:length(P[,1])){logLikevec2[i]<-sum(dbinom(Y,1,P[i,],log=T))} #Likelihood of data given prediction from each iteration
  pDIC <- 2*(logLikelihood1 - mean(logLikevec2))
  dic <- -2*logLikelihood1 + 2*pDIC
  return(dic)
}



