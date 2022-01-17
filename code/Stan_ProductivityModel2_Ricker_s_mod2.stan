

data {
  int Ndays; // number of days
  vector [Ndays] light; // relativized to max value
  vector [Ndays] GPP; // mean estimates from posterior probability distributions
  vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions
  vector [Ndays] tQ; // standardized discharge
}

parameters {
  // Disturbance (persistence) parameters
  real<lower=0> c; // estimate of Qcrit
  real<lower=0> s; // steepness of the transition from P=1 to P=0
  
  // Logistic growth parameters  
  real B [Ndays]; // Biomass
  real r; // growth rate
  real lambda; // r/K
  
  // Error parameters
  real<lower=0> sig_p; // sigma processes error
  real<lower=0> sig_o; // sigma observation error
}


transformed parameters {
  real pred_GPP [Ndays];
  real P [Ndays];
  
  for(i in 1:Ndays){
    P[i] = exp(-exp(s*100*(tQ[i]-c)));
    pred_GPP[i] = light[i]*exp(B[i]);
  }
  
}


model {
  
  // Initial value
  B[1] ~ normal(log(GPP[1]/light[1]), 1);
  
  // Process Model
  for (j in 2:(Ndays)){
    B[j] ~ normal((B[(j-1)] + r + lambda*exp(B[(j-1)]))*P[j], sig_p);
  }
  
  // Observation model
  for (j in 2:(Ndays)) {
    GPP[j] ~ normal(pred_GPP[j], sig_o)T[0,];
  }
  
  // Error priors
  sig_p ~ normal(0,2)T[0,];
  sig_o ~ normal(mean(GPP_sd), sd(GPP_sd))T[0,];
  
  // Param priors
  c ~ normal(0.5,0.25)T[0,];
  s ~ normal(1.5,1)T[0,];
  r ~ normal(0,1);
  lambda ~ normal(0,1)T[,0];
  
}






