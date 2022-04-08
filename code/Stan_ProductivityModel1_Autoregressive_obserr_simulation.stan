
  data {
    int Ndays; // number of days
    vector [Ndays] light; // relativized to max value; unitless
    vector [Ndays] GPP; // estimates from posterior probability distributions; g O2 m-2 d-1
    real prior_sig_o_mean; // mean of daily GPP standard deviations
    real prior_sig_o_sd; // sd of daily GPP standard deviations
    vector [Ndays] tQ; // standardized discharge; unitless
    }

  parameters {
    // Growth & Loss Parameters
    real l_pred_GPP [Ndays]; // predicted GPP (log scale); g O2 m-2 d-1
    real<lower=0> phi; // autoregressive coefficient; unitless
    real<lower=0> alpha; // growth term; g O2 m-2 d-1
    real<upper=0> beta; // discharge to loss conversion term; g O2 m-2 d-1
  
    // Error parameters
    real<lower=0> sig_p; // sigma processes error
    real<lower=0> sig_o; // sigma observation error
  }


  model {
  
  // Initial value
  l_pred_GPP[1] ~ normal(log(GPP[1]), 0.1);
  
  // Process model
  for (j in 2:(Ndays)) {
    l_pred_GPP[j] ~ normal(phi*l_pred_GPP[j-1] + alpha*light[j] + beta*tQ[j], sig_p);
  }
  
  // Observation model
  for (j in 2:(Ndays)) {
    GPP[j] ~ normal(exp(l_pred_GPP[j]), sig_o)T[0,];
  }
  
  // Error priors
  sig_p ~ normal(0,2)T[0,];
  sig_o ~ normal(prior_sig_o_mean, prior_sig_o_sd)T[0,];
  
  // Param priors (weakly informative)  
  phi ~ normal(0,1)T[0,];
  alpha ~ normal(0,1)T[0,];
  beta ~ normal(0,1)T[,0];
  
}

