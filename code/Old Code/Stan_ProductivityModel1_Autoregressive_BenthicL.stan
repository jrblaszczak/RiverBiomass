 
    data {
    int Ndays; // number of days
    vector [Ndays] light; // relativized to max value; unitless
    vector [Ndays] GPP; // estimates from posterior probability distributions; g O2 m-2 d-1
    vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions; g O2 m-2 d-1
    vector [Ndays] tQ; // standardized discharge; unitless
    vector [Ndays] depth; // depth estimate
    vector [Ndays] turb; // mean daily turbidity
    }
    
    parameters {
    // Growth & Loss Parameters
    real l_pred_GPP [Ndays]; // predicted GPP (log scale); g O2 m-2 d-1
    real<lower=0> phi; // autoregressive coefficient; unitless
    real<lower=0> alpha; // growth term; g O2 m-2 d-1
    real<upper=0> beta; // discharge to loss conversion term; g O2 m-2 d-1
    
    // Light adjustment
    real<lower=0> a; // light attenuation coefficient to inform Kd
    
    // Error parameters
    real<lower=0> sig_p; // sigma processes error
    }
    
    transformed parameters {
    real ben_light [Ndays];
    
    for(i in 1:Ndays){
    ben_light[i]=light[i]*exp(-1*a*turb[i]*depth[i]);
    }
    
    }
    
    
    model {
        
    // Initial value
    l_pred_GPP[1] ~ normal(log(GPP[1]), 1e-6);

    // Process model
    for (j in 2:(Ndays)) {
    l_pred_GPP[j] ~ normal(phi*l_pred_GPP[j-1] + alpha*ben_light[j] + beta*tQ[j], sig_p);
    }
    
    // Observation model
    for (j in 2:(Ndays)) {
        GPP[j] ~ normal(exp(l_pred_GPP[j]), GPP_sd[j])T[0,];
    }
    
    // Error priors
    sig_p ~ normal(0,2)T[0,];
    
    // Param priors (weakly informative)  
    phi ~ normal(0,1)T[0,];
    alpha ~ normal(0,1)T[0,];
    beta ~ normal(0,1)T[,0];
    a ~ normal(0,1)T[0,];
    
    }
    
