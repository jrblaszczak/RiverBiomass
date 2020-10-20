data {
    int Ndays; // number of days
    vector [Ndays] light; // relativized to max value
    vector [Ndays] GPP; // estimates from posterior probability distributions
    vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions
    vector [Ndays] tQ; // standardized discharge
    }
    
    parameters {
    // Disturbance (persistence) parameters
    real<lower=0> beta_r; // refugia biomass after Qcrit > Q
    real<lower=0,upper=10> critQ; // estimate of Qcrit
    real<lower=0,upper=20> bQ; // steepness of the transition from P=1 to P=0
    
    // Growth parameters
    real<lower=0> alpha; // light conversion
    real<lower=0> gamma; // multiplied by b it is the intercept of density dependent death
    vector [Ndays] X; // ratio of biomass to carrying capacity (B/K)
    
    // Error parameters
    real<lower=0> sig_p; // sigma processes error
    }
    
    transformed parameters {
    // Persistence dynamics  
    real P [Ndays];
    
    // Growth
    vector [Ndays] b;
    vector [Ndays] ant_b; 
    real pred_GPP [Ndays];
    
    for (i in 1:(Ndays)){
    
    P[i]=exp(-exp(bQ*(tQ[i]-critQ)));
    
    b[i] = alpha*light[i];
    if (i < 21){
    ant_b[i] = mean(b[1:i]);
    } else { ant_b[i] = mean(b[(i-20):i]); // length can be a parameter
    }
    
    pred_GPP[i] = b[i] * X[i];
    
    }
    
    }
    
    model {
    // Process model
    for (j in 2:(Ndays)){
    X[j] ~ normal((X[(j-1)] + (b[j] - ant_b[j]*X[(j-1)]*(gamma+(1-gamma)*X[(j-1)])) - beta_r)*P[j] + beta_r, sig_p);
    }
    
    // Observation model  
    GPP ~ normal(pred_GPP, GPP_sd);
    
    // Error priors
    sig_p ~ normal(0,1);
    
    // Param priors
    bQ ~ normal(0,15);
    critQ ~ normal(0,15);
    alpha ~ normal(0,50);
    gamma ~ beta(1,1);
    beta_r ~ normal(0,1);
    }
    
    
    
