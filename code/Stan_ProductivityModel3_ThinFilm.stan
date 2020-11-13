
    data {
    int Ndays; // number of days
    vector [Ndays] light; // relativized to max value
    vector [Ndays] GPP; // estimates from posterior probability distributions
    vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions
    vector [Ndays] tQ; // standardized discharge
    }
    
    parameters {
    // Disturbance (persistence) parameters
    real<lower=0> c; // estimate of Qcrit
    real<lower=0> s; // steepness of the transition from P=1 to P=0
    
    // Growth parameters
    real<lower=0> alpha; // light conversion
    real<lower=0> gamma; // multiplied by b it is the intercept of density dependent death
    vector [Ndays] N; // ratio of biomass to carrying capacity (B/K)
    
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
    
    P[i]=exp(-exp(s*(tQ[i]-c)));
    
    b[i] = alpha*light[i];
    if (i < 21){
    ant_b[i] = mean(b[1:i]);
    } else { ant_b[i] = mean(b[(i-20):i]); // length can be a parameter but currently fixed
    }
    
    pred_GPP[i] = b[i] * exp(N[i]);
    
    }
    
    }
    
    model {
    // Process model
    for (j in 2:(Ndays)){
    N[j] ~ normal((N[(j-1)] + (b[j] - ant_b[j]*N[(j-1)]*(gamma+(1-gamma)*N[(j-1)])))*P[j], sig_p);
    }
    
    // Observation model
    for (j in 2:(Ndays)) {
        GPP[j] ~ normal(exp(pred_GPP[j]), GPP_sd[j])T[0,];
    }
    
    // Error priors
    sig_p ~ normal(0,2);
    
    // Param priors
    alpha ~ normal(0,50);
    gamma ~ beta(1,1);
    
    
    
    }
    
    
    
