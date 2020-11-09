
    data {
    int Ndays; // number of days
    vector [Ndays] light; // relativized to max value
    vector [Ndays] GPP; // mean estimates from posterior probability distributions
    vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions
    vector [Ndays] tQ; // standardized discharge
    }
    
    parameters {
    // Disturbance (persistence) parameters
    real<lower=0, upper=1> c; // estimate of Qcrit
    real<lower=0, upper=50> s; // steepness of the transition from P=1 to P=0
    
    // Logistic growth parameters  
    real<lower=0> B [Ndays]; // Biomass; g m-2
    real<lower=0> r; // growth rate; d-1
    // real<lower=0> alpha; // light conversion; unitless
    real<lower=0> K; // carrying capacity; g m-2
    
    // Error parameters
    real<lower=0> sig_p; // sigma processes error
    }
    
    transformed parameters {
    real pred_GPP [Ndays];
    real P [Ndays];
    
    for(i in 1:Ndays){
    P[i]=exp(-exp(s*(tQ[i]-c)));
    pred_GPP[i] = light[i]*exp(B[i]);
    }
    
    } 
    
    model {
    
    // Process Model
    for (j in 2:(Ndays)){
    B[j] ~ normal((B[(j-1)]*exp(r*B[(j-1)]*(1-(B[(j-1)]/K))))*P[j], sig_p);
    }
    
    // Observation model
    GPP ~ normal(pred_GPP, GPP_sd);
    
    // Error priors
    sig_p ~ normal(0,2);
    
    // Param priors
    K ~ normal(0,30);
    r ~ normal(0,15);
    c ~ beta(10,1);
    
    }
    
    
    
    
