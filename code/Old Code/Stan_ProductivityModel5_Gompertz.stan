
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
    real B [Ndays]; // Biomass; g m-2
    real beta_0; // r 
    real beta_1; // r/K

    // Error parameters
    real<lower=0> sig_p; // sigma processes error
    }
    
    transformed parameters {
    real pred_GPP [Ndays];
    real P [Ndays];
    
    for(i in 1:Ndays){
    P[i]=exp(-exp(s*(tQ[i]-c)));
    }
    
    } 
    
    model {
    
    // Process Model
    for (j in 2:(Ndays)){
        B[j] ~ normal((beta_0 + beta_1*B[(j-1)])*P[j], sig_p);
      }
 
    // Observation model
    for (j in 2:(Ndays)) {
        GPP[j] ~ normal(light[j]*exp(B[j]), GPP_sd[j])T[0,];
    }
 
    // Error priors
    sig_p ~ normal(0,2)T[0,];
    
    // Param priors
    c ~ rayleigh(0.5);
    s ~ normal(0,100);
    beta_0 ~ normal(0,1);
    beta_1 ~ normal(0,1);
    
    }
    
    
    
    
