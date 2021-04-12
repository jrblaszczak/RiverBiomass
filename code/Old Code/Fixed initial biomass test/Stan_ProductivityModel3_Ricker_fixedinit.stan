
    data {
    int Ndays; // number of days
    vector [Ndays] light; // relativized to max value
    vector [Ndays] GPP; // mean estimates from posterior probability distributions
    vector [Ndays] GPP_sd; // sd estimates from posterior probability distributions
    vector [Ndays] tQ; // standardized discharge
    real B_int; // Initial biomass value where B_int = log(GPP[1]/light[1])
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
    }
    
    transformed parameters {
    real pred_GPP [Ndays];
    real P [Ndays];

    for(i in 1:Ndays){
    P[i]=exp(-exp(s*(tQ[i]-c)));
    pred_GPP[i] =light[i]*exp(B[i]);
    }
    
    }
    
   
    model {

    // Process Model
    for (j in 2:(Ndays)){
      
      if(j == 2){
         B[j] ~ normal((B_int + r + lambda*exp(B_int))*P[j], sig_p);
      }
      else{
        B[j] ~ normal((B[(j-1)] + r + lambda*exp(B[(j-1)]))*P[j], sig_p);
      }
    }
    
    
    // Observation model
    for (j in 2:(Ndays)) {
        GPP[j] ~ normal(pred_GPP[j], GPP_sd[j])T[0,];
    }
 
    // Error priors
    sig_p ~ normal(0,2)T[0,];
    
    // Param priors
    r ~ normal(0,1);
    lambda ~ normal(0,1);
    c ~ rayleigh(0.5);
    s ~ normal(0,50);
 
    }






