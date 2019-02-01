// Inference using three-compartment model of chlamydia testing and data from Natsal.

data {
  
  int N; // number of survey participants
  int N_strata; // number of strata
  int str[N]; // stratum of each participant
  real wt[N]; // survey weights
  int tested[N]; // was the participant tested for chlamydia in the last year?
  int symp[N]; // given the participant was tested, was the test prompted by symptoms?
  int diag[N]; // given the participant was tested, was the test positive?
  
}

parameters {
  
  real<lower=0,upper=1> p_symp; // proportion of infections that develop symptoms
  real<lower=0> lambda_slow; // natural recovery rate of untreated, asymptomatic infections
  
  vector<lower=0>[N_strata] foi; // force of infection
  
  vector<lower=0>[N_strata] scr; // screening rate
  real<lower=0> trt; // treatment rate in symptomatic people
  
}

transformed parameters{
  
  vector<lower=0>[N_strata] alpha_UA;
  vector<lower=0>[N_strata] alpha_AU;
  vector<lower=0>[N_strata] alpha_US;
  vector<lower=0>[N_strata] alpha_SU;
  
  vector<lower=0,upper=1>[N_strata] S; // steady-state proportion symptomatic-infected
  vector<lower=0,upper=1>[N_strata] A; // steady-state proportion symptomatic-infected
  vector<lower=0,upper=1>[N_strata] prev; // steady-state prevalence
  
  alpha_UA = foi * (1 - p_symp);
  alpha_AU = scr + lambda_slow;
  alpha_US = foi * p_symp;
  alpha_SU = trt + scr + lambda_slow;
  
  S = alpha_AU .* alpha_US ./ (alpha_AU .* alpha_US + alpha_SU .* (alpha_AU + alpha_UA));
  A = alpha_SU .* alpha_UA ./ (alpha_AU .* alpha_US + alpha_SU .* (alpha_AU + alpha_UA));
  prev = S + A;
}

model {
  
  // priors

  //////////////////
  // women
  //////////////////
//  p_symp ~ beta(1,1); // uninformative
  p_symp ~ beta(27, 90); // based on Geisler 2008 
//  p_symp ~ beta(18,13); // based on Lewis 2017 
  lambda_slow ~ normal(0.74, (0.89-0.61)/3.919928); // from Price 2013. 2*qnorm(0.975) = 3.919928


  //////////////////
  // same for men and women
  //////////////////
  foi ~ exponential(0.001); // uninformative
  scr ~ exponential(0.001); // uninformative
  trt ~ gamma(14,1); // based on Mercer 2007 and Lewis 2017
  
  // 19 ~ binomial(597, prev);
  //prev ~ beta(20, 579); // steady-state prevalence from Woodhall 2016
   
   for(i in 1:N){
     
     target += wt[i]*bernoulli_lpmf(tested[i] | 1 - exp(poisson_lpmf(0 | scr[str[i]] + S[str[i]] * trt)));
     if(tested[i] == 1){
        target += wt[i]*bernoulli_lpmf(symp[i] | (S[str[i]] * trt)/(scr[str[i]] + S[str[i]] * trt) );
        target += wt[i]*bernoulli_lpmf(diag[i] | (S[str[i]] * (trt+scr[str[i]]) + A[str[i]] * scr[str[i]] ) / (scr[str[i]] + S[str[i]] * trt) );
        
     }
     
   }
   
}

generated quantities{
  
  vector<lower=0,upper=1>[N_strata] p_test;
  vector<lower=0,upper=1>[N_strata] p_symp_given_test;
  vector<lower=0,upper=1>[N_strata] p_diag_given_test;
  
  int n_test_sim[N_strata];
  int n_diag_sim[N_strata];
  
  p_symp_given_test = (S * trt) ./ (scr + S * trt);
  p_diag_given_test = (S .* (trt+scr) + A .* scr ) ./ (scr + S * trt);
  
  for(i in 1:N_strata){
    p_test[i] = 1 - exp(poisson_lpmf(0 | scr[i] + S[i] * trt));
    n_test_sim[i] = poisson_rng(3388842 * p_test[i]);
    n_diag_sim[i] = binomial_rng(n_test_sim[i], p_diag_given_test[i]);
  }
  
}