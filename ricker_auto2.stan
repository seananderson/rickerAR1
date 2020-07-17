data {
  int<lower=0> N;
  vector[N] obs_logR;
  vector[N] obs_S;
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> sigma_obs;
  real<lower=-1, upper=1> rho;
}
transformed parameters {
  real<lower=0> SigAR;
  SigAR = sigma_obs * sqrt(1-rho^2);
}
model {
  vector[N] logRS;
  vector[N] pred_logR;
  vector[N] eps;
  for (i in 1:N) {
    logRS[i] = alpha - beta * obs_S[i];
    if (i == 1) {
      pred_logR[i] = logRS[i] + log(obs_S[i]);
      obs_logR[i] ~ normal(pred_logR[i], sigma_obs);
    } else {
      pred_logR[i] = logRS[i] + log(obs_S[i]) + eps[i-1] * rho;
      obs_logR[i] ~ normal(pred_logR[i], SigAR);
    }
    eps[i] = obs_logR[i] - pred_logR[i];
  }
  // priors:
  sigma_obs ~ normal(0, 10);
  alpha ~ normal(0, 50);
  beta ~ normal(0, 50);
  rho ~ normal(0, 3);
}
