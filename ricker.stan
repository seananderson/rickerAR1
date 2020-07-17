data {
  int<lower=0> N;
  vector[N] obs_logR;
  vector[N] obs_S;
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> sigma_obs;
}
transformed parameters {
  vector[N] L_hat;
  vector[N] eps;
  for (i in 1:N) {
    L_hat[i] = alpha - beta * obs_S[i]; // Eqn. 3 but not log(alpha)
    eps[i] = obs_logR[i] - log(obs_S[i]) - L_hat[i];
  }
}
model {
  eps ~ normal(0, sigma_obs);
  // priors:
  sigma_obs ~ normal(0, 3);
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
}
