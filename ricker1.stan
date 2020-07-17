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
  for (i in 1:N) {
    L_hat[i] = alpha - beta * obs_S[i]; // Eqn. 3 but not log(alpha)
  }
}
model {
  obs_logR ~ normal(log(obs_S) + L_hat, sigma_obs);
  // priors:
  sigma_obs ~ normal(0, 3);
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
}
