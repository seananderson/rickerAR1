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
  vector[N] L_hat;
  vector[N] L_obs;
  vector[N] delta;
  for (i in 1:N) {
    L_hat[i] = alpha - beta * obs_S[i]; // Eqn. 3 but not log(alpha)
    L_obs[i] = obs_logR[i] - log(obs_S[i]);
    if (i == 1) {
      delta[i] = L_obs[i] - L_hat[i];
    } else {
      // Eqn 5 rearranged:
      delta[i] = (L_obs[i] - L_hat[i] - rho * (L_obs[i-1] - L_hat[i-1])) /
                 sqrt(1 - rho^2);
    }
  }
}
model {
  delta ~ normal(0, sigma_obs);
  // priors:
  sigma_obs ~ normal(0, 10);
  alpha ~ normal(0, 50);
  beta ~ normal(0, 50);
  rho ~ normal(0, 3);
}
