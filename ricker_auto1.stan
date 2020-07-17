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
model {
  vector[N] L_hat;
  vector[N] L_obs;
  vector[N] eps;
  vector[N] delta;
  for (i in 1:N) {
    L_hat[i] = alpha - beta * obs_S[i];
    L_obs[i] = obs_logR[i] - log(obs_S[i]);
    eps[i] = L_obs[i] - L_hat[i];
    if (i == 1) {
      delta[i] = eps[i];
    } else {
      delta[i] = (eps[i] - rho * eps[i-1]) / sqrt(1 - rho^2); // Eqn.4
    }
  }
  delta ~ normal(0, sigma_obs);
  // priors:
  sigma_obs ~ normal(0, 3);
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  rho ~ normal(0, 2);
}
generated quantities {
  real umsy;
  real Smsy;

  vector[N] L_hat;
  vector[N] L_obs;
  vector[N] delta;
  vector[N] ppd_obs_logR;
  for (i in 1:N) {
    delta[i] = normal_rng(0, sigma_obs);
    L_hat[i] = alpha - beta * obs_S[i];
    if (i == 1) {
      L_obs[i] = L_hat[i] + delta[i];
    } else {
      L_obs[i] = L_hat[i] + rho * (L_obs[i-1] - L_hat[i-1]) + sqrt(1 - rho^2) * delta[i];
    }
    ppd_obs_logR[i] = L_obs[i] + log(obs_S[i]);
  }

  umsy = 0.5 * alpha - 0.07 * (alpha^2);
  Smsy = alpha/beta * (0.5 - 0.07 * alpha);
}

