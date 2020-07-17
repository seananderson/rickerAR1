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
      obs_logR[i] ~ normal(pred_logR[i], sqrt(1 - rho^2) * sigma_obs);
    }
    eps[i] = obs_logR[i] - pred_logR[i];
  }
  // priors:
  sigma_obs ~ normal(0, 10);
  alpha ~ normal(0, 50);
  beta ~ normal(0, 50);
  rho ~ normal(0, 3);
}
generated quantities {
  real umsy;
  real Smsy;

  vector[N] logRS;
  vector[N] pred_logR;
  vector[N] eps;
  vector[N] ppd_obs_logR;

  for (i in 1:N) {
    logRS[i] = alpha - beta * obs_S[i];
    if (i == 1) {
      pred_logR[i] = logRS[i] + log(obs_S[i]);
      ppd_obs_logR[i] = normal_rng(pred_logR[i], sigma_obs);
    } else {
      pred_logR[i] = logRS[i] + log(obs_S[i]) + eps[i-1] * rho;
      ppd_obs_logR[i] = normal_rng(pred_logR[i], sqrt(1 - rho^2) * sigma_obs);
    }
    eps[i] = ppd_obs_logR[i] - pred_logR[i];
  }

  umsy = 0.5 * alpha - 0.07 * (alpha^2);
  Smsy = alpha/beta * (0.5 - 0.07 * alpha);
}
