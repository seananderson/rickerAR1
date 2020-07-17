data {
  int<lower=0> N;
  vector[N] obs_logR;
  vector[N] obs_S;
}
transformed data {
  vector[N] obs_logRS;
  obs_logRS = obs_logR - log(obs_S);
}
parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  real<lower=0> sigma_obs;
  real<lower=-1, upper=1> rho;
}
model {
  vector[N] pred_logRS;
  vector[N] eps;
  for (i in 1:N) {
    pred_logRS[i] = alpha - beta * obs_S[i];
    if (i == 1) {
      obs_logRS[i] ~ normal(pred_logRS[i], sigma_obs);
    } else {
      pred_logRS[i] = pred_logRS[i] + eps[i-1] * rho;
      obs_logRS[i] ~ normal(pred_logRS[i], sqrt(1 - rho^2) * sigma_obs);
    }
    eps[i] = obs_logRS[i] - pred_logRS[i];
  }
  // priors:
  sigma_obs ~ normal(0, 3);
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  rho ~ normal(0, 2);
}
generated quantities {
  real umsy;
  real Smsy;

  vector[N] ppd_pred_logRS;
  vector[N] ppd_eps;
  vector[N] ppd_obs_logRS;
  vector[N] ppd_obs_logR;

  for (i in 1:N) {
    ppd_pred_logRS[i] = alpha - beta * obs_S[i];
    if (i == 1) {
      ppd_obs_logRS[i] = normal_rng(ppd_pred_logRS[i], sigma_obs);
      ppd_obs_logR[i] = ppd_obs_logRS[i] + log(obs_S[i]);
    } else {
      ppd_obs_logRS[i] = ppd_pred_logRS[i] + ppd_eps[i-1] * rho;
      ppd_obs_logRS[i] = normal_rng(ppd_obs_logRS[i], sqrt(1 - rho^2) * sigma_obs);
      ppd_obs_logR[i] = ppd_obs_logRS[i] + log(obs_S[i]);
    }
    ppd_eps[i] = ppd_obs_logRS[i] - ppd_pred_logRS[i];
  }

  umsy = 0.5 * alpha - 0.07 * (alpha^2);
  Smsy = alpha/beta * (0.5 - 0.07 * alpha);
}
