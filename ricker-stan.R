library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

SR <- read.csv("Harrison_simples_Apr18.csv")
SR$S_adj <- SR$S_adj/1000
SR$R <- SR$R/1000

srm <- lm(log(SR$R / SR$S_adj) ~ SR$S_adj)
(a_srm <- srm$coefficients[1])
(b_srm <- -srm$coefficients[2])
summary(srm)

dat <- list(obs_logR = log(SR$R), obs_S = SR$S_adj, N = length(SR$S_adj))

rick <- stan_model("ricker.stan")
rick1 <- stan_model("ricker1.stan")

m0_map <- optimizing(rick,
  data = dat,
  init = list(alpha = a_srm, beta = b_srm, sigma_obs = 0.8)
)
m0_map$par["alpha"]
m0_map$par["beta"]
m0_map$par["sigma_obs"]

m0 <- sampling(rick, data = dat, control = list(adapt_delta = 0.99), iter = 20000)
print(m0, pars = c("alpha", "beta", "sigma_obs"))

m0.1 <- sampling(rick1, data = dat, control = list(adapt_delta = 0.99), iter = 20000)
print(m0.1, pars = c("alpha", "beta", "sigma_obs"))

m1 <- stan("ricker_auto.stan", data = dat, control = list(adapt_delta = 0.99))
print(m1, pars = c("alpha", "beta", "sigma_obs", "rho"))

m2 <- stan("ricker_auto1.stan", data = dat, control = list(adapt_delta = 0.99))
print(m2, pars = c("alpha", "beta", "sigma_obs", "rho"))

m3 <- stan("ricker_auto2.stan", data = dat, control = list(adapt_delta = 0.99))
print(m3, pars = c("alpha", "beta", "sigma_obs", "rho"))

m4 <- stan("ricker_auto3.stan", data = dat, control = list(adapt_delta = 0.99))
print(m4, pars = c("alpha", "beta", "sigma_obs", "rho", "Smsy", "umsy"))

ricker_auto3 <- stan_model("ricker_auto3.stan")

# MAP estimate agrees with the CH TMB version!
m4_map <- optimizing(ricker_auto3,
  data = dat,
  init = list(alpha = a_srm, beta = b_srm, sigma_obs = 0.8, rho = 0.3)
)
m4_map$par["alpha"]
m4_map$par["beta"]
m4_map$par["sigma_obs"]
m4_map$par["rho"]

bayesplot::mcmc_combo(m4, pars = c("alpha", "beta", "sigma_obs", "rho"))
bayesplot::mcmc_combo(m2, pars = c("alpha", "beta", "sigma_obs", "rho"))

p <- extract(m4)
S <- seq(min(SR$S_adj), max(SR$S_adj), length.out = 100)
R_hat <- matrix(ncol = length(S), nrow = 400)
for (i in seq_len(nrow(R_hat))) {
  for (s in seq_along(S)) {
    R_hat[i, s] <- exp(p$alpha[i] - p$beta[i] * S[s] + log(S[s]))
  }
}
matplot(S, t(R_hat), type = "l", col = "#00000010", lty = 1)
points(SR$S_adj, SR$R)

# PPD:
SR <- mutate(SR, obs = as.factor(seq_len(n())))

library(dplyr)
library(ggplot2)
tidybayes::gather_draws(m4, ppd_obs_logR[obs]) %>%
  mutate(R_ppd = exp(.value)) %>%
  ggplot(aes(as.factor(obs), log(R_ppd))) + geom_violin() +
  geom_point(data = SR, mapping = aes(x = obs, y = log(R))) +
  theme_light()

tidybayes::gather_draws(m2, ppd_obs_logR[obs]) %>%
  mutate(R_ppd = exp(.value)) %>%
  ggplot(aes(as.factor(obs), log(R_ppd))) + geom_violin() +
  geom_point(data = SR, mapping = aes(x = obs, y = log(R))) +
  theme_light()

set.seed(123)
DRAWS <- sample(1:500, 24)
p$rho[DRAWS]

make_ppd_plot <- function(model) {
  tidybayes::gather_draws(model, ppd_obs_logR[obs]) %>%
    mutate(R_ppd = exp(.value), obs = as.factor(obs)) %>%
    left_join(SR) %>%
    filter(.iteration %in% DRAWS, .chain == 1) %>%
    mutate(type = "PPD") %>%
    bind_rows(
      tibble(obs = SR$obs, .iteration = 99999,
        R_ppd = SR$R, S_adj = SR$S_adj, type = "Real")) %>%
    ggplot(aes(S_adj, R_ppd)) + geom_point(aes(shape = type)) + geom_line() +
    facet_wrap(~.iteration, scales = "free_y") +
    theme_light() +
    scale_shape_manual(values = c("PPD" = 21, "Real" = 19))
}
make_ppd_plot(m2)
make_ppd_plot(m4)

tidybayes::spread_draws(m4, rho, sigma_obs) %>%
  filter(.iteration %in% DRAWS, .chain == 1) %>%
  as.data.frame()
