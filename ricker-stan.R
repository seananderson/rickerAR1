library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(dplyr)
library(ggplot2)

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

ctrl <- list(adapt_delta = 0.99)
pars <- c("alpha", "beta", "sigma_obs", "rho")

m0 <- sampling(rick, data = dat, control = ctrl, iter = 4000)
print(m0, pars = c("alpha", "beta", "sigma_obs"))

m0.1 <- sampling(rick1, data = dat, control = ctrl, iter = 4000)
print(m0.1, pars = c("alpha", "beta", "sigma_obs"))

m1 <- stan("ricker_auto.stan", data = dat, control = ctrl)
print(m1, pars = pars)

m2 <- stan("ricker_auto1.stan", data = dat, control = ctrl)
print(m2, pars = pars)

m3 <- stan("ricker_auto2.stan", data = dat, control = ctrl)
print(m3, pars = pars)

m4 <- stan("ricker_auto3.stan", data = dat, control = ctrl, iter = 4000)
print(m4, pars = pars)

# MAP estimate agrees with the CH TMB version!
ricker_auto3 <- stan_model("ricker_auto3.stan")
m4_map <- optimizing(ricker_auto3,
  data = dat,
  init = list(alpha = a_srm, beta = b_srm, sigma_obs = 0.8, rho = 0.3)
)
m4_map$par["alpha"]
m4_map$par["beta"]
m4_map$par["sigma_obs"]
m4_map$par["rho"]

m5 <- stan("ricker_auto4.stan", data = dat, control = ctrl, iter = 4000)
print(m5, pars = pars)

m6 <- stan("ricker_auto5.stan", data = c(dat, use_ar1 = 1), control = ctrl, iter = 4000)
print(m6, pars = pars)

m7 <- stan("ricker_auto5.stan", data = c(dat, use_ar1 = 0), control = ctrl, iter = 4000)
print(m7, pars = pars)

models <- list(thorson_literal = m1, holt_wor = m4, sr_ar1 = m6, sr_iid = m7)

post <- purrr::map_dfr(models,
  ~tidybayes::gather_draws(.x, alpha, beta, sigma_obs, rho), .id = "model") %>%
  filter(!(.variable == "rho" & model == "sr_iid"))

ggplot(post, aes(.value, fill = model, colour = model)) + geom_density(alpha = 0.25) +
  facet_wrap(~.variable, scales = "free") +
  theme_light() +
  scale_fill_brewer(palette = "Dark2") +
  scale_colour_brewer(palette = "Dark2")

p <- extract(m5)
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

tidybayes::gather_draws(m5, ppd_obs_logR[obs]) %>%
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
make_ppd_plot(m5)

tidybayes::spread_draws(m4, rho, sigma_obs) %>%
  filter(.iteration %in% DRAWS, .chain == 1) %>%
  as.data.frame()
