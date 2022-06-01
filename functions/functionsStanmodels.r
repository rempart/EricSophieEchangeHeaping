### compilation of stan models
lapply(c("rstan", "loo"), library, character.only = TRUE)

# For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
# To avoid recompilation of unchanged Stan programs, we recommend calling
rstan::rstan_options(auto_write = TRUE)

# compile available models
stan_nullpoisson <- stan_model(file = paste(getwd(), "code stan", "poisson.stan", sep = "/"))
stan_nullnegbin <- stan_model(file = paste(getwd(), "code stan", "negbin.stan", sep = "/"))
stan_nullpoisson_cov <- stan_model(file = paste(getwd(), "code stan", "poisson_cov.stan", sep = "/"))
stan_nullnegbin_cov <- stan_model(file = paste(getwd(), "code stan", "negbin_cov.stan", sep = "/"))
stan_poisson <- stan_model(file = paste(getwd(), "code stan", "heaping_poisson.stan", sep = "/"))
stan_negbin <- stan_model(file = paste(getwd(), "code stan", "heaping_negbin.stan", sep = "/"))
stan_poisson_cov <- stan_model(file = paste(getwd(), "code stan", "heaping_poisson_cov.stan", sep = "/"))
stan_negbin_cov <- stan_model(file = paste(getwd(), "code stan", "heaping_negbin_cov.stan", sep = "/"))
stan_poisson_Bsplines <- stan_model(file = paste(getwd(), "code stan", "heaping_poisson_Bsplines.stan", sep = "/"))
stan_negbin_Bsplines <- stan_model(file = paste(getwd(), "code stan", "heaping_negbin_Bsplines.stan", sep = "/"))
stan_poisson_cov_multi <- stan_model(file = paste(getwd(), "code stan", "multiple_heaping_poisson_cov.stan", sep = "/"))
stan_negbin_cov_multi <- stan_model(file = paste(getwd(), "code stan", "multiple_heaping_negbin_cov.stan", sep = "/"))
