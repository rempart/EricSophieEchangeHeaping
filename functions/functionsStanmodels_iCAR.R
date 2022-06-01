### compilation of stan models
lapply(c("rstan", "loo"), library, character.only = TRUE)

# For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
# To avoid recompilation of unchanged Stan programs, we recommend calling
rstan::rstan_options(auto_write = TRUE)

stan_negbin_cov_spatial <- stan_model(file = paste(getwd(), "code stan", "heaping_negbin_cov_iCAR.stan", sep = "/"))
stan_nullnegbin_cov_spatial <- stan_model(file = paste(getwd(), "code stan", "negbin_cov_iCAR.stan", sep = "/"))
