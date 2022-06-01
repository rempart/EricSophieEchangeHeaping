### compilation of stan models
lapply(c("rstan", "loo", "bayesplot"), library, character.only = TRUE)

# For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
# To avoid recompilation of unchanged Stan programs, we recommend calling
rstan::rstan_options(auto_write = TRUE)

### poisson likelihood; no heaping
stan_nullpoisson <- stan_model(file = paste(getwd(), "code stan", "poisson.stan", sep = "/"))
### negative binomial likelihood; no heaping
stan_nullnegbin <- stan_model(file = paste(getwd(), "code stan", "negbin.stan", sep = "/"))
### poisson likelihood; heaping
stan_poisson <- stan_model(file = paste(getwd(), "code stan", "heaping_poisson.stan", sep = "/"))
### negative binomial likelihood; heaping
stan_negbin <- stan_model(file = paste(getwd(), "code stan", "heaping_negbin.stan", sep = "/"))
### negative binomial likelihood; heaping + covariates (linear)
stan_negbin_cov <- stan_model(file = paste(getwd(), "code stan", "heaping_negbin_cov.stan", sep = "/"))
### negative binomial likelihood; heaping + covariates (splines)
stan_negbin_Bsplines <- stan_model(file = paste(getwd(), "code stan", "heaping_negbin_Bsplines.stan", sep = "/"))
### negative binomial likelihood; heaping + covariates (linear) + year random effects
stan_negbin_cov_rand <- stan_model(file = paste(getwd(), "code stan", "heaping_negbin_cov_rand.stan", sep = "/"))
### negative binomial likelihood; heaping + covariates (linear) + year random effects + iCAR spatial effect
stan_negbin_cov_rand_spatial <- stan_model(file = paste(getwd(), "code stan", "heaping_negbin_cov_rand_iCAR.stan", sep = "/"))
