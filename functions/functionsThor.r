source('functions/functionsHeapingModel.R')
source('functions/functionsHeapingProbabilityRcpp.R')

### load functions for Stan
source('functions/functionsStan.R')
source('functions/functionsStanmodels.r')

### bayesian estimation (6 models to fit)
# get_bayes <- function(all_data,
#                       n,
#                       distrib,
#                       n_chain = 4,
#                       n_iter = 1e3,
#                       n_warm = 500,
#                       n_thin = 1,
#                       hmc_control = list(adapt_delta = 0.9)
#                       ) {
#   poisson_null <- negbin_null <- vector(mode = 'list', length = length(all_data))
#   poisson_heap <- negbin_heap <- vector(mode = 'list', length = length(all_data))
#   poisson_mheap <- negbin_mheap <- vector(mode = 'list', length = length(all_data))
#   ### fit the models
#   for(i in 1:length(all_data)) {
#     ### null model
#     poisson_null[[i]] <- sampling(stan_nullpoisson_cov,
#                                   data = list(n_obs = length(all_data[[i]]$Y),
#                                               n_y = as.numeric(all_data[[i]]$Y),
#                                               Y = as.numeric(all_data[[i]]$Y),
#                                               prior_location_intercept = log(10),
#                                               prior_scale_intercept = log(5) / 2,
#                                               n_cov = 2,
#                                               prior_location_slope = rep(0, 2),
#                                               prior_scale_slope = rep(log(2) / 2, 2),
#                                               Z = as.matrix(all_data[[i]][, grep("cov", names(all_data[[i]]))])
#                                               ),
#                                   pars = c("intercept", "slope", "log_lik"),
#                                   chains = n_chain,
#                                   iter = n_iter,
#                                   warmup = n_warm,
#                                   thin = n_thin,
#                                   control = hmc_control
#                                   )
#     negbin_null[[i]] <- sampling(stan_nullnegbin_cov,
#                                  data = list(n_obs = length(all_data[[i]]$Y),
#                                              n_y = as.numeric(all_data[[i]]$Y),
#                                              Y = as.numeric(all_data[[i]]$Y),
#                                              prior_location_intercept = log(10),
#                                              prior_scale_intercept = log(5) / 2,
#                                              n_cov = 2,
#                                              prior_location_slope = rep(0, 2),
#                                              prior_scale_slope = rep(log(2), 2),
#                                              Z = as.matrix(all_data[[i]][, grep("cov", names(all_data[[i]]))])
#                                              ),
#                                  pars = c("intercept", "slope", "overdispersion", "log_lik"),
#                                  chains = n_chain,
#                                  iter = n_iter,
#                                  warmup = n_warm,
#                                  thin = n_thin,
#                                  control = hmc_control
#                                  )
#     ### one observer
#     standatalist <- data4stan(countdata = all_data[[i]]$Y,
#                               covar = cbind(rep(1, nrow(all_data[[i]])), as.matrix(all_data[[i]][, grep("cov", names(all_data[[i]]))])),
#                               heaping_levels = heaping_levels,
#                               prior_location_intercept = log(10),
#                               prior_scale_intercept = log(5) / 2,
#                               prior_location = c(log(10), log(0.5), log(0.2), 0, log(10)),
#                               prior_cholmat = diag(c(log(3) ,log(5) /2, log(5) / 2, 5.0, log(5) / 2)) %*%
#                                 diag(5),
#                               threshold = 0.05
#                               )
# 
#     poisson_heap[[i]] <- sampling(stan_poisson_cov,
#                                   data = standatalist,
#                                   pars = c("intercept", "slope", "gamma", "tau_0", "lambda_0", "sup_heaping", "log_lik"),
#                                   chains = n_chain,
#                                   iter = n_iter,
#                                   warmup = n_warm,
#                                   thin = n_thin,
#                                   control = hmc_control
#                                   )
# 
#     negbin_heap[[i]] <- sampling(stan_negbin_cov,
#                                  data = standatalist,
#                                  pars = c("intercept", "slope", "gamma", "tau_0", "lambda_0", "sup_heaping", "overdispersion", "log_lik"),
#                                  chains = n_chain,
#                                  iter = n_iter,
#                                  warmup = n_warm,
#                                  thin = n_thin,
#                                  control = hmc_control
#                                  )
#     ### several observers
#     prior_cholmat <- array(NA, dim = c(length(unique(all_data[[i]]$Observer)), 5, 5))
#     for(k in 1:length(unique(all_data[[i]]$Observer))) {
#       prior_cholmat[k, , ] <- diag(c(log(3) ,log(5) /2, log(5) / 2, 5.0, log(5) / 2)) %*% diag(5)
#     }
#     standatalist <- data4stan(countdata = all_data[[i]]$Y,
#                               covar = cbind(rep(1, nrow(all_data[[i]])), as.matrix(all_data[[i]][, grep("cov", names(all_data[[i]]))])),
#                               heaping_levels = heaping_levels,
#                               prior_location_intercept = log(10),
#                               prior_scale_intercept = log(5) / 2,
#                               prior_location = matrix(rep(c(log(10), log(0.5), log(0.2), 0, log(10)), 
#                                                           length(unique(all_data[[i]]$Observer))
#                                                           ),
#                                                       nrow = 3, byrow = TRUE
#                                                       ),
#                               prior_cholmat = prior_cholmat,
#                               threshold = 0.05
#                               )
#     standatalist$n_observer <- length(unique(all_data[[i]]$Observer))
#     standatalist$OBS <- all_data[[i]]$Observer
#     
#     poisson_mheap[[i]] <- sampling(stan_poisson_cov_multi,
#                                    data = standatalist,
#                                    pars = c("intercept", "slope", "gamma", "tau_0", "lambda_0", "sup_heaping", "log_lik"),
#                                    chains = n_chain,
#                                    iter = n_iter,
#                                    warmup = n_warm,
#                                    thin = n_thin,
#                                    control = hmc_control
#                                    )
#     
#     negbin_mheap[[i]] <- sampling(stan_negbin_cov_multi,
#                                   data = standatalist,
#                                   pars = c("intercept", "slope", "gamma", "tau_0", "lambda_0", "sup_heaping", "overdispersion", "log_lik"),
#                                   chains = n_chain,
#                                   iter = n_iter,
#                                   warmup = n_warm,
#                                   thin = n_thin,
#                                   control = hmc_control
#                                   )
#     
#   }
#   save(file = paste("res/stan_data_", distrib, "_n", n, ".RData", sep = ""),
#        list = c("poisson_null", "negbin_null", 
#                 "poisson_heap", "negbin_heap",
#                 "poisson_mheap", "negbin_mheap"
#                 )
#        )
# }
