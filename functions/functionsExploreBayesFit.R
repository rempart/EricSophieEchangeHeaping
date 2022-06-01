check_conv <- function(fit, param) {
  n_chains <- fit@sim$chains
  n_iter <- (fit@sim$iter -fit@sim$warmup) / fit@sim$thin
  x <- rstan::extract(fit, param)[[1]]
  if(length(dim(x)) == 1) {
    y <- array(NA, dim = c(n_iter, n_chains))
    for(j in 1:n_chains) {
      y[, j] <- x[(j - 1) * n_iter + 1:n_iter]
    }
    out <- rstan::Rhat(sims = y)
  } 
  if(length(dim(x)) == 2) {
    y <- array(NA, dim = c(n_iter, n_chains, dim(x)[-1]))
    for(j in 1:n_chains) {
      y[, j, ] <- x[(j - 1) * n_iter + 1:n_iter, ]
    }
    out <- apply(y, 3, rstan::Rhat)
  }
  if(length(dim(x)) == 3) {
    y <- array(NA, dim = c(n_iter, n_chains, dim(x)[-1]))
    for(j in 1:n_chains) {
      y[, j, , ] <- x[(j - 1) * n_iter + 1:n_iter, , ]
    }
    dd <- apply(y, c(3, 4), rstan::Rhat)
    out <- NULL
    for(k in 1:nrow(dd)) {
      out <- c(out, dd[k, ])
    }
  }
  return(out)
}

get_pointestimate <- function(fit) {
  par_name <- c("intercept", "slope")
  output <- data.frame(mu = mean(exp(rstan::extract(fit, 'intercept')$intercept)),
                       beta_1 = apply(rstan::extract(fit, 'slope')$slope, 2, mean)[1],
                       beta_2 = apply(rstan::extract(fit, 'slope')$slope, 2, mean)[2]
  )
  if(any(fit@model_pars == "tau_0")) {
    par_name <- c(par_name, "tau_0", "lambda_0", "gamma")
    if(length(fit@par_dims$gamma) != 1) {
      output <- cbind(output,
                      data.frame(obs1_tau_0 = apply(rstan::extract(fit, 'tau_0')$tau_0, 2, mean)[1],
                                 obs2_tau_0 = apply(rstan::extract(fit, 'tau_0')$tau_0, 2, mean)[2],
                                 obs3_tau_0 = apply(rstan::extract(fit, 'tau_0')$tau_0, 2, mean)[3],
                                 obs1_lambda_0 = apply(rstan::extract(fit, 'lambda_0')$lambda_0, 2, mean)[1],
                                 obs2_lambda_0 = apply(rstan::extract(fit, 'lambda_0')$lambda_0, 2, mean)[2],
                                 obs3_lambda_0 = apply(rstan::extract(fit, 'lambda_0')$lambda_0, 2, mean)[3],
                                 obs1_gamma_0 = apply(rstan::extract(fit, 'gamma')$gamma, c(2, 3), mean)[1, 1],
                                 obs1_gamma_1 = apply(rstan::extract(fit, 'gamma')$gamma, c(2, 3), mean)[1, 2],
                                 obs1_gamma_2 = apply(rstan::extract(fit, 'gamma')$gamma, c(2, 3), mean)[1, 3],
                                 obs2_gamma_0 = apply(rstan::extract(fit, 'gamma')$gamma, c(2, 3), mean)[2, 1],
                                 obs2_gamma_1 = apply(rstan::extract(fit, 'gamma')$gamma, c(2, 3), mean)[2, 2],
                                 obs2_gamma_2 = apply(rstan::extract(fit, 'gamma')$gamma, c(2, 3), mean)[2, 3],
                                 obs3_gamma_0 = apply(rstan::extract(fit, 'gamma')$gamma, c(2, 3), mean)[3, 1],
                                 obs3_gamma_1 = apply(rstan::extract(fit, 'gamma')$gamma, c(2, 3), mean)[3, 2],
                                 obs3_gamma_2 = apply(rstan::extract(fit, 'gamma')$gamma, c(2, 3), mean)[3, 3]
                      )
      )
    } else {
      output <- cbind(output,
                      data.frame(obs1_tau_0 = mean(rstan::extract(fit, 'tau_0')$tau_0),
                                 obs1_lambda_0 = mean(rstan::extract(fit, 'lambda_0')$lambda_0),
                                 obs1_gamma_0 = apply(rstan::extract(fit, 'gamma')$gamma, 2, mean)[1],
                                 obs1_gamma_1 = apply(rstan::extract(fit, 'gamma')$gamma, 2, mean)[2],
                                 obs1_gamma_2 = apply(rstan::extract(fit, 'gamma')$gamma, 2, mean)[3]
                      )
      )
    }
  }
  if(any(fit@model_pars == "overdispersion")) {
    par_name <- c(par_name, "overdispersion")
    output$overdispersion <- mean(rstan::extract(fit, 'overdispersion')$overdispersion)
  } 
  ### check convergence
  rhat <- do.call('c', lapply(par_name, check_conv, fit = fit))
  output$looic <- loo::loo(loo::extract_log_lik(fit))$estimates[3, 1]
  if(any(rhat > 1.05)) {
    output[1, ] <- rep(NA, ncol(output))
  }
  return(output)
}

get_results <- function(l, j, m) {
  load(paste("D:/Heaping/stan_data_", l, "_n", j, "_", m, ".RData", sep = ""))
  temp <- do.call('rbind', 
                  lapply(get(paste("poisson", m, sep = "_")), 
                         get_pointestimate
                  )
  )
  par_name <- names(temp)
  temp %>% 
    mutate(sim = 1:n(),
           sample_size = j
    ) %>%
    pivot_longer(cols = all_of(par_name),
                 names_to = "param",
                 values_to = "estimate"
    ) %>% 
    mutate(model = paste("poisson", m, sep = "_"),
           truth = l
    ) %>% 
    as.data.frame() %>%
    write.table(file = paste("res/results_", l, "_", j, "_", m, "_poisson.txt", sep = ""),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t"
    )
  ### negbin
  temp <- do.call('rbind', 
                  lapply(get(paste("negbin", m, sep = "_")), 
                         get_pointestimate
                  )
  )
  par_name <- names(temp)
  temp %>% 
    mutate(sim = 1:n(),
           sample_size = j
    ) %>%
    pivot_longer(cols = all_of(par_name),
                 names_to = "param",
                 values_to = "estimate"
    ) %>% 
    mutate(model = paste("negbin", m, sep = "_"),
           truth = l
    ) %>% 
    as.data.frame() %>%
    write.table(file = paste("res/results_", l, "_", j, "_", m, "_negbin.txt", sep = ""),
                quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t"
    )
}