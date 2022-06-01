### some functions
rescale <- function(x) { return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)) }
rescale2 <- function(ynew, y) { return((ynew - mean(y, na.rm = TRUE)) /(sd(y, na.rm = TRUE))) }
unscale <- function(std_y, y) { mean(y, na.rm = TRUE) + sd(y, na.rm = TRUE) * std_y}

get_ci <- function(x, alpha = 0.80){
  c(coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[1], mean(x), coda::HPDinterval(coda::as.mcmc(x), prob = alpha)[2])
}

lower <- function(x) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = 0.80)[1]) }
upper <- function(x) { as.numeric(coda::HPDinterval(coda::as.mcmc(x), prob = 0.80)[2]) }

### prepare splines: B-splines Eilers & Marx (2010)
design_matrix <- function(x, xl, xr, ndx, deg = 3, d = 2) {
  # x = covariate
  # xl = lower bound
  # xr = upper bound
  # ndx = number of knots
  # deg = order of splines
  # d = order difference (1 or 2) for smoothness
  ### from Eilers & Marx
  bbase <- function(x, xl, xr, ndx, deg){
    tpower <- function(x, t, p){
      # Truncated p-th power function
      (x - t) ^ p * (x > t)
    }
    # Construct a B-spline basis of degree 'deg'
    dx <- (xr - xl) / ndx
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    P <- outer(x, knots, tpower, deg)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx ^ deg)
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    return(B)
  }
  B <- bbase(x, xl, xr, ndx, deg)
  # matrix of d-order differences
  D <- diff(diag(dim(B)[2]), diff = d)
  Q <- solve(D %*% t(D), D)
  X <- outer(x, 1:(d-1), "^") # fixed part without the intercept
  Z <- B %*% t(Q) # random part for mixed model formulation
  return(list(B = B, X = X, Z = Z))
}

# function to pass data to stan: no covariates
data4stan <- function(countdata = NULL, 
                      covar = NULL,
                      heaping_levels = NULL,
                      prior_location_intercept = NULL,
                      prior_scale_intercept = NULL,
                      prior_location_slope = 0.0,
                      prior_scale_slope = log(2) / 2,
                      prior_location = NULL,
                      prior_cholmat = NULL,
                      threshold = NULL,
                      splines = FALSE,
                      n_knot = 10,
                      prior_scale_smooth = 0.01
                      ) {
  if(is.null(covar)) {
    writeLines("\t Preparing data for analysis without covariates")
  } else {
    writeLines("\t Preparing data for analysis with covariates")
    if(length(countdata) != nrow(covar)) { stop("\t\t countdata and covar have different number of rows: please check") }
    covar_df <- as.data.frame(covar)
    names(covar_df) <- c("intercept", paste("cov", 1:(ncol(covar) - 1), sep = "_"))
    covar_df$Y <- countdata
    covar_df$ind_Y <- 1:length(countdata)
  }
  if(any(is.null(c(countdata, heaping_levels, prior_location_intercept, prior_scale_intercept, prior_location, prior_cholmat, threshold)))) {
    stop("Must provide values for all arguments")
  }
  # if(length(prior_location) != (length(heaping_levels) + 1)) {
  #   stop("Dimension mismatch between heaping parameters and heaping levels")
  # }
  
  ### table des correspondances
  prepared_data <- prepareData(Y = countdata, 
                               covar = covar, 
                               heaping_levels = heaping_levels
                               )
  
  ### reordonner la table --> normalement match avec covar
  my_rearranged_table <- prepared_data$my_table %>% 
    arrange(ind_Y) %>% 
    as.data.frame()
  
  # nb d'ant?c?dents pour chaque Y
  n_gx <- my_rearranged_table %>% 
    group_by(ind_Y) %>% 
    summarize(n = length(X)) %>% 
    pull() %>% 
    as.numeric()
  
  # ragged arrays
  gx_y <-array(0, dim = c(prepared_data$n, max(n_gx), 3))
  for(i in 1:dim(gx_y)[1]) {
    gx_y[i, 1:n_gx[i], 1] <- my_rearranged_table %>% 
      filter(ind_Y == i) %>% 
      select(X) %>%  
      pull() %>% 
      as.numeric()
    gx_y[i, 1:n_gx[i], 2] <- my_rearranged_table %>% 
      filter(ind_Y == i) %>% 
      select(G_behavior) %>%  
      pull() %>% 
      as.numeric()
    gx_y[i, 1:n_gx[i], 3] <- sapply(gx_y[i, 1:n_gx[i], 1], function(x) { which(prepared_data$unique_X == x) })
  }

  standata <- list(n_obs = prepared_data$n,
                   n_couples = n_gx,
                   n_max = max(n_gx),
                   n_y = as.numeric(prepared_data$repl),
                   Q = length(heaping_levels),
                   n_true = length(prepared_data$unique_X),
                   X = prepared_data$unique_X,
                   n_gx = n_gx,
                   GX_Y = gx_y,
                   prior_location_intercept = prior_location_intercept,
                   prior_scale_intercept = prior_scale_intercept,
                   prior_location = prior_location,
                   prior_cholmat = prior_cholmat,
                   threshold = threshold
                   )
  
  if(!is.null(covar)) {
    ### sanity checks: ensure that ordering is correct
    covar_df <- covar_df %>% 
      arrange(ind_Y) %>% 
      select(starts_with("cov_")) %>% 
      as.data.frame()
    writeLines(paste("\t Nb of covariates: ", ncol(covar_df), sep = ""))
    standata$n_cov <- ncol(covar_df)
    ## standardize
    covar_df <- apply(covar_df, 2, rescale)
    standata$Z <- as.matrix(covar_df)
    standata$prior_location_slope <- rep(prior_location_slope, standata$n_cov)
    standata$prior_scale_slope <- rep(prior_scale_slope, standata$n_cov) 
    if(splines) {
      writeLines("\t\t Beziers-splines model")
      standata$n_knot <- n_knot
      standata$B <- array(NA, dim = c(standata$n_cov, nrow(covar_df), n_knot + 1))
      for(j in 1:standata$n_cov) {
        standata$B[j, , ] <- design_matrix(x = covar_df[, j],
                                           xl = floor(10 * min(covar_df[, j])) / 10,
                                           xr = ceiling(10 * max(covar_df[, j])) / 10,
                                           ndx = n_knot
                                           )$Z
      }
      standata$prior_scale_smooth <- prior_scale_smooth 
    }
  }
  return(standata)
}

### posterior predictive checks
simple_ppc <- function(standata, stanfit, heaping_levels) {
  n_obs <- sum(standata$n_y)
  # extract posterior
  tau_0 <- as.numeric(rstan::extract(stanfit, 'tau_0')$tau_0)
  lambda_0 <- as.numeric(rstan::extract(stanfit, 'lambda_0')$lambda_0)
  gamma <- as.matrix(rstan::extract(stanfit, 'gamma')$gamma)
  intercept <- as.matrix(rstan::extract(stanfit, 'intercept')$intercept)
  
  linpred <- matrix(rep(intercept, n_obs), byrow = FALSE, nrow = length(intercept), ncol = n_obs)
  
  if(any(stanfit@model_pars == "overdispersion")) {
    omega <- as.numeric(rstan::extract(stanfit, 'overdispersion')$overdispersion)
    distrib <- 'nbinom'
  } else {
    omega <- rep(1.01, length(tau_0))
    distrib <- 'pois'
  }
  
  ### handle covariates
  if(any(stanfit@model_pars == "slope")) {
    linpred <- linpred + as.matrix(rstan::extract(stanfit, 'slope')$slope) %*% t(standata$Z)
    ### handle splines
    ##  --> will need a function to generate data simply from the linear predictor
    if(any(stanfit@model_pars == "smooth")) {
      smooth <- as.array(rstan::extract(stanfit, 'smooth')$smooth)
      for(j in 1:standata$n_cov) {
        linpred <- linpred + smooth[, j, ] %*% t(standata$B[j, , ])
      }
    }
  }

  ### handle spatial effect
  if(any(stanfit@model_pars == "spatial") && !is.null(standata$CELL)) {
    linpred <- linpred + as.matrix(rstan::extract(stanfit, 'spatial')$spatial)[, standata$CELL]
  }
  
  ### handle year random effect
  if(any(stanfit@model_pars == "alpha") && !is.null(standata$YEAR)) {
    linpred <- linpred + as.matrix(rstan::extract(stanfit, 'alpha')$alpha)[, standata$YEAR]
  }
  
  make_ppc <- function(iter) {
    out <- simulHeapedData(n = n_obs,
                           param = list(param_distrib = list(transfo_over_dispersion = log(omega[iter] - 1),
                                                             log_mu = log(expm1(linpred[iter, ])) # shifted distribution
                                                             ),
                                        param_heaping = list(lambda_0 = lambda_0[iter],
                                                             tau_0 = tau_0[iter],
                                                             threshold = 0.05,
                                                             gamma = gamma[iter, ]
                                                             )
                                        ),
                           heaping_levels = heaping_levels,
                           model = list(distrib = distrib,
                                        heaping = 'twoStepsRcpp'
                                        ),
                           quiet = TRUE
                           )
    return(out$Y)
  }
  
  y_rep <- t(sapply(1:length(intercept), make_ppc))
  return(y_rep)
}

## regarder l'histogramme des donnees et le comparer a l'histogramme attendu sous le modele
# rootograms
rootogram <- function(countdata, y_rep, min_obs = 1) {
  f_histogram <- function(x, max_obs) { table(factor(x, levels = min_obs:max_obs)) }
  max_obs <- max(countdata, as.numeric(y_rep))
  
  theme_set(theme_bw(base_size = 16))
  g <- ggplot(data.frame(mids = (min_obs:max_obs + (min_obs + 1):(max_obs + 1))/2,
                         y_obs = as.numeric(f_histogram(x = countdata, max_obs = max_obs)),
                         y_rep = apply(apply(y_rep, 1, f_histogram, max_obs = max_obs), 1, mean)
                         ),
              aes(x = mids, y = y_rep)
              ) +
    geom_rect(aes(xmin = mids - 0.5, xmax = mids + 0.5, 
                  ymax = y_obs, ymin = 0),
              fill = "lightgrey", color = grey(0.6), alpha = 0.5
              ) + 
    geom_line(aes(x = mids, y = y_rep), color = "black", size = 1.0) + 
    scale_y_sqrt(name = "Count") +
    scale_x_sqrt(name = "Nb of dolphins", breaks = c(c(1, 5, 10, 50), seq(100, 100 * ceiling(max_obs / 100), 100))) + 
    guides(size = "none") +
    theme(plot.title = element_text(lineheight = 0.8, face = "bold"), 
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          strip.background = element_rect(fill = grey(0.95))
          )
  return(g)
}

### plot linear relationships for glm
plot_linear <- function(data_df, stanfit, cov_name) {
  splines <- FALSE # by default, not a spline model
  X <- data_df[, cov_name]
  x_scale <- apply(X, 2, function(x) {c(mean(x), sd(x))})
  x_range <- apply(apply(X, 2, scale), 2, range)
  stdX <- cbind(rep(1, 1e4),
                apply(x_range, 2, function(vec) {seq(vec[1], vec[2], length.out = 1e4)})
                )
  beta <- cbind(rstan::extract(stanfit, 'intercept')$intercept, 
                rstan::extract(stanfit, 'slope')$slope
                )
  if(ncol(stdX) != ncol(beta)) { stop("Please check model and covariable names:\n\tdimension mismatch between slope parameters and covariables") }
  else {
    if(any(stanfit@model_pars == "smooth")) {
      smooth <- as.array(rstan::extract(stanfit, 'smooth')$smooth)
      splines <- TRUE
      writeLines("\t Splines model")
    } else {
      writeLines("\t Linear model")
    }
    dd <- data.frame(x = numeric(0),
                     y = numeric(0),
                     lower = numeric(0),
                     upper = numeric(0),
                     what = character(0),
                     covariable = character(0)
                     )
    for(j in 1:length(cov_name)) {
      Xmat <- cbind(rep(1, 1e4),
                    matrix(0, nrow = 1e4, ncol = length(cov_name))
                    )
      Xmat[, j+1] <- stdX[, j + 1]
      linpred <- beta %*% t(Xmat)
      if(splines) {
        Z <- design_matrix(x = stdX[, j + 1],
                           xl = floor(10 * min(stdX[, j + 1])) / 10,
                           xr = ceiling(10 * max(stdX[, j + 1])) / 10,
                           ndx = dim(smooth)[3] - 1
                           )$Z
        linpred <- linpred + smooth[, j, ] %*% t(Z)
      }
      dd <- rbind(dd,
                  data.frame(x = Xmat[, j+1] * x_scale[2, cov_name[j]] + x_scale[1, cov_name[j]],
                             y = apply(linpred, 2, mean),
                             lower = apply(linpred, 2, stats::quantile, prob = 0.025),
                             upper = apply(linpred, 2, stats::quantile, prob = 0.975),
                             what = rep("linear predictor", 1e4),
                             covariable = rep(cov_name[j], 1e4)
                             ),
                  data.frame(x = Xmat[, j+1] * x_scale[2, cov_name[j]] + x_scale[1, cov_name[j]],
                             y = apply(exp(linpred), 2, mean),
                             lower = apply(exp(linpred), 2, stats::quantile, prob = 0.025),
                             upper = apply(exp(linpred), 2, stats::quantile, prob = 0.975),
                             what = rep("response", 1e4),
                             covariable = rep(cov_name[j], 1e4)
                             )
                  )
    }
    dd$what <- factor(as.character(dd$what), levels = c("response", "linear predictor"))
    dd$covariable <- factor(as.character(dd$covariable), levels = cov_name)
    
    the_plot <- ggplot(data = dd,
                       aes(x = x, y = y)
                       ) +
      geom_ribbon(aes(x = x, ymin = lower, ymax = upper),
                  fill = "lightgrey", color = grey(0.6), alpha = 0.5
                  ) + 
      geom_line(aes(x = x, y = y), color = "black", size = 1.0) + 
      scale_y_continuous(name = "y") +
      scale_x_continuous(name = "Covariable") +
      facet_grid(what~covariable, scales = "free") +
      theme(plot.title = element_text(lineheight = 0.8, face = "bold"), 
            axis.text = element_text(size = 12),
            strip.text = element_text(size = 12),
            strip.background = element_rect(fill = grey(0.95))
            )
    
    return(the_plot)
  }
}
