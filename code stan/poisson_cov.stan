data {
  int<lower = 1> n_obs;                 // nombre de donnees
  int<lower = 1> n_y[n_obs];            // nombre de repetitions de chaque valeur unique de Y
  int<lower = 1> Y[n_obs];              // comptage reporté par les observateurs
  real prior_location_intercept;
  real<lower = 0.0> prior_scale_intercept;
  int<lower = 1> n_cov;                 // nombre de covariables
  vector[n_cov] prior_location_slope;
  vector<lower = 0.0>[n_cov] prior_scale_slope;
  matrix[n_obs, n_cov] Z;               // design matrix
}

parameters {
  real unscaled_intercept;
  vector[n_cov] unscaled_slope;
}

transformed parameters {
  real intercept;
  vector[n_cov] slope;
  vector[n_obs] linpred;
  vector[n_obs] log_lik;
  intercept = prior_location_intercept + prior_scale_intercept * unscaled_intercept;
  slope = prior_location_slope + prior_scale_slope .* unscaled_slope;
  linpred = rep_vector(intercept, n_obs) + Z * slope;
  // loop over Y
  for(i in 1:n_obs) {
    // likelihood
    log_lik[i] = n_y[i] * poisson_log_lpmf(Y[i]| linpred[i]);
  }
}

model {
  unscaled_intercept ~ normal(0.0, 1.0);
  for(i in 1:n_obs) {
    target += log_lik[i];
  }
}

generated quantities {
  vector[n_obs] y_rep;
  for(i in 1:n_obs) {
    y_rep[i] = poisson_log_rng(linpred[i]);
  }
}
