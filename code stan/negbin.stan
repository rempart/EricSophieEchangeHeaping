data {
  int<lower = 1> n_obs;                 // nombre de donnees uniques
  int<lower = 1> n_y[n_obs];            // nombre de repetitions de chaque valeur unique de Y
  int<lower = 1> Y[n_obs];              // comptage reporté par les observateurs
  real prior_location_intercept;
  real<lower = 0.0> prior_scale_intercept;
}

transformed data {
  int n_tot;
  n_tot = sum(n_y);
}

parameters {
  real unscaled_intercept;
  real<lower = 0.0, upper = 1.0> inv_overdispersion;
}

transformed parameters {
  real overdispersion;
  real intercept;
  vector[n_obs] linpred;
  vector[n_obs] phi;
  vector[n_obs] log_lik;
  overdispersion = 1 / inv_overdispersion;
  intercept = prior_location_intercept + prior_scale_intercept * unscaled_intercept;
  linpred = rep_vector(intercept, n_obs);
  phi = exp(linpred) / (overdispersion - 1);
  // loop over Y
  for(i in 1:n_obs) {
    // likelihood
    log_lik[i] = n_y[i] * neg_binomial_2_log_lpmf(Y[i]| linpred[i], phi[i]);
  }
}

model {
  unscaled_intercept ~ normal(0.0, 1.0);
  for(i in 1:n_obs) {
    target += log_lik[i];
  }
}

generated quantities {
  vector[n_tot] y_rep;
  for(i in 1:n_tot) {
    y_rep[i] = neg_binomial_2_log_rng(intercept, exp(intercept) / (overdispersion - 1));
  }
}
