functions {
  // compute the log probabilities of heaping: matrix output of dimension [length(X), length(gamma) + 1]
  matrix log_prob_heaping(vector X, real tau_0, real lambda_0, vector gamma) { 
    int n = rows(X);
    int Q = rows(gamma);
    matrix[n, Q + 1] res;
    for(i in 1:n) {
     vector[Q - 1] gamma_star = -gamma[2:Q] - gamma[1] * (X[i] - tau_0);
     res[i, 1] = (X[i] <= tau_0) ? 0.0 : -lambda_0 * (X[i] - tau_0);
     res[i, 2] = (X[i] <= tau_0) ? negative_infinity() : log1m_exp(res[i, 1]) + log_inv_logit(gamma_star[1]);
     res[i, Q + 1]  = (X[i] <= tau_0) ? negative_infinity() : log1m_exp(res[i, 1]) + log1m_inv_logit(gamma_star[Q - 1]);
     if(Q > 2) {
      for(k in 3:Q) {
       res[i, k] = (X[i] <= tau_0) ? negative_infinity() : log1m_exp(res[i, 1]) + log_diff_exp(-gamma_star[k - 2], -gamma_star[k - 1]) - log1p_exp(-gamma_star[k - 2]) - log1p_exp(-gamma_star[k - 1]);
      }
     }
    }
    return res;
  }
  // convert delta to gamma parameters
  vector from_delta_to_gamma (vector delta) {
    int Q = rows(delta);
    vector[Q] gamma = delta;
    gamma[1] = exp(delta[1]);
    if(Q > 2) {
      gamma[3:Q] = rep_vector(delta[2], Q - 2) - cumulative_sum(exp(delta[3:Q])); // ensures that gammas are ordered
    }
    return gamma;
  }
}

data {
  int<lower = 1> n_obs;                 // nombre de donnees Y unique
  int<lower = 1> n_couples[n_obs];      // nombre de possibilités en couple (G,X) pour chaque Y
  int<lower = 1> n_max;                 // max de n_gx
  int<lower = 1> n_y[n_obs];            // nombre de repetitions de chaque valeur unique de Y
  int<lower = 3> Q;                     // nombre de comportements d'arrondis
  int<lower = 1> n_true;                // nombre de X unique à considérer
  vector<lower = 1>[n_true] X;          // vrai comptage potentiel
  int<lower = 0> GX_Y[n_obs, n_max, 3]; // couples (G, X)
  real prior_location_intercept;
  real<lower = 0.0> prior_scale_intercept;
  vector[Q + 1] prior_location;
  matrix[Q + 1, Q + 1] prior_cholmat;
  real<lower = 0.0, upper = 1.0> threshold;
  int<lower = 1> n_cov;                  // nombre de covariables
  vector[n_cov] prior_location_slope;
  vector<lower = 0.0>[n_cov] prior_scale_slope;
  matrix[n_obs, n_cov] Z;                // design matrix
}

parameters {
  real unscaled_intercept;
  vector[Q + 1] unscaled_epsilon;
  vector[n_cov] unscaled_slope;
}

transformed parameters {
  real intercept;
  vector[n_cov] slope;
  vector[Q + 1] epsilon;
  real tau_0;
  real lambda_0;
  real sup_heaping;
  vector[n_obs] linpred;
  vector[Q - 1] delta;
  vector[Q - 1] gamma;
  matrix[n_true, Q] log_probaG_givenX;
  vector[n_obs] log_lik;
  intercept = prior_location_intercept + prior_scale_intercept * unscaled_intercept;
  slope = prior_location_slope + prior_scale_slope .* unscaled_slope;
  epsilon = prior_location + prior_cholmat * unscaled_epsilon;
  tau_0 = exp(epsilon[1]);
  lambda_0 = exp(epsilon[2]);
  delta = epsilon[3:(Q + 1)];
  gamma = from_delta_to_gamma(delta);
  sup_heaping  = tau_0 - log(threshold) / lambda_0;
  log_probaG_givenX = log_prob_heaping(X, tau_0, lambda_0, gamma);
  linpred = rep_vector(intercept, n_obs) + Z * slope;
  // loop over Y, implicit in duplet (G, X)
  for(i in 1:n_obs) {
    vector[n_couples[i]] acc;
    // loop over (G, X) duplets
    for(j in 1:n_couples[i]) {
      // log prob of true count given linpred under posited model
      acc[j] = poisson_log_lpmf(GX_Y[i, j, 1]| linpred[i]);// - poisson_log_lpmf(0| linpred[i]);
      // add log(prob of G given X)
      acc[j] += log_probaG_givenX[GX_Y[i, j, 3], GX_Y[i, j, 2]];
    }
    // likelihood
    log_lik[i] = n_y[i] * log_sum_exp(acc);
  }
}

model {
  unscaled_intercept ~ normal(0.0, 1.0);
  unscaled_epsilon ~ normal(0.0, 1.0);
  unscaled_slope ~ normal(0.0, 1.0);
  for(i in 1:n_obs) {
    target += log_lik[i];
  }
}

