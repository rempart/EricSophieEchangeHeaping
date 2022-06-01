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
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param psi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param rho Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of psi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector psi, real tau, real rho, 
  int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_n) {
  row_vector[n] psit_D; // psi' * D
  row_vector[n] psit_W; // psi' * W
  vector[n] ldet_terms;

  psit_D = (psi .* D_sparse)';
  psit_W = rep_row_vector(0, n);
  for (i in 1:W_n) {
    psit_W[W_sparse[i, 1]] = psit_W[W_sparse[i, 1]] + psi[W_sparse[i, 2]];
    psit_W[W_sparse[i, 2]] = psit_W[W_sparse[i, 2]] + psi[W_sparse[i, 1]];
  }

  for (i in 1:n) ldet_terms[i] = log1m(rho * lambda[i]);
    return 0.5 * (n * log(tau) + sum(ldet_terms) - tau * (psit_D * psi - rho * (psit_W * psi)));
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
  // iCAR
  int<lower = 1> n_cell;
  int<lower = 1, upper = n_cell> CELL[n_obs];
  real<lower = 0> prior_scale_sigma;
  matrix[n_cell, n_cell] A;          // adjacency matrix
  int<lower = 1> n_pairs;            // is equal to sum(A)/2
  // random year effect
  int<lower = 1> n_year;
  int<lower = 1, upper = n_year> YEAR[n_obs];
}

transformed data { // could be done outside Stan
 int W_sparse[n_pairs, 2];     // adjacency pairs
 vector[n_cell] D_sparse;      // diagonal of D (number of neigbors for each site)
 vector[n_cell] lambda;        // eigenvalues of invsqrtD * A * invsqrtD

 { // generate sparse representation for A
  int counter;
  counter = 1;
  // loop over upper triangular part of A to identify neighbour pairs
  for (i in 1:(n_cell - 1)) {
    for (j in (i + 1):n_cell) {
      if (A[i, j] == 1) {
        W_sparse[counter, 1] = i;
        W_sparse[counter, 2] = j;
        counter = counter + 1;
      }
    }
  }
 }
 for (i in 1:n_cell) D_sparse[i] = sum(A[i]);
 {
  vector[n_cell] invsqrtD;  
  for (i in 1:n_cell) {
    invsqrtD[i] = 1 / sqrt(D_sparse[i]);
  }
  lambda = eigenvalues_sym(quad_form(A, diag_matrix(invsqrtD)));
 }
}

parameters {
  real unscaled_intercept;
  vector[Q + 1] unscaled_epsilon;
  vector[n_cov] unscaled_slope;
  real<lower = 0.0, upper = 1.0> inv_overdispersion;
  real<lower = 0.0> unscaled_sigma2;
  real<lower = 0.0> tau;
  vector[n_cell] spatial;
  real logit_rho;
  vector[n_year] unscaled_alpha;
  real<lower = 0, upper = 1> prop;
}

transformed parameters {
  real overdispersion;
  real intercept;
  vector[n_cov] slope;
  vector[Q + 1] epsilon;
  real tau_0;
  real lambda_0;
  real sup_heaping;
  vector[n_obs] linpred;
  vector[n_obs] phi;
  vector[Q - 1] delta;
  vector[Q - 1] gamma;
  matrix[n_true, Q] log_probaG_givenX;
  vector[n_obs] log_lik;
  real rho;
  real sill;
  real tau_sill;
  real sigma_year;
  vector[n_year] alpha;
  // iCAR
  sill = prior_scale_sigma * sqrt(prop * unscaled_sigma2 / tau);
  tau_sill = 1 / square(sill);
  rho = exp(1.5 * logit_rho);
  // random year effect
  sigma_year = prior_scale_sigma * sqrt((1 - prop) * unscaled_sigma2 / tau);
  alpha = sigma_year * unscaled_alpha;
  // negbin
  overdispersion = 1 / inv_overdispersion;
  intercept = prior_location_intercept + prior_scale_intercept * unscaled_intercept;
  slope = prior_location_slope + prior_scale_slope .* unscaled_slope;
  // heaping
  epsilon = prior_location + prior_cholmat * unscaled_epsilon;
  tau_0 = exp(epsilon[1]);
  lambda_0 = exp(epsilon[2]);
  delta = epsilon[3:(Q + 1)];
  gamma = from_delta_to_gamma(delta);
  sup_heaping  = tau_0 - log(threshold) / lambda_0;
  log_probaG_givenX = log_prob_heaping(X, tau_0, lambda_0, gamma);
  linpred = rep_vector(intercept, n_obs) + alpha[YEAR] + Z * slope + spatial[CELL];
  phi = exp(linpred) / (overdispersion - 1);
  // loop over Y, implicit in duplet (G, X)
  for(i in 1:n_obs) {
    vector[n_couples[i]] acc;
    // loop over (G, X) duplets
    for(j in 1:n_couples[i]) {
      // log prob of true count given linpred under posited model
      acc[j] = neg_binomial_2_log_lpmf(GX_Y[i, j, 3]| linpred[i], phi[i]);// - neg_binomial_2_log_lpmf(0| linpred[i], phi[i]);
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
  unscaled_alpha ~ normal(0.0, 1.0); 
  // variance
  unscaled_sigma2 ~ gamma(0.5, 1.0);
  tau ~ gamma(1.0, 1.0); 
  // iCAR
  logit_rho ~ normal(0.0, 1.0);
  spatial ~ sparse_car(tau_sill, rho, W_sparse, D_sparse, lambda, n_cell, n_pairs);
  // increment likelihood
  for(i in 1:n_obs) {
    target += log_lik[i];
  }
}

