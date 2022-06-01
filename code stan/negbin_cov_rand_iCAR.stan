functions {
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
  int<lower = 1> n_obs;                 // nombre de donnees uniques
  int<lower = 1> n_y[n_obs];            // nombre de repetitions de chaque valeur unique de Y
  int<lower = 1> Y[n_obs];              // comptage reporté par les observateurs
  real prior_location_intercept;
  real<lower = 0.0> prior_scale_intercept;
  int<lower = 1> n_cov;                 // nombre de covariables
  vector[n_cov] prior_location_slope;
  vector<lower = 0.0>[n_cov] prior_scale_slope;
  matrix[n_obs, n_cov] Z;               // design matrix
  // iCAR
  int<lower = 1> n_cell;
  int<lower = 1, upper = n_cell> CELL[n_obs];
  real<lower = 0> prior_scale_sigma;
  matrix[n_cell, n_cell] A;          // adjacency matrix
  int<lower = 1> n_pairs;            // is equal to sum(A)/2
  // random year effect
  int<lower = 1> n_year;
  int<lower = 1, upper = n_cell> YEAR[n_obs];
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
  vector[n_obs] linpred;
  vector[n_obs] phi;
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
  linpred = rep_vector(intercept, n_obs) + alpha[YEAR] + Z * slope + spatial[CELL];
  phi = exp(linpred) / (overdispersion - 1);
  // loop over Y
  for(i in 1:n_obs) {
    // likelihood
    log_lik[i] = n_y[i] * neg_binomial_2_log_lpmf(Y[i]| linpred[i], phi[i]);
  }
}

model {
  unscaled_intercept ~ normal(0.0, 1.0);
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

generated quantities {
  vector[n_obs] y_rep;
  for(i in 1:n_obs) {
    y_rep[i] = neg_binomial_2_log_rng(linpred[i], phi[i]);
  }
}

