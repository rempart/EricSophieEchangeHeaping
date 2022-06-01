#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;  // Evite d'écrire arma::fonctionArmadillo

// [[Rcpp::export]]   

 
NumericMatrix probHeapingRcppRaw(
    IntegerVector X, 
    double tau0, 
    double lambda0, 
    NumericVector gamma
  ) 
  {
  int  n = X.length(); 
  int  Q = gamma.length();   
  int X_i  = 1;
  
  
  
  double Z_i;
  NumericMatrix res(n,Q + 1); 
  
  for (int i = 0; i < n; ++i ){
   
    X_i  = X[i]; 
    Z_i = gamma[0] * (X_i - tau0); 
    if (X_i <= tau0){ res(i,0)  = 1;} else{ res(i,0) = exp(-lambda0 * (X_i - tau0));} 
    
  
    res(i,1)  = (1-res(i,0))/(1 + exp(gamma[1] + Z_i ));
    res(i,Q)  = (1-res(i,0)) *(1 - 1/(1 + exp(gamma[Q-1] + Z_i)));
    for (int q = 2; q < Q; q++){
      res(i, q)  = (1-res(i,0))*(1/(1 + exp(gamma[q] + Z_i)) - 1/(1 + exp(gamma[q-1] + Z_i)));
    }
    
  }
   return res;
}

// [[Rcpp::export]]   
arma::cube probHeapingVectRcppRaw(
    IntegerVector X, 
    NumericVector tau0, 
    NumericVector lambda0, 
    NumericMatrix gamma
) 
  {
    
    int  n = X.length(); 
    int  Q = gamma.ncol();   
    int X_i  = 1;
    int M = tau0.length();
    
    
    double Z_i;
    arma::cube res(M,n,Q + 1); 
    
    for (int m = 0; m < M; ++m ){
        double lambda0_m  =  lambda0[m];
      
        for (int i = 0; i < n; ++i ){
      
        X_i  = X[i]; 
        Z_i = gamma(m,0) * (X_i - tau0[m]); 
        if (X_i <= tau0[m]){ res(m,i,0)  = 1;} else{ res(m,i,0) = exp(-lambda0_m * (X_i - tau0[m]));} 
      
        res(m,i,1)  = (1-res(m,i,0))/(1 + exp(gamma(m,1) + Z_i ));
        res(m,i,Q)  = (1-res(m,i,0)) *(1 - 1/(1 + exp(gamma(m,Q-1) + Z_i)));
        for (int q = 2; q < Q; q++){
          res(m,i, q)  = (1-res(m,i,0))*(1/(1 + exp(gamma(m,q) + Z_i)) - 1/(1 + exp(gamma(m,q-1) + Z_i)));
        }
      }
    }
    return res;
  }
  
  
// [[Rcpp::export]]   
arma::cube probCondVectRcppRraw(
      IntegerVector X, 
      NumericVector tau0, 
      NumericVector lambda0, 
      NumericMatrix gamma
  ) 
  {
    
    int  n = X.length(); 
    int  Q = gamma.ncol();   
    int X_i  = 1;
    int M = tau0.length();
    
    
    double Z_i;
    arma::cube res(M,n,Q + 1); 
    
    for (int m = 0; m < M; ++m ){
      double lambda_m  =  lambda0[m];
      
      for (int i = 0; i < n; ++i ){
        
        X_i  = X[i]; 
        Z_i = gamma(m,0) * (X_i - tau0[m]); 
        if (X_i <= tau0[m]){ res(m,i,0)  = 1;} else{ res(m,i,0) = exp(-lambda_m * (X_i - tau0[m]));} 
        
        res(m,i,1)  = 1/(1 + exp(gamma(m,1) + Z_i ));
        res(m,i,Q)  = 1 - 1/(1 + exp(gamma(m,Q-1) + Z_i));
        for (int q = 2; q < Q; q++){
          res(m,i, q)  = 1/(1 + exp(gamma(m,q) + Z_i)) - 1/(1 + exp(gamma(m,q-1) + Z_i));
        }
      }
    }
    return res;
  }
    
