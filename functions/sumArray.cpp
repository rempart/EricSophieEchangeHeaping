#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; 
using namespace arma;
 

// [[Rcpp::export]]
NumericVector sumArray(arma::cube U, NumericVector indG) {
  
  
  int M = U.n_rows;
  int N = U.n_cols;
  int place = 0;
  int indGj = 0; 
  NumericVector res(M * N);
  // 
  // 
  for (int m = 0; m < M; m++){
     for (int j = 0; j < N; j++){
       indGj = indG(j) - 1;
       res(place) =  U(m,j,indGj); 
       place = place + 1 ;
     }
  }
  return res;
  
  // for (int m = 0; m < M; ++m ){
  //   U.m
  //   rowSums(U[m,,] * mty  
  // }
  // 
}




 