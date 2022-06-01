library(Rcpp)
library(tidyverse)
source('functions/functionsTransfoParam.R')
#-----------------------------------------------------------------------
Rcpp::sourceCpp('functions/probHeapingRcpp.cpp')
#------------------------------------------------------------------------
probHeapingRcpp = function(X, 
  param_heaping = list(tau_0 = 5,
    lambda_0 = 2,
    gamma =  c(0.27, -4.5, -12))){
  
  if (is.null(param_heaping$delta) & is.null(param_heaping$gamma)) {
    stop('delta or gamma (heaping parameters) have to be specified')
  } else {
    if (is.null(param_heaping$gamma)) { 
      gamma <- fromDeltaToGamma(delta)
    }
    else {
      gamma <- param_heaping$gamma
    }
  }
  
  res <- probHeapingRcppRaw(
    X,
    tau0 = param_heaping$tau_0, 
    lambda0 = param_heaping$lambda_0,
    gamma = param_heaping$gamma
  ) 
  return(res)
  
}

####################################################################
### ---------------------------------------------------------------
####################################################################
probHeapingRcppVect = 
  function(X, 
    param_heaping = list(
      tau_0 = rep(5, 2),
      lambda_0 = rep(2, 2),
      gamma = matrix(rep(c(0.3, -4, -12), 2), nrow = 2, byrow = TRUE),
      delta = NULL#, 
      # cdf = "logit"
    )
  ) {
    # sanity checks
     
    if (is.null(param_heaping$delta) & is.null(param_heaping$gamma)) {
      stop('delta or gamma matrix (heaping parameters) has to be specified')
    } else {
      if (is.null(param_heaping$gamma)) { 
        gamma <- fromDeltaToGamma(delta)
      }
      else {
        gamma <- param_heaping$gamma
      }
    }
    
    M_1 <- length(param_heaping$tau_0)
    M <- nrow(gamma)
    if (M_1 == 1) {M <- 1; gamma <- matrix(gamma,nrow = 1)}
    M_3 <- length(param_heaping$lambda_0)
    
    if (!all(c(M_1,M_3) == rep(M,2))) {
      stop("mismatching sizes of the samples of the parameters")
    }
    
    res <- probHeapingVectRcppRaw(X, 
     tau0 = param_heaping$tau_0, 
     lambda0 = param_heaping$lambda_0,
      gamma = param_heaping$gamma)
    return(res);
  }
    

####################################################################
### ---------------------------------------------------------------
####################################################################
probCondRcppVect = 
  function(X, 
           param_heaping = list(
             tau_0 = rep(5, 2),
             lambda_0 = rep(2, 2),
             gamma = matrix(rep(c(0.3, -4, -12), 2), nrow = 2, byrow = TRUE),
             delta = NULL#, 
             # cdf = "logit"
           )
  ) {
    # sanity checks
    
    if (is.null(param_heaping$delta) & is.null(param_heaping$gamma)) {
      stop('delta or gamma matrix (heaping parameters) has to be specified')
    } else {
      if (is.null(param_heaping$gamma)) { 
        gamma <- fromDeltaToGamma(delta)
      }
      else {
        gamma <- param_heaping$gamma
      }
    }
    
    M_1 <- length(param_heaping$tau_0)
    M <- nrow(gamma)
    if (M_1 == 1) {M <- 1; gamma <- matrix(gamma,nrow = 1)}
    M_3 <- length(param_heaping$lambda_0)
    
    if (!all(c(M_1,M_3) == rep(M,2))) {
      stop("mismatching sizes of the samples of the parameters")
    }
    
    res <- probCondVectRcppRraw(X,tau0 = param_heaping$tau_0, 
                                lambda0 = param_heaping$lambda_0,
                                gamma = param_heaping$gamma)
    return(res);
  }
####################################################################