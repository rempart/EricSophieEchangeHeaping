source('functions/functionsTransfoParam.R')
#-------- Two Steps------------------------------------- 
probHeaping = function(X, 
  param_heaping = list(tau_0 = 5,
    lambda_0 = 2,
    gamma =  c(0.27, -4.5, -12),
    delta = NULL), 
    op.output.cond =  FALSE
  
) {
  # no_heaping : value below which no heaping with certainty
  # sup_heaping : value above which heaping with quasi-certainty
  # gamma : thresholds for heaping behaviour probabilities, with gamma[1] > 0 and gamma[2] > gamma[3] > ...
  # cdf : cumulative density function, e.g. plogis (the inverse logit link)
  # delta : diff(gamma)
  # param_heaping$delta of size Q (if number of regimes = Q + 1)
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
  
  #n <- length(X)
  lambda <- param_heaping$lambda_0; 

  prob_non_heaping <- ifelse(X <= param_heaping$tau_0, 
    1, 
    exp(-lambda * (X - param_heaping$tau_0))
  )
  prob_heaping <-  1 - prob_non_heaping
  
  #####  Heaping behaviour  #####
  Q <- length(gamma)
  Z <- gamma[1] * (X - param_heaping$tau_0)
  # conditional proba
  mat <- matrix(0, nrow = length(X), ncol = Q + 1)
  mat[, Q + 1] <- 1
  mat[, 2:Q] <- sapply(gamma[-1], function(xsi) { plogis(-xsi - Z) })
  res1 <- mat[, 2:(Q + 1)] - mat[, 1:Q]
  # marginal proba
  res2 <- matrix(rep(prob_heaping, Q), ncol = Q, byrow = FALSE) * res1
  
  if(op.output.cond ==  TRUE){
    
    res <- list(prob_heaping = prob_heaping, prob_cond = res1)
    
  } else {
    res <- cbind(matrix(prob_non_heaping, ncol = 1), res2)
  }
  return(res)
}


#-------- Two Steps------------------------------------- 

probHeapingVect = 
  function(X, 
    param_heaping = list(tau_0 = rep(5, 2),
      lambda_0 = rep(2, 2),
      gamma = matrix(rep(c(0.3, -4, -12), 2), nrow = 2, byrow = TRUE),
      delta = NULL#, 
      # cdf = "logit"
    ),
    return_array = TRUE
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
    
    
    n <- length(X)
    M <- nrow(gamma)
    M_1 <- length(param_heaping$tau_0)
    M_3 <- length(param_heaping$lamda_0)
    
    if (!all(c(M_1,M_3)== rep(M,3))) {
     stop("mismatching sizes of the samples of the parameters")
    }
    
    lambda <- param_heaping$lambda_0

    ### matrices of dim(n, Q,M)
    prob_non_heaping <- sapply(1:M, 
      function(i) {
        ifelse(X <= param_heaping$tau_0[i], 1, exp(-lambda[i] * (X - param_heaping$tau_0[i])))
      })
    prob_heaping <-  1 - prob_non_heaping
    
    #####  Heaping behaviour  #####
    Q <- ncol(gamma)
    Z <- outer(X, param_heaping$tau_0,'-') * matrix(gamma[, 1], nrow = n,ncol = p, byrow=TRUE)
    # conditional proba
    foo <- function(i) {
      mat <- matrix(0, nrow = length(X), ncol = Q + 1)
      mat[, Q + 1] <- 1
      mat[, 2:Q] <- sapply(gamma[i, -1], function(xsi) { plogis(-xsi - Z[, i]) })
      res <- mat[, 2:(Q + 1)] - mat[, 1:Q]
      # marginal proba
      res <- matrix(rep(prob_heaping[, i], Q), ncol = Q, byrow = FALSE) * res
      return(cbind(matrix(prob_non_heaping[, i], ncol = 1), res))
    }
    
    if(return_array) {
      Prob_marg <- array(NA, dim = c(M, n, Q + 1))
      for(i in 1:M) { Prob_marg[i, , ] <- foo(i) }
    }
    else {
      Prob_marg <- lapply(1:M, foo)
    }
    return(Prob_marg)
  }

