####################################################################
# Transfo parametres
###################################################################

# From gamma to delta 
fromGammaToDelta <- function(gamma) {
  
  if (is.vector(gamma)) {
    gamma = matrix(gamma, nrow = 1, ncol = length(gamma))
  }
  
  Q <- ncol(gamma)
  delta  <- gamma; 
  delta[,1] <- log(gamma[, 1])
  if (Q > 2) {
    delta[, 3:Q] <- log(gamma[, 2:(Q - 1)] - gamma[, (3:Q)])
  }
  if (nrow(delta) == 1) { delta = c(delta) }
  return(delta)
}


# From delta to gamma 
fromDeltaToGamma <- function(delta) {
  
  if (is.vector(delta)) {
    delta = matrix(delta, nrow = 1, ncol = length(delta))
  }
  M <- nrow(delta)
  Q <- ncol(delta)
  gamma  <- delta; 
  gamma[, 1] <- exp(delta[, 1])
  if (Q > 2) {
    gamma[, 3:Q] <- matrix(rep(delta[, 2], Q - 2), ncol = (Q - 2), byrow = FALSE) - 
      t(apply(matrix(exp(delta[, 3:Q]), nrow = M, ncol = Q - 2), 1, cumsum))
  }
  if (M == 1) { gamma = c(gamma) }
  return(gamma)
}