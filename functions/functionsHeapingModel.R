
source('functions/functionsHeapingProbability.R')


####################################################################
# Arrondir un element --------------------------------------------------
####################################################################
funcHeaping <- function(x, ent){
  # x is a  vector 
  # ent is either a single integer (1,5,10,100,1000... ) or a vector of the same length as x
  if ((length(ent) == 1) | (length(x) == length(ent))) {
    y <- ifelse(x %% ent == 0, 
                x,
                ifelse(x %% ent > ent / 2, 
                       ceiling(x / ent ) * ent,
                       floor(x / ent ) * ent
                       )
              )
    return(y)
  } else {
      stop("ent is\n either a single integer (1,5,10,100,1000... )\n or a vector of the same length as x")
    }
}


####################################################################
# Simulation des données (x, g et y) selon les lois choisis
#
####################################################################
# 
# --------------------- sans covariables
simulHeapedData <- function(n,
                            heaping_levels,
                            param,
                            model = list(distrib = 'nbinom',
                                         heaping = 'twoSteps'
                                         ),
                            oneShifted = TRUE,
                            quiet = TRUE
                            ) {

  #------------Parameters------------
  param_heaping <- param$param_heaping
  param_distrib <- param$param_distrib

  if ((is.null(model$distrib) == FALSE) & (model$distrib %in% c("pois", "nbinom") == FALSE)) {
    stop("Misspecified model\n\tShould be 'pois' or 'nbinom'")
  }
  
  mu  <- exp(param_distrib$log_mu)
  if ((length(mu) == 1) & (n > 1)) { 
    mu <- rep(mu, n)
    sort_X <- TRUE
  } else {
    sort_X <- FALSE
  }

  # about the over dispersion
  over_dispersion <- exp(param_distrib$transfo_over_dispersion) + 1

  if (is.null(over_dispersion) & (model$distrib == 'nbinom')) {
    stop('An over_dispersion parameter must be specified\n\tin the Neg. Binom Model')
  }

  if (is.null(over_dispersion) & (model$distrib != 'nbinom')) {
    over_dispersion <- 1; model <- 'pois'
  }
  if(isFALSE(quiet)) {
   if(model$distrib == 'pois') {
    writeLines("Using a Poisson distribution for data generation")
   } else {
    writeLines("Using a Negative Binomial distribution for data generation")
   }    
  }
  #--------------- simulation of the X (true counts)--------------

  if (model$distrib == 'nbinom') {
    theta_p <- 1 / over_dispersion
    theta_r  <- mu * theta_p / (1 - theta_p)
    X <- rnbinom(n, size = theta_r, prob = theta_p) + oneShifted # X
  } else {
    X <- rpois(n, mu) + oneShifted # X
  }
  # only sort the data when you have no covariates
  # otherwise it destroys the signal
  if(sort_X) { X <- sort(X) }

  #----------------- simulation of the heaped data
  prob_heaping <- switch(model$heaping,
    triangular = prob_heaping_triangular,
    logit = prob_heaping_logit,
    twoSteps = probHeaping,
    twoStepsRcpp = probHeapingRcpp
  )
  n_levels <- length(heaping_levels)# + ifelse(model$heaping == "twoSteps", 1, 0)

  if ((n_levels -1) != max(length(param_heaping$gamma),length(param_heaping$delta))){
    
    stop("Length of gamma or delta not matching with the number of heaping levels. ")
  }

  prob_X <- prob_heaping(X, param_heaping)
  sim_G <- vapply(1:n,
                  function(i) {
                    sample(1:n_levels,
                           1,
                           prob = c(prob_X[i, ])
                           )
                    }, 1
                  )
  Y <- funcHeaping(X, heaping_levels[sim_G])
  return(data.frame(Y, X, sim_G))
}

# --------------------- avec covariables

simulHeapedDataCovar <- function(covar,
                                 param,
                                 heaping_levels,
                                 model = list(distrib = 'nbinom',
                                              heaping = 'twoSteps'
                                              ),
                                 ...
                                 ) {

  # covar : matrix of size n * p (set the first column equal to one if you need an intercept)
  # beta : must be a vector of length beta

  n <- dim(covar)[1];
  p <- dim(covar)[2];
  if (length(param$param_distrib$beta) != p) {
    stop('Mismatched dimensions between covar and beta\n\tPlease check')
  }

  param_2 <- param;
  param_2$param_distrib$log_mu <- covar %*% param$param_distrib$beta
  return(simulHeapedData(n, heaping_levels, param_2, model, ...))
}

 

####################################################################
# Likelihood function
####################################################################


#--------------------------------------------------------------------------------------------------------------
# Get all the possible preimages by heaping functions of a vector of observations  : not public function (if in a package)
#--------------------------------------------------------------------------------------------------------------  
heapingBehaviorTable <- function(dataY, heaping_levels) {
  
  Y_max <- max(dataY) + max(heaping_levels) # to make sure you encompass all possible couples (G, X) that correspond to Y
  G <- length(heaping_levels)
  heaping <- data.frame(X = rep(1:Y_max, each = G), 
                        G_behavior = rep(1:G, time = Y_max)
                        )
  heaping$Y <- funcHeaping(heaping$X, 
                            heaping_levels[heaping$G_behavior]
                           )
  tableau <- cbind(heaping,
                   do.call('rbind', 
                           replicate(Y_max, diag(1, G, G), simplify = FALSE)
                           )
                   )
  
  tableau_4 = c()
  indic_data <- c()
  n  = length(dataY)
  for (i in 1:n) {
    U <- which(tableau$Y == dataY[i])
    indic_data <- c(indic_data, rep(i, length(U)))
    tableau_4 = rbind(tableau_4, tableau[U, ])
  }
  res <- tableau_4
  res$ind_Y <- indic_data
  
  return(res)
}
#---------------------------------------------------------------------------------

# Preparing data
#----------------------------------------------

prepareData <- function(Y,
                        covar = NULL,
                        heaping_levels
                        ) {
  if (is.null(covar)) {
    repl <- table(Y) # ordonne les Y
    Y <- sort(unique(Y))
  } else {
    repl <- rep(1, length(Y))
  }
  data = list()
  data$repl <- repl;
  data$my_table <- heapingBehaviorTable(Y, heaping_levels)
  data$my_table <- data$my_table[order(data$my_table$X), ]
  data$n <- length(Y)
  data$counts_X <- table(data$my_table$X)
  data$my_table$ind_X <- rep(1:length(unique(data$my_table$X)),
                             data$counts_X
                             )
  #### data_mytable construit adhoc. Ne pas changer les ordres des colonnes.
  #data$covar <- covar[unique(data$my_table$ind_Y),]
  data$covar <- covar
  data$unique_X <- unique(data$my_table$X)
  data$n_levels <- length(heaping_levels) - 1
  return(data)
}

 
#################################################################

logLikelihoodFunction <- function(param, 
                                  data, 
                                  model, 
                                  oneShifted = TRUE,
                                  log = TRUE
                                  ) {
  
  #------------- previously prepared data
  n <- data$n
  my_table  <- data$my_table
  repl <- data$repl
  N <- nrow(my_table)
  
  #-------------- heaping model 
  prob_heaping_vect <- switch(model$heaping,
                              triangular = prob_heaping_triangular_vect,
                              logit = prob_heaping_logit_vect,
                              twoSteps = probHeapingVect,
                              twoStepsRcpp = probHeapingRcppVect
                              )
  
  from_gamma_to_delta <- switch(model$heaping,
                                triangular = from_gamma_to_delta_triangular,
                                fromGammaToDelta
                                )
  
  #---------------- parameters 
  param_distrib <- param$param_distrib
  param_heaping <- param$param_heaping
  if (is.null(param_heaping$delta)) {
    if (is.null(param_heaping$gamma)) { stop('Requiring delta or gamma') }
    param_heaping$delta <- t(apply(param_heaping$gamma, 1, from_gamma_to_delta))
  }
  
  if (is.matrix(param_heaping$delta)) {
    Q <- ncol(param_heaping$delta)
    M <- nrow(param_heaping$delta)
  }
  if (is.vector(param_heaping$delta)) {
    M <- 1; 
    Q <- length(param_heaping$delta)
  }
  
  if (is.null(data$covar)) { 
    mu <- matrix(exp(param_distrib$log_mu), nrow = M, ncol = n, byrow = FALSE) 
  } else {
    if (is.null(param_distrib$beta)) { stop('parameter beta required') }
    mu <- t(exp(data$covar %*% t(param_distrib$beta)))
  }
  
  if (model$distrib == "nbinom") {
    over_dispersion <- matrix(exp(param_distrib$transfo_over_dispersion) + 1,
                              nrow = M, ncol = n, byrow = FALSE
                              )
    theta_p <- 1 / over_dispersion
    theta_r  <- mu * theta_p / (1 - theta_p)
  }
  
  #--------------- Computing Prob_heaping(G|X)
  
  U <- prob_heaping_vect(data$unique_X,param_heaping)[, my_table$ind_X, ]
  if (M > 1) {
    ProbG <- sumArray(U, my_table$G_behavior)
    all_prob <- data.frame(ind_M  = rep(1:M, each = N), 
                           ind_Y = rep(my_table$ind_Y, M),
                           X = rep(my_table$X, M),
                           Y = rep(my_table$Y, M),
                           G = rep(my_table$G_behavior, M),
                           ProbG = ProbG
                           )
  } else {
    all_prob <- data.frame(ind_Y = my_table$ind_Y,
                           X = my_table$X,
                           Y = my_table$Y,
                           G = my_table$G_behavior,
                           ProbG = rowSums(U * my_table[, 4 + 0:Q])
                           )
  }
  #--------------- prob(X)
  theta <- switch(model$distrib,
                  nbinom = { cbind(c(t(theta_r[, my_table$ind_Y])), c(t(theta_p[, my_table$ind_Y]))) },
                  pois  = { c(t(mu[, my_table$ind_Y])) }
                  )
 
  all_prob <- all_prob %>% 
    mutate(ProbX = switch(model$distrib,
                          nbinom = dnbinom(all_prob$X - oneShifted, size = theta[, 1], prob = theta[, 2], log = FALSE),
                          pois = dpois(all_prob$X - oneShifted, theta, log = FALSE)
                          ),
           ProbXProbG = ProbX * ProbG
           ) 
  if(M > 1) {
    all_prob %>%  
    group_by(ind_M, ind_Y) %>% 
    summarize(L = sum(ProbXProbG)) %>% 
    select(ind_M, ind_Y, L) -> d
  } else {
    all_prob %>%  
      group_by(ind_Y) %>% 
      summarize(L = sum(ProbXProbG)) %>% 
      select(ind_Y, L) -> d
  }
  
  if (M > 1){
    if (log == FALSE) {
      dd <- data.frame(ind_M = d$ind_M, ll = d$L^rep(repl, M))
      dd %>% group_by(ind_M) %>% summarize(LL = prod(ll)) %>% select(LL) -> res
    } else {
      dd <- data.frame(ind_M = d$ind_M, ll = log(d$L)*rep(repl, M))
      dd %>% group_by(ind_M) %>% summarize(LL = sum(ll)) %>% select(LL) -> res   
    }
    res <- unlist(c(res))
  }
  if (M == 1){
    res <- ifelse(log, sum(repl * log(d$L)), prod(d$L^repl))
  }
  return(res)
}

#################################################################
# for use with optim
nb2steps_logLikelihood <- function(param, # lambda_0, tau_0, gamma, transfo_over_dispersion, beta
                                   data, # output of a call to PrepareData
                                   oneShifted = TRUE
                                   ) {
  
  #------------- previously prepared data
  n <- data$n
  my_table  <- data$my_table
  repl <- data$repl
  N <- nrow(my_table)
  
  #-------------- heaping model 
  prob_heaping <- probHeapingRcpp
  from_gamma_to_delta <- fromGammaToDelta
  
  #---------------- parameters
  param_heaping <- list(lambda_0 = param[1], 
                        tau_0 = param[2],
                        threshold = 0.05,
                        gamma = param[2 + 1:data$n_levels]
                        ) 
  param_distrib <- list(transfo_over_dispersion = param[3 + data$n_levels],
                        beta = param[-c(1:(3 + data$n_levels))]
                        )
  if (is.null(param_distrib$beta)) { stop('parameter beta required') }
  if (is.null(param_heaping$delta)) {
    if (is.null(param_heaping$gamma)) { stop('Requiring delta or gamma') }
    param_heaping$delta <- from_gamma_to_delta(param_heaping$gamma)
  }
  
  M <- 1; 
  Q <- length(param_heaping$delta)
  
  if (is.null(data$covar)) {
    if(length(param_distrib$beta) != M) { stop('wrong dimension of vector param') }
    mu <- matrix(exp(param_distrib$beta), nrow = M, ncol = n, byrow = FALSE) 
  } else {
    mu <- t(exp(data$covar %*% param_distrib$beta))
  }
  
  over_dispersion <- matrix(exp(param_distrib$transfo_over_dispersion) + 1,
                            nrow = M, ncol = n, byrow = FALSE
                            )
  theta_p <- 1 / over_dispersion
  theta_r  <- mu * theta_p / (1 - theta_p)
  
  #--------------- Computing Prob_heaping(G|X)
  U <- prob_heaping(data$unique_X, param_heaping)[my_table$ind_X, ]
  all_prob <- data.frame(ind_Y = my_table$ind_Y,
                         X = my_table$X,
                         Y = my_table$Y,
                         G = my_table$G_behavior,
                         ProbG = rowSums(U * my_table[, 4 + 0:Q])
                         )
  #--------------- prob(X)
  theta <- cbind(c(t(theta_r[, my_table$ind_Y])), c(t(theta_p[, my_table$ind_Y])))
  
  res <- all_prob %>% 
    mutate(ProbX = dnbinom(all_prob$X - oneShifted, size = theta[, 1], prob = theta[, 2], log = FALSE),
           ProbXProbG = ProbX * ProbG
           ) %>%  
    group_by(ind_Y) %>% 
    summarize(L = sum(ProbXProbG)) %>% 
    select(ind_Y, L) %>% 
    mutate(ll = log(L) * repl) %>%
    summarize(LL = sum(ll)) %>% 
    as.numeric()
  
  return(res)
}

# for use with optim
pois2steps_logLikelihood <- function(param, # lambda_0, tau_0, gamma, beta
                                     data, # output of a call to PrepareData
                                     oneShifted = TRUE
                                     ) {
  
  #------------- previously prepared data
  n <- data$n
  my_table  <- data$my_table
  repl <- data$repl
  N <- nrow(my_table)
  
  #-------------- heaping model 
  prob_heaping <- probHeapingRcpp
  from_gamma_to_delta <- fromGammaToDelta
  
  #---------------- parameters
  param_heaping <- list(lambda_0 = param[1], 
                        tau_0 = param[2],
                        threshold = 0.05,
                        gamma = param[2 + 1:data$n_levels]
                        ) 
  param_distrib <- list(beta = param[-c(1:(2 + data$n_levels))])
  if (is.null(param_distrib$beta)) { stop('parameter beta required') }
  if (is.null(param_heaping$delta)) {
    if (is.null(param_heaping$gamma)) { stop('Requiring delta or gamma') }
    param_heaping$delta <- from_gamma_to_delta(param_heaping$gamma)
  }
  
  M <- 1; 
  Q <- length(param_heaping$delta)
  
  if (is.null(data$covar)) {
    if(length(param_distrib$beta) != M) { stop('wrong dimension of vector param') }
    mu <- matrix(exp(param_distrib$beta), nrow = M, ncol = n, byrow = FALSE) 
  } else {
    mu <- t(exp(data$covar %*% param_distrib$beta))
  }
  
  #--------------- Computing Prob_heaping(G|X)
  U <- prob_heaping(data$unique_X, param_heaping)[my_table$ind_X, ]
  all_prob <- data.frame(ind_Y = my_table$ind_Y,
                         X = my_table$X,
                         Y = my_table$Y,
                         G = my_table$G_behavior,
                         ProbG = rowSums(U * my_table[, 4 + 0:Q])
                         )
  #--------------- prob(X)
  theta <- c(t(mu[, my_table$ind_Y]))
  
  res <- all_prob %>% 
    mutate(ProbX = dpois(all_prob$X - oneShifted, theta, log = FALSE),
           ProbXProbG = ProbX * ProbG
           ) %>%  
    group_by(ind_Y) %>% 
    summarize(L = sum(ProbXProbG)) %>% 
    select(ind_Y, L) %>% 
    mutate(ll = log(L) * repl) %>%
    summarize(LL = sum(ll)) %>% 
    as.numeric()
  
  return(res)
}
