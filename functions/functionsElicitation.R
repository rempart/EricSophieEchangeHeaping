##--------------------------------------------------------------------------------------------------------
## SCRIPT : Tester les fonctions arrondies 'Two Steps'
## Authors : Isabelle Albert, Eric Parent, Sophie Donnet & Matthieu Authier
## Last update : 2019-09-13
## R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
## Copyright (C) 2018 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##---


library(ggplot2)
source('functions/functionsHeapingProbabilityRcpp.R')
 




#---------------------------------------------------------- 
# Adjusts Probability  curves to elicited data  
#-----------------------------------------------------------
 

inverseConstraints  <- function(elicitedProb,elicitedValues, method = "Mode", lastCurve = FALSE){
  
  echan <- elicitedValues
  
  if (is.vector(echan)) { echan = as.data.frame(matrix(echan, nrow = 1)) }
  
  N <- ncol(echan)
  Q <- N - 2
  M <- nrow(echan)
  if (length(elicitedProb) != N) { stop("Dimension mismatch between elicitedProb and elicitedValues") }
  
  Names <-   c('tau_0','tau_1','1_1','1_2')
  if (Q > 3) {Names <- c(Names, paste(2:(Q - 2)))}
  if (lastCurve) {
    Names <- c(Names, Q)
  }else{
    Names <- c(Names, Q - 1)
  }
  names(echan) = Names
  
  
  #------------- GET PARAMETERS ----------------------  
  param_heaping_estim = list(tau_0 = echan$tau_0)
  param_heaping_estim$lambda_0 <- -log(elicitedProb[2]) / (echan$tau_1 - echan$tau_0)
  
  #---------------------------------------------
  p11 =  elicitedProb[3]; 
  p12 =  elicitedProb[4]; 
  x11 =  echan[, 3]; 
  x12 =  echan[, 4]; 
  
  param_heaping_estim$gamma = matrix(0,M,Q)
  param_heaping_estim$gamma[,1] <- (logit(1 - p11) - logit(1 - p12)) / (x11 - x12)
  param_heaping_estim$gamma[,2] <- logit(1 - p11) -   param_heaping_estim$gamma[,1]*(x11 - param_heaping_estim$tau_0)
  
  
  q <- 2; 
  
  if (method != 'Mode') {
    if (Q > 3){
      upto = ifelse(lastCurve,Q + 1,Q + 2)
      for (i in 5 : upto) {
        
        pqp1 <- elicitedProb[i]
        
        #---------------------------------------------------
        vq <-  inv.logit(param_heaping_estim$gamma[,1] * (echan[ , i] - param_heaping_estim$tau_0) + param_heaping_estim$gamma[,q]) - pqp1
        w <- which(vq > 0)
        if (length(w) < M) {
          resampl <- sample(w,M,replace = TRUE)
          echan <- echan[resampl,]
          param_heaping_estim$gamma  <- param_heaping_estim$gamma[resampl,]
          param_heaping_estim$tau_0 <- param_heaping_estim$tau_0[resampl]
          param_heaping_estim$lambda_0 <- param_heaping_estim$lambda_0[resampl]
        }
        
        #-------------------------------------------------
        
        vq <-  inv.logit(param_heaping_estim$gamma[,1] * (echan[ , i] - param_heaping_estim$tau_0) + param_heaping_estim$gamma[, q]) - pqp1
        param_heaping_estim$gamma[, q + 1] = logit(vq) - param_heaping_estim$gamma[, 1] * (echan[, i] - param_heaping_estim$tau_0)
        q <- q + 1
      }
    }
    
  }else{
    if ((Q > 3) | (!lastCurve)) {
      upto = ifelse(lastCurve,Q + 1,Q + 2)
      for (i in 5:upto) {
        x_q <- echan[,i]; 
        gamma_0 <- param_heaping_estim$gamma[,1]
        param_heaping_estim$gamma[,q + 1] <- -2 * gamma_0*(x_q- param_heaping_estim$tau_0) - param_heaping_estim$gamma[,q]
        q <- q + 1
      }
    }
  }
  
  if (lastCurve){
    param_heaping_estim$gamma[, Q ] <- logit(elicitedProb[Q + 2]) - param_heaping_estim$gamma[,1]*(echan[,Q + 2] - param_heaping_estim$tau_0)
  }
    
    
   

  
  return(param_heaping_estim)
}


#----------------------------------------------------
# SAMPLING ELICITED VALUE
#---------------------------------------------------------------
sampleElicitedValue <- function(MC, elicited_range) {
  N <- nrow(elicited_range)
  S <- sapply(1:N, function(q) {
    Num <- log(elicited_range[q, 2]) - log(elicited_range[q,1])
    alpha <- (1 - elicited_range[q, 4]) / 2 
    D <- qnorm(1 - alpha) - qnorm(alpha) 
    sigma <- Num/D
    mu <-  log(elicited_range[q, 1]) - sigma * qnorm(alpha)
    echan <- exp(rnorm(MC, mu, sigma))
    return(echan)
  })
  W <- which(S[, 2] > S[, 1])
  S <- S[W, ]
  
  Test <- which(rowSums(t(apply(t(apply(S[, 3:N], 1, diff)), 1, function(x) { as.numeric(x <= 0) }))) == 0)
  S <- S[Test, ]
  colnames(S) <- rownames(elicited_range)
  return(as.data.frame(S)) 
}

#----------------------------------------------------------------
#---------------------------------------------------------------

plotElicited <- function(echanParamElicited,elicitedTable, plotOptions=list() , lastCurve = FALSE){

  
  #------ Size of parameters
  param_elicited <- echanParamElicited
  M <- length(param_elicited$tau_0)
  N <- nrow(elicitedTable)
  n_Q <- Q <- N -  2
  if (ncol(elicitedTable) == 2) { names(elicitedTable) = c('x', 'P') }
  
  
  
  #--------------------------------
  #----- Options
  currentOptions <- list(
    xmax = NULL,
    heaping_levels = 1:Q,
    color_palette_brewer  = "Set2")
  currentOptions[names(plotOptions)] <- plotOptions
  if (is.null(currentOptions$xmax)) {
    currentOptions$xmax  = 1.1 * mean( (logit(0.9999999) - echanParamElicited$gamma[,Q]) / echanParamElicited$gamma[,1] + echanParamElicited$tau_0)
  }
  heaping_levels = as.character(currentOptions$heaping_levels)
  abs = 1:currentOptions$xmax
  palette <- currentOptions$color_palette_brewer
  
  Names <-   c('tau_0','tau_1','1_1','1_2')
  if (Q > 3) {Names <- c(Names, paste(2:(Q - 2)))}
  if (lastCurve) {
    Names <- c(Names, Q)
  }else{
    Names <- c(Names, Q - 1)
  }
  row.names(elicitedTable) = Names
  #-------- Compute Curves-------------------------------- 
  Prestim <- probCondRcppVect(X = abs, param_elicited)

  #----- Formatting data
  
  RES_prob_cond <- Prestim[,, 2:(Q + 1)]
  if (M > 1) {
    mean_cond <- apply(RES_prob_cond, c(2, 3), mean)
  } else{
    mean_cond = RES_prob_cond
  }
  w <- which(mean_cond[, 1] > 1 | mean_cond[, 1] < 0) ### remove negative proba and proba above 1
  if (length(w) > 0) { 
    pr_rm <- max(w)
    mean_cond[1:pr_rm, ] <- NA
  } 
  U <- as.data.frame(mean_cond)
  names(U) <- heaping_levels

  U <- U %>% mutate(xvalue = abs) %>% gather(key = Q, value = value, -xvalue);
  U <- U %>% add_column("type" = rep("cond", length(abs) * Q))

  RES_prob_heaping <- 1 - Prestim[,,1]
  if (M > 1) {mean_heaping <- apply(RES_prob_heaping,2,mean)}else{mean_heaping = RES_prob_heaping}
  V <- data.frame(xvalue = c(abs,abs))
  V$Q <- rep(c('heaping','no heaping'),each = length(abs))
  V$value <- c(mean_heaping,1 - mean_heaping)
  V$type = rep('heaping',2*length(abs))
  
  if (M > 1) {
    max_cond <- apply(RES_prob_cond, c(2, 3), function(x) { quantile(x, 0.95) })
    min_cond <- apply(RES_prob_cond, c(2, 3), function(x) { quantile(x, 0.05) })
    q3_cond  <- apply(RES_prob_cond, c(2, 3), function(x) { quantile(x, 0.75) })
    q1_cond  <- apply(RES_prob_cond, c(2, 3), function(x) { quantile(x, 0.25) })
    U$min <- c(min_cond)
    U$max <- c(max_cond)
    U$q3 <- c(q3_cond)
    U$q1 <- c(q1_cond)
    
    min_heaping <- apply(RES_prob_heaping, 2, function(x) { quantile(x, 0.05) })
    max_heaping <- apply(RES_prob_heaping, 2, function(x) { quantile(x, 0.95) })
    q3_heaping  <- apply(RES_prob_heaping, 2, function(x) { quantile(x, 0.75) })
    q1_heaping  <- apply(RES_prob_heaping, 2, function(x) { quantile(x, 0.25) })
    V$min <- c(min_heaping,1 - min_heaping)
    V$max <- c(max_heaping,1 - max_heaping)
    V$q3 <- c(q3_heaping,1 - q3_heaping)
    V$q1 <- c(q1_heaping,1 - q1_heaping)
  }
  
  W <- rbind(U,V)
  W$Heaping <- factor(W$Q,levels = c("heaping", "no heaping",heaping_levels[1:Q]))
  
  elicitedTable$type = c(rep('heaping', 2), rep('cond', N - 2))
  elicitedTable$P[1] = 1;
  elicitedTable$Q <- c("no heaping","no heaping",heaping_levels[1],heaping_levels[1:(Q - 2)],ifelse(lastCurve,heaping_levels[Q],heaping_levels[Q-1])) 
  elicitedTable$Heaping <- factor(elicitedTable$Q,levels = c("no heaping","heaping",heaping_levels[1:Q]))
  
  # graph part ----
  
  # if first graph with no previous graph
  
  gg <- W %>% 
    mutate(Q = factor(Q, levels = c('no heaping', 'heaping', as.character(1:n_Q)))) %>%
    ggplot() + 
    geom_line(aes(x = xvalue, y = value, color = Heaping, group = Heaping), size = 1) 
  
  if (M > 1) {
    elicitedTable_u <- elicitedTable

    if ((Q > 3) | (!lastCurve)) {
      upto = ifelse(lastCurve,Q + 1,Q + 2)
      u <- c(5:(upto)); elicitedTable_u <- elicitedTable[-u,]
    } 
    gg <- gg + 
      geom_ribbon(aes(ymin = min, ymax = max, x = xvalue, group = Heaping, fill = Heaping), alpha = 0.1) +
      geom_ribbon(aes(ymin = q1, ymax = q3, x = xvalue, group = Heaping, fill = Heaping), alpha = 0.3) +
      geom_vline(data = elicitedTable,
                 mapping = aes(xintercept = xinf, group = type, color = Heaping), 
                 linetype = "dashed"
                 ) +
      geom_vline(data = elicitedTable, 
                 mapping = aes(xintercept = xsup, group = type, color = Heaping), 
                 linetype = "dashed"
                 ) +
      geom_hline(data = elicitedTable_u, mapping = aes(yintercept = P, group = type), color = 'grey')
  }
  if (M == 1) {
    if ((N > 5) | (!lastCurve)) {
      upto = ifelse(lastCurve,Q + 1,Q + 2)
      elicitedTable$P[5:upto] = probCondRcppVect(X = elicitedTable$x[5:upto], param_elicited)[1,,][(2:(Q - 2)) + 1]
    }

    gg <- gg + 
      geom_vline(data = elicitedTable, 
                 mapping = aes(xintercept = x,group = type,color = Heaping),
                 linetype = "dashed") +
      geom_hline(data = elicitedTable,  mapping = aes(yintercept = P,group = type),colour = 'grey') +
      geom_point(data = elicitedTable, mapping = aes(x = x,y =  P, group = type))
  }
  
  gg <- gg + 
    ylim(0.0, 1.0) + scale_x_continuous(name = "True count") + ylab("Probability") +
    scale_fill_brewer(name = "Heaping Behaviour: ", type = "qual", palette = palette) +
    scale_color_brewer(name = "Heaping Behaviour: ", type = "qual", palette = palette) +
    facet_wrap(~ Heaping, ncol = 2) +
    theme(legend.position = "top", 
          plot.title = element_text(lineheight = 0.8, face = "bold"), 
          axis.text = element_text(size = 10),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
  return(gg)
  print(gg)
  
}
#----------------------------------------------------
#------------------------Ajustement des gammas to Q couples (X_q,P_q(X_q))--------------------------------------------
#----------------------------------------------------
computeHyperparametersFromElicited <- function(echanParamElicited){
  
  Q <- ncol(echanParamElicited$gamma)
  Theta <- cbind(log(echanParamElicited$tau_0), log(echanParamElicited$lambda_0));
  Theta <- cbind(Theta,fromGammaToDelta(echanParamElicited$gamma))
  Names <-  c("log_tau_0","log_lambda_0",vapply(0:(Q - 1),function(q){paste('delta',q,sep = '')},'1'))
  colnames(Theta) <-  Names
  mu <- colMeans(Theta)
  Sigma <- cov(Theta)
  return(prior = list(mu = mu,Sigma = Sigma))
}



 

#----------------------------------------------------
#--------- LOGIT ------------------------------------
#----------------------------------------------------

logit <- function(x){log(x / (1 - x))}

inv.logit <- function(x){
  p <- exp(x)/(1 + exp(x))
  p <- ifelse(is.na(p) & !is.na(x), 1, p)
  return(p)
}


