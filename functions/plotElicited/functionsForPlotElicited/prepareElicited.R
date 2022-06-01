prepareElicited <- function(listElicited, plotOptions, M, N, n_Q, Q, lastCurve){
  
  echanParamElicited <- listElicited[[1]]
  elicitedTable <- listElicited[[2]]
  
  #------ Size of parameters
  param_elicited <- echanParamElicited
  
  if (ncol(elicitedTable) == 2) {names(elicitedTable) = c('x','P')}
  
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
  

  
  
  #-------- Compute Curves-------------------------------- 
  Prestim <- probCondRcppVect(X = abs, param_elicited)
  
  #----- Formatting data
  
  RES_prob_cond <- Prestim[,,2:(Q + 1)]
  if (M > 1) {
    mean_cond <- apply(RES_prob_cond,c(2,3),mean)
  }else{
    mean_cond = RES_prob_cond
  }
  w <- which(mean_cond[, 1] > 1 | mean_cond[, 1] < 0)### remove negative proba and proba above 1
  if (length(w) > 0) { 
    pr_rm <- max(w)
    mean_cond[1:pr_rm, ] <- NA
  } 
  U <- as.data.frame(mean_cond)
  names(U) <- heaping_levels
  
  U <- U %>% 
    mutate(xvalue = abs) %>% 
    gather(key = Q, value = value, -xvalue) %>%
    add_column("type" = rep("cond",length(abs)*Q))
  
  
  
  RES_prob_heaping <- 1 - Prestim[,,1]
  if (M > 1) {mean_heaping <- apply(RES_prob_heaping,2,mean)}else{mean_heaping = RES_prob_heaping}
  V <- data.frame(xvalue = c(abs,abs))
  V$Q <- rep(c('heaping','no heaping'),each = length(abs))
  V$value <- c(mean_heaping,1 - mean_heaping)
  V$type = rep('heaping',2*length(abs))
  
  if (M > 1) {
    max_cond <- apply(RES_prob_cond,c(2,3),function(x){quantile(x,0.95)})
    min_cond <- apply(RES_prob_cond,c(2,3),function(x){quantile(x,0.05)})
    q3_cond <- apply(RES_prob_cond, c(2, 3), function(x) { quantile(x, 0.75) })
    q1_cond <- apply(RES_prob_cond, c(2, 3), function(x) { quantile(x, 0.25) })
    U$min <- c(min_cond)
    U$max <- c(max_cond)
    U$q3 <- c(q3_cond)
    U$q1 <- c(q1_cond)
    
    min_heaping <- apply(RES_prob_heaping, 2,  function(x) {quantile(x, 0.05) })
    max_heaping <- apply(RES_prob_heaping, 2, function(x) {quantile(x, 0.95) })
    q3_heaping <- apply(RES_prob_heaping, 2, function(x) { quantile(x, 0.75) })
    q1_heaping <- apply(RES_prob_heaping, 2, function(x) { quantile(x, 0.25) })
    V$min <- c(min_heaping,1 - min_heaping)
    V$max <- c(max_heaping,1 - max_heaping)
    V$q3 <- c(q3_heaping,1 - q3_heaping)
    V$q1 <- c(q1_heaping,1 - q1_heaping)
  }
  
  W <- rbind(U,V)
  W$Heaping <- factor(W$Q,levels = c("heaping","no heaping", heaping_levels[1:Q]))
  
  elicitedTable$type = c(rep('heaping', 2), rep('cond', N - 2))
  elicitedTable$P[1] = 1;
  elicitedTable$Q <- c("no heaping","no heaping",heaping_levels[1],heaping_levels[1:(Q - 2)],ifelse(lastCurve,heaping_levels[Q],heaping_levels[Q-1])) 
  elicitedTable$Heaping <- factor(elicitedTable$Q,levels = c("no heaping","heaping",heaping_levels[1:Q]))
  
  return(list(
    elicitedTable = elicitedTable, 
    W = W
  ))
}
