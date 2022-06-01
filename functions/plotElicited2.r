library(ggrepel)
plotElicited2 <- function(echanParamElicited,
                         elicitedTable, 
                         plotOptions = list(xmax = 300),
                         lastCurve = FALSE
                         ) {
  
  #------ Size of parameters
  param_elicited <- echanParamElicited
  M <- length(param_elicited$tau_0)
  N <- nrow(elicitedTable)
  n_Q <- N -  2
  if (ncol(elicitedTable) == 2) { names(elicitedTable) = c('x', 'P') }
  
  #--------------------------------
  #----- Options
  currentOptions <- list(xmax = NULL, #300,
                         heaping_levels = 1:n_Q,
                         color_palette_brewer = "Set2"
                         )
  if(n_Q == 3) {
    currentOptions$color_palette_manual <- c('#d7191c', '#fdae61', '#6ca872', '#abd9e9', '#2c7bb6')
  }
  if(n_Q == 4) {
    currentOptions$color_palette_manual <- c('#d73027', '#fc8d59', '#fee090', '#e0f3f8', '#91bfdb', '#4575b4')
  }
  if(n_Q == 5) {
    currentOptions$color_palette_manual <- c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628')
  }
  if(n_Q == 6) {
    currentOptions$color_palette_manual <- c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf')
  }
  if(n_Q == 7) {
    currentOptions$color_palette_manual <- c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999')
  }
  # currentOptions <- list(xmax = NULL,
  #   heaping_levels = 1:n_Q,
  #   color_palette_brewer  = "Set2")
  
  currentOptions[names(plotOptions)] <- plotOptions
  if (is.null(currentOptions$xmax)) {
    currentOptions$xmax  = 1.1 * mean( (logit(0.9999999) - echanParamElicited$gamma[, n_Q]) / echanParamElicited$gamma[, 1] + echanParamElicited$tau_0)
  }
  
  heaping_levels <- as.character(currentOptions$heaping_levels)
  abs <- 1:currentOptions$xmax
  palette <- currentOptions$color_palette_brewer
  
  Names <-   c('tau[0]', 'tau[1]', 'x[11]', 'x[12]')
  if (n_Q > 3) { Names <- c(Names, paste0('x[', 2:(n_Q - 2), ']')) }
  if (lastCurve) {
    Names <- c(Names, paste0('x[', n_Q, ']'))
  } else{
    Names <- c(Names, paste0('x[', n_Q - 1, ']'))
  }
  row.names(elicitedTable) <- Names
  elicitedTable <- elicitedTable %>% 
    mutate(x = 0.5 * (xsup + xinf))
  
  #-------- Compute Curves-------------------------------- 
  Prestim <- probCondRcppVect(X = abs, param_elicited)
  
  #----- Formatting data
  ### conditional on heaping
  if (M > 1) {
    mean_cond <- apply(Prestim[, , 2:(n_Q + 1)], c(2, 3), mean)
  } else{
    mean_cond <- Prestim[, , 2:(n_Q + 1)]
  }
  w <- which(mean_cond[, 1] > 1 | mean_cond[, 1] < 0) ### remove negative proba and proba above 1
  if (length(w) > 0) { 
    pr_rm <- max(w)
    mean_cond[1:pr_rm, ] <- NA
  } 
  U <- as.data.frame(mean_cond)
  names(U) <- heaping_levels
  
  U <- U %>% 
    mutate(xvalue = abs) %>% 
    gather(key = Q, value = value, -xvalue) %>% 
    mutate(Q = paste('round to ', Q, '|heaping', sep = "")) %>% 
    add_column("type" = rep("cond", length(abs) * n_Q))
  
  ### decision to heap or not
  if (M > 1) {
    mean_heaping <- apply(1 - Prestim[, , 1], 2, mean)
  } else{
    mean_heaping <- 1 - Prestim[, , 1]
  }
  
  V <- data.frame(xvalue = rep(abs, 2))
  V$Q <- rep(c('heaping', 'no heaping'), each = length(abs))
  V$value <- c(mean_heaping, 1 - mean_heaping)
  V$type = rep('heaping', 2 * length(abs))
  
  if (M > 1) {
    max_cond <- apply(Prestim[, , 2:(n_Q + 1)], c(2, 3), function(x) { quantile(x, 0.95) })
    min_cond <- apply(Prestim[, , 2:(n_Q + 1)], c(2, 3), function(x) { quantile(x, 0.05) })
    q3_cond  <- apply(Prestim[, , 2:(n_Q + 1)], c(2, 3), function(x) { quantile(x, 0.75) })
    q1_cond  <- apply(Prestim[, , 2:(n_Q + 1)], c(2, 3), function(x) { quantile(x, 0.25) })
    U$min <- c(min_cond)
    U$max <- c(max_cond)
    U$q3 <- c(q3_cond)
    U$q1 <- c(q1_cond)
    
    min_heaping <- apply(1 - Prestim[, , 1], 2, function(x) { quantile(x, 0.05) })
    max_heaping <- apply(1 - Prestim[, , 1], 2, function(x) { quantile(x, 0.95) })
    q3_heaping  <- apply(1 - Prestim[, , 1], 2, function(x) { quantile(x, 0.75) })
    q1_heaping  <- apply(1 - Prestim[, , 1], 2, function(x) { quantile(x, 0.25) })
    V$min <- c(min_heaping, 1 - min_heaping)
    V$max <- c(max_heaping, 1 - max_heaping)
    V$q3 <- c(q3_heaping, 1 - q3_heaping)
    V$q1 <- c(q1_heaping, 1 - q1_heaping)
  }
  
  W <- rbind(U, V) %>% 
    mutate(Heaping = factor(Q,
                            levels = c("no heaping", "heaping", 
                                       paste("round to ", heaping_levels[1:n_Q], '|heaping', sep = "")
                                       )
                            ),
           type = ifelse(type == "cond", "P(G | A = 1, x)", "P(A| x)")
           )
  
  elicitedTable$type <- c(rep('heaping', 2), rep('cond', N - 2))
  elicitedTable$type <- ifelse(elicitedTable$type == "cond", "P(G | A = 1, x)", "P(A| x)")
  elicitedTable$P[1] <- 1;
  elicitedTable$Q <- c("no heaping",
                       "no heaping",
                       paste("round to ", heaping_levels[1], '|heaping', sep = ""),
                       paste("round to ", heaping_levels[1:(n_Q - 2)], '|heaping', sep = ""),
                       paste("round to ", ifelse(lastCurve, heaping_levels[n_Q], heaping_levels[n_Q - 1]), '|heaping', sep = "")
                       ) 
  elicitedTable$Heaping <- factor(elicitedTable$Q,
                                  levels = c("no heaping",
                                             "heaping",
                                             paste("round to ", heaping_levels[1:n_Q], '|heaping', sep = "")
                                             )
                                  )
  
  # graph part ----
  
  # if first graph with no previous graph
  theme_set(theme_bw(base_size = 14))
  gg <- ggplot() + 
    geom_line(data = W,
              aes(x = xvalue, y = value, color = Heaping, group = Heaping), size = 1
              ) 
  
  if (M > 1) {
    elicitedTable_u <- elicitedTable
    
    if ((n_Q > 3) | (!lastCurve)) {
      u <- c(5:ifelse(lastCurve, n_Q + 1, n_Q + 2))
      elicitedTable_u <- elicitedTable[-u, ]
    }
    
    gg <- gg + 
      geom_ribbon(data = W,
                  aes(ymin = min, ymax = max, x = xvalue, group = Heaping, fill = Heaping), 
                  alpha = 0.1
                  ) +
      geom_ribbon(data = W,
                  aes(ymin = q1, ymax = q3, x = xvalue, group = Heaping, fill = Heaping), 
                  alpha = 0.3
                  ) +
      geom_vline(data = elicitedTable,
                 mapping = aes(xintercept = xinf, group = type, color = Heaping), 
                 linetype = "dashed"
                 ) +
      geom_vline(data = elicitedTable, 
                 mapping = aes(xintercept = xsup, group = type, color = Heaping), 
                 linetype = "dashed"
                 ) +
      geom_hline(data = elicitedTable_u, 
                 mapping = aes(yintercept = P, group = type), 
                 color = 'grey'
                 )
  }
  if (M == 1) {
    if ((N > 5) | (!lastCurve)) {
      upto = ifelse(lastCurve, n_Q + 1, n_Q + 2)
      elicitedTable$P[5:upto] <- probCondRcppVect(X = elicitedTable$x[5:upto], 
                                                  param_elicited
                                                  )[1,,][(2:(n_Q - 2)) + 1]
    }
    
    gg <- gg + 
      geom_vline(data = elicitedTable, 
                 mapping = aes(xintercept = x,group = type,color = Heaping),
                 linetype = "dashed"
                 ) +
      geom_hline(data = elicitedTable,
                 mapping = aes(yintercept = P, group = type), 
                 colour = 'grey'
                 ) +
      geom_point(data = elicitedTable, 
                 mapping = aes(x = x, y =  P, group = type)
                 )
  }
  
  gg <- gg + 
    geom_label_repel(data = elicitedTable %>% 
                       mutate(what = row.names(.)),
                     mapping = aes(x = x, y = P, label = what),
                     parse = TRUE, box.padding = 0.35, 
                     point.padding = 1,
                     segment.color = 'grey50'
                     ) +
    # geom_text(data = elicitedTable %>% 
    #              mutate(what = row.names(.)),
    #           mapping = aes(x = x, y = P, label = what)
    #           ) +
    ylim(0.0, 1.0) + ylab("Probability") +
    scale_x_continuous(name = "x (true count)") + 
    scale_fill_manual(name = "", values = currentOptions$color_palette_manual) +
    scale_color_manual(name = "", values = currentOptions$color_palette_manual) +
    # scale_fill_brewer(name = "Heaping Behaviour: ", type = "qual", palette = palette) +
    # scale_color_brewer(name = "Heaping Behaviour: ", type = "qual", palette = palette) +
    facet_wrap(~ type, ncol = 1) +
    theme(legend.position = "right", 
          plot.title = element_text(lineheight = 0.8, face = "bold"),  
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
          ) 
  return(gg)
  print(gg)
}


