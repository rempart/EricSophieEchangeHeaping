library(magrittr)
library(gridExtra)
library(RColorBrewer)

# debugage
# Load test parameter ----
# source("functions/functionsElicitation.R")
# source("shiny/R/plotElicited/functionsForPlotElicited/prepareElicited.R")
# source("shiny/R/plotElicited/functionsForPlotElicited/graphFunction.R")

# NEW ----
# 
# Q <- 3
# elicited_range = matrix(1, Q + 2, 4)
# colnames(elicited_range) = c('xinf', 'xsup', 'P', 'Conf')
# elicited_range[1,1:4] = c(5, 10, -1, 0.9);
# elicited_range[2,] = c(35, 50, 0.05, 0.95);
# elicited_range[3,] = c(50, 60, 0.95, 0.9)
# elicited_range[4,] = c(80, 100, 0.05, 0.9)
# elicited_range[5,] = c(200, 250, -1, 0.9)
# 
# # test add more than Q=3 with : 2 -> hl = 200 between 1000 and 2000, 0.95, 0.9
# # elicited_range[6,] = c(500, 700, 0.95, 0.9)
# # elicited_range[7,] = c(1000, 1200, 0.95, 0.9)
# # elicited_range[8,] = c(1500, 1700, 0.95, 0.9)
# 
# lastCurve <- F
# 
# elicitedTable = as.data.frame(elicited_range)
# Names <-   c('tau_0','tau_1','1_1','1_2')
# if (Q > 3) {Names <- c(Names, paste(2:(Q - 2)))}
# 
# if (lastCurve) {
#   Names <- c(Names, Q)
# }else{
#   Names <- c(Names, Q - 1)
# }
# # Names <- c(Names, "Q-1")
# row.names(elicitedTable) = Names
# elicitedProb <- elicitedTable$P
# 
# # sampling elicited values with uncertainty
# echan_elicited <- sampleElicitedValue(100, elicitedTable)
# # transform into parameters
# echanParamElicited <- inverseConstraints(elicitedProb, echan_elicited, lastCurve = F)
# 
# 
# listNewElicited <- list(echanParamElicited = echanParamElicited, elicitedTable = elicitedTable)
# 
# 
# listOldElicited = NULL
# plotOptions = list(
#   xmax = 250 ,
#   heaping_levels = c(5, 10, 50),
#   color_palette_brewer  = 'Set2' # Beware ! this argument can't be NULL
# )
# filter = "all"
# lastCurve = FALSE
# 
# 
# # COMMON PARAMETERS ----
# M <- length(listNewElicited$echanParamElicited$tau_0)
# # N <- nrow(listNewElicited$elicitedTable) 
# N = Q +2
# n_Q <- Q #<- N -  2
# 
# 
# # NEW ----
# prepNew <- prepareElicited(
#   listElicited = listNewElicited,
#   plotOptions = plotOptions,
#   M = M,
#   N = N,
#   n_Q = n_Q,
#   Q = Q,
#   lastCurve = lastCurve
# )
# 
# 
# 
# listElicited = listOldElicited
# plotOptions = plotOptions
# M = M
# N = N
# n_Q = n_Q
# Q = Q
# lastCurve = lastCurve
# 
# 
# W_Old = NULL
# elicitedTableOld = NULL
# filter = "all"
# 
# 
# # # test
# W_New <- prepNew$W
# # W_Old <- prepOld$W
# elicitedTableNew <- prepNew$elicitedTable
# elicitedTableOld <- prepOld$elicitedTable

subGraphFunction <- function(W_New, W_Old = NULL,
                             elicitedTableNew, elicitedTableOld = NULL,
                             M, N, n_Q,
                             filter = "all",
                             lastCurve,
                             Q=Q) {
  
  all_NULL <- is.null(W_Old) & is.null(elicitedTableOld)
  
  W_New_temp <- W_New
  elicitedTableNew_temp <- elicitedTableNew
  
  if(!all_NULL){
    W_Old_temp <- W_Old
    elicitedTableOld_temp <- elicitedTableOld
    colOld <- "grey70"
  }
  
  # filter
  if("all" %in% filter){
    selected_col <- unique(W_New$Heaping)
  } else {
    selected_col <- filter
  }
  
  col_filter <- brewer.pal(length(selected_col),"Set2")
  
  gg_list <- list()
  
  for (g in 1:length(selected_col)) {
    
    palette <- col_filter[g]
    
    # filter
    W_New <- W_New_temp %>% 
      mutate(Q = factor(Q, levels = c('no heaping', 'heaping', as.factor(1:n_Q)))) %>%
      filter(Heaping == selected_col[g])
    
    elicitedTableNew <- elicitedTableNew_temp %>% 
      filter(Heaping == selected_col[g])
    
    if(!all_NULL){
      W_Old <- W_Old_temp %>% 
        mutate(Q = factor(Q, levels = c('no heaping', 'heaping', as.character(1:n_Q)))) %>%
        filter(Heaping == selected_col[g])
      
      elicitedTableNew <- elicitedTableNew_temp %>% 
        filter(Heaping == selected_col[g])
    }
    
    
    # begin graph
    if(!all_NULL){
      gg <- W_Old %>% 
        ggplot() +
        geom_line(data = W_Old, aes(x = xvalue, y = value, group = Heaping), color = colOld, size = 1) +
        geom_line(data = W_New, aes(x = xvalue, y = value, group = Heaping), color = palette, size = 1) +
        labs(title = unique(W_Old$Heaping)) 
        
      
    } else {
      gg <- W_New %>% 
        ggplot() +
        geom_line(aes(x = xvalue, y = value, group = Heaping), color = palette, size = 1) +
        labs(title = unique(W_New$Heaping))
    }
    
    
    # condition on M & N
    if (M > 1) {
      elicitedTableNew_u <- elicitedTableNew
      
      if ((Q > 3) | (!lastCurve)) {
        upto = ifelse(lastCurve,Q + 1,Q + 2)
        u <- c(5:(upto));
        # elicitedTableNew_u <- elicitedTableNew[-u,]
      } 
      
      if(!all_NULL){
        gg <- gg +
          geom_ribbon(data = W_Old, aes(ymin = min, ymax = max, x = xvalue, group = Heaping),
                      alpha = 0.1, fill = colOld) +
          geom_ribbon(data = W_Old, aes(ymin = q1, ymax = q3, x = xvalue, group = Heaping), 
                      alpha = 0.3, fill = colOld) +
          geom_ribbon(data = W_New, aes(ymin = min, ymax = max, x = xvalue, group = Heaping), 
                      fill = palette, alpha = 0.1) +
          geom_ribbon(data = W_New, aes(ymin = q1, ymax = q3, x = xvalue, group = Heaping), 
                      fill = palette, alpha = 0.3) +
          geom_vline(data = elicitedTableNew,
                     mapping = aes(xintercept = xinf, group = type), 
                     linetype = "dashed", color = palette) +
          geom_vline(data = elicitedTableNew, 
                     mapping = aes(xintercept = xsup, group = type), 
                     linetype = "dashed", color = palette) +
          geom_hline(data = elicitedTableNew_u, mapping = aes(yintercept = P, group = type), color = 'grey')
      } else {
        gg <- gg + 
          geom_ribbon(aes(ymin = min, ymax = max, x = xvalue, group = Heaping), 
                      fill = palette, alpha = 0.1) +
          geom_ribbon(aes(ymin = q1, ymax = q3, x = xvalue, group = Heaping), 
                      fill = palette, alpha = 0.3) +
          geom_vline(data = elicitedTableNew,
                     mapping = aes(xintercept = xinf, group = type), 
                     linetype = "dashed", color = palette) +
          geom_vline(data = elicitedTableNew, 
                     mapping = aes(xintercept = xsup, group = type), 
                     linetype = "dashed", color = palette) +
          geom_hline(data = elicitedTableNew_u, mapping = aes(yintercept = P, group = type), color = 'grey')
        
      }
      
      
      
      
    }
    if (M == 1) {
      if ((N > 5) | (!lastCurve)) {
        upto = ifelse(lastCurve,Q + 1,Q + 2)
        elicitedTableNew$P[5:upto] = probCondRcppVect(X = elicitedTableNew$x[5:upto], param_elicited)[1,,][(2:(Q - 2)) + 1]
      }
      
      gg <- gg + 
        geom_vline(data = elicitedTableNew, 
                   mapping = aes(xintercept = x,group = type),
                   linetype = "dashed", color = palette) +
        geom_hline(data = elicitedTableNew,  mapping = aes(yintercept = P,group = type),colour = 'grey') +
        geom_point(data = elicitedTableNew, mapping = aes(x = x, y = P, group = type))
    }
    
    # 2nd part of graph
    gg <- gg + 
      ylim(0.0, 1.0) + 
      # scale_x_continuous(name = "True count",
      #                    breaks = seq(0, max(W_New$xvalue), by = 20)) +
      ylab("Probability") +
      # scale_fill_manual(name = "Heaping Behaviour: ", palette = palette) +
      # scale_color_manual(name = "Heaping Behaviour: ", palette = palette) +
      theme(legend.position = "top", 
            plot.title = element_text(lineheight = 0.8, face = "bold"), 
            axis.text = element_text(size = 10),
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    
    gg_list[[g]] <- gg
  }
  
  gg_arranged <- do.call(grid.arrange, gg_list)
  
  
  return(all_graph = gg_arranged)
}


