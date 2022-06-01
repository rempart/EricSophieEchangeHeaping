# Load test parameter ----
source("functions/functionsElicitation.R")
source("functions/plotElicited/functionsForPlotElicited/prepareElicited.R")
source("functions/plotElicited/functionsForPlotElicited/graphFunction.R")

# NEW ----
# 
# Q <- 3
# elicited_range = matrix(1, Q + 2, 4)
# colnames(elicited_range) = c('xinf', 'xsup', 'P', 'Conf')
# elicited_range[1,1:4] = c(5, 10, -1, 0.9);
# elicited_range[2,] = c(35, 50, 0.05, 0.95);
# elicited_range[3,] = c(50, 60, 0.95, 0.9)
# elicited_range[4,] = c(80, 100, 0.05, 0.9)
# elicited_range[5,] = c(110, 120, 0.95, 0.9)
# 
# # test add more than Q=3 with : 2 -> hl = 200 between 1000 and 2000, 0.95, 0.9
# # elicited_range[6,] = c(500, 700, 0.95, 0.9)
# # elicited_range[7,] = c(1000, 1200, 0.95, 0.9)
# # elicited_range[8,] = c(1500, 1700, 0.95, 0.9)
# 
# lastCurve <- T
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
# echanParamElicited <- inverseConstraints(elicitedProb, echan_elicited, lastCurve = T)
# 
# 
# listNewElicited <- list(echanParamElicited = echanParamElicited, elicitedTable = elicitedTable)




# # OLD ----
# 
# Q <- 3;
# elicited_range = matrix(1, Q + 2, 4)
# colnames(elicited_range) = c('xinf', 'xsup', 'P', 'Conf')
# elicited_range[1,1:4] = c(5, 15, -1, 0.9);
# elicited_range[2,] = c(20, 60, 0.05, 0.9);
# elicited_range[3,] = c(55, 60, 0.95, 0.9)
# elicited_range[4,] = c(100, 150, 0.05, 0.7)
# elicited_range[5,] = c(190, 270, 0.95, 0.7)
# elicitedTable = as.data.frame(elicited_range)
# Names <-   c('tau_0', 'tau_1', '1_1', '1_2')
# if (Q > 3) {Names <- c(Names, paste(2:(Q - 2)))}
# Names <- c(Names, 'Q')
# row.names(elicitedTable) = Names
# elicitedProb <- elicitedTable$P
# 
# # sampling elicited values with uncertainty
# echan_elicited <- sampleElicitedValue(100, elicitedTable)
# # transform into parameters
# echanParamElicited <- inverseConstraints(elicitedProb, echan_elicited)
# 
# 
# listOldElicited <- list(echanParamElicited = echanParamElicited, elicitedTable = elicitedTable)



plotElicitedPrototype <- function(listNewElicited,
                                  listOldElicited = NULL,
                                  plotOptions = list(
                                    xmax = 250 ,
                                    heaping_levels = c(5, 10, 50),
                                    color_palette_brewer  = 'Set2' # Beware ! this argument can't be NULL
                                  ),
                                  filter = "all",
                                  lastCurve,
                                  Q = q) {
  
  
  # COMMON PARAMETERS ----
  M <- length(listNewElicited$echanParamElicited$tau_0)
  # N <- nrow(listNewElicited$elicitedTable)
  N = Q +2
  n_Q <- Q # <- N -  2
  
  
  # NEW ----
  prepNew <- prepareElicited(
    listElicited = listNewElicited,
    plotOptions = plotOptions,
    M = M,
    N = N,
    n_Q = n_Q,
    Q = Q,
    lastCurve = lastCurve
  )
  
  
  # OLD ----
  if(!is.null(listOldElicited)) {
    
    prepOld <- prepareElicited(
      listElicited = listOldElicited,
      plotOptions = plotOptions,
      M = M,
      N = N,
      n_Q = n_Q,
      Q = Q,
      lastCurve = lastCurve
    )
    
    
    
    final_graph <- subGraphFunction(W_New = prepNew$W,
                                    elicitedTableNew = prepNew$elicitedTable,
                                    W_Old = prepOld$W,
                                    elicitedTableOld = prepOld$elicitedTable,
                                    M = M, N = N, n_Q = n_Q,
                                    filter = filter,
                                    lastCurve = lastCurve,
                                    Q=Q)
    
  } else {
    
    final_graph <- subGraphFunction(W_New = prepNew$W,
                                    elicitedTableNew = prepNew$elicitedTable,
                                    M = M, N = N, n_Q = n_Q,
                                    filter = filter,
                                    lastCurve = lastCurve,
                                    Q=Q)
    
  }
  
  return(final_graph)
}

# test <- plotElicitedPrototype(listNewElicited = listNewElicited,
#                               #listOldElicited = listOldElicited
#                               plotOptions = list(heaping_levels = c(5,10,50),
#                                                  xmax = 250),
#                               lastCurve = T,
#                               Q = 3
#                               )

