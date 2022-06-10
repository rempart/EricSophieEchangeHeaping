##--------------------------------------------------------------------------------------------------------
## SCRIPT : Run simulations with Stan
##
## Authors : Isabelle Albert, Matthieu Authier, Daouda Ba, Sophie Donnet, Mathieu Genu & Eric Parent
##
## Last update : 2021-10-02
## R version 4.0.5 (2021-03-31) -- "Shake and Throw"
## Copyright (C) 2021 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##--------------------------------------------------------------------------------------------------------

lapply(c("ggplot2", "ggthemes", "tidyverse"), library, character.only = TRUE)

#rm(list = ls())

source('functions/functionsHeapingModel.R')
source('functions/functionsHeapingProbabilityRcpp.R')

### load functions for Stan
source('functions/functionsStan.R')
source('functions/functionsStanmodels.r')

##############################################################
#   SIMULATION parameters 
##############################################################
n_observer <- 1
n_obs <- c(200, 500, 1000)
n_tot <- n_observer * n_obs

mu <- 90 # mean
omega <- 1000 # overdispersion
heaping_levels <- c(1, 10, 50, 100)

n_cov <- 2
beta_cov <- c(-log(2), log(5) / 2)

make_Xcov <- function(n_obs, n_observer, n_cov) {
  Xcov <- vector(mode = 'list', length = n_observer)
  for(k in 1:n_observer) {
    Xcov[[k]] <- cbind(rep(1, n_obs), replicate(n_cov, runif(n_obs,-1,1)))
  }
  return(Xcov)
}

#=================== Heaping behaviour
# first observer
param_obs1 = list()
param_obs1$param_distrib <- list(transfo_over_dispersion = log(omega - 1),
                                 beta = c(log(mu - 1), beta_cov) # because shifted distribution
                                 )
param_obs1$param_heaping <- list(lambda_0 = 0.1,
                                 tau_0 = 10,
                                 threshold = 0.05,
                                 gamma = c(0.2, -11.5, -35.0)
                                 )

param_sim <- list(param_obs1) #, param_obs2)#, param_obs3)

#-- plots heaping behavior 

data_heap <- do.call("rbind", lapply(1:length(param_sim),function(l){
  true_count <- 0:300
  Heap_behaviour1<- probHeaping(true_count,param_heaping =param_sim[[l]]$param_heaping,op.output.cond =  FALSE)
  data_heap_1 <- as_tibble(c(Heap_behaviour1))
  names(data_heap_1)<- c('Prob')
  data_heap_1$true_count <- rep(true_count,ncol(Heap_behaviour1))
  data_heap_1$heaping_level <- rep(1:ncol(Heap_behaviour1),each = length(true_count))
  data_heap_1$obs <- l
  data_heap_1$obs <- as.factor(data_heap_1$obs)
  data_heap_1$heaping_level <- as.factor(data_heap_1$heaping_level)
  return(data_heap_1)}
))

ggplot(data_heap,aes(x=true_count,y=Prob,col=heaping_level))+geom_line() + facet_grid(obs~.)


##############################################################
#------------ SIMULATING DATA ------------------------
##############################################################

distrib <- c("Poisson") # only Poisson

#---- for reproducibility
set.seed(20211002)
n_sim <- 10
seeds <- sample.int(1e6, size = n_sim, replace = FALSE)

### generate the data
for(l in rev(distrib)) {
  message(paste("Data-Generating Mechanism:", l, sep = " "))
  for(j in rev(n_obs)) {
    message(paste("\t sample size:", j, sep = " "))
    ### generate data
    all_covar <- all_data <- vector(mode = 'list', length = n_sim)
    for(i in 1:n_sim) {
      set.seed(seeds[i])
      all_covar[[i]] <- make_Xcov(n_obs = j, n_observer = n_observer, n_cov = n_cov)
      temp <- vector(mode = 'list', length = n_observer)
      for(k in 1:n_observer) {
        dd <- simulHeapedDataCovar(covar = all_covar[[i]][[k]],
                                   param = param_sim[[k]],
                                   heaping_levels = heaping_levels,
                                   model = list(distrib = ifelse(l == "Poisson", 'pois', 'nbinom'),
                                                heaping = 'twoStepsRcpp'
                                                )
                                   )
        while(any(is.na(dd$Y))) {
          print("NA appaered")
          all_covar[[i]] <- make_Xcov(n_obs = j, n_observer = n_observer, n_cov = n_cov)
          dd <- simulHeapedDataCovar(covar = all_covar[[i]][[k]],
                                     param = param_sim[[k]],
                                     heaping_levels = heaping_levels,
                                     model = list(distrib = ifelse(l == "Poisson", 'pois', 'nbinom'),
                                                  heaping = 'twoStepsRcpp'
                                                  )
                                     )
        }
        dd$Observer <- k
        covX <- as.data.frame(all_covar[[i]][[k]][, -1])
        names(covX) <- paste("cov", 1:ncol(covX), sep = "")
        temp[[k]] <- cbind(dd, covX); rm(dd, covX)
      }
      all_data[[i]] <- do.call('rbind', temp); rm(temp)
    }
    save(list = c("all_data"), file = paste("res_Sophie/simul_data_", l, "_n", j, ".RData", sep = ""))
  }
}
 
###############################################################"
###    check the simulations
############################################################### 


### heaping percentage
percent <- do.call('rbind',
                   lapply(all_data, function(x) { sapply(1:4, function(G) { sum(x$sim_G == G) } ) / nrow(x) })
                   )
percent <- as.data.frame(percent); names(percent) <- c("no_heaping", paste0('heap_at_',heaping_levels[-1]))
skimr::skim(percent)
### summary statistics
sum_stat <- do.call('rbind', lapply(1:n_sim, function(x) {
  c(min(all_data[[x]]$X),mean(all_data[[x]]$X), sd(all_data[[x]]$X), max(all_data[[x]]$X), 
    min(all_data[[x]]$Y), mean(all_data[[x]]$Y), sd(all_data[[x]]$Y), max(all_data[[x]]$Y)
    ) 
  }
)
)
sum_stat <- as.data.frame(sum_stat); names(sum_stat) <- c("min_X","mean_X", "sd_X", "max_X", "min_Y","mean_Y", "sd_Y", "max_Y")
skimr::skim(sum_stat)

# ------------------------

model = 'Poisson'
load(paste0('res_Sophie/simul_data_',model,'_n', max(n_obs),'.RData'))
data_heap <- do.call("rbind", lapply(1:n_sim,function(l){all_data[[l]]$sim=l; return(all_data[[l]])}))
data_heap <- data_heap %>%mutate(Observer =as.factor(Observer)) %>% mutate(Small = as.factor(ifelse(Y<500,'Y < 500','Y >= 500')))
data_heap  %>% ggplot(aes(x=Y,color= Observer)) + geom_histogram(fill="white")+ facet_grid(Small ~ Observer,scales = 'free')

#-------------------  Nombre d'obs par sim_level
data_heap %>% mutate(sim_G = as.factor(sim_G))%>% count(sim_G) %>%mutate(freq = n / sum(n)*100)

#-------------------- 
seuil <- 200
data_heap %>% ggplot(aes(x=X,y=Y,color=Observer)) + geom_point() + scale_x_continuous(trans='log') + scale_y_continuous(trans='log') 
data_heap %>% filter(X>=seuil) %>% ggplot(aes(x=X,y=Y,color=as.factor(sim_G))) + geom_point() + scale_x_continuous(trans='log') + scale_y_continuous(trans='log') + geom_abline(intercept =0, slope=1)
data_heap %>% filter(X <seuil) %>% ggplot(aes(x=X,y=Y,color=as.factor(sim_G))) + geom_point()  + geom_abline(intercept =0, slope=1)
data_heap %>% mutate(RMSE_X_Y= sqrt((X-Y)^2/X^2)*100) %>% filter(sim_G>1) %>% ggplot(aes(x=X,y=RMSE_X_Y,color=as.factor(sim_G))) + geom_point() + facet_grid(.~as.factor(sim_G),scales="free")
data_heap %>% mutate(RMSE_X_Y= sqrt((X-Y)^2/X^2)*100) %>% filter(sim_G == 4) %>%filter(X<100)%>% ggplot(aes(x=X,y=RMSE_X_Y,color=as.factor(sim_G))) + geom_point() + facet_grid(.~as.factor(sim_G),scales="free")
data_heap%>% filter(X>=seuil) %>% ggplot(aes(x=X)) + geom_histogram() 
data_heap%>% filter(X>=seuil)  %>% ggplot(aes(x=Y)) + geom_histogram() 



###############################################################
### bayesian estimation (4 models to fit)
###############################################################"
distrib = 'Poisson'
n_chain = 4
n_iter = 5000
n_warm = 1000
n_thin = 5
hmc_control = list(adapt_delta = 0.9)



get_bayes <- function(all_data,n,distrib,n_chain = 4,n_iter=2000,n_warm=1000,n_thin=2,hmc_control=list(adapt_delta = 0.9)) 
  {
  poisson_null <- negbin_null <- vector(mode = 'list', length = length(all_data))
  poisson_heap <- vector(mode = 'list', length = length(all_data))
  #poisson_mheap <-  negbin_heap<- negbin_mheap <- vector(mode = 'list', length = length(all_data))
  ### fit the models
  for(i in 1:length(all_data)) {
    ### null model
    poisson_null[[i]] <- sampling(stan_nullpoisson_cov,
                                  data = list(n_obs = length(all_data[[i]]$Y),
                                              n_y = as.numeric(all_data[[i]]$Y),
                                              Y = as.numeric(all_data[[i]]$Y),
                                              prior_location_intercept = log(10),
                                              prior_scale_intercept = log(5) / 2,
                                              n_cov = 2,
                                              prior_location_slope = rep(0, 2),
                                              prior_scale_slope = rep(log(2) / 2, 2),
                                              Z = as.matrix(all_data[[i]][, grep("cov", names(all_data[[i]]))])
                                              ),
                                  pars = c("intercept", "slope", "log_lik"),
                                  chains = n_chain,
                                  iter = n_iter,
                                  warmup = n_warm,
                                  thin = n_thin,
                                  control = hmc_control
                                  )
    negbin_null[[i]] <- sampling(stan_nullnegbin_cov,
                                 data = list(n_obs = length(all_data[[i]]$Y),
                                             n_y = as.numeric(all_data[[i]]$Y),
                                             Y = as.numeric(all_data[[i]]$Y),
                                             prior_location_intercept = log(10),
                                             prior_scale_intercept = log(5) / 2,
                                             n_cov = 2,
                                             prior_location_slope = rep(0, 2),
                                             prior_scale_slope = rep(log(2), 2),
                                             Z = as.matrix(all_data[[i]][, grep("cov", names(all_data[[i]]))])
                                             ),
                                 pars = c("intercept", "slope", "overdispersion", "log_lik"),
                                 chains = n_chain,
                                 iter = n_iter,
                                 warmup = n_warm,
                                 thin = n_thin,
                                 control = hmc_control
                                 )
    ### one observer
    standatalist <- data4stan(countdata = all_data[[i]]$Y,
                              covar = cbind(rep(1, nrow(all_data[[i]])), as.matrix(all_data[[i]][, grep("cov", names(all_data[[i]]))])),
                              heaping_levels = heaping_levels,
                              prior_location_intercept = log(10),
                              prior_scale_intercept = log(5) / 2,
                              prior_location = c(log(10), log(0.5), log(0.2), 0, log(10)),
                              prior_cholmat = diag(c(log(3) ,log(5) /2, log(5) / 2, 5.0, log(5) / 2)) %*%
                                diag(5),
                              threshold = 0.05
                              )

    poisson_heap[[i]] <- sampling(stan_poisson_cov,
                                  data = standatalist,
                                  pars = c("intercept", "slope", "gamma", "tau_0", "lambda_0", "sup_heaping", "log_lik"),
                                  chains = n_chain,
                                  iter = n_iter,
                                  warmup = n_warm,
                                  thin = n_thin,
                                  control = hmc_control
                                  )

    
    
  }
  save(file = paste("res_Sophie/stan_data_", distrib, "_n", n, ".RData", sep = ""),
       list = c("poisson_null", "negbin_null", "poisson_heap"))
  #, "poisson_mheap"
  #              )
  #     )
}

for(l in distrib) {
  message(paste("Simulation model:", l, sep = " "))
  for(j in c(max(n_obs))) {
    message(paste("\t sample size:", j, sep = " "))
    load(paste("res_Sophie/simul_data_", l, "_n", j, ".RData", sep = ""))
    get_bayes(all_data = all_data, n = j, distrib = l)
  }
}; rm(i, j, l)



###############################################################
#-------------------------------------------------------------
###############################################################
fit_ss <- extract(poisson_heap[[1]], permuted = TRUE)





