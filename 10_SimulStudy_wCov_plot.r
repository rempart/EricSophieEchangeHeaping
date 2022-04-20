rm(list=ls())
library(ggplot2)
library(tidyverse)

n_obs <- c(100,200,400)
my_res <- c()
for (n in n_obs){
  
  name_file <- paste0('res/results_Poisson_',n,'_null_poisson.txt')
  res_n <- read.table(name_file,header =  TRUE)
  my_res <- rbind(my_res,res_n )
  
  
  name_file <- paste0('res/results_Poisson_',n,'_heap_poisson.txt')
  res_n <- read.table(name_file,header =  TRUE)
  my_res <- rbind(my_res,res_n )
  
  name_file <- paste0('res/results_Poisson_',n,'_mheap_poisson.txt')
  res_n <- read.table(name_file,header =  TRUE)
  my_res <- rbind(my_res,res_n )
  
  
}
my_res%>% group_by(param,model) %>% summarise(v=var(estimate))

lev=levels(factor(my_res$param))
lev[1:4]=c('mu','beta_1','beta_2','looic')

param_true <- data_frame(my_res[1:(3*4),]) %>% rename(value=estimate)
param_true$value<-rep(c(20,-log(2), log(5) / 2,NA),3)
param_true$model=rep(unique(my_res$model),each=4) 


  

#### plot param /  modele 
my_res %>% filter(param %in% c('mu','beta_1','beta_2','looic') )%>% 
  mutate(param=factor(param,levels=lev) )%>% 
  ggplot(aes(x=as.factor(sample_size*2),y=estimate)) +  
  geom_boxplot() + geom_hline(data=param_true %>% 
                                mutate(param=factor(param,levels=lev) ),
                              aes(yintercept = value),col='red')+
facet_grid(param~model,scales="free_y")

# To do:  DONE by ERIC

#### plot zomm sur LOOIC mheap / heap 

my_res %>% filter(param %in% c('looic')) %>% filter(model %in% c('poisson_heap', 'poisson_mheap'))%>% ggplot(aes(x=as.factor(sample_size*2),y=estimate)) +  geom_boxplot() + facet_grid(param~model,scales="free_y")




#