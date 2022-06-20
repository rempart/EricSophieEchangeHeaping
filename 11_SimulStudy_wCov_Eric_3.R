lapply(c("ggplot2", "ggthemes", "tidyverse","rstan","ggmcmc"),
       library, character.only = TRUE)
rm(list = ls())

mu <- 100
beta_cov <- c(-log(2), log(5) / 2)

load("~/Documents/GitHub/DauphinsPapier/EricSophieEchangeHeaping/res_Eric/simul_data_Poisson_n100.RData")
load("~/Documents/GitHub/DauphinsPapier/EricSophieEchangeHeaping/res_Eric/stan_data_Poisson_n100.RData")
StatResume<-tibble()
for(k in (1:100)){
poisson_heap_df<-ggs(poisson_heap[[k]]) %>% 
  mutate(model="Poisson_heap")
poisson_df<-ggs(poisson_null[[k]]) %>% 
  mutate(model="Poisson_null")
negbin_df<-ggs(negbin_null[[k]]) %>% 
  mutate(model="Negbin_null")
d<-rbind(poisson_heap_df,poisson_df,negbin_df) %>% 
group_by(Parameter,model) %>% 
  filter(Parameter %in% c("intercept", "slope[1]","slope[2]")) %>% 
  summarize("m"=mean(value),"s"=sd(value),
  "q05"=quantile(value,probs=c(0.05)),
  "q25"=quantile(value,probs=c(0.25)),
  "q50"=quantile(value,probs=c(0.5)),
  "q75"=quantile(value,probs=c(0.75)),
  "q95"=quantile(value,probs=c(0.95))) %>% 
  mutate(sim=k) 
d$truth=c(rep(log(mu),3),rep(beta_cov[1],3),rep(beta_cov[2],3))

StatResume<-bind_rows(StatResume,d)
}
save(file = "res_Eric/stan_resume_n100.RData",
     list = c("StatResume"))
StatResume %>% ggplot(aes(x=100*abs((m-truth)/truth),color=model))+
  geom_density()+ 
  facet_wrap(~Parameter,scales = "free")

rbind(poisson_heap_df,poisson_df,negbin_df) %>% 
  filter(str_detect(Parameter,"intercept")) %>% 
  ggplot(aes(x=value, col=model))+geom_density()+
  geom_vline(xintercept = log(mu))
rbind(poisson_heap_df,poisson_df,negbin_df) %>% filter(str_detect(Parameter,"slope\\[1\\]")) %>% 
  ggplot(aes(x=value, col=model))+geom_density()+
  geom_vline(xintercept = beta_cov[1])
rbind(poisson_heap_df,poisson_df,negbin_df) %>% filter(str_detect(Parameter,"slope\\[2\\]")) %>% 
  ggplot(aes(x=value, col=model))+geom_density()+
  geom_vline(xintercept = beta_cov[2])

