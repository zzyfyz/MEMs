##################
## Cox
##################
cox_mu <- read.csv("C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Cox/1_1_1_1/effm_results.csv")
cox_sd <- read.csv("C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Cox/1_1_1_1/effsd_results.csv")
cox_prob <- read.csv("C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Cox/1_1_1_1/prob_results.csv")

effm_0 <- cox_mu %>% group_by(cohort,N)%>% summarize(Mean=mean(Mean))%>%
  filter(cohort == 1)  
effm_0$prior <- 'CPHM'
effsd_0 <- cox_sd %>% group_by(cohort,N)%>% summarize(Sd=mean(Sd))%>%
  filter(cohort == 1)  
effsd_0$prior <- 'CPHM'
prob_0 <- cox_prob %>% group_by(cohort,N)%>% summarize(Prob=mean(ind))%>%
  filter(cohort == 1)  
prob_0$prior <- 'CPHM'

##################
## laplace uniform
##################

lp_u_mu <- read.csv("C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Laplace/Uniform/1_1_1_1/effm_results.csv")
lp_u_sd <- read.csv("C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Laplace/Uniform/1_1_1_1/effsd_results.csv")
lp_u_prob<- read.csv("C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Laplace/Uniform/1_1_1_1/prob_results.csv")

effm_1 <- lp_u_mu %>%
  group_by(cohort, N) %>%
  summarize(Mean = mean(Mean), .groups = "drop") %>%
  filter(cohort == 1)  
effm_1$prior <- 'Laplace-Uniform'

effsd_1 <- lp_u_sd %>% group_by(cohort, N) %>%
  summarize(Sd=mean(Sd), .groups = "drop") %>%
  filter(cohort == 1)  
effsd_1$prior <- 'Laplace-Uniform'

prob_1 <- lp_u_prob %>% group_by(cohort,N)%>% group_by(cohort, N) %>%
  summarize(Prob=mean(ind), .groups = "drop") %>%
  filter(cohort == 1)  
prob_1$prior <- 'Laplace-Uniform'

####################
## laplace empirical
####################

lp_e_mu <- read.csv("C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Laplace/Empirical/1_1_1_1/effm_results.csv")
lp_e_sd <- read.csv("C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Laplace/Empirical/1_1_1_1/effsd_results.csv")
lp_e_prob<- read.csv("C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Laplace/Empirical/1_1_1_1/prob_results.csv")

effm_2 <- lp_e_mu %>%
  group_by(cohort, N) %>%
  summarize(Mean = mean(Mean), .groups = "drop") %>%
  filter(cohort == 1)  
effm_2$prior <- 'Laplace-Empirical'

effsd_2 <- lp_e_sd %>% group_by(cohort, N) %>%
  summarize(Sd=mean(Sd), .groups = "drop") %>%
  filter(cohort == 1)  
effsd_2$prior <- 'Laplace-Empirical'

prob_2 <- lp_e_prob %>% group_by(cohort,N)%>% group_by(cohort, N) %>%
  summarize(Prob=mean(ind), .groups = "drop") %>%
  filter(cohort == 1)  
prob_2$prior <- 'Laplace-Empirical'

####################
## BHM
####################


