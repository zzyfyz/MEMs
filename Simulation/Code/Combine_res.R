##################
## Cox
##################
library(dplyr)

file_path <- "C:/Users/feiyi/Desktop/github_MEMs/MEMs/Simulation/Results/Equal/Cox/1_1_1_1/"

cox_mu <- read.csv(paste0(file_path,"effm_results.csv"))
cox_sd <- read.csv(paste0(file_path,"effsd_results.csv"))
cox_prob <- read.csv(paste0(file_path,"prob_results.csv"))

# Define the true value of the parameter for calculating bias and MSE
true_value <- 0  # Modify this based on the true parameter value in your study

# Initialize an empty data frame to store results
results <- data.frame(
  N = integer(),
  Mean = numeric(),
  bias = numeric(),
  mse = numeric(),
  prob = numeric(),
  cover = numeric()
)

# Define the sample sizes
SS <- seq(50, 500, by = 50)

# Loop through each sample size
for (sample_size in SS) {
  # Filter data for the current sample size
  estimators <- cox_mu %>% filter(N == sample_size, cohort==1) %>% pull(Mean)
  sds <- cox_sd %>% filter(N == sample_size, cohort==1) %>% pull(Sd)
  rejections <- cox_prob %>% filter(N == sample_size, cohort==1) %>% pull(ind)
  
  # Calculate metrics
  mean_estimator <- mean(estimators)
  bias <- mean(estimators - true_value)
  mse <- mean((estimators - true_value)^2)
  # Calculate probabilities based on the normal CDF
  p <- pnorm(0, mean = estimators, sd = sds)
  
  # Calculate Type I error rate based on proportion of probs > threshold_prob (0.975)
  reject_rate <- mean(p > 0.975)
  
  # Calculate 95% coverage probability
  lower_bound <- estimators - 1.96 * sds
  upper_bound <- estimators + 1.96 * sds
  coverage <- mean(lower_bound <= true_value & upper_bound >= true_value)
  
  # Append results to data frame
  results <- rbind(results, data.frame(
    N = sample_size,
    Mean = mean_estimator,
    bias = bias,
    mse = mse,
    prob = reject_rate,
    cover = coverage
  ))
}

write.csv(results, file.path(file_path, "cox_null.csv"), row.names = FALSE)


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


