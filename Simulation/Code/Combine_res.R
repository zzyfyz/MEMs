library(tidyr)
library(dplyr)

##################
## Cox
##################


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

lp_u_mu <- read.csv("C:/Users/feiyi/Desktop/MEMs/informative/laplace_unif/1_1_1_1/effm_results.csv")
lp_u_ci <- read.csv("C:/Users/feiyi/Desktop/MEMs/informative/laplace_unif/1_1_1_1/cred_ints_results.csv")
lp_u_prob<- read.csv("C:/Users/feiyi/Desktop/MEMs/informative/laplace_unif/1_1_1_1/prob_results.csv")

# Define the true mean for calculating bias and mse
true_mean <- 0  # Replace with the true mean value if known

# Process mean data
# Process mean data to calculate bias and MSE
effm_1 <- lp_u_mu %>%
  filter(cohort == 1) %>%
  group_by(N) %>%
  summarize(
    mean = mean(Mean),
    bias = mean(Mean - true_mean),
    mse = mean((Mean - true_mean)^2),
    .groups = "drop"
  )
effm_1$prior <- 'Laplace-Uniform'

# Process credible intervals to calculate 95% coverage
effci_1 <- lp_u_ci %>%
  filter(cohort == 1) %>%
  group_by(N) %>%
  summarize(
    cover = mean(true_mean >= CI_value.LowBound & true_mean <= CI_value.UpBound),
    .groups = "drop"
  )
effci_1$prior <- 'Laplace-Uniform'

# Process probability data
prob_1 <- lp_u_prob %>%
  filter(cohort == 1) %>%
  group_by(N) %>%
  summarize(
    prob = mean(ind > 0.975),
    .groups = "drop"
  )
prob_1$prior <- 'Laplace-Uniform'

# Combine all results into a final output
final_output <- effm_1 %>%
  left_join(effci_1 %>% select(N, cover), by = "N") %>%
  left_join(prob_1 %>% select(N, prob), by = "N") %>%
  select(N, mean, bias, mse, prob, cover)

final_output$model <- "Laplace-Uniform"

# Write to CSV
write.csv(final_output, "C:/Users/feiyi/Desktop/MEMs/informative/laplace_unif/laplace_uni_1_1_1_1.csv", row.names = FALSE)

####################
## laplace empirical
####################

lp_e_mu <- read.csv("C:/Users/feiyi/Desktop/MEMs/informative/laplace_empirical/1_1_1_1/effm_results.csv")
lp_e_ci <- read.csv("C:/Users/feiyi/Desktop/MEMs/informative/laplace_empirical/1_1_1_1/cred_ints_results.csv")
lp_e_prob<- read.csv("C:/Users/feiyi/Desktop/MEMs/informative/laplace_empirical/1_1_1_1/prob_results.csv")

# Define the true mean for calculating bias and mse
true_mean <- 0  # Replace with the true mean value if known

# Process mean data
# Process mean data to calculate bias and MSE
effm_1 <- lp_e_mu %>%
  filter(cohort == 1) %>%
  group_by(N) %>%
  summarize(
    mean = mean(Mean),
    bias = mean(Mean - true_mean),
    mse = mean((Mean - true_mean)^2),
    .groups = "drop"
  )
effm_1$prior <- 'Laplace-Empirical'

# Process credible intervals to calculate 95% coverage
effci_1 <- lp_e_ci %>%
  filter(cohort == 1) %>%
  group_by(N) %>%
  summarize(
    cover = mean(true_mean >= CI_value.LowBound & true_mean <= CI_value.UpBound),
    .groups = "drop"
  )
effci_1$prior <- 'Laplace-Empirical'

# Process probability data
prob_1 <- lp_e_prob %>%
  filter(cohort == 1) %>%
  group_by(N) %>%
  summarize(
    prob = mean(ind > 0.975),
    .groups = "drop"
  )
prob_1$prior <- 'Laplace-Empirical'

# Combine all results into a final output
final_output <- effm_1 %>%
  left_join(effci_1 %>% select(N, cover), by = "N") %>%
  left_join(prob_1 %>% select(N, prob), by = "N") %>%
  select(N, mean, bias, mse, prob, cover)

final_output$model <- "Laplace-Empirical"

# Write to CSV
write.csv(final_output, "C:/Users/feiyi/Desktop/MEMs/informative/laplace_empirical/laplace_emp_1_1_1_1.csv", row.names = FALSE)

#############
##PLot
#############



