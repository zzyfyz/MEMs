
# Define the folder names corresponding to the sample sizes
sample_sizes <- seq(50, 500, by=50)
results <- data.frame(N=integer(), Mean=numeric(), Sd=numeric(),
                      Prob=numeric(), bias=numeric(), mse=numeric(),
                      lcl=numeric(), ucl=numeric(), model = character())

# Loop through each sample size folder
for (N in 500) {
  folder_path <- paste0("C:/Users/feiyi/Desktop/github_MEMs/MEMs/Simulation/Results/Equal/MCMC/1_0.8_0.7_0.6/", N)
  
  # Initialize storage for 500 estimates
  posterior_means <- numeric(500)
  posterior_sds <- numeric(500)
  lb_coverage <- logical(500)
  ub_coverage <- logical(500)
  
  # Process each of the 500 files in the folder
  for (i in 1:500) {
    file_path <- file.path(folder_path, paste0("mcmc.result.", i, ".csv"))
    data <- read.csv(file_path)
    
    # Store the Posterior_Mean and Posterior_SD
    posterior_means[i] <- data$Posterior_Mean
    posterior_sds[i] <- data$Posterior_SD
    
    # Check if 0 is within the posterior interval
    lb_coverage[i] <- data$Posterior_LB <= 0
    ub_coverage[i] <- data$Posterior_UB >= 0
  }
  
  # Calculate required metrics
  mean_estimate <- mean(posterior_means)
  sd_estimate <- sqrt(mean(posterior_sds^2))
  bias <- mean_estimate - 0
  mse <- bias^2 + sd_estimate^2
  coverage_probability <- mean(lb_coverage & ub_coverage)
  
  # Calculate 95% credible interval for the mean of posterior means
  ci_lower <- quantile(posterior_means, 0.025)
  ci_upper <- quantile(posterior_means, 0.975)
  
  # Combine results
  results <- rbind(results, data.frame(N=N, Mean=mean_estimate, Sd=sd_estimate,
                                       Prob=coverage_probability, bias=bias,
                                       mse=mse, lcl=ci_lower, ucl=ci_upper, model="MCMC-Uniform"))
}

# Save the results to a CSV file
write.csv(results, "summary_results.csv", row.names=FALSE)