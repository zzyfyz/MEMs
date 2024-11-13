
# Define the folder names corresponding to the sample sizes
sample_sizes <- seq(50, 500, by=50)
results <- data.frame(N=integer(), Mean=numeric(), Sd=numeric(),
                      Prob=numeric(), bias=numeric(), mse=numeric(),
                      lcl=numeric(), ucl=numeric(), model = character())

te <- log(0.8) 

# Loop through each sample size folder
for (N in sample_sizes) {
  folder_path <- paste0("C:/Users/feiyi/Desktop/github_MEMs/MEMs/Simulation/Results/Equal/MCMC_emp/0.8_0.8_0.8_0.8/", N)
  
  # Initialize storage for 500 estimates
  posterior_means <- numeric(500)
  posterior_sds <- numeric(500)
  diff <- numeric(500)
  diff_sq <-  numeric(500)
  reject <- numeric(500)
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
    lb_coverage[i] <- data$Posterior_LB <= log(0.8)
    ub_coverage[i] <- data$Posterior_UB >= log(0.8)
    
    diff[i] <- data$Posterior_Mean - te
    
    diff_sq[i] <- (data$Posterior_Mean - te)^2
    
    reject[i] <- ifelse(data$prob > 0.975,1,0)
  }
  
  # Calculate required metrics
  mean_estimate <- mean(posterior_means)
  bias <- mean(diff)
  mse <- mean(diff_sq)
  rej_rate <- mean(reject)
  coverage_probability <- mean(lb_coverage & ub_coverage)
  
  
  # Combine results
  results <- rbind(results, data.frame(N=N, Mean=mean_estimate, bias=bias,
                                       mse=mse, prob = rej_rate, cover = coverage_probability, model="MCMC-Empirical"))
}

# Save the results to a CSV file
write.csv(results, "C:/Users/feiyi/Desktop/github_MEMs/MEMs/Simulation/Results/Equal/MCMC_emp/0.8_0.8_0.8_0.8/mcmc_emp_0.8_0.8_0.8_0.8.csv", row.names=FALSE)