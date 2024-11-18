# Define the function to generate a single dataset
f <- function(n1, n2, h0, G, P, tau, censor_rate) {
  # Set time intervals
  DT <- diff(tau)
  
  # Create a helping matrix
  LD <- matrix(0, nrow=length(tau), ncol=length(G))
  LD[lower.tri(LD)] <- 1
  
  # Function to generate time for each group
  generate_times <- function(n, GP) {
    LS <- log(1 - runif(n))
    LSM <- -h0 * as.vector(LD %*% (GP * DT))
    time <- rep(NA, n)
    for (i in 1:(length(tau) - 1)) {
      time <- ifelse(LSM[i] >= LS & LS > LSM[i + 1], tau[i] + (LSM[i] - LS) / (h0 * GP[i]), time)
    }
    time
  }
  
  # Generate event times for treatment and placebo groups
  GP1 <- P[1] * G
  t1 <- generate_times(n1, GP1)
  trt1 <- rep(1, n1)
  
  GP2 <- P[2] * G
  t2 <- generate_times(n2, GP2)
  trt2 <- rep(0, n2)
  
  # Combine event times and treatment indicators
  time <- c(t1, t2)
  trt <- c(trt1, trt2)
  
  # Generate censoring times
  censoring_time <- rexp(n1 + n2, rate = censor_rate)
  
  # Determine observed times and censoring indicators
  observed_time <- pmin(time, censoring_time)
  censor <- ifelse(time <= censoring_time, 1, 0)
  
  # Handle subjects with time exceeding the last interval
  observed_time <- ifelse(is.na(observed_time), tau[length(tau)], observed_time)
  censor <- ifelse(is.na(censor), 0, censor)
  
  # Assign IDs and prepare dataset
  id <- 1:(n1 + n2)
  data <- data.frame(time = observed_time, trt = trt, censor = censor, id = id)
  data
}

set.seed(123)

# Define tau points and hazard ratios
tau <- c(0, 4, 8, 12, 24)
h1 <- 1
h2 <- 1.25
h3 <- 0.83
h4 <- 1
censor_rate = 0.008


# Values of SS to iterate over
sample_sizes <- seq(50, 500, by = 50)

# Base directory
base_dir <- "C:/Users/feiyi/Desktop/github_MEMs/MEMs/Simulation/Data/Equal/1_0.8_1.2_1"

# Generate and save datasets for each sample size
for (i in 1:500) {
  for (SS in sample_sizes) {
    # Define the output directory for the current sample size
    output_dir <- file.path(base_dir, as.character(SS))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Generate and save 500 datasets for the current sample size
    
    PC <- f(SS, SS, 0.03, c(1, 1.01, 1.3, 1.6), c(1, h1), tau, censor_rate)
    PC$cohort <- 1
    S1 <- f(SS, SS, 0.03, c(1, 1.01, 1.3, 1.6), c(1, h2), tau, censor_rate)
    S1$cohort <- 2
    S2 <- f(SS, SS, 0.03, c(1, 1.01, 1.3, 1.6), c(1, h3), tau, censor_rate)
    S2$cohort <- 3
    S3 <- f(SS, SS, 0.03, c(1, 1.01, 1.3, 1.6), c(1, h4), tau, censor_rate)
    S3$cohort <- 4
    
    # Combine cohorts into a list
    dataset <- list(PC = PC, S1 = S1, S2 = S2, S3 = S3)
    
    # Save each list of datasets as a separate .RData file
    save(dataset, file = file.path(output_dir, paste0("dat_", i, ".RData")))
  }
  
}