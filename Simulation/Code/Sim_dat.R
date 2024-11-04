# Define the function to generate a single dataset
f <- function(n1, n2, h0, G, P, tau) {
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
  
  # Treatment group
  GP1 <- P[1] * G
  t1 <- generate_times(n1, GP1)
  trt1 <- rep(1, n1)
  
  # Placebo group
  GP2 <- P[2] * G
  t2 <- generate_times(n2, GP2)
  trt2 <- rep(0, n2)
  
  # Combine and prepare dataset
  time <- c(t1, t2)
  trt <- c(trt1, trt2)
  censor <- ifelse(is.na(time), 0, 1)
  time <- ifelse(is.na(time), tau[length(tau)], time)
  id <- 1:(n1 + n2)
  data <- data.frame(time = time, trt = trt, censor = censor, id = id)
  data
}

set.seed(123)

# Define tau points and hazard ratios
tau <- c(0, 1, 4.33, 26, 52)
h1 <- 1
h2 <- 1
h3 <- 1
h4 <- 1

# Values of SS to iterate over
sample_sizes <- seq(50, 500, by = 50)

# Base directory
base_dir <- "C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Data/Unequal/1_1_1_1"

# Generate and save datasets for each sample size
for (SS in sample_sizes) {
  # Define the output directory for the current sample size
  output_dir <- file.path(base_dir, as.character(SS))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Generate and save 500 datasets for the current sample size
  for (i in 1:500) {
    PC <- f(SS, SS, exp(-4.01), c(1, 1.01, 1.3, 1.6), c(1, h1), tau)
    PC$cohort <- 1
    S1 <- f(SS/2, SS/2, exp(-4.01), c(1, 1.01, 1.3, 1.6), c(1, h2), tau)
    S1$cohort <- 2
    S2 <- f(SS/2, SS/2, exp(-4.01), c(1, 1.01, 1.3, 1.6), c(1, h3), tau)
    S2$cohort <- 3
    S3 <- f(SS/2, SS/2, exp(-4.01), c(1, 1.01, 1.3, 1.6), c(1, h4), tau)
    S3$cohort <- 4
    
    # Combine cohorts into a list
    dataset <- list(PC = PC, S1 = S1, S2 = S2, S3 = S3)
    
    # Save each list of datasets as a separate .RData file
    save(dataset, file = file.path(output_dir, paste0("dat_", i, ".RData")))
  }
  
}
