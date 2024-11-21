#!/usr/bin/env Rscript
library(dplyr)
library(survival)
library(rstan)
library(bridgesampling)


as.matrix(load(list.files(pattern="dat_")))

# Extract datasets
PC <- dataset$PC
S1 <- dataset$S1
S2 <- dataset$S2
S3 <- dataset$S3


# Update the Stan code to handle cohort-specific thetas
stan_code <- "
data {
  int<lower=0> N;                     // Number of observations
  int<lower=0> K;                     // Number of intervals
  int<lower=0> C;                     // Number of cohorts in this model
  int num_subjects[C]; 
  int<lower=1> max_subjects;
  int<lower=0, upper=1> x[C, max_subjects];                        // Covariate for each observation
  array[C, max_subjects, K] real delta;    // Delta matrices by cohort
  array[C, max_subjects] int nu;        // Indicator if survival time is observed (with event)
  int<lower=1, upper=C> cohort[N];    // Cohort identifier for each observation
  matrix[C, K] eta_tilde;             // Precomputed eta_tilde values for each cohort and interval
  matrix[C, K] phi0;                  // Base phi values for each cohort and interval
  array[C, max_subjects, K] real phi_tilde_summand; // Summation component for each cohort, subject, and interval
}

parameters {
  matrix[1, C] theta;                 // Cohort-specific theta effects
}

model {
  // Prior for theta
  theta[1] ~ normal(0, sqrt(0.5));

  // Marginal likelihood using dynamically calculated phi_tilde
  for (c in 1:C) {  // For each cohort
    for (k in 1:K) {  // For each interval
      // Initialize phi_tilde with the base phi0 for this cohort and interval
      real phi_tilde = phi0[c, k];

      // Sum over actual subjects in the cohort, applying the exp(x * theta) weighting
      for (i in 1:num_subjects[c]) {
          phi_tilde += phi_tilde_summand[c, i, k] * exp(x[c][i] * theta[1, c]);
      }

      // Convert delta and nu to vectors for element-wise operations
      vector[num_subjects[c]] delta_k = to_vector(delta[c, 1:num_subjects[c], k]);
      vector[num_subjects[c]] nu_vec = to_vector(nu[c, 1:num_subjects[c]]);
      vector[num_subjects[c]] x_vec = to_vector(x[c, 1:num_subjects[c]]);

      // Marginal likelihood contribution for this cohort and interval
      target += lgamma(eta_tilde[c, k]) - eta_tilde[c, k] * log(phi_tilde) 
                + sum((delta_k .* nu_vec) .* (x_vec * theta[1, c]));
    }
  }
}
"

# Function to calculate marginal likelihood for each model
calculate_marginal_likelihood <- function(model_data, stan_model, C, K, gamma0) {
  # Prepare data for Stan model
  gamma0 <- gamma0
  N <- nrow(model_data)      # Number of observations
  K <- K                        # Number of intervals
  C <- length(unique(model_data$cohort)) # Number of cohorts in this model
  delta <- matrix(0, N, K)      # Event or censor indicators for each interval
  nu <- model_data$censor     # Observed survival indicator (1 = observed with event, 0 = censored)
  
  #int.cuts <- c( 0, quantile( model_data$time[which(model_data$censor == 1)],
  #probs = seq(1/K, 1, length = K)[-K] ) )
  #model_data$interval <- apply( cbind(model_data$time), 1, function(x){ sum(x > int.cuts) } )
  
  tau_p <- c(0, 4, 8, 12, 24)
  intervals <- c(0, 4, 8, 12, 24)
  int.cuts <- tau_p[1:length(tau_p)-1]
  model_data$interval <- apply( cbind(model_data$time), 1, function(x){ sum(x > int.cuts) } )
  
  # Populate delta based on the event/censoring time and intervals
  for (i in 1:N) {
    for (k in 1:K) {
      if (nu[i] == 1 && model_data$time[i] <= intervals[k + 1] && model_data$time[i] > intervals[k]) {
        delta[i, k] <- 1  # Observed event in interval k
      } else if (nu[i] == 0 && model_data$time[i] <= intervals[k + 1] && model_data$time[i] > intervals[k]) {
        delta[i, k] <- 1  # Censored in interval k
      }
    }
  }
  
  # Create lists to hold `delta` and `nu` matrices for each cohort
  delta_list <- vector("list", C)
  nu_list <- vector("list", C)
  
  for (i in 1:C) {
    cohort_indices <- which(model_data$cohort == i)  # Indices for subjects in this cohort
    delta_list[[i]] <- delta[cohort_indices, ]       # Subset delta for this cohort
    nu_list[[i]] <- nu[cohort_indices]               # Subset nu for this cohort
  }
  
  eta0 <- matrix( 0.01, nrow = C, ncol = K )      
  phi0 <- matrix( 0.01, nrow = C, ncol = K ) 
  
  phi.tilde.summand.list <- list()
  for(i in 1:C){
    
    dat.i <- model_data[ which(model_data$cohort == i), ]    # data for only cohort i
    
    # matrix to store summand for all subjects in cohort i for each time interval k, k=1,...,K
    summand.mat.i <- matrix( 0, nrow = nrow(dat.i), ncol = K )
    for(k in 1:K){
      summand.val <- ifelse( dat.i$interval == k, dat.i$time - int.cuts[k], 0 )
      summand.val <- ifelse( dat.i$interval > k, int.cuts[k+1] - int.cuts[k], summand.val )
      summand.mat.i[,k] <- summand.val
    }
    phi.tilde.summand.list[[i]] <- summand.mat.i
    
  }
  
  eta.tilde <- matrix( 0, nrow = C, ncol = K )
  for(i in 1:C){
    for(k in 1:K){
      d.ijk.times.nu.ij <- ifelse( model_data$cohort == i & model_data$censor == 1 &
                                     model_data$interval == k, 1, 0 )
      eta.tilde[i,k] <- sum( d.ijk.times.nu.ij ) + eta0[i,k]
    }
  }
  
  num_subjects <- as.integer(sapply(phi.tilde.summand.list, nrow))
  if (length(num_subjects) == 1) {
    num_subjects <- array(num_subjects, dim = C)  # Explicitly set it as an array if C = 1
  }
  
  # Convert `phi.tilde.summand.list` into a 3D array for Stan, only using the required number of rows per cohort
  phi_tilde_summand <- array(0, dim = c(C, max(num_subjects), K))
  
  for (i in 1:C) {
    phi_tilde_summand[i, 1:num_subjects[i], ] <- phi.tilde.summand.list[[i]]
  }
  
  
  # Create covariate matrix, e.g., using trt and cohort as predictors
  x <- model_data$trt
  
  x_list <- vector("list", C)
  for (i in 1:C) {
    cohort_indices <- which(model_data$cohort == i)
    x_list[[i]] <- x[cohort_indices]
  }
  
  max_subjects <- max(num_subjects)
  
  # Initialize arrays with dimensions based on max_subjects, K, and C
  x_padded <- array(0, dim = c(C, max_subjects))            
  delta_padded <- array(0, dim = c(C, max_subjects, K))     
  nu_padded <- array(0, dim = c(C, max_subjects))           
  
  # Populate padded arrays with cohort-specific data
  for (i in 1:C) {
    x_padded[i, 1:num_subjects[i]] <- x_list[[i]]               
    delta_padded[i, 1:num_subjects[i], ] <- delta_list[[i]]     
    nu_padded[i, 1:num_subjects[i]] <- nu_list[[i]]             
  }
  
  
  # Define Stan data, with cohort indicator
  stan_data <- list(
    N = N,
    K = K,
    C = C,                                      
    num_subjects = num_subjects, 
    max_subjects = max_subjects,
    x = x_padded,                        
    delta = delta_padded,                
    nu = nu_padded,                                           
    cohort = model_data$cohort,
    eta_tilde = eta.tilde,
    phi0 = phi0,
    phi_tilde_summand = phi_tilde_summand
  )
  
  
  init_fn <- function() {
    list(theta=matrix(0, nrow = 1, ncol = C))
  }
  
  # Run the model with example data
  fit <- sampling(stan_model, data = stan_data, init = init_fn, iter = 2000, warmup = 500, chains = 4, refresh = 100, cores = 4)
  
  
  # Calculate posterior mean values for theta
  posterior_mean_theta <- summary(fit, pars = "theta")$summary[, "mean"]
  posterior_sd_theta <- summary(fit, pars = "theta")$summary[, "sd"]
  posterior_lb_theta <- summary(fit, pars = "theta")$summary[, "2.5%"]
  posterior_ub_theta <- summary(fit, pars = "theta")$summary[, "97.5%"]
  
  posterior_samples <- extract(fit, pars = "theta")$theta
  posterior_samples_pc <- posterior_samples[, 1, 1]
  prob_less_gamma0 <- mean(posterior_samples_pc < gamma0)
  
  
  bridge_result <- bridge_sampler(fit)
  log_marginal_likelihood <- bridge_result$logml
  
  return(list(
    mean_theta = posterior_mean_theta,
    sd_theta = posterior_sd_theta,
    lb_theta = posterior_lb_theta,
    ub_theta = posterior_ub_theta,
    log_marginal_likelihood = log_marginal_likelihood,
    samples = posterior_samples_pc,
    prob_less_gamma0 = prob_less_gamma0
  ))
}


# Define the 8 models
models <- list(
  c(1, 1, 1, 1),
  c(1, 1, 1, 2),
  c(1, 2, 1, 1),
  c(1, 1, 2, 1),
  c(1, 1, 2, 3),
  c(1, 2, 1, 3),
  c(1, 2, 3, 1),
  c(1, 2, 3, 4)
)


# Function to combine datasets based on model indicator
combine_datasets <- function(model, datasets) {
  # Initialize a list to hold selected datasets with assigned cohort indicators
  selected_datasets <- list()
  
  # Loop over each element in the model vector to select datasets based on indicators
  for (i in seq_along(model)) {
    # Check if the dataset should be included based on the model indicator
    if (model[i] > 0) {
      # Copy the dataset and add a 'cohort' column based on the model indicator
      dataset_copy <- dataset[[i]]
      dataset_copy$cohort <- ifelse(i == 1, 1, model[i])
      selected_datasets[[length(selected_datasets) + 1]] <- dataset_copy
    }
  }
  
  # Combine all selected datasets by row
  combined_data <- do.call(rbind, selected_datasets)
  return(combined_data)
}

# Apply the function to each model to create combined datasets
combined_datasets <- lapply(models, combine_datasets, datasets = dataset)

# Naming the list for easier identification
names(combined_datasets) <- paste0("Model_", seq_along(models))


stan_model <- stan_model(model_code = stan_code)

# Apply the function to each model to get theta and marginal likelihoods
results <- lapply(1:length(models), function(i) {
  model_data <- combined_datasets[[paste0("Model_", i)]]
  C <- length(unique(model_data$cohort))
  K <- 4  # Number of intervals
  gamma0 <- 0
  calculate_marginal_likelihood(model_data, stan_model, C, K, gamma0)
})

# Extract theta and marginal likelihoods from results
posterior_mean <- sapply(results, function(res) res$mean_theta)
posterior_sd <- sapply(results, function(res) res$sd_theta)
log_marginal_likelihoods <- sapply(results, function(res) res$log_marginal_likelihood)
prob_less_gamma0s <- sapply(results, function(res) res$prob_less_gamma0)
post_sample <- sapply(results, function(res) res$samples)


# Define constraint for optimal source
constraint <- 0.67

# Identify the optimal model(s) based on log marginal likelihoods
max_logml <- max(log_marginal_likelihoods)
optimal_indices <- which(log_marginal_likelihoods == max_logml)

# Initialize a vector to store exchangeability probabilities for each supplemental cohort
num_supp_cohorts <- length(models[[1]]) - 1  # Exclude primary cohort
prior_weights <- rep(0, num_supp_cohorts)    # Start with zero probabilities

# Apply exchangeability probability based on optimal model(s)
for (i in optimal_indices) {
  model_config <- models[[i]][-1]  # Exclude primary cohort (first element)
  prior_weights <- prior_weights | (model_config == 1)  # Set to 1 if included in any optimal model
}

# Multiply by constraint to set the exchangeability probabilities to either 0 or 0.67
prior_weights <- prior_weights * constraint



# Calculate empirical prior for each model configuration based on p1, p2, p3
mod_prior <- sapply(models, function(model) {
  # Determine probability for this configuration
  prob_config <- (ifelse(model[2] == 1, prior_weights[1], 1 - prior_weights[1])) *
    (ifelse(model[3] == 1, prior_weights[2], 1 - prior_weights[2])) *
    (ifelse(model[4] == 1, prior_weights[3], 1 - prior_weights[3]))
  prob_config
})

# Calculate posterior model weights
max_log <- max(log_marginal_likelihoods)
ml_prop <- exp(log_marginal_likelihoods - max_log)
pmp <- (ml_prop * mod_prior)/sum(ml_prop * mod_prior)

est_theta <- sum(sapply(seq_along(posterior_mean), function(i) {
  posterior_mean[[i]][1] * pmp[i]
}))

est_sd <- sum(sapply(seq_along(posterior_sd), function(i) {
  posterior_sd[[i]][1] * pmp[i]
}))

prob <- sum(sapply(seq_along(prob_less_gamma0s), function(i) {
  prob_less_gamma0s[[i]] * pmp[i]
}))

# get the 95% credible interval
n.samp <- 100000  # number of samples to draw
cred.ints <- matrix(0, nrow = 1, ncol = 2)  
colnames(cred.ints) <- c("LowBound", "UpBound")

# Sample from models based on pmp
which.mods <- sample(seq_along(posterior_mean), n.samp, replace = TRUE, prob = pmp)

# Generate samples for global effect
glob.samp <- sapply(1:n.samp, function(i) {
  sample(post_sample[, which.mods[i]], 1)
})

# Compute 95% credible interval
cred.ints[1, ] <- quantile(glob.samp, c(0.025, 0.975))


results_df <- data.frame(
  Posterior_Mean = est_theta,
  Posterior_SD = est_sd,
  Posterior_LB = cred.ints[1, 1],
  Posterior_UB = cred.ints[1, 2],
  prob = prob
)

text <- list.files(pattern="dat_")
num <- unlist(lapply(text, function(x) {
  sub("dat_(\\d+)\\.RData", "\\1", x)
}))
write.csv(results_df, paste0("mcmc.result.",num,".csv"), row.names=FALSE)

