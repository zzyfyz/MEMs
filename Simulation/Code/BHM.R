# Load required libraries
library(dplyr)
library(survival)
library(ggplot2)
library(ggpubr)
library(rstan)

# Set up parameters and progress bar
set.seed(123)
nsim <- 500
gamma0 <- 0
xind <- seq(50, 500, by = 50) # Assuming sample_sizes is defined and contains sample sizes
base_dir <- "C:/Users/feiyi/Desktop/github_MEMs/MEMs/Simulation/Data/Equal/1_0.8_0.7_0.6"
pb <- txtProgressBar(0, nsim, style = 3)

results <- data.frame(
  sample_size = rep(xind, each = nsim),
  posterior_mean_theta = NA,
  posterior_sd_theta = NA,
  posterior_lb_theta = NA,
  posterior_ub_theta = NA,
  prob_less_gamma0 = NA
)

# Loop over simulations and sample sizes
for (i in 1:nsim) {
  for (j in 1:length(xind)) {
    
    data_file <- file.path(base_dir, as.character(xind[j]), paste0("dat_", i, ".RData"))
    load(data_file) 
    
    # Extract datasets
    PC <- dataset$PC
    S1 <- dataset$S1
    S2 <- dataset$S2
    S3 <- dataset$S3
    
    dat <- rbind(PC, S1, S2, S3)
    dat1 <- rbind(PC, S1)
    dat1$cohort <- ifelse(dat1$cohort!=1, 2, 1)
    dat2 <- rbind(PC, S2)
    dat2$cohort <- ifelse(dat2$cohort!=1, 2, 1)
    dat3 <- rbind(PC, S3)
    dat3$cohort <- ifelse(dat3$cohort!=1, 2, 1)
    
    mean_1 <- log(sum(PC$censor[which(PC$trt==1)])/
                                 sum(PC$time[which(PC$trt==1)])/
                                 (sum(PC$censor[which(PC$trt==0)])/
                                    sum(PC$time[which(PC$trt==0)])))
    mean_2 <- log(sum(S1$censor[which(S1$trt==1)])/
                                 sum(S1$time[which(S1$trt==1)])/
                                 (sum(S1$censor[which(S1$trt==0)])/
                                    sum(S1$time[which(S1$trt==0)])))
    mean_3 <- log(sum(S2$censor[which(S2$trt==1)])/
                                 sum(S2$time[which(S2$trt==1)])/
                                 (sum(S2$censor[which(S2$trt==0)])/
                                    sum(S2$time[which(S2$trt==0)])))
    mean_4 <- log(sum(S3$censor[which(S3$trt==1)])/
                                 sum(S3$time[which(S3$trt==1)])/
                                 (sum(S3$censor[which(S3$trt==0)])/
                                    sum(S3$time[which(S3$trt==0)])))
    
    sd_1 <- 1/(sum(PC$censor[PC$trt==0])) + 1/(sum(PC$censor[PC$trt==1]))
    sd_2 <- 1/(sum(PC$censor[S1$trt==0])) + 1/(sum(S1$censor[S1$trt==1]))
    sd_3 <- 1/(sum(PC$censor[S2$trt==0])) + 1/(sum(S2$censor[S2$trt==1]))
    sd_4 <- 1/(sum(PC$censor[S3$trt==0])) + 1/(sum(S3$censor[S3$trt==1]))
    
    dat_BHM <- data.frame(
      cohort = 1:4,
      mean = c(mean_1, mean_2, mean_3, mean_4),
      sd = c(sd_1, sd_2, sd_3, sd_4)
    )
    
    stan_data <- list(
      cohort = nrow(dat_BHM),           # Number of cohorts/regions
      y = dat_BHM$mean,                     # Means from each cohort
      s2i = dat_BHM$sd                    # Variances (sd^2) for each cohort
    )
    
    stan_code <- "
data {
  int<lower=1> cohort;     
  real y[cohort];         
  real s2i[cohort];       
}

parameters {
  real mu;                     // Overall mean
  real<lower=0> tau2;       
  real mui[cohort];        // Cohort-level means
}

model {
  // Likelihood
  for (i in 1:cohort) {
    y[i] ~ normal(mui[i], s2i[i]);
  }
  
  // Prior distributions
  mu ~ normal(0, 10);
  tau2 ~ inv_gamma(0.001, 0.001);  
  for (i in 1:cohort) {
    mui[i] ~ normal(mu, sqrt(tau2));
  }
}

"
stan_model <- stan_model(model_code = stan_code)

# Compile the model
fit <- sampling(stan_model, data = stan_data, iter = 1000, warmup = 500, chains = 4, cores=4)

# Calculate posterior mean values for theta
posterior_mean_theta <- summary(fit, pars = "mui")$summary[1, "mean"]
posterior_sd_theta <- summary(fit, pars = "mui")$summary[1, "sd"]
posterior_lb_theta <- summary(fit, pars = "mui")$summary[1, "2.5%"]
posterior_ub_theta <- summary(fit, pars = "mui")$summary[1, "97.5%"]

posterior_samples <- extract(fit)$mui[, 1]
prob_less_gamma0 <- mean(posterior_samples < gamma0)

# Store results
results$posterior_mean_theta[results$sample_size == xind[j] & is.na(results$posterior_mean_theta)][1] <- posterior_mean_theta
results$posterior_sd_theta[results$sample_size == xind[j] & is.na(results$posterior_sd_theta)][1] <- posterior_sd_theta
results$posterior_lb_theta[results$sample_size == xind[j] & is.na(results$posterior_lb_theta)][1] <- posterior_lb_theta
results$posterior_ub_theta[results$sample_size == xind[j] & is.na(results$posterior_ub_theta)][1] <- posterior_ub_theta
results$prob_less_gamma0[results$sample_size == xind[j] & is.na(results$prob_less_gamma0)][1] <- prob_less_gamma0

  }
setTxtProgressBar(pb,i)
}

summary_results <- results %>%
  group_by(sample_size) %>%
  summarize(
    mean_posterior_mean_theta = mean(posterior_mean_theta, na.rm = TRUE),
    bias = mean(posterior_mean_theta - 0, na.rm = TRUE),
    mse = mean((posterior_mean_theta - 0)^2, na.rm = TRUE),
    coverage_prob = mean((posterior_lb_theta <= 0 & posterior_ub_theta >= 0), na.rm = TRUE),
    rej_rate = mean(ifelse(prob_less_gamma0 > 0.975, 1, 0), na.rm = TRUE)
  )
