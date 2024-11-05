# Load required libraries
library(dplyr)
library(survival)
library(ggplot2)
library(ggpubr)

# Set up parameters and progress bar
set.seed(123)
nsim <- 500

xind <- seq(50, 500, by = 50) # Assuming sample_sizes is defined and contains sample sizes
base_dir <- "C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Data/Equal/0.8_0.7_0.6_0.5"
pb <- txtProgressBar(0, nsim, style = 3)

# Initialize result arrays
result_effm <- array(NA_real_, dim = c(nsim, 4, length(xind)),
                     dimnames = list("simulation" = 1:nsim, "cohort" = 1:4, "N" = xind))

result_effsd <- array(NA_real_, dim = c(nsim, 4, length(xind)),
                      dimnames = list("simulation" = 1:nsim, "cohort" = 1:4, "N" = xind))

result_prob <- array(NA_real_, dim = c(nsim, 4, length(xind)),
                     dimnames = list("simulation" = 1:nsim, "cohort" = 1:4, "N" = xind))

# Loop over simulations and sample sizes
for (i in 1:nsim) {
  for (j in 1:length(xind)) {
    
    data_file <- file.path(base_dir, as.character(xind[j]), paste0("dat_", i, ".RData"))
    load(data_file)  # This loads the `dataset` list with elements PC, S1, S2, S3
    
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
    
    cox <- coxph(Surv(time, censor) ~ trt, data = PC)
    cox1 <- coxph(Surv(time, censor) ~ trt, data = S1)
    cox2 <- coxph(Surv(time, censor) ~ trt, data = S2)
    cox3 <- coxph(Surv(time, censor) ~ trt, data = S3)
    
    result_effm[i,1,j] <- summary(cox)$coeff[1]
    result_effm[i,2,j] <- summary(cox1)$coeff[1]
    result_effm[i,3,j] <- summary(cox2)$coeff[1]
    result_effm[i,4,j] <- summary(cox3)$coeff[1]
    
    result_effsd[i,1,j] <- summary(cox)$coeff[3]
    result_effsd[i,2,j] <- summary(cox1)$coeff[3]
    result_effsd[i,3,j] <- summary(cox2)$coeff[3]
    result_effsd[i,4,j] <- summary(cox3)$coeff[3]
    
    result_prob[i,1,j] <- summary(cox)$coeff[5]
    result_prob[i,2,j] <- summary(cox1)$coeff[5]
    result_prob[i,3,j] <- summary(cox2)$coeff[5]
    result_prob[i,4,j] <- summary(cox3)$coeff[5]
    
  }
  
  setTxtProgressBar(pb,i)
  
}

effm <- as.data.frame.table(result_effm, responseName = "Mean", stringsAsFactors = F)
effsd <- as.data.frame.table(result_effsd, responseName = "Sd", stringsAsFactors = F)
prob <- as.data.frame.table(result_prob, responseName = "Prob", stringsAsFactors = F)
prob$ind <- ifelse(prob$Prob < 0.05, 1, 0)

output_dir <- "C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Cox/0.8_0.7_0.6_0.5"
write.csv(effm, file.path(output_dir, "effm_results.csv"), row.names = FALSE)
write.csv(effsd, file.path(output_dir, "effsd_results.csv"), row.names = FALSE)
write.csv(prob, file.path(output_dir, "prob_results.csv"), row.names = FALSE)



