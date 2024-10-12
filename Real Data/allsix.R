library(dplyr)
library(survival)
library(ggplot2)
library(ggpubr)
library(changepoint)

options(scipen = 999)

result_pmp <- array(NA_real_, dim=32,
                    dimnames=list("model"=1:32))

result_prior <- array(NA_real_, dim=32,
                      dimnames=list("model"=1:32))

result_ss<- array(NA_real_, dim=1,)

result_effm <- array(NA_real_, dim=6,
                     dimnames=list("cohort"=1:6))

result_effsd <- array(NA_real_, dim=6,
                      dimnames=list("cohort"=1:6))

result_effci <- array(NA_real_, dim=c(6,2),
                      dimnames=list("cohort"=1:6,
                                    "CI"=1:2))

result_prob <- array(NA_real_, dim=6,
                     dimnames=list("cohort"=1:6))

dat <- read.csv("./combine.csv")
prim <- subset(dat,cohort==1)

fit <- survfit(Surv(time, censor) ~ 1, data = prim)
cumhaz <- -log(fit$surv)
timepoints <- fit$time
cp <- cpt.mean(cumhaz, method = "PELT")
change_points <- cpts(cp)
change_points

mod <- main(6,rbind(
  c(1, 1, 1, 1, 1, 1),
  c(1, 2, 1, 1, 1, 1),
  c(1, 1, 2, 1, 1, 1),
  c(1, 1, 1, 2, 1, 1),
  c(1, 1, 1, 1, 2, 1),
  c(1, 1, 1, 1, 1, 2),
  c(1, 2, 3, 1, 1, 1),
  c(1, 2, 1, 3, 1, 1),
  c(1, 2, 1, 1, 3, 1),
  c(1, 2, 1, 1, 1, 3),
  c(1, 1, 2, 3, 1, 1),
  c(1, 1, 2, 1, 3, 1),
  c(1, 1, 2, 1, 1, 3),
  c(1, 1, 1, 2, 3, 1),
  c(1, 1, 1, 2, 1, 3),
  c(1, 1, 1, 1, 2, 3),
  c(1, 2, 3, 4, 1, 1),
  c(1, 2, 3, 1, 4, 1),
  c(1, 2, 3, 1, 1, 4),
  c(1, 2, 1, 3, 4, 1),
  c(1, 2, 1, 3, 1, 4),
  c(1, 2, 1, 1, 3, 4),
  c(1, 1, 2, 3, 4, 1),
  c(1, 1, 2, 3, 1, 4),
  c(1, 1, 2, 1, 3, 4),
  c(1, 1, 1, 2, 3, 4),
  c(1, 2, 3, 4, 5, 1),
  c(1, 2, 1, 3, 4, 5),
  c(1, 2, 3, 1, 4, 5),
  c(1, 2, 3, 4, 1, 5),
  c(1, 1, 2, 3, 4, 5),
  c(1, 2, 3, 4, 5, 6)
),2,c(0, 24, 90),cuttype="primary",priortype=1,W=NA,alpha0=0)


for (j in 1:32){
  result_pmp[j] <- mod$pmp[j]
}

result_ss <- mod$ss

for (j in 1:32){
  result_prior[j] <- mod$mod_prior[j]
}

for (j in 1:6){
  result_effm[j] <- exp(mod$mean[j])
}

for (j in 1:6){
  result_effsd[j] <- mod$sd[j]
}

for (j in 1:6){
  result_effci[j,1] <- exp(mod$mean[j]+1.96*mod$sd[j])
  result_effci[j,2] <- exp(mod$mean[j]-1.96*mod$sd[j])
}

for (j in 1:6){
  result_prob[j] <- exp(mod$test[j])
}


pmp <- as.data.frame.table(result_pmp, responseName = "percentage", stringsAsFactors = F)  
ss <- as.data.frame.table(result_ss, responseName = "N", stringsAsFactors = F)
prior <- as.data.frame.table(result_prior, responseName = "Weight", stringsAsFactors = F) 
effm <- as.data.frame.table(result_effm, responseName = "Mean", stringsAsFactors = F)
effsd <- as.data.frame.table(result_effsd, responseName = "Sd", stringsAsFactors = F) 
effci <- as.data.frame.table(result_effci, responseName = "CI", stringsAsFactors = F)
prob <- as.data.frame.table(result_prob, responseName = "Prob", stringsAsFactors = F) 
prob$ind <- ifelse(prob$Prob > 0.95, 1, 0)


write.csv(effm, "./effm_6.csv", row.names=FALSE)
write.csv(pmp, "./pmp_6.csv", row.names=FALSE)
write.csv(ss, "./ss_6.csv", row.names=FALSE)
write.csv(effci, "./ci_6.csv", row.names=FALSE)
