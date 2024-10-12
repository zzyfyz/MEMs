library(dplyr)
library(survival)
library(ggplot2)
library(ggpubr)
library(changepoint)

set.seed(123)

result_pmp <- array(NA_real_, dim=2,
                    dimnames=list("model"=1:2))

result_prior <- array(NA_real_, dim=2,
                      dimnames=list("model"=1:2))

result_ss<- array(NA_real_, dim=1,)

result_effm <- array(NA_real_, dim=2,
                     dimnames=list("cohort"=1:2))

result_effsd <- array(NA_real_, dim=2,
                      dimnames=list("cohort"=1:2))

result_effci <- array(NA_real_, dim=c(2,2),
                      dimnames=list("cohort"=1:2,
                                    "CI"=1:2))

result_prob <- array(NA_real_, dim=2,
                     dimnames=list("cohort"=1:2))

tmp <- read.csv("./combine.csv")
dat <- subset(tmp,cohort==1|cohort==5)
dat$cohort <- ifelse(dat$cohort!=1,2,1)
prim <- subset(dat,cohort==1)

fit <- survfit(Surv(time, censor) ~ 1, data = prim)
cumhaz <- -log(fit$surv)
timepoints <- fit$time
cp <- cpt.mean(cumhaz, method = "PELT")
change_points <- cpts(cp)
change_points


mod <- main(2,rbind(c(1,1),c(1,2)),2,c(0, 24, 90),cuttype="primary",priortype=1,W=NA,alpha0=0)

result_pmp[1] <- mod$pmp[1]
result_pmp[2] <- mod$pmp[2]

result_ss <- mod$ss

result_prior[1] <- mod$mod_prior[1]
result_prior[2] <- mod$mod_prior[2]

result_effm[1] <- exp(mod$mean[1])
result_effm[2] <- exp(mod$mean[2])

result_effsd[1] <- mod$sd[1]
result_effsd[2] <- mod$sd[2]

result_effci[1,1] <- exp(mod$mean[1]+1.96*mod$sd[1])
result_effci[1,2] <- exp(mod$mean[1]-1.96*mod$sd[1])
result_effci[2,1] <- exp(mod$mean[2]+1.96*mod$sd[2])
result_effci[2,2] <- exp(mod$mean[2]-1.96*mod$sd[2])

result_prob[1] <- mod$test[1]
result_prob[2] <- mod$test[2]



pmp <- as.data.frame.table(result_pmp, responseName = "percentage", stringsAsFactors = F)  
ss <- as.data.frame.table(result_ss, responseName = "N", stringsAsFactors = F)
prior <- as.data.frame.table(result_prior, responseName = "Weight", stringsAsFactors = F) 
effm <- as.data.frame.table(result_effm, responseName = "Mean", stringsAsFactors = F)
effsd <- as.data.frame.table(result_effsd, responseName = "sd", stringsAsFactors = F)
effci <- as.data.frame.table(result_effci, responseName = "CI", stringsAsFactors = F) 
prob <- as.data.frame.table(result_prob, responseName = "Prob", stringsAsFactors = F) 
prob$ind <- ifelse(prob$Prob > 0.95, 1, 0)


write.csv(effm, "./effm_2.csv", row.names=FALSE)
write.csv(pmp, "./pmp_2.csv", row.names=FALSE)
write.csv(ss, "./ss_2.csv", row.names=FALSE)
write.csv(effci, "./ci_2.csv", row.names=FALSE)

