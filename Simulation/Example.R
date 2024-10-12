library(dplyr)
library(survival)
library(ggplot2)
library(ggpubr)

set.seed(123)
nsim <- 500
xind <- seq(50,500,len=10)
pb <- txtProgressBar(0,nsim,style=3)

result_pmp <- array(NA_real_, dim=c(nsim,8,length(xind)),
                    dimnames=list("simulation"=1:nsim,
                                  "model"=1:8,
                                  "N"=xind))

result_pmp1 <- array(NA_real_, dim=c(nsim,8,length(xind)),
                     dimnames=list("simulation"=1:nsim,
                                   "model"=1:8,
                                   "N"=xind))

result_pmp2 <- array(NA_real_, dim=c(nsim,8,length(xind)),
                     dimnames=list("simulation"=1:nsim,
                                   "model"=1:8,
                                   "N"=xind))

result_pmp3 <- array(NA_real_, dim=c(nsim,8,length(xind)),
                     dimnames=list("simulation"=1:nsim,
                                   "model"=1:8,
                                   "N"=xind))

result_prior <- array(NA_real_, dim=c(nsim,8,length(xind)),
                      dimnames=list("simulation"=1:nsim,
                                    "model"=1:8,
                                    "N"=xind))

result_prior1 <- array(NA_real_, dim=c(nsim,8,length(xind)),
                       dimnames=list("simulation"=1:nsim,
                                     "model"=1:8,
                                     "N"=xind))

result_prior2 <- array(NA_real_, dim=c(nsim,8,length(xind)),
                       dimnames=list("simulation"=1:nsim,
                                     "model"=1:8,
                                     "N"=xind))

result_prior3 <- array(NA_real_, dim=c(nsim,8,length(xind)),
                       dimnames=list("simulation"=1:nsim,
                                     "model"=1:8,
                                     "N"=xind))

result_ss<- array(NA_real_, dim=c(nsim,length(xind)),
                  dimnames=list("simulation"=1:nsim,
                                "N"=xind))

result_ss1<- array(NA_real_, dim=c(nsim,length(xind)),
                   dimnames=list("simulation"=1:nsim,
                                 "N"=xind))

result_ss2<- array(NA_real_, dim=c(nsim,length(xind)),
                   dimnames=list("simulation"=1:nsim,
                                 "N"=xind))

result_ss3<- array(NA_real_, dim=c(nsim,length(xind)),
                   dimnames=list("simulation"=1:nsim,
                                 "N"=xind))

result_effm <- array(NA_real_, dim=c(nsim,4,length(xind)),
                     dimnames=list("simulation"=1:nsim,
                                   "cohort"=1:4,
                                   "N"=xind))

result_effm1 <- array(NA_real_, dim=c(nsim,4,length(xind)),
                      dimnames=list("simulation"=1:nsim,
                                    "cohort"=1:4,
                                    "N"=xind))

result_effm2 <- array(NA_real_, dim=c(nsim,4,length(xind)),
                      dimnames=list("simulation"=1:nsim,
                                    "cohort"=1:4,
                                    "N"=xind))

result_effm3 <- array(NA_real_, dim=c(nsim,4,length(xind)),
                      dimnames=list("simulation"=1:nsim,
                                    "cohort"=1:4,
                                    "N"=xind))

result_effsd <- array(NA_real_, dim=c(nsim,4,length(xind)),
                      dimnames=list("simulation"=1:nsim,
                                    "cohort"=1:4,
                                    "N"=xind))

result_effsd1 <- array(NA_real_, dim=c(nsim,4,length(xind)),
                       dimnames=list("simulation"=1:nsim,
                                     "cohort"=1:4,
                                     "N"=xind))

result_effsd2 <- array(NA_real_, dim=c(nsim,4,length(xind)),
                       dimnames=list("simulation"=1:nsim,
                                     "cohort"=1:4,
                                     "N"=xind))

result_effsd3 <- array(NA_real_, dim=c(nsim,4,length(xind)),
                       dimnames=list("simulation"=1:nsim,
                                     "cohort"=1:4,
                                     "N"=xind))

result_prob <- array(NA_real_, dim=c(nsim,4,length(xind)),
                     dimnames=list("simulation"=1:nsim,
                                   "cohort"=1:4,
                                   "N"=xind))

result_prob1 <- array(NA_real_, dim=c(nsim,4,length(xind)),
                      dimnames=list("simulation"=1:nsim,
                                    "cohort"=1:4,
                                    "N"=xind))

result_prob2 <- array(NA_real_, dim=c(nsim,4,length(xind)),
                      dimnames=list("simulation"=1:nsim,
                                    "cohort"=1:4,
                                    "N"=xind))

result_prob3 <- array(NA_real_, dim=c(nsim,4,length(xind)),
                      dimnames=list("simulation"=1:nsim,
                                    "cohort"=1:4,
                                    "N"=xind))

for(i in 1:nsim){
  for(j in 1:length(xind)){
    
    tau_p <- c(0, 1, 4.33, 26, 52)
    tau_s1 <- c(0, 1, 4.33, 26, 52)
    tau_s2 <- c(0, 1, 4.33, 26, 52)
    tau_s3 <- c(0, 1, 4.33, 26, 52)
    
    PC <- f(n1=xind[j], n2=xind[j], h0=exp(-4.01), G=c(1, 1.01, 1.3, 1.6), P=c(1, 1), tau=tau_p)
    table(PC$trt,PC$censor)
    PC$cohort <- 1
    S1 <- f(n1=xind[j]/2, n2=xind[j]/2, h0=exp(-4.01), G=c(1, 1.01, 1.3, 1.6), P=c(1, 1), tau=tau_s1)
    S1$cohort <- 2
    S2 <- f(n1=xind[j]/2, n2=xind[j]/2, h0=exp(-4.01), G=c(1, 1.01, 1.3, 1.6), P=c(1, 1), tau=tau_s2)
    S2$cohort <- 3
    S3 <- f(n1=xind[j]/2, n2=xind[j]/2, h0=exp(-4.01), G=c(1, 1.01, 1.3, 1.6), P=c(1, 1), tau=tau_s3)
    S3$cohort <- 4
    
    dat <- rbind(PC, S1, S2, S3)
    dat1 <- rbind(PC, S1)
    dat1$cohort <- ifelse(dat1$cohort!=1, 2, 1)
    dat2 <- rbind(PC, S2)
    dat2$cohort <- ifelse(dat2$cohort!=1, 2, 1)
    dat3 <- rbind(PC, S3)
    dat3$cohort <- ifelse(dat3$cohort!=1, 2, 1)
    
    #test <- coxph(Surv(time,censor) ~ trt,data=S1)
    #summary(test)
    
    ##get empirical prior
    
    p_1 <- prior_post(dat1,2,rbind(c(1,1),c(1,2)),4,c(0, 1, 4.33, 26, 52),cuttype="primary")
    p_2 <- prior_post(dat2,2,rbind(c(1,1),c(1,2)),4,c(0, 1, 4.33, 26, 52),cuttype="primary")
    p_3 <- prior_post(dat3,2,rbind(c(1,1),c(1,2)),4,c(0, 1, 4.33, 26, 52),cuttype="primary")
    
    mod <- main(4,rbind(c(1,1,1,1),c(1,2,1,1),c(1,1,1,2),c(1,1,2,1),c(1,2,3,1),c(1,2,1,3),c(1,1,2,3),c(1,2,3,4)),
                4,c(0, 1, 4.33, 26, 52),cuttype="primary",priortype=1,W=NA,alpha0=3)
    
    result_pmp[i,1,j] <- mod$pmp[1]
    result_pmp[i,2,j] <- mod$pmp[2]
    result_pmp[i,3,j] <- mod$pmp[3]
    result_pmp[i,4,j] <- mod$pmp[4]
    result_pmp[i,5,j] <- mod$pmp[5]
    result_pmp[i,6,j] <- mod$pmp[6]
    result_pmp[i,7,j] <- mod$pmp[7]
    result_pmp[i,8,j] <- mod$pmp[8]
    
    result_ss[i,j] <- mod$ss
    
    result_prior[i,1,j] <- mod$mod_prior[1]
    result_prior[i,2,j] <- mod$mod_prior[2]
    result_prior[i,3,j] <- mod$mod_prior[3]
    result_prior[i,4,j] <- mod$mod_prior[4]
    result_prior[i,5,j] <- mod$mod_prior[5]
    result_prior[i,6,j] <- mod$mod_prior[6]
    result_prior[i,7,j] <- mod$mod_prior[7]
    result_prior[i,8,j] <- mod$mod_prior[8]
    
    result_effm[i,1,j] <- mod$mean[1]
    result_effm[i,2,j] <- mod$mean[2]
    result_effm[i,3,j] <- mod$mean[3]
    result_effm[i,4,j] <- mod$mean[4]
    
    result_effsd[i,1,j] <- mod$sd[1]
    result_effsd[i,2,j] <- mod$sd[2]
    result_effsd[i,3,j] <- mod$sd[3]
    result_effsd[i,4,j] <- mod$sd[4]
    
    result_prob[i,1,j] <- mod$test[1]
    result_prob[i,2,j] <- mod$test[2]
    result_prob[i,3,j] <- mod$test[3]
    result_prob[i,4,j] <- mod$test[4]
    
    
    
    mod1 <- main(4,rbind(c(1,1,1,1),c(1,2,1,1),c(1,1,1,2),c(1,1,2,1),c(1,2,3,1),c(1,2,1,3),c(1,1,2,3),c(1,2,3,4)),
                 4,c(0, 1, 4.33, 26, 52),cuttype="primary",priortype=3,W=NA,alpha0=NA)
    
    result_pmp1[i,1,j] <- mod1$pmp[1]
    result_pmp1[i,2,j] <- mod1$pmp[2]
    result_pmp1[i,3,j] <- mod1$pmp[3]
    result_pmp1[i,4,j] <- mod1$pmp[4]
    result_pmp1[i,5,j] <- mod1$pmp[5]
    result_pmp1[i,6,j] <- mod1$pmp[6]
    result_pmp1[i,7,j] <- mod1$pmp[7]
    result_pmp1[i,8,j] <- mod1$pmp[8]
    
    result_ss2[i,j] <- mod1$ss
    
    result_prior1[i,1,j] <- mod1$mod_prior[1]
    result_prior1[i,2,j] <- mod1$mod_prior[2]
    result_prior1[i,3,j] <- mod1$mod_prior[3]
    result_prior1[i,4,j] <- mod1$mod_prior[4]
    result_prior1[i,5,j] <- mod1$mod_prior[5]
    result_prior1[i,6,j] <- mod1$mod_prior[6]
    result_prior1[i,7,j] <- mod1$mod_prior[7]
    result_prior1[i,8,j] <- mod1$mod_prior[8]
    
    result_effm1[i,1,j] <- mod1$mean[1]
    result_effm1[i,2,j] <- mod1$mean[2]
    result_effm1[i,3,j] <- mod1$mean[3]
    result_effm1[i,4,j] <- mod1$mean[4]
    
    result_effsd1[i,1,j] <- mod1$sd[1]
    result_effsd1[i,2,j] <- mod1$sd[2]
    result_effsd1[i,3,j] <- mod1$sd[3]
    result_effsd1[i,4,j] <- mod1$sd[4]
    
    result_prob1[i,1,j] <- mod1$test[1]
    result_prob1[i,2,j] <- mod1$test[2]
    result_prob1[i,3,j] <- mod1$test[3]
    result_prob1[i,4,j] <- mod1$test[4]
    
    
    ## Cox model
    cox <- coxph(Surv(time, censor) ~ trt, data = PC)
    cox1 <- coxph(Surv(time, censor) ~ trt, data = S1)
    cox2 <- coxph(Surv(time, censor) ~ trt, data = S2)
    cox3 <- coxph(Surv(time, censor) ~ trt, data = S3)
    
    result_effm2[i,1,j] <- summary(cox)$coeff[1]
    result_effm2[i,2,j] <- summary(cox1)$coeff[1]
    result_effm2[i,3,j] <- summary(cox2)$coeff[1]
    result_effm2[i,4,j] <- summary(cox3)$coeff[1]
    
    result_effsd2[i,1,j] <- summary(cox)$coeff[3]
    result_effsd2[i,2,j] <- summary(cox1)$coeff[3]
    result_effsd2[i,3,j] <- summary(cox2)$coeff[3]
    result_effsd2[i,4,j] <- summary(cox3)$coeff[3]
    
    result_prob2[i,1,j] <- summary(cox)$coeff[5]
    result_prob2[i,2,j] <- summary(cox1)$coeff[5]
    result_prob2[i,3,j] <- summary(cox2)$coeff[5]
    result_prob2[i,4,j] <- summary(cox3)$coeff[5]
    
    ## BHM
    result_effm3[i,1,j] <- log(sum(PC$censor[which(PC$trt==1)])/
                                 sum(PC$time[which(PC$trt==1)])/
                                 (sum(PC$censor[which(PC$trt==0)])/
                                    sum(PC$time[which(PC$trt==0)])))
    result_effm3[i,2,j] <- log(sum(S1$censor[which(S1$trt==1)])/
                                 sum(S1$time[which(S1$trt==1)])/
                                 (sum(S1$censor[which(S1$trt==0)])/
                                    sum(S1$time[which(S1$trt==0)])))
    result_effm3[i,3,j] <- log(sum(S2$censor[which(S2$trt==1)])/
                                 sum(S2$time[which(S2$trt==1)])/
                                 (sum(S2$censor[which(S2$trt==0)])/
                                    sum(S2$time[which(S2$trt==0)])))
    result_effm3[i,4,j] <- log(sum(S3$censor[which(S3$trt==1)])/
                                 sum(S3$time[which(S3$trt==1)])/
                                 (sum(S3$censor[which(S3$trt==0)])/
                                    sum(S3$time[which(S3$trt==0)])))
    
    result_effsd3[i,1,j] <- 1/(sum(PC$censor[PC$trt==0])) + 1/(sum(PC$censor[PC$trt==1]))
    result_effsd3[i,2,j] <- 1/(sum(PC$censor[S1$trt==0])) + 1/(sum(S1$censor[S1$trt==1]))
    result_effsd3[i,3,j] <- 1/(sum(PC$censor[S2$trt==0])) + 1/(sum(S2$censor[S2$trt==1]))
    result_effsd3[i,4,j] <- 1/(sum(PC$censor[S3$trt==0])) + 1/(sum(S3$censor[S3$trt==1]))
  }
  
  setTxtProgressBar(pb,i)
  
}

pmp <- as.data.frame.table(result_pmp, responseName = "percentage", stringsAsFactors = F)  
ss <- as.data.frame.table(result_ss, responseName = "N", stringsAsFactors = F)
prior <- as.data.frame.table(result_prior, responseName = "Weight", stringsAsFactors = F) 
effm <- as.data.frame.table(result_effm, responseName = "Mean", stringsAsFactors = F)
effsd <- as.data.frame.table(result_effsd, responseName = "Sd", stringsAsFactors = F) 
prob <- as.data.frame.table(result_prob, responseName = "Prob", stringsAsFactors = F) 
prob$ind <- ifelse(prob$Prob > 0.95, 1, 0)

pmp1 <- as.data.frame.table(result_pmp1, responseName = "percentage", stringsAsFactors = F)  
ss1 <- as.data.frame.table(result_ss1, responseName = "N", stringsAsFactors = F)
prior1 <- as.data.frame.table(result_prior1, responseName = "Weight", stringsAsFactors = F)
effm1 <- as.data.frame.table(result_effm1, responseName = "Mean", stringsAsFactors = F)
effsd1 <- as.data.frame.table(result_effsd1, responseName = "Sd", stringsAsFactors = F) 
prob1 <- as.data.frame.table(result_prob1, responseName = "Prob", stringsAsFactors = F)
prob1$ind <- ifelse(prob1$Prob > 0.95, 1, 0)

effm2 <- as.data.frame.table(result_effm2, responseName = "Mean", stringsAsFactors = F)
effsd2 <- as.data.frame.table(result_effsd2, responseName = "Sd", stringsAsFactors = F)
prob2 <- as.data.frame.table(result_prob2, responseName = "Prob", stringsAsFactors = F)
prob2$ind <- ifelse(prob2$Prob < 0.05, 1, 0)

effm3 <- as.data.frame.table(result_effm3, responseName = "Mean", stringsAsFactors = F)
effsd3 <- as.data.frame.table(result_effsd3, responseName = "Sd", stringsAsFactors = F)
eff3_all <- merge(effm3,effsd3,by=c('simulation','cohort','N'))

write.csv(eff3_all, "./cox_unequal_ss_1.csv", row.names=FALSE)
bhm <- read.csv("./bhm_unequal_ss_1.csv")
bhm <- subset(bhm, Parameter=='mui_1')
bhm$p <- pnorm(0, bhm$Mean, bhm$StdDev)
bhm$ind <- ifelse(bhm$p > 0.95, 1, 0)

#pmp <- as.data.frame.table(result_pmp, responseName = "percentage", stringsAsFactors = F)
pmp_1 <- pmp %>% group_by(model,N)%>% summarize(percentage=mean(percentage))
ss_1 <- ss %>% group_by(N)%>% summarize(SS=mean(N.1))
ss_1$prior <- 'MEMs-Exponential'
prior_1 <- prior %>% group_by(model,N)%>% summarize(Weight=mean(Weight))
effm_1 <- effm %>% group_by(cohort,N)%>% summarize(Mean=mean(Mean))
effm_1$prior <- 'MEMs-Exponential'
effsd_1 <- effsd %>% group_by(cohort,N)%>% summarize(Sd=mean(Sd))
effsd_1$prior <- 'MEMs-Exponential'
prob_1 <- prob %>% group_by(cohort,N)%>% summarize(Prob=mean(ind))
prob_1$prior <- 'MEMs-Exponential'

pmp_2 <- pmp1 %>% group_by(model,N)%>% summarize(percentage=mean(percentage))
ss_2 <- ss1 %>% group_by(N)%>% summarize(SS=mean(N.1))
ss_2$prior <- 'Empirical Beta-Binomial'
prior_2 <- prior1 %>% group_by(model,N)%>% summarize(Weight=mean(Weight))
effm_2 <- effm1 %>% group_by(cohort,N)%>% summarize(Mean=mean(Mean))
effm_2$prior <- 'MEMs-Empirical'
effsd_2 <- effsd1 %>% group_by(cohort,N)%>% summarize(Sd=mean(Sd))
effsd_2$prior <- 'MEMs-Empirical'
prob_2 <- prob1 %>% group_by(cohort,N)%>% summarize(Prob=mean(ind))
prob_2$prior <- 'MEMs-Empirical'

effm_3 <- effm2 %>% group_by(cohort,N)%>% summarize(Mean=mean(Mean))
effm_3$prior <- 'CPHM'
effsd_3 <- effsd2 %>% group_by(cohort,N)%>% summarize(Sd=mean(Sd))
effsd_3$prior <- 'CPHM'
prob_3 <- prob2 %>% group_by(cohort,N)%>% summarize(Prob=mean(ind))
prob_3$prior <- 'CPHM'

effm_4 <- bhm %>% group_by(Parameter,N)%>% summarize(Mean=mean(Mean))
effm_4$N <- as.character(effm_4$N)
effm_4$prior <- 'BHM'
effm_4$cohort <- '1'
effm_4 <- subset(effm_4,select=c("cohort","N","Mean","prior"))
effsd_4 <- bhm %>% group_by(Parameter,N)%>% summarize(Sd=mean(StdDev))
effsd_4$cohort <- '1'
effsd_4$N <- as.character(effsd_4$N)
effsd_4$prior <- 'BHM'
effsd_4 <- subset(effsd_4,select=c("cohort","N","Sd","prior"))
prob_4 <- bhm %>% group_by(Parameter,N)%>% summarize(Prob=mean(ind))
prob_4$cohort <- '1'
prob_4$N <- as.character(effsd_4$N)
prob_4$prior <- 'BHM'
prob_4 <- subset(prob_4,select=c("cohort","N","Prob","prior"))

effm_all <- rbind(effm_1,effm_2,effm_3,effm_4)
effsd_all <- rbind(effsd_1,effsd_2,effsd_3,effsd_4)
prob_all <- rbind(prob_1,prob_2,prob_3,prob_4)
eff_all <- merge(effm_all,effsd_all,by=c("cohort","N","prior"))
eff_all <- merge(eff_all,prob_all,by=c("cohort","N","prior"))
eff_all <- subset(eff_all,cohort==1)
eff_all$bias <- eff_all$Mean - log(1)
eff_all$mse <- (eff_all$Mean - log(1))^2
eff_all$lcl <-  eff_all$Mean - 1.96*eff_all$Sd
eff_all$ucl <-  eff_all$Mean + 1.96*eff_all$Sd

options(repr.plot.width =3, repr.plot.height =3) 

bias_plot <- eff_all %>% mutate(N=as.numeric(N)) %>%  
  group_by(prior) %>%
  ggplot(aes(N,bias,group=prior,color=prior)) + geom_line(aes(linetype=prior, color=prior)) + 
  scale_y_continuous(breaks = seq(-1,1,by=0.01)) + theme(text = element_text(size = 9),element_line(size =0.2))+
  labs(color  = "Model", linetype = "Model")
p_bias <- bias_plot + labs(title="hr0=1, hr_s1=1, hr_s2=1, hr_s3=1",
                           x ="Sample Size", y = "Bias") +
  theme(plot.title = element_text(hjust = 0.5))

mse_plot <- eff_all %>% mutate(N=as.numeric(N)) %>%  
  group_by(prior) %>%
  ggplot(aes(N,mse,group=prior,color=prior)) + geom_line(aes(linetype=prior, color=prior)) + 
  scale_y_continuous(breaks = seq(0,1,by=0.0005)) + theme(text = element_text(size = 9),element_line(size =0.2))+
  labs(color  = "Model", linetype = "Model")
p_mse <- mse_plot + ggtitle("hr0=1, hr_s1=1, hr_s2=1, hr_s3=1") +
  xlab("Sample Size") + ylab("Mean Squared Error") +
  theme(plot.title = element_text(hjust = 0.5)) 

CI_plot <- eff_all %>% mutate(N=as.numeric(N)) %>% 
  group_by(prior) %>%
  ggplot(aes(N,Mean,group=prior,color=prior)) + geom_errorbar(aes(ymin=lcl, ymax=ucl,group=prior,linetype=prior),position="dodge") + 
  geom_point(position=position_dodge(45)) +
  geom_hline(yintercept=log(1)) + theme(text = element_text(size = 9),element_line(size =0.2))+
  labs(color  = "Model", linetype = "Model")
p_CI <- CI_plot + ggtitle("hr0=1, hr_s1=1, hr_s2=1, hr_s3=1") +
  xlab("Sample Size") + ylab("95% Coverage Probability") +
  theme(plot.title = element_text(hjust = 0.5)) 

prob_plot <- eff_all %>% mutate(N=as.numeric(N)) %>% 
  group_by(prior) %>%
  ggplot() + aes(N,Prob,group=prior,color=prior,linetype=prior) + geom_line() +
  scale_y_continuous(breaks = seq(0,1,by=0.02)) + theme(text = element_text(size = 9),element_line(size =0.2))+
  labs(color  = "Model", linetype = "Model")
p_prob <- prob_plot + ggtitle("hr0=1, hr_s1=1, hr_s2=1, hr_s3=1") +
  xlab("Sample Size") + ylab("Rejection rate") +
  theme(plot.title = element_text(hjust = 0.5)) 

#write.csv(eff_all, "./unequal_ss_1.csv", row.names=FALSE)

#pdf("unequal_ss_1.pdf")
#par(mfrow=c(2,2))
#plot(p_bias)
#plot(p_mse)
#plot(p_CI)
#plot(p_prob)
#dev.off()