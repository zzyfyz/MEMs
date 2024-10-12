library(tidyr)
library(dplyr)

dat1 <- read.csv("./equal_ss_1.csv")
dat2 <- read.csv("./equal_ss_2.csv")
dat3 <- read.csv("./equal_ss_3.csv")
dat4 <- read.csv("./equal_ss_4.csv")
dat5 <- read.csv("./equal_ss_5.csv")

dat1$mod <- "Scenario 1"
dat2$mod <- "Scenario 2"
dat3$mod <- "Scenario 3"
dat4$mod <- "Scenario 4"
dat5$mod <- "Scenario 5"

dat_all <- rbind(dat1,dat2,dat3,dat4,dat5)
dat_all$Model <- factor(dat_all$prior,levels=c('MEMs-Empirical','MEMs-Exponential','BHM','CPHM'))
dat_all$bias <- dat_all$Mean - log(1)
dat_all$mse <- (dat_all$Mean - log(1))^2 + (dat_all$Sd)^2
dat_all$N <- dat_all$N*2

subset <- subset(dat_all,select=c(N,Model,mse,mod))
data_wide <- subset %>%
  pivot_wider(names_from = Model, values_from = mse)
data_wide$pct_emp_cph <- (1 - data_wide$`MEMs-Empirical`/data_wide$CPHM)*100
data_wide$pct_exp_cph <- (1 - data_wide$`MEMs-Exponential`/data_wide$CPHM)*100
data_wide$pct_emp_bhm <- (1 - data_wide$`MEMs-Empirical`/data_wide$BHM)*100
data_wide$pct_exp_bhm <- (1 - data_wide$`MEMs-Exponential`/data_wide$BHM)*100

subset <- subset(dat_all,select=c(N,Model,Prob,mod))
data_wide <- subset %>%
  pivot_wider(names_from = Model, values_from = Prob)

myTheme2 <- theme(
  axis.title = element_text(size = rel(1.4),vjust = -1),
  axis.text = element_text(size = rel(1)),
  axis.text.x = element_text(angle=30),
  plot.title = element_text(size = rel(1.6)),
  plot.subtitle = element_text(size = rel(1.5)),
  plot.caption = element_text(size = rel(1.5)),
  legend.position = "right",
  legend.text=element_text(size=10)
)

png('./bias_equal_null.png', height = 8, width = 14, units = "in", res = 480)
p <- ggplot(dat_all, aes(x = N, y = bias, group = Model)) +
  geom_line(aes(linetype=Model)) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  geom_hline(yintercept=0, col = "orange") +
  scale_x_continuous(breaks = seq(100, 1000, by = 100)) +
  theme_bw() +  
  xlab("Sample Size") +
  ylab("Bias") + 
  myTheme2

p + facet_wrap(~factor(mod), ncol=5) +
  theme(strip.text.x = element_text(size = 12))
invisible(dev.off())


png('./MSE_equal_null.png', height = 8, width = 14, units = "in", res = 480)
p <- ggplot(dat_all, aes(x = N, y = mse, group = Model)) +
  geom_line(aes(linetype=Model)) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  geom_hline(yintercept=0, col = "orange") +
  scale_x_continuous(breaks = seq(100, 1000, by = 100)) +
  theme_bw() +  
  xlab("Sample Size") +
  ylab("MSE") + 
  myTheme2

p + facet_wrap(~factor(mod), ncol=5) +
  theme(strip.text.x = element_text(size = 12))
invisible(dev.off())

png('./CI_equal_null.png', height = 8, width = 14, units = "in", res = 480)
p <- ggplot(dat_all, aes(x = N, y = Mean, group = Model)) +
  geom_point(size = 2, aes(shape = Model), position = position_dodge(30)) +
  geom_pointrange(aes(ymin = lcl, ymax = ucl, shape = Model,linetype = Model), 
                  position = position_dodge(30), size = .2) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  geom_hline(yintercept=0, col = "orange") +
  scale_x_continuous(breaks = seq(100, 1000, by = 100)) +
  theme_bw() +  
  xlab("Sample Size") +
  ylab("95% Coverage Probability") + 
  myTheme2

p + facet_wrap(~factor(mod), ncol=2) +
  theme(strip.text.x = element_text(size = 12),
        legend.position = c(0.6, 0.2))
invisible(dev.off())


png('./type1_equal_null.png', height = 8, width = 14, units = "in", res = 480)
p <- ggplot(dat_all, aes(x = N, y = Prob, group = Model)) +
  geom_line(aes(linetype=Model)) +
  scale_linetype_manual(values=c("solid", "twodash", "dotted", "longdash")) +
  geom_hline(yintercept=0.05, col = "orange") +
  scale_x_continuous(breaks = seq(100, 1000, by = 100)) +
  theme_bw() +  
  xlab("Sample Size") +
  ylab("False Rejection Rate") + 
  myTheme2

p + facet_wrap(~factor(mod), ncol=5) +
  theme(strip.text.x = element_text(size = 12))
invisible(dev.off())
