#Simulate data 

f <- function(n1,n2,h0,G,P,tau) {
  # set the time points
  DT = tau[2:length(tau)] - tau[1:length(tau)-1]
  
  # create a helping matrix
  LD <- matrix(0,nrow=length(tau), ncol=length(G))
  LD[lower.tri(LD)]<-1
  
  # start with population treatment group
  LS <- log(1-runif(n1))
  GP <- P[1]*G
  # determine the ln(S) for all TAU
  LSM <- -h0 * as.vector(LD %*% (GP * DT))
  # find the appropriate time interval
  t1 <- rep(NA,n1)
  for (i in 1:4) {
    t1 <- ifelse(LSM[i]>=LS & LS>LSM[i+1], tau[i] + (LSM[i] - LS)/h0/GP[i], t1)
  }
  # end of population treatment group
  
  # start with population placebo group
  LS <- log(1-runif(n2))
  GP <- P[2]*G
  # determine the ln(S) for all TAU
  LSM <- -h0 * as.vector(LD %*% (GP * DT))
  # find the appropriate time interval
  t2 <- rep(NA,n2)
  for (i in 1:length(tau)-1) {
    t2 <- ifelse(LSM[i]>=LS & LS>LSM[i+1], tau[i] + (LSM[i] - LS)/h0/GP[i], t2)
  }
  # end of population placebo group
  
  # combine both populations
  sim.data <- data.frame( rbind(cbind(time=t1,trt=1),cbind(t2,0)))
  sim.data$censor <- ifelse(is.na(sim.data$time), 0, 1)
  sim.data$time <- ifelse(is.na(sim.data$time), tau[length(tau)], sim.data$time)
  sim.data$id <- row(sim.data)[,1]
  # the data set is ready
  return(sim.data)
}