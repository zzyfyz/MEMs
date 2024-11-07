# Load required libraries
library(dplyr)
library(survival)
library(ggplot2)
library(ggpubr)

# Set up parameters and progress bar
set.seed(123)
nsim <- 500

xind <- seq(50, 500, by = 50) # Assuming sample_sizes is defined and contains sample sizes
base_dir <- "C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Data/Equal/1_1_1_1"
pb <- txtProgressBar(0, nsim, style = 3)

# Initialize result arrays
result_pmp <- array(NA_real_, dim = c(nsim, 8, length(xind)),
                    dimnames = list("simulation" = 1:nsim, "model" = 1:8, "N" = xind))

result_ss <- array(NA_real_, dim = c(nsim, length(xind)),
                   dimnames = list("simulation" = 1:nsim, "N" = xind))

result_effm <- array(NA_real_, dim = c(nsim, 4, length(xind)),
                     dimnames = list("simulation" = 1:nsim, "cohort" = 1:4, "N" = xind))

result_effsd <- array(NA_real_, dim = c(nsim, 4, length(xind)),
                      dimnames = list("simulation" = 1:nsim, "cohort" = 1:4, "N" = xind))

result_cred_ints <- array(NA_real_, dim = c(nsim, 4, 2, length(xind)),
                          dimnames = list("simulation" = 1:nsim, "cohort" = 1:4, 
                                          "CI" = c("LowBound", "UpBound"), "N" = xind))

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
    
    #test <- coxph(Surv(time,censor) ~ trt,data=S1)
    #summary(test)
    
    #main functions
    #S: umber of total cohorts(primary + supplemental)
    #mod.mat: matrix of model-specific cohort classifications
    #K: number of intervals for piecewise exponential model
    #tau: change points in the primary cohort
    #cuttype: primary-intervals follow change points from primary cohort
    #         equal-create intervals with equal lengths
    #         events-interval bounds such that approximately same number of events are in each interval    
    
    main <- function(S,mod.mat,K,tau_p,cuttype,priortype,W,alpha0){
      ### Derive values from data
      S <- S                              
      T0 <- S                             # max number of distinct cohorts to consider
      N <- nrow(dat)                      # total number of subjects
      X.mat <- NULL                       # N x p matrix of p optional covariates (excluding treatment)
      
      num_models <- 2^(S-1)               # total number of models
      mod.mat <- mod.mat 
      gamma0 = 0
      W = W                               # Subjective probability of exchangeability
      
      ### Create inter
      K <- K                                           
      
      if (cuttype == "primary") {
        
        int.cuts <- tau_p[1:length(tau_p)-1]
        dat$interval <- apply( cbind(dat$time), 1, function(x){ sum(x > int.cuts) } )
      }
      else if (cuttype == "equal") {
        
        int.cuts <- c( 0, max(dat$time)/K, max(dat$time)/K*2, max(dat$time)/4*3)
        dat$interval <- apply( cbind(dat$time), 1, function(x){ sum(x > int.cuts) } )
      }
      else if (cuttype == "events") {
        
        int.cuts <- c( 0, quantile( dat$time[which(dat$censor == 1)],
                                    probs = seq(1/K, 1, length = K)[-K] ) )
        dat$interval <- apply( cbind(dat$time), 1, function(x){ sum(x > int.cuts) } )
      }     
      
      
      ### Construct design matrix
      W.trt.mat <- matrix(0, nrow = N, ncol = S)
      for(i in 1:S){
        W.trt.mat[,i] <- ifelse( dat$trt == 1 & dat$cohort == unique(dat$cohort)[i], 1, 0 )
      }
      W.mat <- cbind(W.trt.mat, X.mat)
      
      
      ### Prior elicitation
      cov.preds <- NULL                                # predictions of true covariate effects (NULL if no covariates)
      mu0 <- rep(0, S)    # mean vector for prior on treatment and covariate effects
      Sig0 <- diag( 1000, S)  # covariance matrix for prior on treatment effects
      
      eta0 <- matrix( 0.01, nrow = S, ncol = K )       # shape hyperparameters for prior on baseline hazards
      phi0 <- matrix( 0.01, nrow = S, ncol = K )       # rate hyperparameters for prior on baseline hazards
      alpha0 <- alpha0                                      # pre-specified value for prior model probabilities
      
      if (priortype == 1) {
        
        ##Model prior from bmabasket approach
        construct_modPriors <- function(S, T0, mod_mat, alpha = 0){
          
          
          # Determine number of models
          numMods <- 2^(S-1) 
          ml.prop <- rep(0, numMods)
          
          # Initialize model priors
          mod_priors <- numeric(numMods)
          num_distinct <- 0
          for(l in 1:numMods){
            num_distinct <- length(unique(mod_mat[l,]))
            mod_priors[l] <- exp(num_distinct * alpha)
          }
          mod_priors <- mod_priors / sum(mod_priors)
          
          return(mod_priors)
          
        }
        
        mod.priors <- construct_modPriors(S, T0, mod_mat=mod.mat,alpha = alpha0)
      }
      else if (priortype == 2) {
        
        construct_modPriors <- function(S, W, mod_mat){
          
          # Determine number of models
          numMods <- 2^(S-1) 
          ml.prop <- rep(0, numMods)
          
          # Initialize model priors
          mod_priors <- numeric(numMods)
          num_distinct <- 0
          p <- matrix(NA_real_, nrow = numMods, ncol = S-1)
          for(l in 1:numMods){
            num_distinct <- length(unique(mod_mat[l,]))
            for(j in 2:S){
              p[l,j-1] <- ifelse(mod_mat[l,j]>1, 1-W[j-1],W[j-1])
            }
            mod_priors[l] <- prod(p[l,])
          }
          
          return(mod_priors)
          
        }
        
        mod.priors <- construct_modPriors(S, W=W,mod_mat=mod.mat)
      }
      else if (priortype == 3) {
        
        W <- cbind(p_1$pmp[1],p_2$pmp[1],p_3$pmp[1])
        
        construct_modPriors <- function(S, W, mod_mat){
          
          # Determine number of models
          numMods <- 2^(S-1) 
          ml.prop <- rep(0, numMods)
          
          # Initialize model priors
          mod_priors <- numeric(numMods)
          num_distinct <- 0
          p <- matrix(NA_real_, nrow = numMods, ncol = S-1)
          for(l in 1:numMods){
            num_distinct <- length(unique(mod_mat[l,]))
            for(j in 2:S){
              p[l,j-1] <- ifelse(mod_mat[l,j]>1, 1-W[j-1],W[j-1])
            }
            mod_priors[l] <- prod(p[l,])
          }
          
          return(mod_priors)
          
        }
        mod.priors <- construct_modPriors(S, W=W,mod_mat=mod.mat)
      }     
      
      
      # Calculate delta_ijk * (y_ij - m_{k-1}) + sum_(g=k+1)^K( delta_ijg * (m_k - m_{k-1}) )
      # for all subjects in cohort i for each time interval k. Save (n.i x K) matrix for each
      # cohort in a list. These summands are later used to construct phi.tilde.ik.
      phi.tilde.summand.list <- list()
      for(i in 1:S){
        
        dat.i <- dat[ which(dat$cohort == i), ]    # data for only cohort i
        
        # matrix to store summand for all subjects in cohort i for each time interval k, k=1,...,K
        summand.mat.i <- matrix( 0, nrow = nrow(dat.i), ncol = K )
        for(k in 1:K){
          summand.val <- ifelse( dat.i$interval == k, dat.i$time - int.cuts[k], 0 )
          summand.val <- ifelse( dat.i$interval > k, int.cuts[k+1] - int.cuts[k], summand.val )
          summand.mat.i[,k] <- summand.val
        }
        phi.tilde.summand.list[[i]] <- summand.mat.i
        
      }
      
      eta.tilde <- matrix( 0, nrow = S, ncol = K )
      for(i in 1:S){
        for(k in 1:K){
          d.ijk.times.nu.ij <- ifelse( dat$cohort == i & dat$censor == 1 &
                                         dat$interval == k, 1, 0 )
          eta.tilde[i,k] <- sum( d.ijk.times.nu.ij ) + eta0[i,k]
        }
      }
      
      ### Functions for marginal posterior distributions of regression effects
      
      ## Approximate h(theta_l) function, where h(theta_l) = log( p(theta_l|D, M_l) )
      ##   i.e., the log of the full conditional distribution of theta_l, l = 1,...,L
      # x: place holder for (T_l + p)-dimensional vector theta_l
      # S: number of cohorts
      # mu.l: prior mean vector for theta_l
      # Sig.l: prior covariance matrix for theta_l
      # phi.0: (S x K) matrix of phi.ik values, i=1,...,S, k=1,...,K
      # eta.tilde: (S x K) matrix with values of eta.tilde.ik
      # phi.tilde.summand.list: list with all S (n.i x K) matrices of summands for phi.tilde.ik
      # dat.all.regs: full dataset for all subjects from all cohorts
      # W.l: (N x (T_l + p)) design matrix
      h.theta.l <- function( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                             dat.all.regs, W.l ){
        
        # Verify that x is saved as a (T_l + p) x 1 matrix
        x <- cbind(x)
        
        # Calculate phi.tilde.ik for all i and k
        # Also calculate log of exp( sum_j^n.i delta_ijk * nu_ij * w_ijl' theta_l ) for all i
        K <- ncol(phi.0)
        phi.tilde <- matrix( 0, nrow = nrow(phi.0), ncol = ncol(phi.0) )
        log.exp.piece <- matrix( 0, nrow = nrow(phi.0), ncol = ncol(phi.0) )
        for(i in 1:S){
          
          # phi.tilde.ik
          reg.i <- sort( unique(dat.all.regs$cohort) )[i]
          W.il <- as.matrix( W.l[ which(dat.all.regs$cohort == reg.i), ] )
          phi.tilde[i,] <- phi.0[i,] + colSums( phi.tilde.summand.list[[i]] * as.numeric(exp(W.il %*% x)) )
          
          # log of exp( sum_j^n.i delta_ijk * nu_ij * w_ijl' theta_l )
          dat.i <- dat.all.regs[ which(dat.all.regs$cohort == reg.i), ]
          for(k in 1:K){
            log.exp.piece[i,k] <- sum( dat.i$censor[ which(dat.i$interval == k) ] *
                                         as.numeric(W.il[ which(dat.i$interval == k), ] %*% x) )
          }
          
        }
        
        # Log of full conditional of theta_l
        log.val <- sum( eta.0 * log(phi.0) - lgamma(eta.0) + lgamma(eta.tilde) -
                          eta.tilde * log(phi.tilde) + log.exp.piece ) -
          .5*length(x)*log(2*pi) - .5*log(det(Sig.l)) - .5*t(x - mu.l) %*% solve(Sig.l) %*% (x - mu.l)
        
        return(log.val)
        
      }
      
      
      ## Negative of h(theta_l) function
      # See function "h.theta.l" for function input descriptions
      neg.h.theta.l <- function( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                                 dat.all.regs, W.l ){
        val <- h.theta.l( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                          dat.all.regs, W.l )
        return( -1 * val )
      }
      
      ## Gradient of h(theta_l) function
      # See function "h.theta.l" for function input descriptions
      h.p.theta.l <- function( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                               dat.all.regs, W.l ){
        
        # Verify that x is saved as a (T_l + p) x 1 matrix
        x <- cbind(x)
        
        # First piece
        p1 <- -solve(Sig.l) %*% (x - mu.l)
        
        # Second piece
        p2.i <- matrix(0, nrow = length(x), ncol = S)
        for(i in 1:S){
          reg.i <- sort( unique(dat.all.regs$cohort) )[i]
          c.hat.i <- phi.tilde.summand.list[[i]]
          W.il <- as.matrix( W.l[ which(dat.all.regs$cohort == reg.i), ] )
          p2.i[,i] <- t( rowSums( t( as.numeric( colSums( as.numeric(exp(W.il %*% x)) * c.hat.i ) +
                                                   phi.0[i,] )^-1 * as.numeric(eta.tilde[i,]) *
                                       t( t(W.il) %*% ( as.numeric(exp(W.il %*% x)) * c.hat.i ) ) ) ) )
        }
        p2 <- rowSums(p2.i)
        
        # Third piece
        p3 <- colSums( W.l[ which(dat.all.regs$censor == 1), ] )
        
        # Gradient of h(gamma_l) function
        grad <- p1 - p2 + p3
        return(grad)
        
      }
      
      
      ## Hessian of h(theta_l) function
      # See function "h.theta.l" for function input descriptions
      h.pp.theta.l <- function( x, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde, phi.tilde.summand.list,
                                dat.all.regs, W.l ){
        
        # Verify that x is saved as a (T_l + p) x 1 matrix
        x <- cbind(x)
        
        # First piece
        p1 <- -solve(Sig.l)
        
        # Second piece
        p2 <- matrix(0, nrow = ncol(W.l), ncol = ncol(W.l))
        K <- ncol(phi.0)
        for(i in 1:S){
          reg.i <- sort( unique(dat.all.regs$cohort) )[i]
          W.il <- as.matrix( W.l[ which(dat.all.regs$cohort == reg.i), ] )
          for(k in 1:K){
            c.hat.ik <- phi.tilde.summand.list[[i]][,k]
            p2 <- p2 + eta.tilde[i,k] *
              ( (t( c.hat.ik * as.numeric(exp(W.il %*% x)) * W.il ) %*% W.il) *
                  (phi.0[i,k] + sum(exp(W.il %*% x) * c.hat.ik) )^-1 -
                  cbind( colSums( c.hat.ik * as.numeric(exp(W.il %*% x)) * W.il ) ) %*%
                  rbind( colSums( c.hat.ik * as.numeric(exp(W.il %*% x)) * W.il ) ) *
                  (phi.0[i,k] + sum( exp(W.il %*% x) * c.hat.ik ) )^-2 )
          }
        }
        
        # Hessian of h(gamma_(l,t)) function
        hess <- p1 - p2
        return(hess)
        
      }
      
      
      
      ### Function to construct W.l (updated design matrix)
      update_W_l <- function(W, modMat.l){
        
        # Extract necessary values
        N <- nrow(W)                               # number of subjects
        S <- length(modMat.l)                      # number of cohorts
        unq.parms.l <- sort(unique(modMat.l))
        T.l <- length(unq.parms.l)                 # number of distinct cohort labels
        
        # Construct W.l (updated design matrix)
        W.lt <- matrix(0, nrow = N, ncol = T.l)    # columns indicate distinct cohortal trtmt assigments
        
        for(t in 1:T.l){
          for(i in 1:S){
            for(k in 1:N){
              if(modMat.l[i] == unq.parms.l[t] & W[k,i] == 1){
                W.lt[k,t] <- 1
              }
            }
          }
        }
        
        return(W.lt)
        
      }
      
      ## Approximated log of marginal distribution of data given model M_l
      # theta.mode: posterior mode of pr(theta|D,M_l)
      # See function "h.theta.l" for function remaining input descriptions
      log.p.D.Ml <- function( theta.mode, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde,
                              phi.tilde.summand.list, dat.all.regs, W.l ){
        
        # Negative inverse of Hessian matrix
        psi.mat <- -solve( h.pp.theta.l( theta.mode, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde,
                                         phi.tilde.summand.list, dat.all.regs, W.l ) )
        
        # h(theta_l) evaluated at posterior mode of theta
        log.post.theta <- h.theta.l( theta.mode, S, mu.l, Sig.l, phi.0, eta.0, eta.tilde,
                                     phi.tilde.summand.list, dat.all.regs, W.l )
        
        # Approximated log of marginal likelihood
        log.val <- .5 * length(theta.mode) * log(2*pi) + .5 * log(det(psi.mat)) + log.post.theta
        return(log.val)
        
      }
      
      # Derive additional values based on function inputs
      n.i <- table(dat$cohort)              # cohortal sample size
      num.covs <- ncol(W.mat) - S           # number of covariates
      
      
      # Vector to store log marginal likelihoods for each model
      L <- nrow(mod.mat)        # number of models in model space
      log.p.D <- numeric(L)
      
      
      # Matrices to hold model-specific posterior summary statistics for each cohort
      rte.mean.l.mat <- matrix( 0, nrow = L, ncol = S )
      rte.sd.l.mat <- matrix( 0, nrow = L, ncol = S )
      rte.gamma0.prob.l.mat <- matrix( 0, nrow = L, ncol = S )
      n.l.mat <- matrix( 0, nrow = L, ncol = 1 )
      
      
      # Matrices to hold model-specific posterior summary statistics for covariate effects
      beta.mean.l.mat <- matrix( 0, nrow = L, ncol = num.covs )
      beta.sd.l.mat <- matrix( 0, nrow = L, ncol = num.covs )
      
      
      ml.prop <- rep(0, 2^(S-1))
      for(l in 1:L){
        # Model 1
        # Details for model M_1
        mod.l <- mod.mat[l,]               # cohort groupings for model l
        T.l <- length( unique(mod.l) )     # number of distinct treatment effects T_l
        n.l <- 0                           # vector to store combined sample sizes
        num.regs.l <- numeric(T.l)         # vector to store number of cohorts in each group
        rte.means.l <- numeric(T.l)        # vector to store distinct gamma_(l,t) means
        rte.var.l <- numeric(T.l)          # vector to store distinct gamma_(l,t) variances
        log.ml.cntrbtn <- numeric(T.l)     # contribution of each cohort to marginal likelihood
        
        # Update dimensions of mean vector and covariance matrix hyperparameters in prior
        # distributions, and update design matrix W.l
        mu0.l <-  mu0[1:T.l]
        Sig0.l <-  diag( 1000, T.l)
        W.l <- update_W_l( W.mat, mod.mat[l,] )
        
        # Calculate posterior mode of p(theta_l|D, M_l)
        # Nelder-Mead method (default) is robust and nearly as fast as L-BFGS-B method
        if( ncol(W.l) > 1 ){
          theta.mode.l <- optim( par = numeric(ncol(W.l)), fn = neg.h.theta.l,
                                 method = "Nelder-Mead", S = S, mu.l = mu0.l, Sig.l = Sig0.l,
                                 phi.0 = phi0, eta.0 = eta0, eta.tilde = eta.tilde,
                                 phi.tilde.summand.list = phi.tilde.summand.list,
                                 dat.all.regs = dat, W.l = W.l )$par
        } else{
          theta.mode.l <- optim( par = numeric(ncol(W.l)), fn = neg.h.theta.l, method = "Brent",
                                 lower = -100, upper = 100, S = S, mu.l = mu0.l, Sig.l = Sig0.l,
                                 phi.0 = phi0, eta.0 = eta0, eta.tilde = eta.tilde,
                                 phi.tilde.summand.list = phi.tilde.summand.list,
                                 dat.all.regs = dat, W.l = W.l )$par
        }
        
        # Calculate posterior covariance matrix of theta_l|D,M_l
        theta.cov.l <- -h.pp.theta.l( x = theta.mode.l, S = S, mu.l = mu0.l, Sig.l = Sig0.l,
                                      phi.0 = phi0, eta.0 = eta0, eta.tilde = eta.tilde,
                                      phi.tilde.summand.list = phi.tilde.summand.list,
                                      dat.all.regs = dat, W.l = W.l )^(-1)
        
        
        # Calculate log of marginal likelihood for model M_l
        log.p.D <- log.p.D.Ml( theta.mode = theta.mode.l, S = S, mu.l = mu0.l, Sig.l = Sig0.l,
                               phi.0 = phi0, eta.0 = eta0, eta.tilde = eta.tilde,
                               phi.tilde.summand.list = phi.tilde.summand.list,
                               dat.all.regs = dat, W.l = W.l )
        
        
        ml.prop[l] <- log.p.D
        
        
        # Calculate combined sample size, number of events, and sum of event times (y)
        # for cohorts that share the t^th distinct trtmt effect
        
        for(i in 1:S){
          if(mod.l[i] == 1){
            
            # Combined cohortal values based on groupings for model M_l
            n.l = n.l + n.i[i]
          }else{
            n.l = n.l
          }
          n.l.mat[l] <- n.l
        }
        
        
        
        # Calculate summary stats for posterior distribution of gamma_(l,t) and store stats
        # with cohorts that share the treatment effect gamma_(l,t), t=1,...,T_l, l=1,...,L
        for(t in 1:T.l){
          
          # Save approximated posterior mean and variance for cohortal treatment effects
          rte.means.l[t] <- theta.mode.l[t]
          rte.var.l[t] <- theta.cov.l[t,t]
          
          # Calculate Pr(gamma_(l,t) < gamma0|D,M_l)
          rte.prob.less.gamma0 <- pnorm(gamma0, rte.means.l[t], sqrt(rte.var.l[t]))
          
          # Store posterior values for each cohort
          for(k in 1:S){
            if(mod.l[k] == t){
              
              # Posterior summary stats
              rte.mean.l.mat[l,k] <- rte.means.l[t]
              rte.sd.l.mat[l,k] <- sqrt(rte.var.l[t])
              rte.gamma0.prob.l.mat[l,k] <- rte.prob.less.gamma0
              
            }
          }
          
        }
      }
      # Calculate posterior model probability for each model
      max.log.p.D <- max(ml.prop)
      ml.prop_ <- exp( ml.prop - max.log.p.D )
      pmp.vec <- ( ml.prop_ * mod.priors ) / sum( ml.prop_ * mod.priors )
      #pmp.vec
      
      
      rte.means <- rbind(pmp.vec) %*% rte.mean.l.mat
      rte.sds <- rbind(pmp.vec) %*% rte.sd.l.mat
      rte.less.gamma0 <- rbind(pmp.vec) %*% rte.gamma0.prob.l.mat
      
      ## Sample size
      ss <- rbind(pmp.vec) %*% n.l.mat
      #ss
      cred.ints <- matrix(0, nrow = S, ncol = 2)
      
      for(i in 1:S){
      # get the 95% credible interval
      n.samp <- 100000  # number of samples to draw
  
      colnames(cred.ints) <- c("LowBound", "UpBound")
      
      # Sample from models based on pmp
      which.mods <- sample(1:L, n.samp, replace = TRUE, prob = pmp.vec)
      
      # Generate samples for global effect
      glob.samp <- rnorm(n.samp, rte.mean.l.mat[which.mods, i], rte.sd.l.mat[which.mods, i])
      
      # Compute 95% credible interval
      cred.ints[i,] <- quantile(glob.samp, c(.025, .975))
      }
      
      my_list<- list("pmp"=pmp.vec,"mean"=rte.means,"sd"=rte.sds, "CI"=cred.ints, "ss"=ss,"test"=rte.less.gamma0,"mod_prior"=mod.priors)
      return(my_list)
    }
    
    mod <- main(4,rbind(c(1,1,1,1),c(1,2,1,1),c(1,1,1,2),c(1,1,2,1),c(1,2,3,1),c(1,2,1,3),c(1,1,2,3),c(1,2,3,4)),
                4,c(0, 1, 4.33, 26, 52),cuttype="primary",priortype=1,W=NA,alpha0=0)
    
    result_pmp[i,1,j] <- mod$pmp[1]
    result_pmp[i,2,j] <- mod$pmp[2]
    result_pmp[i,3,j] <- mod$pmp[3]
    result_pmp[i,4,j] <- mod$pmp[4]
    result_pmp[i,5,j] <- mod$pmp[5]
    result_pmp[i,6,j] <- mod$pmp[6]
    result_pmp[i,7,j] <- mod$pmp[7]
    result_pmp[i,8,j] <- mod$pmp[8]
    
    result_ss[i,j] <- mod$ss
    
    
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
    
    result_cred_ints[i, 1, , j] <- mod$CI[1, ]
    result_cred_ints[i, 2, , j] <- mod$CI[2, ]
    result_cred_ints[i, 3, , j] <- mod$CI[3, ]
    result_cred_ints[i, 4, , j] <- mod$CI[4, ]
  }
  
  setTxtProgressBar(pb,i)
  
}

pmp <- as.data.frame.table(result_pmp, responseName = "percentage", stringsAsFactors = F)  
ss <- as.data.frame.table(result_ss, responseName = "N", stringsAsFactors = F)
effm <- as.data.frame.table(result_effm, responseName = "Mean", stringsAsFactors = F)
effsd <- as.data.frame.table(result_effsd, responseName = "Sd", stringsAsFactors = F) 
prob <- as.data.frame.table(result_prob, responseName = "Prob", stringsAsFactors = F) 
prob$ind <- ifelse(prob$Prob > 0.95, 1, 0)
cred_ints <- as.data.frame.table(result_cred_ints, responseName = "CI_value", stringsAsFactors = FALSE)
cred_ints <- reshape(cred_ints, idvar = c("simulation", "cohort", "N"), timevar = "CI", direction = "wide")

output_dir <- "C:/Users/feiy/OneDrive - The University of Colorado Denver/Documents 1/MEMs/Simulation/Results/Equal/Laplace/Uniform/1_1_1_1"
write.csv(pmp, file.path(output_dir, "pmp_results.csv"), row.names = FALSE)
write.csv(ss, file.path(output_dir, "ss_results.csv"), row.names = FALSE)
write.csv(effm, file.path(output_dir, "effm_results.csv"), row.names = FALSE)
write.csv(effsd, file.path(output_dir, "effsd_results.csv"), row.names = FALSE)
write.csv(prob, file.path(output_dir, "prob_results.csv"), row.names = FALSE)
write.csv(cred_ints, file.path(output_dir, "cred_ints_results.csv"), row.names = FALSE)


