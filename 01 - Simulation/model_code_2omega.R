# Simulation Study
library(tidyverse)
library(nimble)
# Loading graphics packages
library(ggplot2)
library(ggmcmc) # more detail: http://xavier-fim.net/post/using_ggmcmc/

# Ground Truth Parameters
#############################################
true_lambda = 50
true_gamma = 5
true_omega1 = 0.25
true_omega2 = 0.5
true_p_rec = 0.4
true_p_mor = 0.1
true_p_det = 0.70

par_true = data.frame(
  lambda=true_lambda,
  gamma=true_gamma,
  omega1=true_omega1,
  omega2=true_omega2,
  p_det=true_p_det,
  p_rec=true_p_rec,
  p_mor=true_p_mor)

# MCMC parameters
iterations <- 1000000
burnin <- 100000
chains <- 1
thinning <- 200


nimModel <- nimbleCode({ 
  ## Priors and constraints
  mean.lambda ~  dgamma(shape=5,scale=20)  # Prior for mean lambda
  mean.gamma ~ dunif(0,30)                 # Prior for mean gamma
  mean.omega1 ~ dunif(0,5)                 # Prior for mean omega (unobserved cases)
  mean.omega2 ~ dunif(0,5)                 # Prior for mean omega (observed cases)
  mean.p_det ~ dunif(0, 1)                 # Prior for mean probability of detection
  mean.p_rec ~ dunif(0, 1)                 # Prior for mean probability of recovery
  mean.p_mor ~ dunif(0, mean.p_rec)        # Prior for mean probability of death
  
  mean.p_rec2 <- mean.p_rec/(1-mean.p_mor) # probability of recovery among all cases after accounting for all deaths
  # define conjugate probabilities
  pd_p     <- mean.p_mor/mean.p_det # probability of death within observed cases (all deaths are observed), alpha = 1/mean.p_det
  pr_1pd_p <- mean.p_rec/(1-pd_p)   # probability of recovery among observed cases after accounting for observed deaths

  ## Likelihood
  for(site in 1:SITES) {
    for(time in 1:(TIMES-1)) {
      # Latent ARD State Process  
      D[site, time] ~ dbinom(mean.p_mor, N[site, time]) 
      R[site, time] ~ dbinom(mean.p_rec2, N[site, time] - D[site, time]) 
      A[site, time] <- N[site, time] - R[site, time] - D[site, time]
      
      # Detection Process (split active cases by remain active, recover, die)
      d[site, time] ~ dbinom(pd_p, m[site,time+1] + a[site, time])
      r[site, time] ~ dbinom(pr_1pd_p, m[site,time+1] + a[site, time] - d[site, time])    
    }  # 1:(TIMES-1)
    
    # Abundance and Detection
    N[site, 1] ~ dpois(mean.lambda)
    for (time in 2:TIMES){
      # Latent NSG Abundance Process
      N[site, time] <- A[site, time-1] + S[site, time-1] + G[site, time-1]
      S[site, time-1] ~ dpois(mean.omega1*(N[site, time-1]-a[site, time-1]) + mean.omega2*a[site, time-1])
      G[site, time-1] ~ dpois(mean.gamma)
      
      # Detection Process
      n[site, time] ~ dbinom(mean.p_det, N[site, time]-a[site, time-1]) 
    } # TIMES
    n[site, 1] ~ dbinom(mean.p_det, N[site, 1])  # args: p, N 
  }
})


################ BEGIN SIMULATION ###############################
for(iter in 1:50) {
  
  # Create Data
  #############################################
  # 2 site, 30 sampling occassions
  impossible = TRUE
  counter = 0
  while(impossible) {
    counter = counter+1
    if(counter > 1000000) { stop("ERROR: COULD NOT GENERATE DATA AFTER 1e6 ATTEMPTS") }
    Nit = matrix(0, nrow=2, ncol=30)
    Ait = matrix(0, nrow=2, ncol=29)
    Rit = matrix(0, nrow=2, ncol=29)
    Dit = matrix(0, nrow=2, ncol=30)
    Sit = matrix(0, nrow=2, ncol=29)
    Git = matrix(0, nrow=2, ncol=29)
    ait = matrix(0, nrow=2, ncol=30)
    nit = matrix(0, nrow=2, ncol=30)
    rit = matrix(0, nrow=2, ncol=30)
    
    Nit[,1] = rpois(2,true_lambda)
    nit[,1] = rbinom(2,Nit[,1], true_p_det)
    ait[,1] = nit[,1]
    for(site in 1:2) {
      for(time in 2:30) {
        
        Dit[site,time-1] = rbinom(1,Nit[site,time-1], true_p_mor)
        Rit[site,time-1] = rbinom(1,Nit[site,time-1]-Dit[site,time-1], true_p_rec/(1-true_p_mor)) 
        Ait[site,time-1] = Nit[site,time-1]-Dit[site,time-1]-Rit[site,time-1]
        
        Sit[site,time-1] = rpois(1,true_omega1*(Nit[site, time-1]-ait[site, time-1]) + true_omega2*ait[site, time-1])
        Git[site,time-1] = rpois(1,true_gamma)
        Nit[site,time]   = Ait[site,time-1] + Sit[site,time-1] + Git[site,time-1]
        
        nit[site,time]   = rbinom(1,Nit[site,time]-ait[site,time-1],true_p_det)
        rit[site,time-1] = rbinom(1,ait[site,time-1]+nit[site,time]-Dit[site,time-1],true_p_rec/(1-true_p_mor/true_p_det))
        ait[site,time]   = ait[site,time-1]+nit[site,time]-Dit[site,time-1]-rit[site,time-1]
      }
    }
    if(!any(is.na(nit))) { impossible=FALSE }
  }
  
  #############################################
  maxT = 30
  dit = Dit
  n_dim <- dim(nit) # SITES = n_dim[1], TIMES = n_dim[2]
  #############################################
  # initialize model parameters
  Ninit <- 2*(ait+50)
  Rinit <- rit+20
  Ginit <- Rinit
  Sinit <- Ninit-Ginit
  Ninit[,-1] <- NA
  Rinit[,-1] <- NA
  Ginit[,-1] <- NA
  Sinit[,-1] <- NA
  
  model.inits <- list(mean.lambda = 100, 
                      mean.gamma  = 1, 
                      mean.omega1 = 0.45, 
                      mean.omega2 = 0.45, 
                      mean.p_det  = 0.70, 
                      mean.p_rec  = 0.34, 
                      mean.p_mor  = 0.005,
                      N = Ninit,
                      R = Rinit,
                      S = Sinit,
                      G = Ginit)
  
  
  ######################################################
  # create the model graph
  model.fit.nim <- nimbleModel(code = nimModel, 
                               name = "nimble model", 
                               data = list(n=nit, m=nit, r=rit, D=Dit, d=Dit, a=ait),
                               constants = list(SITES=n_dim[1], TIMES=n_dim[2]),
                               inits = model.inits)
  
  # compile the model in C++
  Cmodel.fit.nim <- compileNimble(model.fit.nim)
  
  # Run MCMC
  parameters <- c("mean.lambda", "mean.gamma", "mean.omega1", "mean.omega2", "mean.p_det", "mean.p_rec", "mean.p_mor")
  model.fit.nim.MCMC <- buildMCMC(Cmodel.fit.nim, 
                                  monitors = parameters)
  Cmodel.fit.nim.MCMC <- compileNimble(model.fit.nim.MCMC)
  
  model.samples.nim <- runMCMC(Cmodel.fit.nim.MCMC, 
                               nburnin = burnin, 
                               niter   = iterations, 
                               nchains = chains, 
                               thin    = thinning, 
                               samplesAsCodaMCMC = TRUE)
  
  saveRDS(model.samples.nim, paste0("./MCMC.model.samples.sim_",iter,".RDS"))
  
  model_data = list(Nit=Nit, nit=nit, rit=rit, Dit=Dit, Ait=Ait, Rit=Rit, Sit=Sit, Git=Git, ait=ait)
  saveRDS(model_data, paste0("./model.data.sim_",iter,".RDS"))
  
  par_medians = data.frame(
    median.lambda=median(model.samples.nim[, 'mean.lambda']),
    median.gamma=median(model.samples.nim[, 'mean.gamma']),
    median.omega1=median(model.samples.nim[, 'mean.omega1']),
    median.omega2=median(model.samples.nim[, 'mean.omega2']),
    median.p_det=median(model.samples.nim[, 'mean.p_det']),
    median.p_rec=median(model.samples.nim[, 'mean.p_rec']),
    median.p_mor=median(model.samples.nim[, 'mean.p_mor']))
  
  par_means = data.frame(
    mean.lambda=mean(model.samples.nim[, 'mean.lambda']),
    mean.gamma=mean(model.samples.nim[, 'mean.gamma']),
    mean.omega1=mean(model.samples.nim[, 'mean.omega1']),
    mean.omega2=mean(model.samples.nim[, 'mean.omega2']),
    mean.p_det=mean(model.samples.nim[, 'mean.p_det']),
    mean.p_rec=mean(model.samples.nim[, 'mean.p_rec']),
    mean.p_mor=mean(model.samples.nim[, 'mean.p_mor']))
  
  par_lambda = data.frame(lambda_hat=mean(nit[,1]/par_means$mean.p_det))
  
  save_dat = cbind(data.frame(iter=iter), par_lambda, par_medians, par_means, par_true)

  bool = file.exists("./simulation_data.csv")
  if(bool) {
    write.table(x = save_dat, row.names = F, col.names = FALSE,
                file = "./simulation_data.csv",
                append = bool)
  } else {
    write.table(x = save_dat, row.names = F, col.names = colnames(save_dat),
                file = "./simulation_data.csv",
                append = bool)
  }
  
}

