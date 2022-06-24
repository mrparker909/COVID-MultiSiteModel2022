# In this code we fit the single site COVID model for Northern Health only
library(tidyverse)
library(nimble)

# Load Covariates
#############################################
cov_dat = read.csv("Data/BC Covariates by Health Authority Region.csv")

# Load Data
#############################################
# load data:
dat = read.csv("Data/BC Case Counts by Health Authority Region.csv")

pop_dat = data.frame(HA = dat$HA %>% unique(), Npop=c(1.8,0.74,0.84,0.3,1.25)) # approximate population size in millions for each Health Authority Region

dat$HA %>% unique()
# [1] "Fraser"            "Interior"          "Vancouver Island"  "Northern"          "Vancouver Coastal"
USED_SITES = 4

temp1 <- reshape(dat %>% select(HA, Date, NewCases), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
nit_temp = data.matrix(temp2)

temp1 <- reshape(dat %>% select(HA, Date, NewRecoveries), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
rit_temp = data.matrix(temp2)

temp1 <- reshape(dat %>% select(HA, Date, NewDeaths), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
Dit_temp = data.matrix(temp2)

temp1 <- reshape(dat %>% select(HA, Date, ActiveCases), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
ait_temp = data.matrix(temp2)

maxT = 32
nit = matrix(matrix(nit_temp[,1:maxT],nrow=1)[,-1],nrow=1)
rit = matrix(matrix(rit_temp[,1:maxT],nrow=1)[,-c(1,2)],nrow=1)
Dit = matrix(matrix(Dit_temp[,1:maxT],nrow=1)[,-c(1,2)],nrow=1)
dit = Dit
ait = matrix(matrix(ait_temp[,1:maxT],nrow=1)[,-1],nrow=1)
n_dim <- dim(nit) # SITES = n_dim[1], TIMES = n_dim[2]

temp1 <- reshape(cov_dat %>% select(HA, Date, New.Tests), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
testvol_temp = data.matrix(temp2)

testvolit = matrix(matrix(testvol_temp[,1:maxT],nrow=1)[,-1],nrow=1)
testvolnorm = t(apply(testvolit, 1, function(x)(x)/(sd(x)))) # here we standardize testing volumes

temp1 <- reshape(cov_dat %>% select(HA, Date, BC_Phase1), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
BCphase1 = data.matrix(temp2[,-1])

temp1 <- reshape(cov_dat %>% select(HA, Date, BC_Phase2), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
BCphase2 = data.matrix(temp2[,-1])

temp1 <- reshape(cov_dat %>% select(HA, Date, BC_Phase3a), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
BCphase3a = data.matrix(temp2[,-1])

temp1 <- reshape(cov_dat %>% select(HA, Date, BC_Phase3b), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
BCphase3b = data.matrix(temp2[,-1])


# MCMC parameters
iterations <- 1200000
burnin <- 1000000
chains <- 1
thinning <- 1

#############################################
# initialize model parameters
Ninit <- 2*(ait+50)
Rinit <- cbind(rit,0)+20
Ginit <- Rinit
Sinit <- Ninit[,-maxT]-Ginit
Ninit[,-1] <- NA
Rinit[,-1] <- NA
Ginit[,-1] <- NA
Sinit[,-1] <- NA

model.inits <- list(p_det.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    mean.p_rec  = 0.4, 
                    mean.p_mor  = 0.005,
                    lambda = 500,
                    gamma = 2,
                    site.beta.pdet = rep(0, times=n_dim[1]),
                    site.beta.omega1 = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    site.beta.omega2 = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    beta.test.vol = rep(0.1, times=n_dim[1]),
                    beta.phase1.omega1 = 1.05,
                    beta.phase2.omega1 = 1.05,
                    beta.phase3a.omega1 = 1.05,
                    beta.phase3b.omega1 = 1.05,
                    beta.phase1.omega2 = 0,
                    beta.phase2.omega2 = 0,
                    beta.phase3a.omega2 = 0,
                    beta.phase3b.omega2 = 0,
                    N = Ninit,
                    R = Rinit,
                    S = Sinit,
                    G = Ginit
                    )

nimModel <- nimbleCode({ 
  ## Priors and constraints
  mean.p_rec ~ dunif(0, 1)              
  mean.p_mor ~ dunif(0, mean.p_rec)     
  
  lambda ~ dgamma(shape=5,scale=20)
  gamma  ~ dunif(0,50)
  
  beta.phase1.omega1 ~ dunif(0,15)
  beta.phase2.omega1 ~ dunif(0,15)
  beta.phase3a.omega1 ~ dunif(0,15)
  beta.phase3b.omega1 ~ dunif(0,15)
  beta.phase1.omega2 ~ dunif(0,15)
  beta.phase2.omega2 ~ dunif(0,15)
  beta.phase3a.omega2 ~ dunif(0,15)
  beta.phase3b.omega2 ~ dunif(0,15)
  
  
  for(site in 1:SITES) {
    beta.test.vol[site] ~ dunif(0,5)
    site.beta.pdet[site] ~ dunif(0,1)
    
    for(time in 1:TIMES) {
      site.beta.omega1[site, time] <-
        beta.phase1.omega1 * phase1.cov[site, time] + 
        beta.phase2.omega1 * phase2.cov[site, time] + 
        beta.phase3a.omega1 * phase3a.cov[site, time] +
        beta.phase3b.omega1 * phase3b.cov[site, time]
      site.beta.omega2[site, time] <- 
        beta.phase1.omega2 * phase1.cov[site, time] + 
        beta.phase2.omega2 * phase2.cov[site, time] + 
        beta.phase3a.omega2 * phase3a.cov[site, time] +
        beta.phase3b.omega2 * phase3b.cov[site, time]
      
      p_det.matrix[site, time] <-  site.beta.pdet[site] + beta.test.vol[site] * test.volume[site, time]
      
      # define conjugate probabilities
      pd_p[site, time]     <- mean.p_mor/p_det.matrix[site, time] # probability of death within observed cases (all deaths are observed), alpha = 1/p
      pr_1pd_p[site, time] <- mean.p_rec/(1-pd_p[site, time])     # probability of recovery among observed cases after accounting for observed deaths
    }
  }
  
  mean.p_rec2 <- mean.p_rec/(1-mean.p_mor) # probability of recovery among all cases after accounting for all deaths
  
  ## Likelihood
  for(site in 1:SITES) {
    for(time in 1:(TIMES-1)) {
      # Latent ARD State Process  
      D[site, time] ~ dbinom(mean.p_mor, N[site, time]) 
      R[site, time] ~ dbinom(mean.p_rec2, N[site, time] - D[site, time]) 
      A[site, time] <- N[site, time] - R[site, time] - D[site, time]
      
      # Detection Process (split active cases by remain active, recover, die)
      d[site, time] ~ dbinom(pd_p[site, time], m[site,time+1] + a[site, time])
      r[site, time] ~ dbinom(pr_1pd_p[site, time], m[site,time+1] + a[site, time] - d[site, time])    
    }  # 1:(TIMES-1)
    
    # Abundance and Detection
    N[site, 1] ~ dpois(lambda)
    for (time in 2:TIMES){
      # Latent Abundance Process
      N[site, time] <- A[site, time-1] + S[site, time-1] + G[site, time-1]
      S[site, time-1] ~ dpois(site.beta.omega1[site, time]*(N[site, time-1]-a[site, time-1])*
                                (Npop - N[site, time-1])/Npop + # damping effect for population saturation
                                site.beta.omega2[site, time]*a[site, time-1])
      G[site, time-1] ~ dpois(gamma)
      
      # Detection Process
      n[site, time] ~ dbinom(p_det.matrix[site, time], N[site, time]-a[site, time-1]) 
    } # TIMES
    n[site, 1] ~ dbinom(p_det.matrix[site, 1], N[site, 1])
  }
})

# create the model graph
model.fit.nim <- nimbleModel(code = nimModel,
                             name = "nimble model", 
                             data = list(n=nit, m=nit, r=rit, D=Dit, d=Dit, a=ait),
                             constants = list(SITES=n_dim[1], 
                                              TIMES=n_dim[2],
                                              test.volume=testvolnorm,
                                              phase1.cov=BCphase1,
                                              phase2.cov=BCphase2,
                                              phase3a.cov=BCphase3a,
                                              phase3b.cov=BCphase3b,
                                              Npop=pop_dat$Npop[USED_SITES]*1e6
                                              ),
                             inits = model.inits)

# compile the model in C++
Cmodel.fit.nim <- compileNimble(model.fit.nim)

# Run MCMC
WAICparameters <- c("N", "G", "d", "S", "D", "R",
                    "site.beta.pdet",
                    "beta.test.vol",
                    "beta.phase1.omega1",
                    "beta.phase2.omega1",
                    "beta.phase3a.omega1",
                    "beta.phase3b.omega1",
                    "beta.phase1.omega2",
                    "beta.phase2.omega2",
                    "beta.phase3a.omega2",
                    "beta.phase3b.omega2",
                    "lambda",  
                    "gamma",   
                    "mean.p_rec", 
                    "mean.p_mor")
parameters <- c("beta.test.vol",
                "site.beta.pdet",
                "beta.phase1.omega1",
                "beta.phase2.omega1",
                "beta.phase3a.omega1",
                "beta.phase3b.omega1",
                "beta.phase1.omega2",
                "beta.phase2.omega2",
                "beta.phase3a.omega2",
                "beta.phase3b.omega2",
                "lambda",  
                "gamma", 
                "mean.p_rec", 
                "mean.p_mor")
model.fit.nim.MCMC <- buildMCMC(conf = configureMCMC(Cmodel.fit.nim,
                                                     enableWAIC = TRUE,
                                                     monitors = WAICparameters, 
                                                     monitors2 = parameters))
Cmodel.fit.nim.MCMC <- compileNimble(model.fit.nim.MCMC)

model.samples.nim <- runMCMC(Cmodel.fit.nim.MCMC, 
                             nburnin = burnin, 
                             niter   = iterations, 
                             nchains = chains, 
                             thin    = thinning, 
                             samplesAsCodaMCMC = TRUE,
                             WAIC = TRUE)

WAIC_dat = data.frame(model = "single_site_NHA", 
                      WAIC = model.samples.nim$WAIC, 
                      parameters=paste(parameters, collapse = "/"))

bool = file.exists("./WAIC_dat.csv")
if(bool) {
  write.table(x = WAIC_dat, row.names = F, col.names = FALSE,
              file = "./WAIC_dat.csv", 
              sep = ",",
              append = bool)
} else {
  write.table(x = WAIC_dat, row.names = F, col.names = colnames(WAIC_dat),
              file = "./WAIC_dat.csv",
              sep = ",",
              append = bool)
}

saveRDS(model.samples.nim, "./MCMC.model.samples.RDS")

