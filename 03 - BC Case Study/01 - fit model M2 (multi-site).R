# In this code we fit the COVID multi site model to the BC HA data
# lam~region, gam~pha+region, w1~pha, w2~pha, pr, pd, p~vol
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
USED_SITES = 1:5

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
nit = nit_temp[,1:maxT][,-1]
rit = rit_temp[,1:maxT][,-c(1,2)]
Dit = Dit_temp[,1:maxT][,-c(1,2)]
dit = Dit
ait = ait_temp[,1:maxT][,-1]
n_dim <- dim(nit) # SITES = n_dim[1], TIMES = n_dim[2]

temp1 <- reshape(cov_dat %>% select(HA, Date, New.Tests), idvar = "HA", timevar = "Date", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",0:(ncol(temp2)-1))
testvol_temp = data.matrix(temp2)

testvolit = testvol_temp[,1:maxT][,-1]
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

model.inits <- list(p_rec = 0.35,
                    p_mor = 0.01,
                    vol.pdet = 0.15,
                    site.lambda = 2*(ait[,1]+50),
                    site.gamma  = c(10,10,10,10,10),
                    phase1.gamma = 0,
                    phase2.gamma = 0,
                    phase3a.gamma = 0,
                    phase3b.gamma = 0,
                    phase1.omega1 = 0.5,
                    phase2.omega1 = 0.5,
                    phase3a.omega1 = 0.5,
                    phase3b.omega1 = 0.5,
                    phase1.omega2 = 0.5,
                    phase2.omega2 = 0.5,
                    phase3a.omega2 = 0.5,
                    phase3b.omega2 = 0.5,
                    p_det.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    gamma.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    omega1.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    omega2.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    N = Ninit,
                    R = Rinit,
                    S = Sinit,
                    G = Ginit
)

nimModel <- nimbleCode({ 
  ## Priors and constraints
  p_rec ~ dunif(0, 1)          
  p_mor ~ dunif(0, p_rec)     
  
  phase1.gamma ~ dunif(0,50)
  phase2.gamma ~ dunif(0,50)
  phase3a.gamma ~ dunif(0,50)
  phase3b.gamma ~ dunif(0,50)
  
  phase1.omega1 ~ dunif(0,5)
  phase2.omega1 ~ dunif(0,5)
  phase3a.omega1 ~ dunif(0,5)
  phase3b.omega1 ~ dunif(0,5)
  
  phase1.omega2 ~ dunif(0,5)
  phase2.omega2 ~ dunif(0,5)
  phase3a.omega2 ~ dunif(0,5)
  phase3b.omega2 ~ dunif(0,5)
  
  vol.pdet ~ dunif(0,1)
  
  for(site in 1:SITES) {
    site.lambda[site] ~ dgamma(shape=5,scale=30)
    site.gamma[site]  ~ dunif(0,50)
    
    for(time in 1:TIMES) {
      gamma.matrix[site, time] <- site.gamma[site] +
        phase1.gamma * phase1.cov[site, time] + 
        phase2.gamma * phase2.cov[site, time] + 
        phase3a.gamma * phase3a.cov[site, time] +
        phase3b.gamma * phase3b.cov[site, time]
      
      omega1.matrix[site, time] <- 
        phase1.omega1 * phase1.cov[site, time] + 
        phase2.omega1 * phase2.cov[site, time] + 
        phase3a.omega1 * phase3a.cov[site, time] +
        phase3b.omega1 * phase3b.cov[site, time]
      
      omega2.matrix[site, time] <- 
        phase1.omega2 * phase1.cov[site, time] + 
        phase2.omega2 * phase2.cov[site, time] + 
        phase3a.omega2 * phase3a.cov[site, time] +
        phase3b.omega2 * phase3b.cov[site, time]
      
      p_det.matrix[site, time] <- vol.pdet * vol.cov[site, time]
      
      # define conjugate probabilities
      pd_p[site, time]     <- p_mor/p_det.matrix[site, time] # probability of death within observed cases (all deaths are observed), alpha = 1/p
      pr_1pd_p[site, time] <- p_rec/(1-pd_p[site, time])     # probability of recovery among observed cases after accounting for observed deaths
    }
  }
  
  p_rec2 <- p_rec/(1-p_mor) # probability of recovery among all cases after accounting for all deaths
  
  ## Likelihood
  for(site in 1:SITES) {
    for(time in 1:(TIMES-1)) {
      # Latent ARD State Process  
      D[site, time] ~ dbinom(p_mor, N[site, time]) 
      R[site, time] ~ dbinom(p_rec2, N[site, time] - D[site, time]) 
      A[site, time] <- N[site, time] - R[site, time] - D[site, time]
      
      # Detection Process (split active cases by remain active, recover, die)
      d[site, time] ~ dbinom(pd_p[site, time], m[site,time+1] + a[site, time])
      r[site, time] ~ dbinom(pr_1pd_p[site, time], m[site,time+1] + a[site, time] - d[site, time])    
    }  # 1:(TIMES-1)
    
    # Abundance and Detection
    N[site, 1] ~ dpois(site.lambda[site])
    for (time in 2:TIMES){
      # Latent Abundance Process
      N[site, time] <- A[site, time-1] + S[site, time-1] + G[site, time-1]
      S[site, time-1] ~ dpois(omega1.matrix[site, time]*(N[site, time-1]-a[site, time-1])*
                                (Npop[site] - N[site, time-1])/Npop[site] + # damping effect for population saturation
                                omega2.matrix[site, time]*a[site, time-1])
      G[site, time-1] ~ dpois(gamma.matrix[site, time-1])
      
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
                                              phase1.cov=BCphase1,
                                              phase2.cov=BCphase2,
                                              phase3a.cov=BCphase3a,
                                              phase3b.cov=BCphase3b,
                                              vol.cov=testvolnorm,
                                              Npop=pop_dat$Npop*1e6
                             ),
                             inits = model.inits)

# compile the model in C++
Cmodel.fit.nim <- compileNimble(model.fit.nim)

# Run MCMC
WAICparameters <- c("N", "G", "d", "S", "D", "R",
                    "vol.pdet",
                    "site.lambda", 
                    "site.gamma",
                    "phase1.gamma",
                    "phase2.gamma",
                    "phase3a.gamma",
                    "phase3b.gamma",
                    "phase1.omega1",
                    "phase2.omega1",
                    "phase3a.omega1",
                    "phase3b.omega1",
                    "phase1.omega2",
                    "phase2.omega2",
                    "phase3a.omega2",
                    "phase3b.omega2",
                    "p_rec", 
                    "p_mor")
parameters <- c("vol.pdet",
                "site.lambda", 
                "site.gamma",
                "phase1.gamma",
                "phase2.gamma",
                "phase3a.gamma",
                "phase3b.gamma", 
                "phase1.omega1",
                "phase2.omega1",
                "phase3a.omega1",
                "phase3b.omega1",
                "phase1.omega2",
                "phase2.omega2",
                "phase3a.omega2",
                "phase3b.omega2",
                "p_rec", 
                "p_mor")
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

WAIC_dat = data.frame(model = "M2", 
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

