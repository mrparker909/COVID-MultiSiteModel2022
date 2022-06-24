# In this code we fit the COVID multi site model to the Canada Wide data
library(tidyverse)
library(nimble)

# Load Data
#############################################
# load data:
dat1 = read.csv("Data/covid-case-counts-canada.csv")
dat1$WeekNum = dat1$WeekNum - 4 # shift to starting week 1 = 2020-04-23
dat <- dat1 %>% filter(Province %in% c("Alberta",
                                       "British Columbia",
                                       "Manitoba",
                                       "New Brunswick",
                                       "Newfoundland and Labrador",
                                       "Nova Scotia",
                                       "Ontario",
                                       "Prince Edward Island",
                                       "Quebec",
                                       "Saskatchewan",
                                       "Northwest Territories",
                                       "Nunavut",
                                       "Yukon"))
unique(dat$Province)
# [1] "Alberta"                   "British Columbia"          "Manitoba"                  "New Brunswick"            
# [5] "Newfoundland and Labrador" "Northwest Territories"     "Nova Scotia"               "Nunavut"                  
# [9] "Ontario"                   "Prince Edward Island"      "Quebec"                    "Saskatchewan"             
# [13] "Yukon"
pop_dat = data.frame(Province = dat$Province %>% unique(), Npop=c(4464170,
                                                                  5249635,
                                                                  1386333,
                                                                  794300,
                                                                  521758,
                                                                  45515,
                                                                  998832,
                                                                  39589,
                                                                  14915270,
                                                                  165936,
                                                                  8631147,
                                                                  1180867,
                                                                  43095
                                                                  )) # approximate population size for each province (2021)

# > pop_dat
# Province     Npop
# 1                    Alberta  4464170
# 2           British Columbia  5249635
# 3                   Manitoba  1386333
# 4              New Brunswick   794300
# 5  Newfoundland and Labrador   521758
# 6      Northwest Territories    45515
# 7                Nova Scotia   998832
# 8                    Nunavut    39589
# 9                    Ontario 14915270
# 10      Prince Edward Island   165936
# 11                    Quebec  8631147
# 12              Saskatchewan  1180867
# 13                     Yukon    43095
USED_SITES = 1:13

temp1 <- reshape(dat %>% select(Province, WeekStart, nit), idvar = "Province", timevar = "WeekStart", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",1:(ncol(temp2)))
nit_temp = data.matrix(temp2)

temp1 <- reshape(dat %>% select(Province, WeekStart, rit), idvar = "Province", timevar = "WeekStart", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",1:(ncol(temp2)))
rit_temp = data.matrix(temp2)

temp1 <- reshape(dat %>% select(Province, WeekStart, dit), idvar = "Province", timevar = "WeekStart", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",1:(ncol(temp2)))
Dit_temp = data.matrix(temp2)

temp1 <- reshape(dat %>% select(Province, WeekStart, ait), idvar = "Province", timevar = "WeekStart", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",1:(ncol(temp2)))
ait_temp = data.matrix(temp2)

maxT = 90
nit = nit_temp[,1:maxT]
rit = rit_temp[,1:maxT][,-c(1)]
Dit = Dit_temp[,1:maxT][,-c(1)]
dit = Dit
ait = ait_temp[,1:maxT]
n_dim <- dim(nit) # SITES = n_dim[1], TIMES = n_dim[2]


##############################################################
# correcting data discrepancies:
cbind(which(nit+ait-cbind(0,rit)-cbind(0,dit) < 0, arr.ind = T), discrepancy=(nit+ait-cbind(0,rit)-cbind(0,dit))[which(nit+ait-cbind(0,rit)-cbind(0,dit) < 0, arr.ind = T)])
#                           row col discrepancy
# Yukon                      13   2          -1
# Nova Scotia                 7   4          -5
# Newfoundland and Labrador   5   5          -1
# Newfoundland and Labrador   5   9          -1
# Yukon                      13  17          -2
# Prince Edward Island       10  28          -1
# Yukon                      13  39          -4
# Nunavut                     8  48          -9
# Northwest Territories       6  49          -1
# Yukon                      13  53          -1
# Northwest Territories       6  57          -2
# Prince Edward Island       10  60          -1
# Nunavut                     8  79          -3

# Here we are adding the discrepancy to nit, making the assumption that the extra recoveries were not properly reported as active cases
nit[which(nit+ait-cbind(0,rit)-cbind(0,dit) < 0, arr.ind = T)] = nit[which(nit+ait-cbind(0,rit)-cbind(0,dit) < 0, arr.ind = T)] - (nit+ait-cbind(0,rit)-cbind(0,dit))[which(nit+ait-cbind(0,rit)-cbind(0,dit) < 0, arr.ind = T)]
##############################################################

## First Omicron cases confirmed in Canada: November 28, 2021 (Week 84)
omicron = matrix(0, nrow=n_dim[1], ncol=n_dim[2])
omicron[,84:n_dim[2]] = 1

## First Delta cases confirmed in Canada: April 4, 2021 (Week 54)
delta = matrix(0, nrow=n_dim[1], ncol=n_dim[2])
delta[,54:83] = 1

temp1 <- reshape(dat %>% select(Province, WeekStart, TestVol), idvar = "Province", timevar = "WeekStart", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",1:(ncol(temp2)))
testvol_temp = data.matrix(temp2)

testvolit = testvol_temp[,1:maxT]
testvolnorm = t(apply(testvolit, 1, function(x)(x)/(sd(x, na.rm = T)))) # here we standardize testing volumes


temp1 <- reshape(dat %>% select(Province, WeekStart, Vaccine1Dose), idvar = "Province", timevar = "WeekStart", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",1:(ncol(temp2)))
vacc1_temp = data.matrix(temp2)

vacc1it = vacc1_temp[,1:maxT]
vacc1pct = vacc1it/pop_dat$Npop


temp1 <- reshape(dat %>% select(Province, WeekStart, Vaccine2Dose), idvar = "Province", timevar = "WeekStart", direction = "wide")
temp2 <- temp1[USED_SITES,-1]
rownames(temp2) <- temp1[USED_SITES,1]
colnames(temp2) <- paste0("W",1:(ncol(temp2)))
vacc2_temp = data.matrix(temp2)

vacc2it = vacc2_temp[,1:maxT]
vacc2pct = vacc2it/pop_dat$Npop


# MCMC parameters
iterations <- 300000
burnin <- 200000
chains <- 1
thinning <- 1

#############################################
# initialize model parameters
Ninit <- (ait+nit+cbind(rit,0))
Rinit <- 2*cbind(rit,0)
Ginit <- Rinit
Sinit <- Ninit-Ginit
Ninit[,-1] <- NA
Rinit[,-1] <- NA
Ginit[,-1] <- NA
Sinit[,-1] <- NA

any(Ninit < 0, na.rm = T)
any(Sinit < 0, na.rm = T)
any(Rinit < 0, na.rm = T)
any(Ginit < 0, na.rm = T)

model.inits <- list(prec      = 0.24, 
                    prec.vac1 = 0,
                    prec.vac2 = 0,
                    prec.omicron = 0,
                    prec.delta = 0,
                    pmor      = 0.005,
                    pmor.vac1 = 0,
                    pmor.vac2 = 0,
                    pmor.omicron = 0,
                    pmor.delta = 0,
                    pdet = rep(0.34, times=n_dim[1]),
                    pdet.vol   = 0,
                    lambda     = pop_dat$Npop/1e4,
                    omega1 = rep(1, times=n_dim[1]),
                    omega2 = rep(1, times=n_dim[1]),
                    omega1.omicron = 0,
                    omega1.delta = 0,
                    omega1.vac1 = 0,
                    omega1.vac2 = 0,
                    omega2.omicron = 0,
                    omega2.delta = 0,
                    omega2.vac1 = 0,
                    omega2.vac2 = 0,
                    gamma  = rep(1, times=n_dim[1]),
                    gamma.omicron = 0,
                    gamma.delta = 0,
                    gamma.vac1 = 0,
                    gamma.vac2 = 0,
                    prec.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    pmor.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    p_det.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    gamma.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    omega1.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    omega1.matrix = matrix(0, nrow=n_dim[1], ncol=n_dim[2]),
                    N = Ninit,
                    R = Rinit,
                    S = Sinit,
                    G = Ginit
)

nimModel <- nimbleCode({ 
  ## Priors and constraints
  prec ~ dunif(0, 1)               
  prec.vac1 ~ dunif(-1, 1)
  prec.vac2 ~ dunif(-1, 1)
  prec.delta ~ dunif(-1, 1)
  prec.omicron ~ dunif(-1, 1)
  
  pmor ~ dunif(0, prec)        
  pmor.vac1 ~ dunif(-1, 1)
  pmor.vac2 ~ dunif(-1, 1)
  pmor.delta ~ dunif(-1, 1)
  pmor.omicron ~ dunif(-1, 1)
  
  gamma.vac1 ~ dunif(-1,1)
  gamma.vac2 ~ dunif(-1,1)
  gamma.delta   ~ dunif(-10000,10000)
  gamma.omicron ~ dunif(-10000,10000)
  
  omega1.vac1 ~ dunif(-5,5)
  omega1.vac2 ~ dunif(-5,5)
  omega1.delta   ~ dunif(-5,5)
  omega1.omicron ~ dunif(-5,5)
  
  omega2.vac1 ~ dunif(-5,5)
  omega2.vac2 ~ dunif(-5,5)
  omega2.delta   ~ dunif(-5,5)
  omega2.omicron ~ dunif(-5,5)
  
  pdet.vol ~ dunif(0,5)
  
  for(site in 1:SITES) {
    lambda[site] ~ dunif(0,1e6)
    pdet[site]   ~ dunif(0,1)
    
    gamma[site] ~ dunif(0,1e5)
    
    omega1[site] ~ dunif(0,5)
    omega2[site] ~ dunif(0,5)
    
    for(time in 1:TIMES) {
      temp.prec[site, time] <- prec + prec.vac1 * vac1[site, time] + prec.vac2 * vac2[site, time] + prec.omicron * omicron[site, time] + prec.delta * delta[site, time]
      prec.matrix[site, time] <- min(1, max(0, temp.prec[site, time]))
      
      temp.pmor[site, time] <- pmor + pmor.vac1 * vac1[site, time] + pmor.vac2 * vac2[site, time] + pmor.omicron * omicron[site, time] + pmor.delta * delta[site, time]
      pmor.matrix[site, time] <- min(max(0, temp.pmor[site,time]),1-prec.matrix[site, time])
      
      gamma.matrix[site, time] <- (gamma[site] + gamma.vac1 * vac1[site, time] + gamma.vac2 * vac2[site, time] + gamma.omicron * omicron[site, time] + gamma.delta * delta[site, time])
      
      temp.omega1[site, time] <- omega1[site] + omega1.vac1 * vac1[site, time] + omega1.vac2 * vac2[site, time] + omega1.omicron * omicron[site, time] + omega1.delta * delta[site, time]
      omega1.matrix[site, time] <- max(0, temp.omega1[site,time])
      
      temp.omega2[site, time] <- omega2[site] + omega2.vac1 * vac1[site, time] + omega2.vac2 * vac2[site, time] + omega2.omicron * omicron[site, time] + omega2.delta * delta[site, time]
      omega2.matrix[site, time] <- max(0, temp.omega2[site,time])
      
      # p = B0 + B1 * vol
      p_det.matrix[site, time] <- pdet[site] + pdet.vol * vol[site,time]
      
      # define conjugate probabilities
      pd_p[site, time]     <- pmor.matrix[site, time]/p_det.matrix[site, time]   # probability of death within observed cases (all deaths are observed), alpha = 1/p
      pr_1pd_p[site, time] <- prec.matrix[site, time]/(1-pd_p[site, time])       # probability of recovery among observed cases after accounting for observed deaths
      prec2.matrix[site, time]   <- prec.matrix[site, time]/(1-pmor.matrix[site, time]) # probability of recovery among all cases after accounting for all deaths
    }
  }
  
  ## Likelihood
  for(site in 1:SITES) {
    for(time in 1:(TIMES-1)) {
      # Latent ARD State Process  
      D[site, time] ~ dbinom(pmor.matrix[site, time], N[site, time]) 
      R[site, time] ~ dbinom(prec2.matrix[site, time], N[site, time] - D[site, time]) 
      A[site, time] <- N[site, time] - R[site, time] - D[site, time]
      
      # Detection Process (split active cases by remain active, recover, die)
      d[site, time] ~ dbinom(pd_p[site, time], m[site,time+1] + a[site, time])
      # NOTE: for the canada data, ait data is t+1 compared to BC data, which was t
      r[site, time] ~ dbinom(pr_1pd_p[site, time], m[site,time+1] + a[site, time+1] - d[site, time])    
    }  # 1:(TIMES-1)
    
    # Abundance and Detection
    N[site, 1] ~ dpois(lambda[site])
    for (time in 2:TIMES){
      # Latent Abundance Process
      N[site, time] <- A[site, time-1] + S[site, time-1] + G[site, time-1]
      S[site, time-1] ~ dpois(omega1.matrix[site, time]*(N[site, time-1]-a[site, time-1])*
                                (Npop[site] - N[site, time-1])/Npop[site] + # damping effect for population saturation
                                omega2.matrix[site, time]*a[site, time-1])
      G[site, time-1] ~ dpois(gamma.matrix[site, time])
      
      # Detection Process
      n[site, time] ~ dbinom(p_det.matrix[site, time], N[site, time]-a[site, time-1]) 
    } # TIMES
    n[site, 1] ~ dbinom(p_det.matrix[site, 1], N[site, 1])  # args: p, N 
  }
})

# create the model graph
model.fit.nim <- nimbleModel(code = nimModel,
                             name = "nimble model", 
                             data = list(n=nit, m=nit, r=rit, D=Dit, d=Dit, a=ait),
                             constants = list(SITES=n_dim[1], 
                                              TIMES=n_dim[2],
                                              Npop=pop_dat$Npop,
                                              vol=testvolnorm,
                                              vac1=vacc1pct,
                                              vac2=vacc2pct,
                                              omicron=omicron,
                                              delta=delta
                             ),
                             inits = model.inits)

# compile the model in C++
Cmodel.fit.nim <- compileNimble(model.fit.nim)

# Run MCMC
WAICparameters <- c("N", "G", "d", "S", "D", "R", "r",
                    "prec.matrix",
                    "pmor.matrix",
                    "p_det.matrix",
                    "gamma.matrix",
                    "omega1.matrix",
                    "omega2.matrix",
                    "prec", 
                    "prec.vac1",
                    "prec.vac2", 
                    "prec.omicron",
                    "prec.delta",
                    "pmor",
                    "pmor.vac1",
                    "pmor.vac2", 
                    "pmor.omicron",
                    "pmor.delta",
                    "pdet",
                    "pdet.vol",
                    "lambda",
                    "omega1",
                    "omega2",
                    "omega1.omicron",
                    "omega1.delta",
                    "omega1.vac1",
                    "omega1.vac2",
                    "omega2.omicron",
                    "omega2.delta",
                    "omega2.vac1",
                    "omega2.vac2",
                    "gamma",
                    "gamma.omicron",
                    "gamma.delta",
                    "gamma.vac1",
                    "gamma.vac2")
parameters <- c("prec", 
                "prec.vac1",
                "prec.vac2", 
                "prec.omicron",
                "prec.delta",
                "pmor",
                "pmor.vac1",
                "pmor.vac2", 
                "pmor.omicron",
                "pmor.delta",
                "pdet",
                "pdet.vol",
                "lambda",
                "omega1",
                "omega2",
                "omega1.omicron",
                "omega1.delta",
                "omega1.vac1",
                "omega1.vac2",
                "omega2.omicron",
                "omega2.delta",
                "omega2.vac1",
                "omega2.vac2",
                "gamma",
                "gamma.omicron",
                "gamma.delta",
                "gamma.vac1",
                "gamma.vac2")

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

WAIC_dat = data.frame(model = "canada_full_model", 
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

