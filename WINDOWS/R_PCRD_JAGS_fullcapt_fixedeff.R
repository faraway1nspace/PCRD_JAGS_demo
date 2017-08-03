# this is a demostration of using the Bayesian PCRD in R and JAGS, as used in: "Rankin RW, Nicholson K, Allen S, Kr\'{u}tzen M, Bejder L, Pollock KH. 2016. A full-capture Hierarchical Bayesian model of Pollock's Closed Robust Design and application to dolphins. Frontiers of Marine Science, in Press". There are 3 demonstrations. THIS FILE contains the "fixed effect" using the full-capture histories (i.e., recruitment model); see also the  "fixed effect" model conditioning on first-capture, and the  Hierarchical Bayesian model, using Gelman half-Student-t hyper-priors that try to explicitly shrink time-varying parameters to their time-invariant, plus random-effects for individual heterogeneity.

library(rjags)
source("R_PCRD_JAGS_SOURCE.R") # load some handy functions

# load the Useless Loop (Western gulf Shark Bay) Tursiops aduncus captures histories from "Nicholson, Bejder, Allen,Krützen, Pollock. 2012. Abundance, survival and temporary emigration of bottlenose dolphins (Tursiops sp.) off Useless Loop in the western gulf of Shark Bay, Western Australia. Marine and Freshwater Research 63:1059–1068."
# CH are stored as MARK file (.inp) -> convert to 3D array
MARK.file.name <- "mark_capture_histories.inp"
T2 <- c(5, 5, 10, 5,3) # number of secondary periods per primary period
T <- length(T2) # number of primary periods
capture.histories = import.mark.inp(MARK.file.name,T2,header=FALSE) # import MARK inp
Y.tt <- capture.histories[["Y.tt"]] # get the 3D array N x T2 x T
# Y.t <- capture.histories[["Y.t"]] # if you only want total counts per primary period (per individual)
N = nrow(Y.tt) # number of (observed) individuals 

# Jags code "fixed effect model"  phi(.) gamma'(.) gamma''(t) p(t,s)
# time-invariant phi, time-invariant gamma', time-varying gamma'', time- and session-varying detection probabilities.
# model's the full-capture history's using PXDA techniques (see parameters PSI and LAMBDA).
# jags model text
jags.txt.1 <-"model{
  # priors
  phi ~ dbeta(pr.phi[1],pr.phi[2]) # apparent survival probability ADJUST FOR ECOLOGICAL CONTEXT
  gamma1 ~ dbeta(pr.gamma1[1],pr.gamma1[2]) #  temp migration: probabilty stays out conditional on being out
  # loop through primary periods: parameters for detection probability pd_t, and gamma2_t
  for(t_ in 1:T){ # T primary periods...
     pd_a[t_] ~ dgamma(pr.pd[1],pr.pd[2]) # hyperprior on detect probs
     pd_b[t_] ~ dgamma(pr.pd[1],pr.pd[2]) # hyperprior on detect probs
     for(tt_ in 1:T2[t_]){ # loop through secondary periods...
        pd[tt_,t_] ~ dbeta(pd_a[t_],pd_b[t_]) # secondard-period-level detection probs
     }
     gamma2[t_] ~ dbeta(pr.gamma2[1],pr.gamma2[2]) # temp migration: prob migrates outside conditional on being inside
     # recruitment process from eigenvector decomposition
     lambda[1,t_] <- (1-gamma1)/(gamma2[t_]-gamma1+1) # recruitment ratio, or long-term prob of being inside
     lambda[2,t_] <- 1-lambda[1,t_] # 
     psi[t_] ~ dunif(0,1) # inclusion probability
     # trmat: transition matrix for Markovian latent-states
     # 1 =not yet in population;2=dead;3=offsite;4=onsite (only observable state)
     # transition are from the column --> rows
     # trmat[row,column,time] = [state at time=t_; state at time t_-1; time=t_]
     trmat[1,1,t_] <- 1-psi[t_] # excluded from pop
     trmat[2,1,t_] <- 0 # dead
     trmat[3,1,t_] <- psi[t_]*lambda[2,t_] # inclusion into pop, outside study are
     trmat[4,1,t_] <- psi[t_]*lambda[1,t_] # inclusion into pop, inside study area
     trmat[1,2,t_]<- 0
     trmat[2,2,t_]<- 1 # stay dead
     trmat[3,2,t_]<- 0
     trmat[4,2,t_]<- 0
     trmat[1,3,t_]<- 0
     trmat[2,3,t_]<- 1-phi # dies outside
     trmat[3,3,t_]<- gamma1*phi # stays outside | outside
     trmat[4,3,t_]<- (1-gamma1)*phi # reenters study area | outside
     trmat[1,4,t_]<- 0 # 
     trmat[2,4,t_]<- 1-phi # dies inside
     trmat[3,4,t_]<- gamma2[t_]*phi # leaves study area | inside
     trmat[4,4,t_]<- (1-gamma2[t_])*phi # stays inside | inside
  } # t_ 
  # likelihood: loop through M individuals, both real and pseudo-individuals
  for (i in 1:M){ 
    #draw latent state at primary period 1: 
    # ... by definition, everyone starts in z=1 (not-in-population) at time=0 
    z[i,1]~ dcat(trmat[1:4,1,1]) # first z strictly excluded from pop
    # likelihood for first primary period
    for(tt_ in 1:T2[1]){ # loop through secondary periods
         # Bernouilli process, conditional on z=4, otherwise no observation
         y[i,tt_,1] ~ dbern(pd[tt_,1]*equals(z[i,1],4)) 
    }
    alive_i[i,1] <- step(z[i,1]-3) # count if i is alive or not
    Nin_i[i,1] <- equals(z[i,1],4) # count if i is within study area
    # loop through primary periods after 1st primary periods
    for(t_ in 2:T){
      # state process: draw z(t_) conditional on z(t_-1)
      z[i,t_] ~ dcat(trmat[1:4, z[i,t_-1] , t_])
      # likelihood: loop through secondary period observations
      for(tt_ in 1:T2[t_]){
         # Bernouilli observation error, conditional on z=4
         y[i,tt_,t_] ~ dbern(pd[tt_,t_]*equals(z[i,t_],4)) 
      }
      # tally population size 
      alive_i[i,t_] <- step(z[i,t_]-3) # check alive or not
      Nin_i[i,t_] <- equals(z[i,t_],4) # count if i is within study area
    } # t_
  } # i
  # estimate population size per primary periods
  for(t_ in 1:T){
     alive[t_] <- sum(alive_i[,t_]) # number alive
     Nin[t_] <- sum(Nin_i[,t_]) # number in study area
  } # t_
}"
# save jags file to local disk
modname1 <- "JAGS_fullcapt_fixedeff.JAG"
sink(file=modname1)
cat(jags.txt.1,fill=TRUE)
sink()

# use Parameter-Expansion Data Augmentation method to model unseen individuals and full-capture histories
n.aug <- N # a N
Y.aug <- array(NA,c(N+n.aug, max(T2),T))
mm <- nrow(Y.aug)
dimnames(Y.aug)[[1]] <- c(row.names(Y.tt),paste0("aug",1:n.aug))
Y.aug[1:N,,]<-Y.tt # insert observed capture-histories
Y.aug[(N+1):mm,,]<-0*Y.tt[1:n.aug,,] # all-zero capture histories

# Data for jags: Capture-histories and Prior parameters
jags.data=list(y=Y.aug, T=T, T2=T2,M=mm,
              pr.pd=c(3.00,2.00), # Gamma hyperprior on hierarchical detection probibililities
              pr.gamma1=c(1.3,1.3), #  Beta prior on gamma'
              pr.gamma2=c(1,1), # Beta prior on gamma''
              pr.phi = c(1,1)  # Beta prior on phi, apparent surival                            
#              pr.phi = c(6,1) # useful for dolphins 
              )

# MCMC parameters
  nchains <- 3 # number of MCMC chains
  nadapt <-10000 # adaption phase
  nburn <- 40000 # discard these draws
  niter <- 150000 # length of MCMC chains
  nsamp <- 1000 # number of samples to take from each MCMC chains
  thin_ <- round(niter/nsamp) # only

# initialize start values for each MCMC chain
# user is expected to be able to initialize the univariate parameters;
# I provide 'generate.z.psi' function which imputs values of latent states 'z' (forwards-backwards algorithm) as well as values of 'Psi' (full-capture conditioning)
init.func = function(){
    RET=list(
        phi = rbeta(1,60,10), # user to specify
        gamma1=rbeta(1,20,20), # user to specify
        gamma2=rbeta(T,20,20),# user to specify
        pd_a =rgamma(T,30,20),# user to specify
        pd_b = rgamma(T,30,20),# user to specify
        pd = do.call(function(T2,T){ # user to specify
            pd = matrix(NA,max(T2),T);
            for(t_ in 1:T){ pd[1:T2[t_],t_]<-rbeta(T2[t_],12,65)}
            return(pd)}, args=list(T2=T2,T=T))
    )
    RET=c(RET, generate.z.psi(y=Y.aug,T2=T2,first.capture=FALSE,z.priors = list(phi.beta=c(shape1=30,shape2=5),g1.beta=c(shape1=20,shape2=20),g2.beta=c(shape1=20,shape2=20), pd.beta=c(shape1=12,shape2=65)))) # handy function to initialize latent states z
    return(RET)
}

# initialize random variables 
if(nchains == 1){
    jags.inits = init.func()
} else {
    jags.inits = lapply(1:nchains, function(x) init.func())
}

# RUN JAGS MODEL AND SAMPLE FROM POSTERIOR
  m1 <- jags.model(file=modname1,data=jags.data,inits=jags.inits,n.chains=nchains,n.adapt=nadapt)
# burn-in phase
  update(m1,nburn) # discard nburn iterations (ensure that chains reach equilibium states)
# which variables to summarize
  variable.names <- c("phi","gamma1",paste0("gamma2[",1:T,"]"),unlist(mapply(as.list(1:T),sapply(X=T2,FUN=seq,from=1,by=1),FUN=function(x,y) paste0("pd[",y,",",x,"]"))),paste0("psi[",1:T,"]"),paste0("alive[",1:T,"]"),paste0("Nin[",1:T,"]"))
# sample from posterior
  samp <- coda.samples(m1,variable.names=variable.names,n.iter=niter,thin=thin_) 

  summary(samp) # summarize output
  plot(samp, ask =TRUE) # inspect chains for adequate mixing and convergence
  gelman.diag(samp) # inspect chains for convergence (statistics ~ 1 and < 1.1)
# DONE PART 1
    
