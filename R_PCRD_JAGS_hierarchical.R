# this is a demostration of using the Bayesian PCRD in R and JAGS, as used in: "Rankin RW, Nicholson K, Allen S, Kr\'{u}tzen M, Bejder L, Pollock KH. 2016. A full-capture Hierarchical Bayesian model of Pollock's Closed Robust Design and application to dolphins. Frontiers of Marine Science, in Press". There are 3 demonstrations. THIS FILE contains the the Hierarchical Bayesian model, using Gelman half-Student-t hyper-priors that try to explicitly shrink time-varying parameters to their time-invariant, plus random-effects for individual heterogeneity. See other files for the "fixed effect" version using the full-capture histories (i.e., recruitment model); and the  "fixed effect" version which conditions on first-capture.
# Users should focus on the hyperpriors pr.taug2,pr.taug1,pr.tauphi,pr.taupd, prtaupd2nd, prtaupdi as these control the amount of shrinkage between a fully-time-varying specification to a time-invariant Theta(.) specification. I.e., its only the dispersion parameters (sigma.pdmu,sigma.g1,sigma.g2,sigma.phi,sigma.eps) that receive strongly informative (shrinkage-inducing) hyperpriors. The priors on the 'location' of parameters (pd.mu,g1.mu,g2.mu,phi.mu) are non-informative. 

library(rjags)
library(boot)
source("R_PCRD_JAGS_SOURCE.R") # load some handy functions

# load the Useless Loop (Western gulf Shark Bay) Tursiops aduncus captures histories from "Nicholson, Bejder, Allen,Kr\'{u}tzen, Pollock. 2012. Abundance, survival and temporary emigration of bottlenose dolphins (Tursiops sp.) off Useless Loop in the western gulf of Shark Bay, Western Australia. Marine and Freshwater Research 63:1059-1068."
# CH are stored as MARK file (.inp) -> convert to 3D array
MARK.file.name <- "mark_capture_histories.inp"
T2 <- c(5, 5, 10, 5,3) # number of secondary periods per primary period
T <- length(T2) # number of primary periods
capture.histories = import.mark.inp(MARK.file.name,T2,header=FALSE) # import MARK inp
Y.tt <- capture.histories[["Y.tt"]] # get the 3D array N x T2 x T
# Y.t <- capture.histories[["Y.t"]] # if you only want total counts per primary period (per individual)
N = nrow(Y.tt) # number of (observed) individuals 

# Jags code hierarchical Bayesian PCRD (including individual heterogeneity)
# time-vary phi, time-invary gamma', time-varying gamma'', time- and session-varying detection probabilities: all shrunk towards time-invariant theta(dot) using hyperpriors on the dispersion parameters (sigma)
# model's the full-capture history's using PXDA techniques (see parameters PSI and LAMBDA).
# jags model text
jags.txt.3 <-"model{
  #hyperpriors
  phi.mu ~ dunif(pr.phiunif[1],pr.phiunif[2]) #mean survival with a Uniform prior
  sigma.phi ~ dt(0,pr.tauphi[1],pr.tauphi[2]) T(0,) #mean survival dispersion, half-t hyperprior
  g1.mu ~ dnorm(0,pr.g1mu) #mean gamma1, temp migration out | out
  sigma.g1~dt(0,pr.taug1[1],pr.taug1[2]) T(0,) #mean gamma1 dispersion, half-t hyperprior
  g2.mu ~ dnorm(0,pr.g2mu) #mean gamma2,  temp migration out | in
  sigma.g2~dt(0,pr.taug2[1],pr.taug2[2]) T(0,) #mean gamma2 dispersion, half-t hyperprior
  pd.mu ~ dnorm(0,pr.pdmu) #mean detection prob, overall 
  sigma.pdmu~dt(0,pr.taupdmu[1],pr.taupdmu[2]) T(0,) #primary period detection prob dispersion
  sigma.pd2nd~dt(0,pr.taupd2nd[1],pr.taupd2nd[2]) T(0,) #secondary periods detection prob dispersion 
  sigma.eps ~ dt(0,pr.taueps[1],pr.taueps[2]) T(0,) #individual detection prob dispersion
  #time-variant parameters 
  for(t_ in 1:T){ #loop through primary periods
    pd_mu[t_]~dnorm(pd.mu,pow(sigma.pdmu,-2)) #primary period mean detaction prob (logit)
    lgamma1[t_]~dnorm(g1.mu,pow(sigma.g1,-2)) #prob migrate out|out (logit)
    gamma1[t_] <- ilogit(lgamma1[t_]) #prob migrate out|out (probability)
    lgamma2[t_]~dnorm(g2.mu,pow(sigma.g2,-2)) #prob migrate out|in (logit)
    gamma2[t_] <- ilogit(lgamma2[t_]) #prob migrate out|in (probability)
    #RECRUITMENT: psi is the 'inclusion probability' and lambda is the 'recruitment ratio'
    psi[t_]~dunif(0,1) #inclusion probability
    lambda[t_] <- (1-gamma1[t_])/(gamma2[t_]-gamma1[t_]+1) #long-term probability inside study area
    #NOTE, lambda could also be a random variable with a beta prior
    #secondary-period detection probabilities
    for(tt_ in 1:T2[t_]){ #loop through secondary periods
      pd[tt_,t_] ~ dnorm(pd_mu[t_],pow(sigma.pd2nd,-2))
    } #tt_
  } #t_
  #first state transition (pure nusance; strictly from outside-pop to part of marked-pop)
  trmat0[1] <- (1-psi[1]) #remains not-yet-in-pop
  trmat0[2] <- 0
  trmat0[3] <- psi[1]*(1-lambda[1]) #inclusion into pop, goes outside study are
  trmat0[4] <- psi[1]*lambda[1] #inclusion into pop, goes inside study
  #state transitions (2:T)
  for(t_ in 1:(T-1)){
    lphi[t_]~dnorm(log(phi.mu/(1-phi.mu)), pow(sigma.phi,-2)) #survival prob (logit)
    phi[t_]<-ilogit(lphi[t_])
    #state transitions 
     #trmat: transition matrix for Markovian latent-states
     #1 =not yet in population; 2=dead;3=offsite;4=onsite (only observable state)
     #transition are from the column --> rows
     #trmat[row,column,time] = [state at time=t_; state at time t_-1; time=t_]
     #notice that the primary periods are offset by 1 (because we already dealt with T=1)
     trmat[1,1,t_]<- 1-psi[t_+1] #excluded from pop
     trmat[2,1,t_] <- 0 #dead
     trmat[3,1,t_] <- psi[t_+1]*(1-lambda[t_+1]) #inclusion into pop,outside study
     trmat[4,1,t_] <- psi[t_+1]*lambda[t_+1] #inclusion into pop,inside study
     trmat[1,2,t_]<- 0
     trmat[2,2,t_]<- 1 #stay dead
     trmat[3,2,t_]<- 0
     trmat[4,2,t_]<- 0
     trmat[1,3,t_]<- 0
     trmat[2,3,t_]<- 1-phi[t_] #dies outside
     trmat[3,3,t_]<- gamma1[t_+1]*phi[t_] #stays outside | outside
     trmat[4,3,t_]<- (1-gamma1[t_+1])*phi[t_] #reenters study area | outside
     trmat[1,4,t_]<- 0 #
     trmat[2,4,t_]<- 1-phi[t_] #dies inside
     trmat[3,4,t_]<- gamma2[t_+1]*phi[t_] #leaves study area | inside
     trmat[4,4,t_]<- (1-gamma2[t_+1])*phi[t_] #stays inside | inside
  } #t_
  #loop through M individuals
  for (i in 1:M){
    #state transitions and likelihood for the first primary period
    z[i,1]~ dcat(trmat0) #z at time 0 is strictly 'not-yet-in-pop'
    alive_i[i,1] <- step(z[i,1]-3) #count if i is alive or not
    Nin_i[i,1] <- equals(z[i,1],4) #count if i is within study area
    eps_i[i] ~ dnorm(0,pow(sigma.eps,-2)) #random effects at individual levels
    #Observation error y[i,tt_,t_] ~ Bernoulli conditional on being inside z=4
    for(tt_ in 1:T2[1]){  #loop through secondary periods
       y[i,tt_,1] ~ dbern(equals(z[i,1],4)/(1+exp(-pd[tt_,1]-eps_i[i]))) #inverse-logit transform
    }
    #state transition and likelihood for primary periods 2:T
    for(t_ in 2:T){ 
      #State process: draw z(t_) conditional on  z(t_-1)
      z[i,t_] ~ dcat(trmat[1:4, z[i,t_-1] , t_-1])
      #Observation error y[i,tt_,t_] ~ Bernoulli condition on being inside z=4
      for(tt_ in 1:T2[t_]){ #loop through secondary periods
        y[i,tt_,t_] ~ dbern(equals(z[i,t_],4)/(1+exp(-pd[tt_,t_]-eps_i[i]))) #inverse-logit transform
      }
      #check whether i individual is alive and inside
      alive_i[i,t_] <- step(z[i,t_]-3) #check alive
      Nin_i[i,t_] <- equals(z[i,t_],4) #count if i is within study area
    } #t_
  } #i
  #tally population size
  for(t_ in 1:T){
     alive[t_] <- sum(alive_i[,t_]) #number alive
     Nin[t_] <- sum(Nin_i[,t_]) #number in study area
  } #t_
} #end model
"
# save jags text/code in object jags.txt.3 to local disk as a file to import into external program JAGS
modname3 <- "JAGS_HierBayes.JAG"
sink(file=modname3) # this creates a new file called 'JAGS_HierBayes.JAG'
cat(jags.txt.3,fill=TRUE) # this places the text from object jags.txt.3 into file JAGS_HierBayes.JAG
sink() # this cloes the connection to file 'JAGS_HierBayes.JAG'

# use Parameter-Expansion Data Augmentation method to model unseen individuals and full-capture histories
n.aug <- 2*N # a N
Y.aug <- array(NA,c(N+n.aug, max(T2),T))
mm <- nrow(Y.aug)
dimnames(Y.aug)[[1]] <- c(row.names(Y.tt),paste0("aug",1:n.aug))
Y.aug[1:N,,]<-Y.tt # insert observed capture-histories
Y.aug[(N+1):mm,,] <- 0*Y.tt[(1+(0:(n.aug-1)) %% N),,]
# the previous line may seem odd: really it is just increasing the data with an augmented amount of data. The modulus operation part is just to ensure that any NAs in the original data (for non-sampled secondary periods) are included in the augmented array as well

# Data for jags: Capture-histories and Hyper-Prior parameters
# Users should familiarize themselves with the half-student-t distribution for controlling the dispersion hyperparameters on sigma.g1, sigma.g2, sigma.phi, sigma.pdmu,sigma.pd2mu,sigma.eps. In particular, the 'degrees of freedeom' hyperparameter 'nu' in half-t(0,tau,nu). nu=1 is the half-cauchy, and likely too much variation for a logit-Normal distribution. nu=2 is recommended. Tau should be << 1 to ensure that values of theta don't cluster on 0 and 1 on the probability scale. We suggest tau = 0.3. 
jags.data=list(y=Y.aug, T=T, T2=T2,M=mm,
#    pr.phiunif= c(0.8,0.999999999), # used in dolphin study
    pr.phiunif= c(0.0001,0.999999999), # uniform prior on mean surival phi               
    pr.tauphi =c(tau=1/(0.2)^2,nu=13), # half-t hyperprior on sigma.phi, controlling the dispersion of logit(phi_t) around phi.mu
    pr.g1mu=c(tau=1/(1.1)^2), # tau prior on mean logit(gamma1): constrain to be uniform on the prob scale
    pr.taug1=c(tau=1/(0.175^2),nu=4), # half t hyperprior on sigma.g1, the dispersion of logit(gamma1_t) around g1.mu 
    pr.g2mu=c(tau=1/2.41), # prior on mean logit(gamma2): constraining to be uniform on prob scale
    pr.taug2=c(tau=1/(0.3^2),nu=2), # half t hyperprior on sigma.g2, dispersion of logit(gamma2_t) around g2.mu
    pr.pdmu=c(tau=1/2.41),# expected variance of pd.mu; flat on prob scale
    pr.taupdmu=c(tau=1/(0.3^2),nu=2), # half t hyperprior on sigma.pdmu, dispersion of logit(pd_mu) (mean detection prob per primary period) around pd.mu
    pr.taupd2nd=c(tau=1/(0.11^2),nu=2), # half t hyperprior on sigma.pd2nd, dispersion of logit(pd_{t,s}) at secondary periods level around their primary-period mean pd_mu[t]
   pr.taueps=c(tau=1/(0.11^2),nu=2) # half t prior on sigma.eps, the dispersion of eps_i around logit(pd_{t,s}), representing individual heterogeneity in detection probabilities
  )
               
# MCMC parameters
  nchains <- 3 # number of MCMC chains
  nadapt <-40000 # adaption phase
  nburn <- 40000 # discard these draws
  niter <- 300000 # length of MCMC chains
  nsamp <- 1000 # number of samples to take from each MCMC chains
  thin_ <- round(niter/nsamp) # only

# initialize start values for each MCMC chain
# user is expected to be able to initialize the univariate parameters;
# I provide 'generate.z.psi' function which imputs values of latent states 'z' (forwards-backwards algorithm) as well as values of 'Psi' (full-capture conditioning)
init.func = function(){
    RET=list(
     phi.mu=runif(1,0.87,0.96), 
     sigma.phi=runif(1,0.0001,0.001)^0.5,
     g1.mu=logit(runif(1,0.2,0.8)),
     sigma.g1=runif(1,0.0001,0.001)^0.5,
     g2.mu=logit(runif(1,0.2,0.8)),
     sigma.g2=runif(1,0.01,0.2)^0.5,
     pd.mu =logit(runif(1,0.1,0.5)),
     sigma.pdmu=runif(1,0.01,0.025)^0.5,
     sigma.pd2nd=runif(1,0.005,0.015)^0.5,
     sigma.eps = runif(1,0.0005,0.0023)^0.5
    );
    RET=c(RET,list(
      lphi=rnorm(T-1,RET$phi.mu,sd=(RET$sigma.phi)),
      lgamma1=rnorm(T,RET$g1.mu,sd=(RET$sigma.g1)),
      lgamma2=rnorm(T,RET$g2.mu,sd=(RET$sigma.g2)),
      pd_mu=rnorm(T,RET$pd.mu,sd=(RET$sigma.pdmu)),
      eps_i=rnorm(mm,0,sd=RET$sigma.eps))
      );
    pd = matrix(NA,max(T2),T);    
    for(t_ in 1:T){ pd[1:T2[t_],t_]<-rnorm(T2[t_],RET$pd_mu[t_],(RET$sigma.pd2nd))};
    RET=c(RET, list(pd=pd));
    RET=c(RET, generate.z.psi(y=Y.aug,T2=T2,first.capture=FALSE,z.priors = list(phi.beta=c(shape1=30,shape2=5),g1.beta=c(shape1=20,shape2=20),g2.beta=c(shape1=20,shape2=20), pd.beta=c(shape1=12,shape2=65)))); # handy function to initialize latent state
    return(RET)
}

# initialize random variables 
if(nchains == 1){
    jags.inits = init.func()
} else {
    jags.inits = lapply(1:nchains, function(x) init.func())
}

# RUN JAGS MODEL AND SAMPLE FROM POSTERIOR
  m3 <- jags.model(file=modname3,data=jags.data,inits=jags.inits,n.chains=nchains,n.adapt=nadapt)
# burn-in phase
  update(m3,nburn) # discard nburn iterations (ensure that chains reach equilibium states)
# which variables to summarize
  variable.names <- c("sigma.phi","sigma.g1","sigma.g2","sigma.pdmu","sigma.pd2nd","sigma.eps",paste0("phi[",1:(T-1),"]"),paste0("gamma1[",1:T,"]"),paste0("gamma2[",1:T,"]"),unlist(mapply(as.list(1:T),sapply(X=T2,FUN=seq,from=1,by=1),FUN=function(x,y) paste0("pd[",y,",",x,"]"))),paste0("psi[",1:T,"]"),paste0("alive[",1:T,"]"),paste0("Nin[",1:T,"]"))
# sample from posterior
  samp <- coda.samples(m3,variable.names=variable.names,n.iter=niter,thin=thin_) 

  summary(samp) # summarize output
  plot(samp, ask =TRUE) # inspect chains for adequate mixing and convergence
  gelman.diag(samp) # inspect chains for convergence (statistics ~ 1 and < 1.1)
# DONE PART 1
    
