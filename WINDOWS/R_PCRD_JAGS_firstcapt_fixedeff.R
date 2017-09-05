# this is a demostration of using the Bayesian PCRD in R and JAGS, as used in: "Rankin RW, Nicholson K, Allen S, Kr\'{u}tzen M, Bejder L, Pollock KH. 2016. A full-capture Hierarchical Bayesian model of Pollock's Closed Robust Design and application to dolphins. Frontiers of Marine Science, in Press". There are 3 demonstrations. THIS FILE contains the "fixed effect" while conditiong on first-capture (i.e., non-recruitment model); see also the  "fixed effect" model with full-capturehistory modelling, and the  Hierarchical Bayesian model, using Gelman half-Student-t hyper-priors that try to explicitly shrink time-varying parameters to their time-invariant, plus random-effects for individual heterogeneity.

# WARNING: this JAGS script uses the WinBUGS 'zero's trick' to pass a non-standard likelihood to the JAGS sampler. Please read up about the zero's trick, and importantly, set the parameter C so that the Poisson density never sees a negative number

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
Y.t <- capture.histories[["Y.t"]] # if you only want total counts per primary period (per individual)
N = nrow(Y.tt) # number of (observed) individuals 

# Jags code "fixed effect model"  phi(.) gamma'(.) gamma''(t) p(t,s), conditioning on first capture
# time-invariant phi, time-invariant gamma', time-varying gamma'', time- and session-varying detection probabilities.
jags.txt.2 <-"model{
  # priors
  phi ~ dbeta(pr.phi[1],pr.phi[2]) #apparent survival probability 
  gamma1 ~ dbeta(pr.gamma1[1],pr.gamma1[2]) #temporary migration: probabilty stays out conditional on being out
  for(t_ in 1:T){ #T primary periods...
     pd_a[t_] ~ dgamma(pr.pd[1],pr.pd[2]) #hyperprior on detect probs
     pd_b[t_] ~ dgamma(pr.pd[1],pr.pd[2]) #hyperprior on  detect probs
     for(tt_ in 1:T2[t_]){ #loop through secondary periods...
        pd[tt_,t_] ~ dbeta(pd_a[t_],pd_b[t_]) #secondard-period-level detection probs
     }
     p.eff[t_] <- 1-prod(1-pd[1:T2[t_],t_]) #effective detection prob per primary period
  }
  #loop through (T-1) primary periods
  #NOTE: trmat's are offset -1 in time, eg. t_=1 implies a transition between period 1 to period 2.
  for(t_ in 1:(T-1)){ 
     gamma2[t_] ~ dbeta(pr.gamma2[1],pr.gamma2[2]) #temp migration: prob migrate out conditional on being inside
     #trmat: transition matrix for Markovian latent-states
     #1=dead;2=offsite;3=onsite
     #transition are from the column --> rows
     #trmat[row,column,time] = [state at time=t_; state at time t_-1; time=t_]
     trmat[1,1,t_]<- 1 #stay dead
     trmat[2,1,t_]<- 0
     trmat[3,1,t_]<- 0
     trmat[1,2,t_]<- 1-phi #dies outside
     trmat[2,2,t_]<- gamma1*phi #stays outside | outside
     trmat[3,2,t_]<- (1-gamma1)*phi #reenters study area | outside
     trmat[1,3,t_]<- 1-phi #dies inside
     trmat[2,3,t_]<- gamma2[t_]*phi #leaves study area | inside
     trmat[3,3,t_]<- (1-gamma2[t_])*phi #stays inside | inside
  } #t_ state process
  #PART1: likelihood for first-capture (all individuals)
  for(i in 1:N){
    #Observation error during 1st capture: condition on at least one capture: 
    #the following formula is the (conditional) multinomial log-likelihood of
    #...the sequence of observations in a primary period, conditional on that we 
    #...know they were seen at least once
        LL[i]<-sum(y[i,1:T2[first[i]],first[i]]*log(pd[1:T2[first[i]],first[i]])+ (1-y[i,1:T2[first[i]],first[i]])*log(1-pd[1:T2[first[i]],first[i]])) - log(p.eff[first[i]]) #multinomial log-likelihood
	zeros[i] ~ dpois(-1*LL[i]+C) #Winbugs zeros trick, likelihood passed to JAGS as ZIP
  } #i 
  mintrick <- max(LL[1:N]) #strictly for monitoring the first-capture likelihood trick
  #PART2: loop through individuals potentially available for 2 or more primary periods
  for(i in 1:length(N.ix2)){
    #state process for latent state after their first primary period
    #draw z conditional on being seen during previous primary period
    z[N.ix2[i],first[N.ix2[i]]+1]~ dcat(trmat[1:3,3,first[N.ix2[i]]]) 
    #loop through secondary periods
    for(tt_ in 1:T2[first[N.ix2[i]]+1]){
        #likelihood of secondary periods observation, conditional on z=3
        y[N.ix2[i],tt_,first[N.ix2[i]]+1] ~ dbern(pd[tt_,first[N.ix2[i]]+1]* equals(z[N.ix2[i],first[N.ix2[i]]+1],3))
    } #tt_
 } #N.ix2
 #PART3: loop through individuals potentially available for 3 or more primary periods
 for(i in 1:length(N.ix3)){
   #loop through remain primary periods after first and second primary periods 
   for(t_ in (first[N.ix3[i]]+2):T){ 
       #state process: draw z(t_) conditional on  z(t_-1)
       z[N.ix3[i],t_]~ dcat(trmat[1:3, z[N.ix3[i],t_-1], t_-1]) #
       #Observation error: Bernoulli
       for(tt_ in 1:T2[t_]){
           #likelihood of secondary periods observation, conditional on z=3
           y[N.ix3[i],tt_,t_] ~ dbern(pd[tt_,t_]*equals(z[N.ix3[i],t_],3))
      } #tt_
    } #t_
  } #N.ix3
# estimate number of individuals available for capture (inside) per primary period
for(t_ in 1:T){
  Nin[t_] <- n.vector[t_]/p.eff[t_] 
} #t_
}"
# save jags file to local disk
modname2 <- "JAGS_firstcapt_fixedeff.JAG"
sink(file=modname2)
cat(jags.txt.2,fill=TRUE)
sink()

# Special data for conditioning on first capture
first = apply(Y.t,1,function(x) min(which(x>0))) # which primary period were individuals first encountered
seen.potentially.twice <- which(first <T) # potential to be seen at least twice
seen.potentially.thrice<- which(first < (T-1)) # potential to be seen at least twice
n.vector = apply(Y.t, 2,function(x) sum(x>0)) # number of unique animals per primary period, seen
# Data for jags: Capture-histories and Prior parameters
jags.data=list(y=Y.tt, # 3D array of capture histories
               T=length(T2), # number of primary periods
               T2=T2, # number of secondary periods per period period
               N= nrow(Y.tt), # number of observed individuals
               first = first, # vector (N), which primary period each individual was first encountered;
               N.ix2 = seen.potentially.twice,# vector of indices pointing to those individuals who were potentially available for 2 or more primary periods, i.e., their first primary period was at least T-1 or earlier.
               N.ix3 =seen.potentially.thrice, # vector of indices pointing to those individuals who were potentially available for 3 or more primary periods, i.e., their first primary period was at least T-2 or earlier, which may include some of the same individuals in N.ix2
               zeros = rep(0,nrow(Y.tt)),# vector of zeros (N), and is used for the WinBUGS ``zeros trick'' to manually calculate likelihood of observations during secondary period in which we know there was at least one observation during the entire primary period;
               C = 20, # for the zero's trick, to ensure that the poisson density never goes negative
               n.vector=n.vector,# number of unique animals seen per primary period (useful for estimating unseen animals)
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
# user is expected to be able to initialize the univariate parameters;I provide 'generate.z.psi' function which imputs values of latent states 'z' (forwards-backwards algorithm) 
init.func = function(){
    RET=list(
        phi = rbeta(1,60,10), # user to specify
        gamma1=rbeta(1,20,20), # user to specify
        gamma2=rbeta(T-1 ,20,20),# user to specify (notice it is T-1)
        pd_a =rgamma(T,30,20),# user to specify
        pd_b = rgamma(T,30,20),# user to specify
        pd = do.call(function(T2,T){ # user to specify
            pd = matrix(NA,max(T2),T);
            for(t_ in 1:T){ pd[1:T2[t_],t_]<-rbeta(T2[t_],12,65)};
            return(pd)}, args=list(T2=T2,T=T))
    );
    RET=c(RET, generate.z.psi(y=Y.tt,T2=T2,first.capture=TRUE,z.priors = list(phi.beta=c(shape1=30,shape2=5),g1.beta=c(shape1=20,shape2=20),g2.beta=c(shape1=20,shape2=20), pd.beta=c(shape1=12,shape2=65)))); # handy function to initialize latent states z
    return(RET)
}

# initialize random variables 
if(nchains == 1){
    jags.inits = init.func()
} else {
    jags.inits = lapply(1:nchains, function(x) init.func())
}

# RUN JAGS MODEL AND SAMPLE FROM POSTERIOR
  m2 <- jags.model(file=modname2,data=jags.data,inits=jags.inits,n.chains=nchains,n.adapt=nadapt)
# burn-in phase
  update(m2,nburn) # discard nburn iterations (ensure that chains reach equilibium states)
# which variables to summarize
  variable.names <- c("phi","gamma1",paste0("gamma2[",1:(T-1),"]"),unlist(mapply(as.list(1:T),sapply(X=T2,FUN=seq,from=1,by=1),FUN=function(x,y) paste0("pd[",y,",",x,"]"))),paste0("Nin[",1:T,"]"))
# sample from posterior
  samp <- coda.samples(m2,variable.names=variable.names,n.iter=niter,thin=thin_) 

  summary(samp) # summarize output
  plot(samp, ask =TRUE) # inspect chains for adequate mixing and convergence
  gelman.diag(samp) # inspect chains for convergence (statistics ~ 1 and < 1.1)
# DONE PART 1
    
