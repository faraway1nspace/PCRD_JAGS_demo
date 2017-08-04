# THIS IS A SOURCE FILE FOR ARTICLE
# see article Rankin RW, Nicholson K, Allen S, Kr\'{u}tzen M, Bejder L, Pollock KH. 2016 "A full-capture Hierarchical Bayesian model of Pollock's Closed Robust Design and application to dolphins" Frontiers of Marine Science, in Press.

# function to import MARK inp capture histories ONLY the capture histories.
import.mark.inp <- function(
                            MARK.file.name, # mark file name
                            T2, # vector of number of secondary periods per primary period
                            header=FALSE # is there a header?
                            ){
    # load the capture histories, saved as a MARK '.inp' file.
    warning("this function only imports the capture histories -- not other individual covariates. You'll have to do that manually")
    T = length(T2);
    rawinp <- paste(scan(MARK.file.name,what="character"),collapse=""); # raw import 
    rawinp2 <- strsplit(gsub("\\*\\/",";",gsub(pattern="\\/\\*","",rawinp)),split=";"); # split by row names and CH
    if(header){rawinp2=rawinp2[-1]};
    rawinp3 <- rawinp2[[1]][seq(2,length(rawinp2[[1]]),by=2)]; # extract CH
    id <- rawinp2[[1]][seq(1,length(rawinp2[[1]]),by=2)]; # extract ID's of animals
    T = length(T2); # number of primary periods
    N = length(rawinp3); # number of (observed) individuals
    # expand the CH into a (ragged) 3D array: [id, secondary periods, primary periods]
    Y.tt <- array(NA,c(N,max(T2),T)); # N x max(T2) x T
    for(x in 1:T){
        splitx=substring(rawinp3,first=sum(c(0,T2)[1:x])+1,last=sum(c(0,T2)[1:(x+1)]));        
        Y.tt[,,x] <- as.matrix(t(as.data.frame(lapply(strsplit(splitx,""),function(q){as.numeric(q)[1:max(T2)]}))));
    };
    dimnames(Y.tt)[[1]] <- id;
    # collect other, non chapture history data
    other.data=matrix(substring(rawinp3,first=sum(T2)+1),N,1); row.names(other.data)=id;
    # collapse 3D array just into NxT, summing counts per primary periods
    Y.t <-apply(Y.tt,c(1,3),function(x) sum(x,na.rm=TRUE)); # sum captures per primary period
    return(
        list(Y.t = Y.t, # ch:returns N x T matrix of counts per primary period
           Y.tt = Y.tt, # ch: returns 3D array, N x T2 x T, shows which 2nd-ary periods had captures
           other.data=other.data # other data (not processed by this function)
           ));
};

# backwards forward sampler to draw latent states (consistent with the data)
generate.z.psi <- function(y,T2,first.capture=FALSE,z.priors = list(phi.beta=c(shape1=30,shape2=5),g1.beta=c(shape1=20,shape2=20),g2.beta=c(shape1=20,shape2=20), pd.beta=c(shape1=12,shape2=65))){
    exclude_=1; dead_=2;out_=3;in_=4; nstates = 4;           
     T=length(T2); # number of primary periods
     Y.t <-apply(y,c(1,3),function(x) sum(x,na.rm=TRUE)); # sum captures per primary period
     mm<-nrow(Y.t);
     first.capt <- apply(Y.t,1,function(x){ if(any(x>0)){ min(which(x>0))} else {T+1}});
     priors.str = strsplit(names(z.priors),split=".",fixed=TRUE); # get name of parameter and its density
     random.var = lapply(priors.str, function(x) x[1]); # get names of random variables
     for(par_ in 1:length(z.priors)){
         rdraw=do.call(get(paste0("r",priors.str[[par_]][2])),args=c(list(n=T),as.list(z.priors[[par_]])));
         assign(priors.str[[par_]][1],value=rdraw);
     };
     # estimate psi values
     psi.counts = tabulate(first.capt,(T+1))[1:T]; # estimate of number of entries per primary periods
     psi= rbeta(T,10+psi.counts,5+(mm + mm*first.capture)-psi.counts); # random draws from beta
     # construct transition matrix
     t.mat = array(0,c(4,4,T));
     if(!any(random.var=="lambda")){ # esimate lambda, only if its not a random variable
         lambda = (1-g1)/(g2-g1+1); # eigenvector decomposition: probability of entry into 'onsite' state
     };
     for(t_ in 1:T){
         t.mat[,1,t_]<-c(1-psi[t_],0,psi[t_]*(1-lambda[t_]),psi[t_]*lambda[t_]);
         t.mat[,2,t_]<-c(0,1,0,0);
         t.mat[,3,t_]<-c(0,1-phi[t_],g1[t_]*phi[t_],(1-g1[t_])*phi[t_]);
         t.mat[,4,t_]<-c(0,1-phi[t_],g2[t_]*phi[t_],(1-g2[t_])*phi[t_]);
     };
     z = matrix(0,mm,T); # latent states 
     a = matrix(0,nstates,T+1);a[,1]=1*((1:4)==exclude_); # backwards pass-values
     z.vec = numeric(T); # latent states
     for(i in 1:mm){
         for(t_ in 1:T){
             alph = dbinom(x = Y.t[i,t_],size=T2[t_],prob = pd[t_]*((1:nstates) == in_)) * (t.mat[,,t_])%*%a[,t_]; # notice a has T+1 elements; here, a is really indexed t_-1, but because it has 1 more column, t_-1+1 equals t_
             a[,t_+1] = alph/sum(alph);
         };
         # forward sampling
         z.vec[T] <- sample(1:nstates,1,replace=FALSE,prob=a[,T+1]);
         for(t_ in (T-1):1){
             p_z = dbinom(x = Y.t[i,t_+1],size=T2[t_+1],prob = pd[t_+1]*(z.vec[t_+1] == in_)) * t.mat[z.vec[t_+1],,t_+1]*a[,t_+1]/a[z.vec[t_+1],t_+2];
             z.vec[t_]<-sample(1:nstates,1,replace=FALSE,prob=p_z/sum(p_z));
         };
         z[i,] <- z.vec;
     }; # i
     # if conditioning on first capture 
     if(first.capture){ # replace z values as NA, if only modelling fullcapture
         z <- z-1; # need to remove unseen 'not-yet-entered' state
         for(i in 1:mm){
             z[i,1:first.capt[i]]<-rep(NA,first.capt[i]); 
         };
         RET = list(z=z);
     } else { # 
     # if modelling full capture histories, need psi estimates too
         RET = list(z = z, psi = psi); # return latents states and psi
     };
     return(RET); # inner function
  };



             

         
    
    



