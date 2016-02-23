# this is stricly to generat the image on the tutorial website
# various scenarios:
# sparse data, medium data, high data
# nu = 1; nu =2; nu =4; nu=10
# sig = 0.2
library(rjags)
library(boot)
# start with the JAGS model
jagstxt = 'model{
 # hyperprior
 sigma ~ dt(0,pr.sigma.tau, pr.sigma.nu) T(0,) # scaled half student t
 logit.mu ~ dnorm(0,pow(1.55,-2)) # location
 p.mu = 1/(1+exp(-logit.mu))
 for( i in 1:length(y)){
    logit.p[i] ~ dnorm(logit.mu, pow(sigma,-2)) # logit(probability)
    p[i] <- 1/(1+exp(-logit.p[i]))
    y[i] ~ dbin(p[i], N) # binomial distribution
 }
}'
sink(file="make_studentT_graph.JAG")
cat(jagstxt,fill=TRUE)
sink()

#(seeed=round(runif(1)*10000))
#set.seed(seeed)
#(y.ratio.true = sample(data.mles, size=N.groups, replace=TRUE,prob=dbinom(data.mles*(length(data.mles)-1),size=min(data.size.vec),prob=mu.true))) # 
#table(y.ratio.true)

set.seed(9216)
# scenario 1: low data
data.size.vec = c(5,10,20,30,50,200)
nu.vec <- c(1,2,4,10)
sigma.vec <- c(0.01, 0.1, 0.2)
scenarios <- expand.grid(data.size = data.size.vec, nu=nu.vec, sigma=sigma.vec)
# generate data at particular MLE ratios; which will given different amounts of data
N.groups = 10
data.mles = seq(0,1,length.out=min(data.size.vec)+1) # pre-define MLE's
mu.true = 0.5
# each group is simulating a binomial response, 
y.ratio.true = sample(data.mles, size=N.groups, replace=TRUE,prob=dbinom(data.mles*(length(data.mles)-1),size=min(data.size.vec),prob=mu.true)) # 
dhalft.prop <- function(x,sigma,nu) (((1 + (1/nu)*(x/sigma)^2)^(-(nu+1)/2))*(x>=0))
dhalft <- function(x,sigma,nu) dhalft.prop(x,sigma,nu)/integrate(f=function(x) dhalft.prop(x,sigma=sigma,nu=nu),lower=0,upper=20)[[1]]
for(s_ in 1:nrow(scenarios)){
    png(filename = paste0("img_",ifelse(nchar(s_)==1,paste0("0",s_),as.character(s_)),".png"),width=1000,height=600) #quality=80
    data.size = scenarios[s_,"data.size"]
    nu = scenarios[s_,"nu"]
    sigma = scenarios[s_,"sigma"]

    # get the data
    N = data.size; plotname = paste0(c("5"="SPARSE   ","10"="LOW      ","20"="MODERATE","30"="HIGH     ","50"="VERY HIGH ","200"="EXTREME! ")[as.character(N)],"\n(",N," observations per group)")#c("5"="sparse","10"="low","20"="moderate","200"="high")[as.character(N)]
    y = as.integer(y.ratio.true*(N))

    # first estimate the mle's
    m.glm = glm(y~group-1,family="binomial",data=data.frame(y=as.numeric(sapply(y,function(y_) c(rep(1,y_),rep(0,N-y_)))), group=factor(rep(1:N.groups,each=N))))
    m.se = (summary(m.glm))[["coefficients"]][paste0("group",1:N.groups),"Std. Error"]
    m.hici=inv.logit(logit(y.ratio.true) + qnorm(0.975)*m.se)
    m.loci=inv.logit(logit(y.ratio.true) + qnorm(0.025)*m.se)    
    
    nchains=4;nadapt=2000;nburn=1000;niter=20000;nsamp=300;thin_=round(niter/nsamp)
    jagsdata = list(y=y,N=N,pr.sigma.tau = 1/(sigma^2), pr.sigma.nu = nu)
    jagsinitsf = function(x){ list(sigma=runif(1,0,0.1), logit.mu=runif(1,0,1), logit.p=logit(runif(N.groups,0.2,0.8)))}
    jagsinits = lapply(1:nchains,jagsinitsf)

    # run the model
    m <- jags.model(file="make_studentT_graph.JAG",data=jagsdata,inits=jagsinits,n.chains=nchains,n.adapt=nadapt)
    update(m,nburn)
    s <- coda.samples(m,variable.names=c("p","p.mu"),n.iter=niter,thin=thin_)
    s2 = do.call("rbind",s);rm(s)

    # make a plot of the hyperprior curve
    par(mfrow=c(1,2),mgp=c(1.8,0.6,0),bty="l",family="Times",mar=c(3,3,3,0))
    plot(c(0,0.8),c(0,8),xlab = expression(sigma), ylab = expression("density f("*sigma*")"),cex.lab=1.5,type="n",main="Scaled Half Student-t: Hyperprior",cex.main=1.8)
    # make a polygon of the curent density
    mutecol = "grey20"; curcol="purple" #
    xseq = seq(0.0001,1.5,0.01)
    densy = dhalft(xseq,sigma=sigma,nu=nu)
    polygon(c(xseq,rev(xseq)), c(densy,rep(0,length(densy))),col="grey80",border=NA)
#    lines(xseq,densy,col=curcol,lwd=3)    
    curves_ = unique(scenarios[,c("sigma","nu")])
    for(c_ in 1:nrow(curves_)){
        lines(xseq,dhalft(xseq,sigma=curves_[c_,"sigma"],nu=curves_[c_,"nu"]),col=mutecol,lwd=1)
    }
    lines(xseq,dhalft(xseq,sigma=sigma,nu=nu),col=curcol,lwd=4) # current nu
#    eval(parse(text =paste0('text(x=1.2,y=6, labels=expression("f("*sigma*")" %prop% "T(0, s="*',sigma,'*","~nu~"=',nu,')I["~sigma>0~"]"),cex=2,pos=2,font=2)')))
    eval(parse(text =paste0('text(x=0.55,y=6, labels=expression("f("*sigma*")" %prop% "T("*mu*",s,"*nu*")I["~sigma>0~"]"),cex=2,font=2)')))
    eval(parse(text =paste0('text(x=0.55,y=5.5, labels=expression(mu*"=0"),cex=2,font=2)')))
    # three s values
#    eval(parse(text =paste0('text(x=0.55,y=5.0, labels=expression("s="),cex=2,font=2,col="purple",pos=1)')))
    space=c(0.08,0.055)
    eval(parse(text =paste0('text(x=0.52-space[1]*1,y=5.2, labels="s=",cex=2,font=2,col="purple",pos=1)')))    
    eval(parse(text =paste0('text(x=0.52+space[1]*0,y=5.2, labels="',sigma.vec[1],'",cex=2,font=2,col=c("grey80",curcol)[1+1*(sigma==sigma.vec[1])],pos=1)')))
    eval(parse(text =paste0('text(x=0.52+space[1]*1,y=5.2, labels="',sigma.vec[2],'",cex=2,font=2,col=c("grey80",curcol)[1+1*(sigma==sigma.vec[2])],pos=1)')))
    eval(parse(text =paste0('text(x=0.52+space[1]*2,y=5.2, labels="',sigma.vec[3],'",cex=2,font=2,col=c("grey80",curcol)[1+1*(sigma==sigma.vec[3])],pos=1)')))
    # nu
    eval(parse(text =paste0('text(x=0.52-space[1]*1,y=4.7, labels=expression(nu*"="),cex=2,font=2,col="purple",pos=1)')))    
    eval(parse(text =paste0('text(x=0.52+space[2]*0,y=4.7, labels="',nu.vec[1],'",cex=2,font=2,col=c("grey80",curcol)[1+1*(nu==nu.vec[1])],pos=1)')))
    eval(parse(text =paste0('text(x=0.52+space[2]*1,y=4.7, labels="',nu.vec[2],'",cex=2,font=2,col=c("grey80",curcol)[1+1*(nu==nu.vec[2])],pos=1)')))
    eval(parse(text =paste0('text(x=0.52+space[2]*2,y=4.7, labels="',nu.vec[3],'",cex=2,font=2,col=c("grey80",curcol)[1+1*(nu==nu.vec[3])],pos=1)')))
    eval(parse(text =paste0('text(x=0.52+space[2]*3,y=4.7, labels="',nu.vec[4],'",cex=2,font=2,col=c("grey80",curcol)[1+1*(nu==nu.vec[4])],pos=1)')))            
    #eval(parse(text =paste0('text(x=0.55,y=4.5, labels=expression(nu*"=',nu,'"),cex=2,font=2,col="purple")')))    
    # PLOTs: MLEs and posteriors,
    plot(c(0,1),c(1,N.groups),xlab = "Capture probabilities p", ylab = "Group",cex.lab=1.5,type="n",main = paste("data quantity:",plotname),cex.main=1.8)
    lines(x=rep(mean(s2[,"p.mu"]),2),y=c(1,N.groups),col="purple",lwd=1,lty=2)    
    for(k in 1:N.groups){
        # draws posterior 95%CI
        lines(x = HPDinterval(as.mcmc(s2[,paste0("p[",k,"]")]))[1:2],y=rep(k,2),col="grey40",lwd=1)
        lines(x = HPDinterval(as.mcmc(s2[,paste0("p[",k,"]")]),prob=0.68)[1:2],y=rep(k,2),col="purple",lwd=3)
        # draw confidence intervals
        lines(x = c(m.loci[k],m.hici[k]),y=rep(k,2)-0.1,col="red",lwd=1,lty=2)
        # point estimates
        points(c(mean(s2[,paste0("p[",k,"]")]), y.ratio.true[k]), y=c(k,k-0.1), pch = c(19,8), col=c("purple","red"), cex=2.5)
    }
    legend(x="bottomleft",legend=c("MLE and 95%CI","posterior mean\n and 95%CI"),pch = c(8,19), lwd=c(1,2), lty=c(2,1), col=c("red","purple"), cex=1.5,bty="n")
    text(x=0,y=10,pos=4,labels=expression("logit("*p[k]*")" %~%  "N("*mu*","*sigma*")"),cex=1.3)
    text(x=0,y=9.6,pos=4,labels=expression(y %~% "Bin(N,"*p[k]*")"),cex=1.3)    
    dev.off()
} # scenarios
# make animated gif from png (on unix-like environments)
# convert -delay 120 -loop 0 +map *.png HalfTdemo.gif  # 2.5
# convert -delay 100 -loop 0 +dither -map *.png animation.gif # 2.6
# convert -delay 120 -loop 0 -layers optimize *.png HalfTdemo.gif  # 681kb
