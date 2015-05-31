#################################################################################################
#
#  FIRST ORDER DLM 
#
#  AUXILIARY PARTICLE FILTER + PARAMETER LEARNING (Liu and West, 2001)
#
#  y(t) = x(t)   + e(t)   e(t) ~ N(0,sig2)
#  x(t) = x(t-1) + w(t)   w(t) ~ N(0,tau2)
#  
#  x(0) ~ N(m0,V0)
#  sig2 ~ IG(a0,b0)
#
#################################################################################################
#
# Author : Hedibert Freitas Lopes
#          The University of Chicago Booth School of Business
#          5807 South Woodlawn Avenue
#          Chicago, Illinois, 60637
#          Email : hlopes@Chicagobooth.edu
#
#################################################################################################
q025 = function(x){quantile(x,0.025)}

q975 = function(x){quantile(x,0.975)}

dinvgamma = function(x,a,b){dgamma(1/x,a,b)*(1/x^2)}

drawnorm = function(x){rmvnorm(1,x[1:2],matrix(x[3:6],2,2))}

likelihood = function(x,mu,sig){dnorm(x,mu,sig)}

ffbslike=function(y,sig2,tau2,d0,V0){
  n=length(y)
  like = 0.0
  for (t in 1:n){
    if (t==1){
      a = d0
      R = V0 + tau2
    }else{
      a = m
      R = C + tau2
    }
    f = a
    Q = R + sig2
    A = R/Q
    m = a+A*(y[t]-f)
    C = R - Q*A**2
    like = like + dnorm(y[t],f,sqrt(Q),log=T)
  }
  return(like)
}

liuwest1 = function(y,tau2,N,a0,b0,m0,V0,delta){
  n       = length(y)
  ws      = NULL
  thetass = NULL
  sig2s   = NULL
  sig2    = 1/rgamma(N,a0,b0)
  thetas  = rnorm(N,m0,sqrt(V0))
  h2      = 1-((3*delta-1)/(2*delta))^2
  a       = sqrt(1-h2)
  for (t in 1:n){
    lsig2   = log(sig2)
    mlsig2  = mean(lsig2)
    vl      = var(lsig2)
    
    mulsig2 = a*lsig2+(1-a)*mlsig2
    mus     = thetas
    weight  = dnorm(y[t],mus,exp(mulsig2/2))
    
    k       = sample(1:N,size=N,replace=T,prob=weight)
    sig2    = exp(rnorm(N,mulsig2[k],sqrt(h2*vl)))
    thetas1 = rnorm(N,thetas[k],sqrt(tau2))
    weight  = dnorm(y[t],thetas1,sqrt(sig2))/dnorm(y[t],mus[k],exp(mulsig2[k]/2))

    ind     = sample(1:N,size=N,replace=T,prob=weight)
    thetas  = thetas1[ind]
    sig2    = sig2[ind]
    ws      = rbind(ws,weight)
    sig2s   = rbind(sig2s,sig2)
    thetass = rbind(thetass,thetas)
  }
  return(sig2s=sig2s)
}

liuwest2 = function(y,tau2,N,a0,b0,m0,V0,delta){
  n       = length(y)
  ws      = NULL
  thetass = NULL
  sig2s   = NULL
  sig2    = 1/rgamma(N,a0,b0)
  thetas  = rnorm(N,m0,sqrt(V0)) # this is sampling of x
  h2      = 1-((3*delta-1)/(2*delta))^2
  a       = sqrt(1-h2)
  for (t in 1:n){
    # calculating  weight
    msig2   = mean(sig2)
    vl      = h2*var(sig2)
    musig2  = a*sig2+(1-a)*msig2
    as      = musig2^2/vl+2
    bs      = musig2*(as-1)
    mus     = thetas 
    weight  = dnorm(y[t],mus,sqrt(musig2)) # q(x_t,theta_t|y_(t+1))
    # resampling
    k       = sample(1:N,size=N,replace=T,prob=weight)
    # sampling twice
    sig2    = 1/rgamma(N,as[k],bs[k])
    thetas1 = rnorm(N,thetas[k],sqrt(tau2)) #x_t+1 only depend on tau not sigma
    weight  = dnorm(y[t],thetas1,sqrt(sig2))/dnorm(y[t],mus[k],sqrt(musig2[k]))

    ind     = sample(1:N,size=N,replace=T,prob=weight)
    thetas  = thetas1[ind]
    sig2    = sig2[ind]
    # resample so that we don't have to keep weight for the next round
    # but can still pass on distribution information
    # keep in mind that we replace weight at line (RP)    
    ws      = rbind(ws,weight)
    sig2s   = rbind(sig2s,sig2)
    thetass = rbind(thetass,thetas)
  }
  return(sig2s=sig2s)
}

liuwest3 = function(y,tau2,N,a0,b0,m0,V0,delta){
  n       = length(y)
  ws      = NULL
  thetass = NULL
  sig2s   = NULL
  sig2    = 1/rgamma(N,a0,b0)
  thetas  = rnorm(N,m0,sqrt(V0))
  h2      = 1-((3*delta-1)/(2*delta))^2
  a       = sqrt(1-h2)
  for (t in 1:n){
    lsig2   = log(sig2)
    mlsig2  = mean(lsig2)
    vl      = var(lsig2)
    mulsig2 = a*lsig2+(1-a)*mlsig2
    mus     = thetas 
    weight  = dnorm(y[t],mus,exp(mulsig2/2)) # (RP)
    k       = sample(1:N,size=N,replace=T,prob=weight)
    sig2    = exp(rnorm(N,mulsig2[k],sqrt(h2*vl)))
    #optimal propogation
    var     = 1/(1/sig2+1/tau2)
    mean    = var*(y[t]/sig2+thetas[k]/tau2)
    #
    thetas1 = rnorm(N,mean,sqrt(var))
    weight  = dnorm(y[t],thetas1,sqrt(sig2))/dnorm(y[t],mus[k],exp(mulsig2[k]/2))
    w1      = dnorm(thetas1,thetas[k],sqrt(tau2))/dnorm(thetas1,mean,sqrt(var))
    weight  = weight*w1
    
    # To be used next loop
    ind     = sample(1:N,size=N,replace=T,prob=weight)
    thetas  = thetas1[ind]
    sig2    = sig2[ind]

    ws      = rbind(ws,weight)
    sig2s   = rbind(sig2s,sig2)
    thetass = rbind(thetass,thetas)
  }
  return(sig2s=sig2s)
}

liuwest4 = function(y,tau2,N,a0,b0,m0,V0,delta){
  n       = length(y)
  ws      = NULL
  thetass = NULL
  sig2s   = NULL
  sig2    = 1/rgamma(N,a0,b0)
  thetas  = rnorm(N,m0,sqrt(V0))
  h2      = 1-((3*delta-1)/(2*delta))^2
  a       = sqrt(1-h2)
  for (t in 1:n){
    msig2   = mean(sig2)
    vl      = h2*var(sig2)
    musig2  = a*sig2+(1-a)*msig2
    as      = musig2^2/vl+2
    bs      = musig2*(as-1)
    mus     = thetas 
    weight  = dnorm(y[t],mus,sqrt(musig2))
    k       = sample(1:N,size=N,replace=T,prob=weight)
    sig2    = 1/rgamma(N,as[k],bs[k])
    
    var     = 1/(1/sig2+1/tau2)
    mean    = var*(y[t]/sig2+thetas[k]/tau2)
    thetas1 = rnorm(N,mean,sqrt(var))
    weight  = dnorm(y[t],thetas1,sqrt(sig2))/dnorm(y[t],mus[k],sqrt(musig2[k]))
    w1      = dnorm(thetas1,thetas[k],sqrt(tau2))/dnorm(thetas1,mean,sqrt(var))
    weight  = weight*w1
    ind     = sample(1:N,size=N,replace=T,prob=weight)
    thetas  = thetas1[ind]
    sig2    = sig2[ind]
    ws      = rbind(ws,weight)
    sig2s   = rbind(sig2s,sig2)
    thetass = rbind(thetass,thetas)
  }
  return(sig2s=sig2s)
}

#######################################################################################################
# Data simulation
set.seed(12345)
n     = 200
tau2  = 0.05
sig2  = 0.1
sqrt(tau2/sig2)
x     = rep(0,n+1)
y     = rep(0,n+1)
x[1]  = 25
y[1]  = rnorm(1,x[1],sqrt(sig2))
for (t in 2:(n+1)){
  x[t] = rnorm(1,x[t-1],sqrt(tau2))
  y[t] = rnorm(1,x[t],sqrt(sig2))
}
x0 = x[1]
x  = x[2:(n+1)]
y  = y[2:(n+1)]

# True values
x0true   = x0
xtrue    = x
sig2true = sig2
tau2true = tau2

# Prior hyperparameters
a0 = 5
b0 = (a0-1)*sig2
m0 = x0true
V0 = 100

# Posterior for sig2
# ------------------
L      = 0.0001
U      = 5*sig2true
N      = 100
sig2s  = seq(L,U,length=N)
hsig2  = sig2s[2]-sig2s[1]
ns     = seq(1,n,by=4)
qs     = NULL
for (i in ns){
  logpost = rep(0,N)
  for (j in 1:N)
    logpost[j] = ffbslike(y[1:i],sig2s[j],tau2true,m0,V0)+log(dinvgamma(sig2s[j],a0,b0))
  postsig2 = exp(logpost-max(logpost))
  postsig2 = postsig2/sum(postsig2)
  postsig2 = postsig2/hsig2
  q        = cumsum(postsig2*hsig2)
  quant    = c(max(sig2s[q<0.025]),max(sig2s[q<0.5]),max(sig2s[q<0.975]))
  qs       = rbind(qs,quant)
}

# APF+LW
# ------
set.seed(54321)
M             = 2000
delta1        = 0.75
delta2        = 0.95
deltas        = c(rep(delta1,4),rep(delta2,4))
filters       = array(0,c(n,M,10))
filters[,,1]  = liuwest1(y,tau2true,M,a0,b0,m0,V0,delta1)
filters[,,2]  = liuwest2(y,tau2,N,a0,b0,m0,V0,delta1)
filters[,,3]  = liuwest3(y,tau2true,M,a0,b0,m0,V0,delta1)
filters[,,4]  = liuwest4(y,tau2true,M,a0,b0,m0,V0,delta1)
filters[,,5]  = liuwest1(y,tau2true,M,a0,b0,m0,V0,delta2)
filters[,,6]  = liuwest2(y,tau2,N,a0,b0,m0,V0,delta2)
filters[,,7]  = liuwest2(y,tau2true,M,a0,b0,m0,V0,delta2)
filters[,,8]  = liuwest4(y,tau2true,M,a0,b0,m0,V0,delta2)

mfilters = matrix(0,n,8)
Lfilters = matrix(0,n,8)
Ufilters = matrix(0,n,8)
for (i in 1:8){
  Lfilters[,i] = apply(filters[,,i],1,q025)
  mfilters[,i] = apply(filters[,,i],1,mean)
  Ufilters[,i] = apply(filters[,,i],1,q975)
}

ind = ns
L = min(qs[,1],Lfilters[ind,])
U = max(qs[,3],Ufilters[ind,])

pdf(file="tausig3.pdf",width=10,height=7.5)
inds=c(1,2,3,4,1,2,3,4)
par(mfrow=c(2,4))
for (i in 1:8){ 
  plot(ind,mfilters[ind,i],type="l",xlab="time",ylab="",lwd=1,ylim=c(L,U),col=2)
  title(paste("LW",inds[i]," with delta=",deltas[i],sep=""))
  lines(ind,Lfilters[ind,i],col=2,lwd=1,lty=1)
  lines(ind,Ufilters[ind,i],col=2,lwd=1,lty=1)
  lines(ind,qs[,1],col=4)
  lines(ind,qs[,2],col=4)
  lines(ind,qs[,3],col=4)
  abline(h=sig2true,lwd=1,col=3)
}
dev.off()






