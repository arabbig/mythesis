#####################################################################################
#
# Application 2: LOG STOCHASTIC VOLATILITY AR(1) MODEL
#
# Lopes and Tsay (2011)
# Particle Filters and Bayesian Inference in Financial Econometrics.
# Journal of Forecasting.
# 
#####################################################################################
#
# For t=1,...,n
#
#    y[t]|x[t]   ~ N(0,exp(x[t]/2))
#    x[t]|x[t-1] ~ N(alpha+beta*x[t-1],tau2)
#
# and
#
#    x[0]       ~ N(m0,C0)
#    alpha|tau2 ~ N(b0[1],tau2*B0[1])
#    beta|tau2  ~ N(b0[2],tau2*B0[2])
#    tau2       ~ IG(nu0/2,nu0*tau20/2)
#
# with known hyperparameters m0,C0,b0,B0,nu0 and tau20.
#
#####################################################################################
#
# Data: Monthly log returns of GE stock
# Period: January 1926 to December 1999 (888 observations)
# Source: Tsay (2005), Chapter 12, Example 12.6, page 591.
# http://faculty.chicagobooth.edu/ruey.tsay/teaching/fts2/m-geln.txt
#
######################################################################################
#
# Author : Hedibert Freitas Lopes
#          Associate Professor of Econometrics and Statistics
#          The University of Chicago Booth School of Business
#          5807 South Woodlawn Avenue
#          Chicago, Illinois, 60637
#          Email : hlopes@ChicagoBooth.edu
#          URL: http://faculty.chicagobooth.edu/hedibert.lopes
#
#####################################################################################
rm(list=ls())

getindex = function(x){c(x[1+x[1]],x[8+x[1]])}

ldt = function(x,nu,m,sd){
  dt((x-m)/sd,nu,log=TRUE)-log(sd)
}

LW = function(y,alphas,betas,tau2s,xs,delta){
  n  = length(y)
  N  = length(xs)
  quants = array(0,c(n,4,3))
  h2 = 1-((3*delta-1)/(2*delta))^2
  a  = sqrt(1-h2)
  pars = cbind(alphas,betas,log(tau2s))
  like = rep(0,n)
  for (t in 1:n){
    like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
    # Resampling
    mus     = pars[,1]+pars[,2]*xs
    mpar    = apply(pars,2,mean)
    vpar    = var(pars)
    ms      = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    weight  = dnorm(y[t],0.0,exp(mus/2),log=TRUE)
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    # Propagating
    if (delta<1){
      ms1 = ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    }else{
      ms1 = ms
    }
    xt   = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=T,prob=w)
    xs   = xt[ind]
    pars = ms1[ind,]  	 
    # Storing quantiles
    quants[t,1,] = quantile(pars[,1],c(0.5,0.025,0.975))
    quants[t,2,] = quantile(pars[,2],c(0.5,0.025,0.975))
    quants[t,3,] = quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    quants[t,4,] = quantile(exp(xs/2),c(0.5,0.025,0.975))  
  }
  return(list(like=like,quants=quants))
}


LWt = function(y,nu,alphas,betas,tau2s,xs,delta){
  n  = length(y)
  N  = length(xs)
  quants = array(0,c(n,4,3))
  h2 = 1-((3*delta-1)/(2*delta))^2
  a  = sqrt(1-h2)
  pars = cbind(alphas,betas,log(tau2s))
  like = rep(0,n)
  for (t in 1:n){
    like[t] = sum(exp(ldt(y[t],nu,0.0,exp(xs/2))))
    # Resampling
    mus     = pars[,1]+pars[,2]*xs
    mpar    = apply(pars,2,mean)
    vpar    = var(pars)
    ms      = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    weight  = ldt(y[t],nu,0.0,exp(mus/2))
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    # Propagating
    if (delta<1){
      ms1 = ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    }else{
      ms1 = ms
    }
    xt   = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    w    = ldt(y[t],nu,0.0,exp(xt/2))-weight[k]
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=T,prob=w)
    xs   = xt[ind]
    pars = ms1[ind,]  	 
    # Storing quantiles
    quants[t,1,] = quantile(pars[,1],c(0.5,0.025,0.975))
    quants[t,2,] = quantile(pars[,2],c(0.5,0.025,0.975))
    quants[t,3,] = quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    quants[t,4,] = quantile(exp(xs/2),c(0.5,0.025,0.975))  
  }
  return(list(like=like,quants=quants))
}

# Dataset
y = c(4.204,-0.59,-9.376,5.45,-0.894,8.763,4.712,0.971,-1.953,-7.062,6.92,-3.768,-1.048,3.114,1.737,10.435,
      8.614,4.674,16.528,9.496,-2.336,-11.128,8.531,2.781,-3.944,-1.35,21.73,2.791,-2.321,-3.341,-1.32,10.403,
      -0.752,0.152,19.362,10.984,14.865,-4.578,-2.684,2.49,9.102,20.101,15.051,4.658,-9.614,-35.18,-15.647,12.83,
      17.477,4.879,10.91,2.484,-2.928,-17.457,0.533,4.33,-18.327,-19.138,-3.793,-9.681,0.856,17.662,-7.873,-9.246,
      -12.745,12.801,-8.388,1.551,-39.475,12.516,-9.409,-9.298,-26.136,2.564,-10.677,-26.881,-43.872,10.71,42.286,35.839,
      -10.269,-19.177,-8.48,12.503,-0.823,-23.144,4.879,44.469,13.193,7.987,-8.135,10.705,-21.556,-10.677,11.935,-2.02,
      14.864,-8.651,8.208,-3.99,-7.859,0.752,-9.217,4.715,-1.191,-4.812,15.616,7.652,5.466,-1.609,-4.307,8.179,
      2.584,5.073,13.231,6.213,7.949,6.15,6.125,1.507,2.899,-1.6,-0.324,-7.085,4.779,0.665,15.415,7.696,
      -4.325,8.516,6.661,5.877,14.826,-3.607,-4.338,-8.04,-1.17,0.282,9.896,-10.37,-15.598,-3.129,0,-2.699,
      -4.348,1.887,-36.248,15.092,-3.565,26.299,3.97,-0.299,2.256,9.557,-7.804,1.324,-7.765,3.359,-18.774,-0.733,
      6.062,-9.809,15.635,-8.224,16.733,-2.774,-2.467,5.081,-6.392,-0.995,4.815,-4.256,-21.98,5.984,6.528,1.109,
      4.247,0.709,-4.925,0.757,-0.378,-0.38,-0.841,-11.054,-1.308,13.441,-1.562,3.101,-1.616,-14.31,-3.312,2.25,
      1.861,-4.235,-8.108,-2.682,8.338,4.784,1.439,0.947,4.969,6.596,-0.427,5.326,10.863,4.668,5.396,-2.377,
      4.37,5.379,-8.53,1.361,5.521,-5.981,-4.914,6.88,-1.022,-4.197,3.44,-0.349,1.044,5.977,-1.659,1.658,
      -1.059,3.625,1.605,1.517,-1.917,9.237,-4.451,8.338,1.137,-2.807,2.032,8.269,3.177,-2.353,-1.869,3.756,
      5.599,-5.338,-0.997,-2.423,8.36,-4.041,-1.592,-8.361,-8.111,-5.903,-5.557,3.578,8.672,-3.909,-4.693,-2.136,
      -1.085,4.34,7.133,-3.674,0.407,0,-4.197,3.233,-3.922,-5.225,10.949,-3.175,18.59,-1.137,-7.247,3.216,
      0,8.004,-12.595,5.548,-1.621,-6.062,6.389,-2.337,-5.201,0.355,6.612,0.336,1.002,2,5.459,7.523,
      4.639,4.704,1.289,4.493,1.793,-6.881,1.368,3.733,0.73,2.34,0,5.564,5.381,6.231,-0.899,1.589,
      -3.436,-0.233,5.755,6.495,4.107,-8.719,-0.445,7.315,-2.123,-5.514,8.475,-3.871,4.502,6.494,0.199,0.198,
      1.376,-1.594,12.09,4.522,-5.293,0,0.361,2.7,2.457,-0.521,3.12,-3.473,4.492,10.825,10.322,-1.123,
      6.096,7.631,7.206,14.256,-4.884,17.054,-4.445,-5.85,6.19,-3.785,8.255,3.014,6.22,4.664,-4.363,6.359,
      -1.914,3.372,-3.593,3.358,-5.03,-6.721,13.953,8.04,-6.252,7.114,9.609,-3.193,-7.364,7.161,4.605,-3.993
      -8.276,6.481,-0.419,2.079,-11.415,2.978,5.916,8.432,5.188,3.315,1.809,-6.67,-7.352,0.415,7.573,-4.917,
      2.211,-4.265,0.207,-1.049,0.84,1.247,4.681,1.186,4.42,3.169,2.357,12.294,-0.962,2.701,2.63,2.727,
      -2.574,-0.924,1.238,0.766,-3.101,2.349,12.634,8.371,-14.786,4.572,0.973,-1.401,0.842,3.706,-8.628,-2.846,
      -10.738,-0.512,2.703,0,-4.988,-7.696,-0.19,-5.942,5.557,-1.165,4.428,6.56,6.322,-1.342,6.221,-3.884,
      -0.333,-0.501,3.619,-11.353,-4.279,-9.151,10.935,3.152,-6.211,5.89,9.832,3.289,1.935,-7.288,3.213,5.381,
      6.454,-5.506,-1.114,4.229,-2.014,5.365,-1.504,6.084,0.857,-2.593,2.508,-5.749,-1.683,-0.868,4.445,0.747,
      6.816,-2.12,5.692,1.261,6.361,0.502,0.175,4.664,-0.843,-6.96,4.587,4.147,11.733,1.062,-3.879,4.218,
      -2.793,-5.141,1.049,4.792,-9.825,2.513,-10.549,-11.648,0.762,11.517,2.715,-9.465,-0.709,-1.866,-0.845,11.138,
      -9.821,1.757,21.766,1.15,3.176,-5.618,-1.785,-7.45,-8.844,-0.399,-0.578,10.435,-7.694,-2.151,-2.198,-0.296,
      2.088,10.829,4.753,-5.255,-1.205,-6.228,6.435,5.791,-3.647,-3.895,-4.401,-2.797,2.244,-1.034,-4.872,-2.496,
      -6.842,-2.625,5.384,-3.449,-7.087,2.639,13.005,2.092,7.985,1.908,2.725,6.728,5.568,8.565,4.1,9.267,
      -2.472,1.921,-12.943,14.979,-0.242,-5.43,2.542,5.256,-1.003,-2.243,6.929,4.537,1.65,-3.931,-4.273,5.407,
      0.337,-3.839,6.588,7.107,-4.025,-4.045,-3.616,-8.49,0.421,-3.028,8.955,-5.951,6.981,3.701,-2.871,-0.198,
      -5.087,-4.877,-4.284,-6.178,-5.543,2.861,-13.571,-12.691,-17.956,20.212,-2.985,-9.457,14.609,14.215,4.161,0,
      -1.368,15.581,-10.774,-2.139,-4.082,8.961,0.26,-3.64,19.399,-4.801,-0.894,2.353,-5.497,12.285,-4.704,-2.791,
      2.929,-1.396,-5.289,10.224,-5.07,-4.596,-1.849,9.23,0.694,5.032,-3.834,-1.856,-2.995,-2.219,1.091,-0.5,
      -8.633,-2.21,5.277,11.39,1.189,-3.316,6.004,1.159,-1.345,-10.72,0.577,0,2.359,-3.372,4.256,2.317,
      1.415,1.765,3.44,2.858,-3.438,-4.793,-3.256,9.309,8.289,-7.649,-4.321,-1.308,5.628,3.961,8.822,-3.16,
      -2.086,2.81,13.382,0.202,0,8.599,1.487,-2.073,0.19,-4.79,-1.424,-8.776,-0.359,-1.146,10.697,-3.712,
      8.555,0.598,2.047,1.176,-3.572,4.116,3.285,13.496,0.464,13.874,8.347,2.536,8.701,4.718,-2.472,6.343,
      -7.902,7.072,-9.781,2.475,4.026,-1.914,10.754,2.592,-7.525,-4.466,6.288,1.133,-4.373,-0.235,0,7.564,
      -0.664,3.091,-2.777,2.232,12.048,0.469,-7.336,0,3.408,2.041,3.571,-4.996,-4.957,0.217,13,11.101,
      -2.611,9.682,1.282,-0.479,2.059,2.269,-10.399,7.423,-8.172,5.909,8.482,4.28,15.332,3.069,2.039,-1.32,
      0.482,5.071,8.587,5.34,-1.078,-26.094,-11.154,4.837,2.241,0.828,-11.088,0,3.946,5.162,-2.597,-5.716,
      7.971,0.574,3.103,0.355,7.79,-5.579,-1.853,9.377,11.123,-4.857,12.502,-1.507,-1.249,-2.011,11.301,4.679,
      -3.55,-0.048,3.969,-0.39,8.022,0.543,3.372,-14.431,-12.035,-4.696,5.153,5.568,10.928,6.81,2.362,1.603,
      9.447,-4.298,-1.019,2.194,-6.896,-0.542,-6.357,17.392,-1.647,4.387,-3.002,1.148,-0.327,2.489,-1.621,-3.323,
      6.336,-1.936,8.129,3.401,0.728,-2.349,6.478,1.669,2.318,3.839,2.832,-0.254,-1.792,1.166,1.408,7.082,
      2.704,-2.229,-4.519,-4.867,4.365,-5.718,7.736,-1.249,-2.576,1.546,-6.062,11.119,0.975,6.12,-0.623,3.637,
      3.509,-2.117,4.551,-0.212,8.597,-0.787,5.946,7.648,6.389,-1.642,3.686,-0.806,6.878,4.721,-4.769,1.058,
      9.556,6.127,7.226,-4.529,4.572,-0.606,-3.064,11.189,8.419,7.381,7.96,-11.411,8.807,-5.183,13.377,-0.271,
      5.47,0.322,10.65,-1.167,-2.151,8.614,-1.26,-11.151,-0.172,9.509,3.233,12.443,2.78,-4.447,10.101,-4.862,
      -3.562,10.549,-3.283,2.994,5.71,13.353,-4.048,17.596)
n = length(y)

# Prior hyperparameters
# ---------------------
m0    = 0
C0    = 10
nu0   = 3
tau20 = 0.01
b0    = c(0,1)
B0    = c(10,10)/tau20
sC0   = sqrt(C0)
sB0   = sqrt(B0)

set.seed(246521)
N        = 100000
delta    = 0.975
xs       = rnorm(N,m0,sC0)
tau2s    = 1/rgamma(N,nu0/2,nu0*tau20/2)
taus     = sqrt(tau2s)
alphas   = rnorm(N,b0[1],taus*sB0[1])
betas    = rnorm(N,b0[2],taus*sB0[2])
nus      = 2:20
nnu      = length(nus)
pf       = array(0,c(1+nnu,n,4,3))
like     = matrix(0,n,1+nnu)
run      = LW(y,alphas,betas,tau2s,xs,delta)
pf[1,,,] = run$quants
like[,1] = run$like
for (i in 1:nnu){
  print(nus[i])
  run      = LWt(y,nus[i],alphas,betas,tau2s,xs,delta)
  pf[1+i,,,] = run$quants
  like[,i+1] = run$like
}

pred = apply(log(like),2,cumsum)
pmp  = matrix(0,n,1+nnu)
for (t in 1:n){
  prob = exp(pred[t,]-max(pred[t,]))
  pmp[t,] = prob/sum(prob)
}

pdf(file="pmp.pdf",height=10,width=10)
par(mfrow=c(2,2))
tt = c(1,n/2,3*n/4,n)
for (t in tt){
  plot(1:20,pmp[t,],type="h",axes=FALSE,xlab="Degrees of freedom",ylab="PMP",lwd=3)
  axis(2);box();axis(1,at=1:20,label=c("N",nus))
  title(paste("t=",t,sep=""))
}
dev.off()  

pdf(file="sv-par.pdf",height=10,width=12)
par(mfrow=c(3,4))
ind = 200:n
alphas1 = alphas[abs(alphas)<20]
hist(alphas1,breaks=seq(min(alphas1),max(alphas1),length=30),prob=TRUE,xlab="",main=expression(p(alpha)));box()
plot(ind,pf[1,ind,1,1],ylim=range(pf[c(1,12,18),ind,1,]),type="l",xlab="Months",main="Normal",ylab=expression(alpha))
for (i in 2:3) lines(ind,pf[1,ind,1,i])
plot(ind,pf[12,ind,1,1],ylim=range(pf[c(1,12,18),ind,1,]),type="l",xlab="Months",main=expression(t[12]),ylab=expression(alpha))
for (i in 2:3) lines(ind,pf[12,ind,1,i])
plot(ind,pf[18,ind,1,1],ylim=range(pf[c(1,12,18),ind,1,]),type="l",xlab="Months",main=expression(t[18]),ylab=expression(alpha))
for (i in 2:3) lines(ind,pf[18,ind,1,i])

betas1 = betas[abs(betas)<20]
hist(betas1,breaks=seq(min(betas1),max(betas1),length=30),prob=TRUE,xlab="",main=expression(p(beta)));box()
plot(ind,pf[1,ind,2,1],ylim=range(pf[c(1,12,18),ind,2,]),type="l",xlab="Months",main="Normal",ylab=expression(beta))
for (i in 2:3) lines(ind,pf[1,ind,2,i])
plot(ind,pf[12,ind,2,1],ylim=range(pf[c(1,12,18),ind,2,]),type="l",xlab="Months",main=expression(t[12]),ylab=expression(beta))
for (i in 2:3) lines(ind,pf[12,ind,2,i])
plot(ind,pf[18,ind,2,1],ylim=range(pf[c(1,12,18),ind,2,]),type="l",xlab="Months",main=expression(t[18]),ylab=expression(beta))
for (i in 2:3) lines(ind,pf[12,ind,2,i])

tau2s1 = tau2s[abs(tau2s)<0.2]
hist(tau2s1,breaks=seq(min(tau2s1),max(tau2s1),length=30),prob=TRUE,xlab="",main=expression(p(tau^2)));box()
plot(ind,pf[1,ind,3,1],ylim=range(pf[c(1,12,18),ind,3,]),type="l",xlab="Months",main="Normal",ylab=expression(tau^2))
for (i in 2:3) lines(ind,pf[1,ind,3,i])
plot(ind,pf[12,ind,3,1],ylim=range(pf[c(1,12,18),ind,3,]),type="l",xlab="Months",main=expression(t[12]),ylab=expression(tau^2))
for (i in 2:3) lines(ind,pf[12,ind,3,i])
plot(ind,pf[18,ind,3,1],ylim=range(pf[c(1,12,18),ind,3,]),type="l",xlab="Months",main=expression(t[18]),ylab=expression(tau^2))
for (i in 2:3) lines(ind,pf[12,ind,3,i])
dev.off()

svall = matrix(0,n,3)
for (t in 1:n)
  svall[t,] = apply(pf[,t,4,]*matrix(pmp[t,],20,3),2,sum)

pdf(file="sv-vars.pdf",height=8,width=10)
par(mfrow=c(2,2))
ts.plot(y,xlab="Months",ylab="Returns",main="(a)",type="l",lwd=2)
ts.plot(pf[12,,4,],xlab="Months",ylab="Standard deviation",main="(b)",type="l",lwd=2,ylim=range(pf[,,4,]))
ts.plot(pf[18,,4,],xlab="Months",ylab="Standard deviation",main="(c)",type="l",lwd=2,ylim=range(pf[,,4,]))
ts.plot(svall,xlab="Months",ylab="Standard deviation",main="(d)",type="l",lwd=2,ylim=range(pf[,,4,]))
dev.off()





