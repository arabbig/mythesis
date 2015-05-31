#################################################################################################
#
#  LOG STOCHASTIC VOLATILITY AR(1) MODEL
#
#  AUXILIARY PARTICLE FILTER + PARAMETER LEARNING 
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
pri        = function(x,m0,sC0){dnorm(x,m0,sC0)}
likelihood = function(y,x){dnorm(y,0,exp(x/2))}
rlike      = function(x){rnorm(1,0,exp(x/2))}
post       = function(y,x,m0,C0){pri(x,m0,sC0)*likelihood(y,x)}
quant025   = function(x){quantile(x,0.025)}
quant975   = function(x){quantile(x,0.975)}

# Simulating the data
# -------------------
set.seed(12345)
n     =  500
alpha = -0.0031
beta  =  0.9951
tau2  =  0.0074
tau   = sqrt(tau2)
y     = rep(0,n)
x     = rep(0,n)
x[1]  = alpha/(1-beta)
y[1]  = rlike(x[1])
for (t in 2:n){
  x[t] = rnorm(1,alpha+beta*x[t-1],tau)
  y[t] = rlike(x[t])
}
alpha.true = alpha
beta.true  = beta
tau2.true  = tau2

# Data and prior hyperparameters 
# ------------------------------
m0      = 0.0
C0      = 0.1
sC0     = sqrt(C0)
ealpha  = alpha
valpha  = 0.01
ephi    = beta
vphi    = 0.01
nu      = 3
lambda  = tau2

# Liu and West filter
# -------------------
set.seed(8642)
N      = 5000
xs     = rnorm(N,m0,sC0)
pars   = cbind(rnorm(N,ealpha,sqrt(valpha)),rnorm(N,ephi,sqrt(vphi)),log(1/rgamma(N,nu/2,nu*lambda/2)))
delta  = 0.75
h2     = 1-((3*delta-1)/(2*delta))^2
a      = sqrt(1-h2)
parss  = array(0,c(N,3,n))
xss    = NULL
ws     = NULL
ESS    = NULL
par(mfrow=c(1,1))
for (t in 1:n){
  mpar        = apply(pars,2,mean)
  vpar        = var(pars)
  ms          = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
  mus         = pars[,1]+pars[,2]*xs # mean x_(t+1) by auto regressive
  #mus         = ms[,1] + ms[,2]*xs # try it myself
  weight      = likelihood(y[t],mus) + 0.5 # for resampling ... with p(y_(t+1)|g,m)
  k           = sample(1:N,size=N,replace=T,prob=weight)
  ms1         = ms[k,] + matrix(rnorm(3*N),N,3)%*%chol(h2*vpar) # sampling param_(t+1)
  xt          = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2)) # sampling x_(t+1)
  w           = likelihood(y[t],xt)/(likelihood(y[t],mus[k])+0.5)
  w           = w/sum(w)
  # Resampling to embbed weight into new sample
  ind         = sample(1:N,size=N,replace=T,prob=w)
  xs          = xt[ind]
  pars        = ms1[ind,]
  xss         = rbind(xss,xs)
  parss[,,t]  = pars 
  ws          = rbind(ws,w)
  cv2         = var(w)/(mean(w)^2)
  ESS         = c(ESS,N/(1+cv2))
  # ts.plot(ESS,xlim=c(1,n))
}

# Posterior summary
# -----------------
mvol   = apply(exp(xss),1,mean)
lvol   = apply(exp(xss),1,quant025)
uvol   = apply(exp(xss),1,quant975)
malpha = apply(parss[,1,],2,mean)
lalpha = apply(parss[,1,],2,quant025)
ualpha = apply(parss[,1,],2,quant975)
mbeta  = apply(parss[,2,],2,mean)
lbeta  = apply(parss[,2,],2,quant025)
ubeta  = apply(parss[,2,],2,quant975)
mtau2  = apply(exp(parss[,3,]),2,mean)
ltau2  = apply(exp(parss[,3,]),2,quant025)
utau2  = apply(exp(parss[,3,]),2,quant975)

# Graphical analysis
# ------------------
pdf(file="vol.pdf",width=10,height=10)
par(mfrow=c(2,1))
ts.plot(y,xlab="time",ylab="",main=expression(y[t]))
lines(mvol,lwd=1.25)
lines(lvol,lwd=1.25)
lines(uvol,lwd=1.25)
lines(exp(x),lwd=1.25,col=2)
dev.off()


pdf(file="par.pdf",width=10,height=4)
par(mfrow=c(1,3))
ts.plot(malpha,ylim=range(lalpha,ualpha),ylab="",main=expression(alpha))
lines(lalpha,lwd=1.25)
lines(malpha,lwd=1.25)
lines(ualpha,lwd=1.25)
abline(h=alpha.true,col=2,lwd=2)

ts.plot(mbeta,ylim=range(lbeta,ubeta),ylab="",main=expression(beta))
lines(lbeta,lwd=1.25)
lines(mbeta,lwd=1.25)
lines(ubeta,lwd=1.25)
abline(h=beta.true,col=2,lwd=2)

ts.plot(mtau2,ylim=range(ltau2,utau2),ylab="",main=expression(tau^2))
lines(ltau2,lwd=1.25)
lines(mtau2,lwd=1.25)
lines(utau2,lwd=1.25)
abline(h=tau2.true,col=2,lwd=2)
dev.off()
