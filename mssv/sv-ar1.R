############################################################################################
#
#  STOCHASTIC VOLATILITY MODEL
#
############################################################################################
#
# HEDIBERT FREITAS LOPES
# Associate Professor of Econometrics and Statistics
# The University of Chicago Booth School of Business
# 5807 South Woodlawn Avenue
# Chicago, Illinois, 60637
# Email : hlopes@ChicagoGSB.edu
# URL: http://faculty.chicagobooth.edu/hedibert.lopes/research/
#
############################################################################################
rm(list=ls())
source("sv-ar1-routines.R")

# Simulating the data
# -------------------
set.seed(1243)
mu   = -0.00645
phi  =  0.99
tau2 =  0.15^2
tau  = sqrt(tau2)
h0   = 0.0
n    = 500
h    = rep(0,n)
h[1] = rnorm(1,mu+phi*h0,tau)
for (t in 2:n) 
  h[t] = rnorm(1,mu+phi*h[t-1],tau)
vol = exp(h)
sd  = sqrt(vol)
y   = rnorm(n,0,sd)

# True values
# -----------
ptr = c(mu,phi,tau2)
h0tr = h0
htr = h
vtr = vol
sdtr= sd

pdf(file="simulateddata.pdf",width=12,height=12)
par(mfrow=c(2,1))
plot(y,type="l",xlab="time",ylab="")
plot(vol,type="l",xlab="time",ylab="")
dev.off()

# Prior hyperparameters
# ---------------------
theta0 = c(mu,phi)
V0     = diag(100,2)
nu0    = 10.0
s02    = (nu0-2)/nu0*tau2
m0     = h0
C0     = 100.0
c(sqrt((nu0*s02/2)/(nu0/2-1)),sqrt(((nu0*s02/2)/(nu0/2-1))^2/(nu0/2-2)))
iV0    = solve(V0)
iC0    = 1/C0
iC0m0  = iC0*m0

# General MCMC setup
# ------------------
M0 = 1000
M  = 3000

#####################################################################################
# STOCHASTIC VOLATILITY: INDIVIDUAL MOVES VIA RANDOM WALK METROPOLIS                #            
#####################################################################################
# Initial values and specific MCMC setup
# --------------------------------------
h0   = h0tr
h    = htr
mu   = ptr[1]
phi  = ptr[2]
tau2 = ptr[3]
vh   = 0.1

set.seed(12578)
hs = NULL
ps = NULL
for (iter in 1:(M0+M)){
  h = svol.rw(y,mu,phi,tau2,h0,h,vh)
  var   = 1/(iC0+phi^2/tau2)
  mean  = var*(iC0m0+phi*(h[1]-mu)/tau2)
  h0    = rnorm(1,mean,sqrt(var))
  X     = cbind(1,c(h0,h[1:(n-1)]))
  par   = fixedpar(h[1:n],X,theta0,iV0,nu0,s02)
  mu    = par[1]
  phi   = par[2]
  tau2  = par[3]
  ps    = rbind(ps,par)
  hs    = rbind(hs,h)
}
hs.rw = hs
ps.rw = ps
vol   = exp(hs[(M0+1):(M0+M),])
mvol  = apply(vol,2,median)
lvol  = apply(vol,2,quant05)
uvol  = apply(vol,2,quant95)

names = c("mu","phi","tau2")
pdf(file="par-rw.pdf",width=10,height=7)
par(mfrow=c(3,3))
for (i in 1:3){
  plot(ps[(M0+1):(M0+M),i],type="l",xlab="iteration",ylab="",main=names[i])
  acf(ps[(M0+1):(M0+M),i],main="")
  hist(ps[(M0+1):(M0+M),i],xlab="",ylab="",main="")
}
dev.off()

pdf(file="vol-rw.pdf",width=10,height=7)
par(mfrow=c(1,1))
plot(vtr,type="l",ylim=range(lvol,uvol,vtr),xlab="time",ylab="")
lines(lvol,col=2)
lines(uvol,col=2)
lines(mvol,col=4)
lines(vtr,col=1)
legend(250,7,legend=c("TRUE","POSTERIOR MEAN","90%CI"),col=c(1,4,2),lty=rep(1,3),lwd=rep(1,3),bty="n")
dev.off()

#####################################################################################
# STOCHASTIC VOLATILITY: INDIVIDUAL MOVES VIA INDEPENDENT METROPOLIS-HASTINGS       #
#####################################################################################
# Initial values and specific MCMC setup
# --------------------------------------
h0   = h0tr
h    = htr
mu   = ptr[1]
phi  = ptr[2]
tau2 = ptr[3]
vh   = 0.1

set.seed(12578)
hs = NULL
ps = NULL
for (iter in 1:(M0+M)){
  h = svolsingle(y,mu,phi,tau2,h0,h)
  var   = 1/(iC0+phi^2/tau2)
  mean  = var*(iC0m0+phi*(h[1]-mu)/tau2)
  h0    = rnorm(1,mean,sqrt(var))
  X     = cbind(1,c(h0,h[1:(n-1)]))
  par   = fixedpar(h[1:n],X,theta0,iV0,nu0,s02)
  mu    = par[1]
  phi   = par[2]
  tau2  = par[3]
  ps    = rbind(ps,par)
  hs    = rbind(hs,h)
}
hs.ind = hs
ps.ind = ps
vol    = exp(hs[(M0+1):(M0+M),])
mvol   = apply(vol,2,median)
lvol   = apply(vol,2,quant05)
uvol   = apply(vol,2,quant95)

names = c("mu","phi","tau2")
pdf(file="par-ind.pdf",width=10,height=7)
par(mfrow=c(3,3))
for (i in 1:3){
  plot(ps[(M0+1):(M0+M),i],type="l",xlab="iteration",ylab="",main=names[i])
  acf(ps[(M0+1):(M0+M),i],main="")
  hist(ps[(M0+1):(M0+M),i],xlab="",ylab="",main="")
}
dev.off()

pdf(file="vol-ind.pdf",width=10,height=7)
par(mfrow=c(1,1))
plot(vtr,type="l",ylim=range(lvol,uvol,vtr),xlab="time",ylab="")
lines(lvol,col=2)
lines(uvol,col=2)
lines(mvol,col=4)
lines(vtr,col=1)
legend(250,7,legend=c("TRUE","POSTERIOR MEAN","90%CI"),col=c(1,4,2),lty=rep(1,3),lwd=rep(1,3),bty="n")
dev.off()

acfind=NULL
for (t in 1:n)
  acfind = rbind(acfind,acf(hs.ind[(M0+1):(M0+M),t],plot=FALSE)$acf)
pdf(file="acf-ind.pdf",width=10,height=7)
ts.plot(t(acfind),xlab="lag",ylab="ACF",col=gray(0.8))
dev.off()

pdf(file="acf-rw-ind.pdf",width=10,height=7)
par(mfrow=c(1,2))
ts.plot(t(acfrw),xlab="lag",ylab="ACF",ylim=range(acfrw,acfind),main="RANDOM WALK",col=gray(0.8))
ts.plot(t(acfind),xlab="lag",ylab="ACF",ylim=range(acfrw,acfind),main="INDEPENDENT",col=gray(0.8))
dev.off()

#####################################################################################
# STOCHASTIC VOLATILITY: NORMAL APPROXIMATION + FORWARD FILTERING BACKWARD SAMPLING #            
#####################################################################################
x      = seq(-20,10,length=10000)
den    = exp(-(x)/2)*exp(-0.5*exp(x))*exp(x)/sqrt(2*pi)
norm   = dnorm(x,-1.270399,sqrt(4.934854))
pdf(file="normalapproximation.pdf",width=10,height=7)
par(mfrow=c(1,1))
plot(x,den,ylab="density",main="",type='n')
lines(x,den,col=1,lty=1,lwd=2)
lines(x,norm,col=2,lty=2,lwd=2)
legend(-15,0.2,legend=c("log chi^2_1","normal"),lty=1:2,col=1:2,lwd=c(2,2),bty="n")
dev.off()

# Transformations
# ---------------
y1   = log(y^2)+1.27
sig2 = pi^2/2

# Initial values and specific MCMC setup
# --------------------------------------
h0   = h0tr
h    = htr
mu   = ptr[1]
phi  = ptr[2]
tau2 = ptr[3]
vh   = 0.1

set.seed(12578)
hs = NULL
ps = NULL
for (iter in 1:(M0+M)){
  h = ffbsu(y1,1.0,0.0,sig2,mu,phi,tau2,m0,C0)  
  var   = 1/(iC0+phi^2/tau2)
  mean  = var*(iC0m0+phi*(h[1]-mu)/tau2)
  h0    = rnorm(1,mean,sqrt(var))
  X     = cbind(1,c(h0,h[1:(n-1)]))
  par   = fixedpar(h[1:n],X,theta0,iV0,nu0,s02)
  mu    = par[1]
  phi   = par[2]
  tau2  = par[3]
  ps    = rbind(ps,par)
  hs    = rbind(hs,h)
}
hs.nor = hs
ps.nor = ps
vol    = exp(hs[(M0+1):(M0+M),])
mvol   = apply(vol,2,median)
lvol   = apply(vol,2,quant05)
uvol   = apply(vol,2,quant95)

names = c("mu","phi","tau2")
pdf(file="par-normal.pdf",width=10,height=7)
par(mfrow=c(3,3))
for (i in 1:3){
  plot(ps[(M0+1):(M0+M),i],type="l",xlab="iteration",ylab="",main=names[i])
  acf(ps[(M0+1):(M0+M),i],main="")
  hist(ps[(M0+1):(M0+M),i],xlab="",ylab="",main="")
}
dev.off()

pdf(file="vol-normal.pdf",width=10,height=7)
par(mfrow=c(1,1))
plot(vtr,type="l",ylim=range(lvol,uvol,vtr),xlab="time",ylab="")
lines(lvol,col=2)
lines(uvol,col=2)
lines(mvol,col=4)
lines(vtr,col=1)
legend(250,7,legend=c("TRUE","POSTERIOR MEAN","90%CI"),col=c(1,4,2),lty=rep(1,3),lwd=rep(1,3),bty="n")
dev.off()

acfnor=NULL
for (t in 1:n)
  acfnor = rbind(acfnor,acf(hs.nor[(M0+1):(M0+M),t],plot=FALSE)$acf)
pdf(file="acf-normal.pdf",width=10,height=7)
ts.plot(t(acfnor),xlab="lag",ylab="ACF",col=gray(0.8))
dev.off()

pdf(file="acf-rw-ind-normal.pdf",width=14,height=7)
par(mfrow=c(1,3))
ts.plot(t(acfrw),xlab="lag",ylab="ACF",ylim=range(acfrw,acfind,acfnor),main="RANDOM WALK",col=gray(0.8))
ts.plot(t(acfind),xlab="lag",ylab="ACF",ylim=range(acfrw,acfind,acfnor),main="INDEPENDENT",col=gray(0.8))
ts.plot(t(acfnor),xlab="lag",ylab="ACF",ylim=range(acfrw,acfind,acfnor),main="NORMAL",col=gray(0.8))
dev.off()

#################################################################################################
# STOCHASTIC VOLATILITY: MIXTURE OF NORMALS APPROXIMATION + FORWARD FILTERING BACKWARD SAMPLING #            
#################################################################################################
# Defining the mixture of 7 normals that approximate the log-chi-square with one degree of freedom
# ------------------------------------------------------------------------------------------------
set.seed(1576)
ms  = c(-11.40039,-5.24321,-9.83726,1.50746,-0.65098,0.52478,-2.35859)
s2s = c(  5.79596, 2.61369, 5.17950,0.16735, 0.64009,0.34023, 1.26261)
qs  = c(  0.00730, 0.10556, 0.00002,0.04395, 0.34001,0.24566, 0.25750)
mm  = sum(qs*ms)
vv  = sum(qs*(s2s+ms^2))-(mm)^2
MM  = 10000
x  = seq(-20,10,length=MM)
den = exp(-(x)/2)*exp(-0.5*exp(x))*exp(x)/sqrt(2*pi)
mix = rep(0,MM)
for (i in 1:MM) 
  mix[i] = sum(qs*dnorm(x[i],ms,sqrt(s2s)))
norm   = dnorm(x,mm,sqrt(vv))

pdf(file="logchi2.pdf",width=10,height=7)
par(mfrow=c(1,1))
plot(x,den,ylab="density",main="",type='n')
lines(x,den,col=1,lty=1,lwd=2)
lines(x,mix,col=2,lty=2,lwd=2)
legend(-15,0.2,legend=c("log chi^2_1","mixture of 7 normals"),lty=1:2,col=1:2,lwd=c(2,2),bty="n")
dev.off()

# Initial values and specific MCMC setup
# --------------------------------------
y1   = log(y^2)
h0   = h0tr
h    = htr
mu   = ptr[1]
phi  = ptr[2]
tau2 = ptr[3]
vh   = 0.1

set.seed(12578)
hs = NULL
ps = NULL
for (iter in 1:(M0+M)){
  h     = svu(y1,h,mu,phi,tau2,m0,C0)
  var   = 1/(iC0+phi^2/tau2)
  mean  = var*(iC0m0+phi*(h[1]-mu)/tau2)
  h0    = rnorm(1,mean,sqrt(var))
  X     = cbind(1,c(h0,h[1:(n-1)]))
  par   = fixedpar(h[1:n],X,theta0,iV0,nu0,s02)
  mu    = par[1]
  phi   = par[2]
  tau2  = par[3]
  ps    = rbind(ps,par)
  hs    = rbind(hs,h)
}
hs.mix = hs
ps.mix = ps
vol    = exp(hs[(M0+1):(M0+M),])
mvol   = apply(vol,2,median)
lvol   = apply(vol,2,quant05)
uvol   = apply(vol,2,quant95)

names = c("mu","phi","tau2")
pdf(file="par-mixture.pdf",width=10,height=7)
par(mfrow=c(3,3))
for (i in 1:3){
  plot(ps[(M0+1):(M0+M),i],type="l",xlab="iteration",ylab="",main=names[i])
  acf(ps[(M0+1):(M0+M),i],main="")
  hist(ps[(M0+1):(M0+M),i],xlab="",ylab="",main="")
}
dev.off()

pdf(file="vol-mixture.pdf",width=10,height=7)
par(mfrow=c(1,1))
plot(vtr,type="l",ylim=range(lvol,uvol,vtr),xlab="time",ylab="")
lines(lvol,col=2)
lines(uvol,col=2)
lines(mvol,col=4)
lines(vtr,col=1)
legend(250,7,legend=c("TRUE","POSTERIOR MEAN","90%CI"),col=c(1,4,2),lty=rep(1,3),lwd=rep(1,3),bty="n")
dev.off()

acfmix=NULL
for (t in 1:n)
  acfmix = rbind(acfmix,acf(hs.nor[(M0+1):(M0+M),t],plot=FALSE)$acf)
pdf(file="acf-mixture.pdf",width=10,height=7)
ts.plot(t(acfmix),xlab="lag",ylab="ACF",col=gray(0.8))
dev.off()

pdf(file="acf-rw-ind-normal-mixture.pdf",width=14,height=9)
L = min(acfrw,acfind,acfnor,acfmix)
U = max(acfrw,acfind,acfnor,acfmix)
par(mfrow=c(2,2))
ts.plot(t(acfrw),xlab="lag",ylab="ACF",ylim=c(L,U),main="RANDOM WALK",col=gray(0.8))
ts.plot(t(acfind),xlab="lag",ylab="ACF",ylim=c(L,U),main="INDEPENDENT",col=gray(0.8))
ts.plot(t(acfnor),xlab="lag",ylab="ACF",ylim=c(L,U),main="NORMAL",col=gray(0.8))
ts.plot(t(acfmix),xlab="lag",ylab="ACF",ylim=c(L,U),main="MIXTURE",col=gray(0.8))
dev.off()


# Comparing the four schemes
# --------------------------
mvol.rw   = apply(exp(hs.rw[(M0+1):(M0+M),]),2,mean)
mvol.ind = apply(exp(hs.ind[(M0+1):(M0+M),]),2,mean)
mvol.nor  = apply(exp(hs.nor[(M0+1):(M0+M),]),2,mean)
mvol.mix  = apply(exp(hs.mix[(M0+1):(M0+M),]),2,mean)
L = min(mvol.rw,mvol.ind,mvol.nor,mvol.mix,vtr)
U = max(mvol.rw,mvol.ind,mvol.nor,mvol.mix,vtr)
pdf(file="volatilities.pdf",width=14,height=9)
ts.plot(mvol.rw,xlab="time",ylab="",ylim=c(L,U),col=2)
lines(mvol.ind,col=3)
lines(mvol.nor,col=4)
lines(mvol.mix,col=5)
lines(vtr,col=1)
legend(300,4,legend=c("TRUE","RANDOM WALK","INDEPENDENT","NORMAL","MIXTURE"),lty=rep(1,5),bty="n",col=1:5)
dev.off()

pdf(file="parameters.pdf",width=12,height=6)
par(mfrow=c(1,3))
plot(density(ps.rw[,1]),xlab="",ylab="Density",main=expression(mu),col=2,ylim=c(0,70))
lines(density(ps.ind[,1]),col=3)
lines(density(ps.nor[,1]),col=4)
lines(density(ps.mix[,1]),col=5)
abline(v=ptr[1],lwd=2)

plot(density(ps.rw[,2]),xlab="",ylab="Density",main=expression(phi),col=2,ylim=c(0,70))
lines(density(ps.ind[,2]),col=3)
lines(density(ps.nor[,2]),col=4)
lines(density(ps.mix[,2]),col=5)
abline(v=ptr[2],lwd=2)

plot(density(ps.rw[,3]),xlab="",ylab="Density",main=expression(tau^2),col=2,ylim=c(0,70))
lines(density(ps.ind[,3]),col=3)
lines(density(ps.nor[,3]),col=4)
lines(density(ps.mix[,3]),col=5)
abline(v=ptr[3],lwd=2)
dev.off()
