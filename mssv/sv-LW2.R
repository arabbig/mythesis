require(msm)
require(gtools)
pri        = function(x,m0,sC0){dnorm(x,m0,sC0)}
likelihood = function(y,x){dnorm(y,0,exp(x/2))}
rlike      = function(x){rnorm(1,0,exp(x/2))}
post       = function(y,x,m0,C0){pri(x,m0,sC0)*likelihood(y,x)}
quant025   = function(x){quantile(x,0.025)}
quant975   = function(x){quantile(x,0.975)}
tranf  = function(v) {
  c(v[1],
    log(v[2]/(1-v[2])), # p00_tf = log(p00/(1-p00))
    log(v[3]/(1-v[3])), # p10_tf = log(p10/(1-p10))
    log(v[4]), # tau2_tf = log(tau2)
    v[5], # gam1_tf = gam1 
    log(v[6])) # gam2_tf = log(gam2)
}
itranf  = function(v) {
  c(v[1],
    exp(v[2])/(1+exp(v[2])), # p00_tf = log(p00/(1-p00))
    exp(v[3])/(1+exp(v[3])), # p10_tf = log(p10/(1-p10))
    exp(v[4]), # tau2_tf = log(tau2)
    v[5], # gam1_tf = gam1 
    exp(v[6])) # gam2_tf = log(gam2)
}   
# Simulating the data
# -------------------
set.seed(12345)
n     =  500
alpha = c(-2.5,-1)
beta  =  0.5
tau2  =  0.0074
tau   = sqrt(tau2)
y     = rep(0,n)
xtrue = rep(0,n)
strue     = rep(0,n)
xtrue[1]  = alpha[0]/(1-beta[0])
y[1]  = rlike(xtrue[1])
strue[1]  = 0

for (t in 2:n){
  if(strue[t-1]){
    strue[t] = rbinom(1,1,0.985)
  } else {
    strue[t] = rbinom(1,1,0.01)
  }
  xtrue[t] = rnorm(1,alpha[strue[t]+1]+beta*xtrue[t-1],tau)
  y[t] = rlike(xtrue[t])
}
alpha.true = alpha
beta.true  = beta
tau2.true  = tau2

# Data and prior hyperparameters 
# ------------------------------
m0      = 0.0
a0      = -4
b0      = 0.0
A0      = 3
B0      = 3
C0      = 1
nu      = 4.01
tau0   = 0.1
g0      = 3
G0      = 3
# Liu and West filter
# -------------------
set.seed(3210)
N      = 5000
x     = rnorm(N,m0,sqrt(C0))
s     = rbinom(N,1,0.5)
#initialize param
tau2  = 1/rgamma(N,nu/2,nu/2*tau0^2)
phi   = rtnorm(N,b0,sqrt(B0),lower = -1,upper = 1)
p00   = rdirichlet(N,c(0.5,0.5))[,1]
p10   = rdirichlet(N,c(0.5,0.5))[,1]
gam1  = rtnorm(N,a0,sqrt(A0), upper = -1)
gam2  = rtnorm(N,g0,sqrt(G0), lower = 2)
para   = cbind(phi,p00,p10,tau2,gam1,gam2)


delta  = 0.75
h2     = 1-((3*delta-1)/(2*delta))^2
a      = sqrt(1-h2)
parss  = array(0,c(N,6,n))
xss    = NULL
ss     = NULL
ws     = NULL
ESS    = NULL
par(mfrow=c(1,1))
for (t in 1:n){
  if(!(t%%10)) print(t)
  # par   = cbind(phi,p00,p10,tau2,gam1,gam2)
  # calculate most likely next state s 
  p           = apply(cbind(para[,'p00'],para[,'p10'],s),1,function(v){ifelse(v[3],v[2],v[1])})
  s_argm      = as.integer(p < 0.5)
    
  # calulate mu (t+1)
  alp       = para[,'gam1'] + s_argm*para[,'gam2'];
  mx          = alp + para[,'phi']*x
  # resample k, calculate mu (state required)
  weight      = likelihood(y[t],mx) # y only depend on x, not alpha,beta,tau
  k           = sample(1:N,size=N,replace=T,prob=weight)
  
  # sample of theta(t+1)

  par_tf      = t(apply(para,1,tranf))
  mpar        = apply(par_tf,2,mean)
  mix         = a*par_tf+(1-a)*matrix(mpar,N,6,byrow=T)
  V           = var(par_tf)
  par_tfnext  = mix[k,] + matrix(rnorm(6*N),N,6)%*%chol(h2*V)
  par_next    = t(apply(par_tfnext,1,itranf))

  # sample of s(t+1), x(t+1)
  p           = apply(cbind(par_next[,'p00'],par_next[,'p10'],s[k]),1,function(v){ifelse(v[3],v[2],v[1])})
  #p           = p[k]
  s_next      = rbinom(N,1,1-p)
  alp_next  = par_next[,'gam1'] + s_next*par_next[,'gam2']
  x_next      = rnorm(N,alp_next+par_next[,'phi']*x[k],exp(par_next[,'tau2']/2))
  w           = likelihood(y[t],x_next)/likelihood(y[t],mx[k])
  w           = w/sum(w)
  
  # resample
  ind         = sample(1:N,size=N,replace=T,prob=w)
  x           = x_next[ind]
  s           = s_next[ind]
  para         = par_next[ind,]
  xss         = rbind(xss,x)
  ss          = rbind(ss,s)
  parss[,,t]  = as.matrix(para) 
  ws          = rbind(ws,w)
  cv2         = var(w)/(mean(w)^2)
  ESS         = c(ESS,N/(1+cv2))
  #ts.plot(ESS,xlim=c(1,n))
}

# Posterior summary
# -----------------
mlogvol   = apply(xss,1,mean)
mvol   = apply(exp(xss/2),1,mean)
mstate = apply(ss,1,mean)
# Graphical analysis
par(mfrow=c(3,1));
ts.plot(log(y^2), col='blue', main='log-variance');
lines(xtrue,col='black')
lines(mlogvol,col='blue')
ts.plot(y,col='black',main='return')
lines(mvol,col='red')
ts.plot(mstate,col='red', main='state',ylim=c(-0.1,1.1));
lines(strue,col='black');