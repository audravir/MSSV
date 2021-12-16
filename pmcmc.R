rm(list=ls(all=TRUE))
#set.seed(1)
q025  = function(x) quantile(x,0.025)
q975  = function(x) quantile(x,0.975)
bst2b = function(x) 1/(1+exp(-x))
b2bst = function(x) log(x)-log(1-x)

nn       = 1000
alpha.tr = 0.01
beta.tr  = 0.95
tau.tr   = 0.5
sig.tr   = 1
x.tr     = rep(NA,nn)
signaltonoise=tau.tr/sig.tr

x0.tr    = rnorm(1,alpha.tr/(1-beta.tr),sqrt(tau.tr^2/(1-beta.tr^2)))
x.tr[1]  = rnorm(1,alpha.tr+beta.tr*x0.tr,tau.tr)
for(t in 2:nn) x.tr[t] = rnorm(1,alpha.tr+beta.tr*x.tr[t-1],tau.tr)
y        = rnorm(nn,x.tr,sig.tr)

par(mfrow=c(1,1))
plot(y,type='l',lwd=2,col='gray80')
lines(x.tr,col=2,lwd=2)

N   = 100000 #number of particles
# ONLY filtering
# resample-propagate-resample
xss = matrix(NA,ncol=3,nrow=nn)
x0  = rnorm(N,alpha.tr/(1-beta.tr),sqrt(tau.tr^2/(1-beta.tr^2)))
for(t in 1:nn){
    ws  = dnorm(y[t],alpha.tr+beta.tr*x0,sig.tr)
    ind = sample(N,prob = ws,replace = TRUE)
    xs  = rnorm(N,alpha.tr+beta.tr*x0[ind],tau.tr)
    ws2 = dnorm(y[t],xs,sig.tr)/ws[ind]
    ind = sample(N,prob = ws2,replace = TRUE)
    x0  = xs[ind]
    xss[t,] = c(q025(x0),median(x0),q975(x0))
}
