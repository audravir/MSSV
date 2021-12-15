rm(list=ls(all=TRUE))
#set.seed(1)
q025  = function(x) quantile(x,0.025)
q975  = function(x) quantile(x,0.975)
bst2b = function(x) 1/(1+exp(-x))
b2bst = function(x) log(x)-log(1-x)

nn       = 500
alpha.tr = 0.01
beta.tr  = 0.95
tau.tr   = 0.5
sig.tr   = 0.5
x.tr     = rep(NA,nn)

x0.tr    = rnorm(1,alpha.tr/(1-beta.tr),sqrt(tau.tr^2/(1-beta.tr^2)))
x.tr[1]  = rnorm(1,alpha.tr+beta.tr*x0.tr,tau.tr)
for(t in 2:nn) x.tr[t] = rnorm(1,alpha.tr+beta.tr*x.tr[t-1],tau.tr)
y        = rnorm(nn,x.tr,sig.tr)

par(mfrow=c(1,1))
plot(y,type='l',lwd=2,col='gray80')
lines(x.tr,col=2,lwd=2)

N   = 20000 #number of particles
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

lines(xss[,1],col='gray80')
lines(xss[,3],col='gray80')
lines(xss[,2],col=3,lwd=2)


###################
## parameter estimation
###################

delta = 0.99 #0.95-0.99
b2    = 1-((3*delta-1)/(2*delta))^2
a     = sqrt(1-b2)

th0  = cbind(rnorm(N,0,0.1),rnorm(N,1,1))
x0   = rnorm(N,th0[,1]/(1-bst2b(th0[,2])),sqrt(tau.tr^2/(1-bst2b(th0[,2])^2)))

npar = ncol(th0)
xss  = matrix(NA,ncol=3,nrow=nn)
thss  = array(NA,c(nn,npar,3))

for(t in 1:nn){

    ws  = dnorm(y[t],th0[,1]+bst2b(th0[,2])*x0,sig.tr)
    ind = sample(N,prob = ws,replace = TRUE)

    V    = b2*cov(th0[ind,])*(N-1)/N
    L    = chol(V)

    m    = a*th0[ind,]+(1-a)*apply(th0[ind,],2,mean)[col(th0)]
    rmat = matrix(rnorm(N*ncol(th0)), nrow = N)
    ths  = rmat%*%L+m
    xs   = rnorm(N,ths[,1]+bst2b(ths[,2])*x0[ind],tau.tr)

    ws2  = dnorm(y[t],xs,sig.tr)/ws[ind]
    ind  = sample(N,prob = ws2,replace = TRUE)
    x0   = xs[ind]

    th0  = ths[ind,]
    xss[t,] = c(q025(x0),median(x0),q975(x0))
    thss[t,,1] = apply(th0,2,q025)
    thss[t,,2] = apply(th0,2,median)
    thss[t,,3] = apply(th0,2,q975)
}



par(mfrow=c(3,1))
plot(xss[,2],type='l')
lines(x.tr,col=2)
plot(thss[,1,2],type='l',ylim = c(min(thss[,1,]),max(thss[,1,])))
lines(thss[,1,1],col='gray80')
lines(thss[,1,3],col='gray80')
abline(h=alpha.tr,col=2)
plot(thss[,2,2],type='l',ylim = c(min(thss[,2,]),max(thss[,2,])))
lines(thss[,2,1],col='gray80')
lines(thss[,2,3],col='gray80')
abline(h=b2bst(beta.tr),col=2)

