rm(list=ls(all=TRUE))
set.seed(1)
q025 = function(x) quantile(x,0.025)
q975 = function(x) quantile(x,0.975)

nn       = 1000
alpha.tr = 0.01
beta.tr  = 0.5
tau.tr   = 0.5
sig.tr   = 0.5
x.tr     = rep(NA,nn)

x0       = rnorm(1,alpha.tr/(1-beta.tr),sqrt(tau.tr^2/(1-beta.tr^2)))
x.tr[1]  = rnorm(1,alpha.tr+beta.tr*x0,tau.tr)
for(t in 2:nn) x.tr[t] = rnorm(1,alpha.tr+beta.tr*x.tr[t-1],tau.tr)
y        = rnorm(nn,x.tr,sig.tr)

par(mfrow=c(1,1))
plot(y,type='l',lwd=2,col='gray80')
lines(x.tr,col=2,lwd=2)

N   = 5000 #number of particles
xss1 = matrix(NA,ncol=4,nrow=nn)

x0 = rnorm(N,alpha.tr/(1-beta.tr),sqrt(tau.tr^2/(1-beta.tr^2)))


# just resample-propagae
for(t in 1:nn){
    ws  = dnorm(y[t],x0,sig.tr)
    ind = sample(N,prob = ws,replace = TRUE)
    xs  = rnorm(N,alpha.tr+beta.tr*x0[ind],tau.tr)
    # ws2 = dnorm(y[t],xs,sig.tr)/ws
    # ind = sample(N,prob = ws2,replace = TRUE)
    # x0  = xs[ind]
    x0=xs
    xss1[t,] = c(q025(x0),mean(x0),median(x0),q975(x0))
}



# resample-propagate-resample
xss2 = matrix(NA,ncol=4,nrow=nn)
x0 = rnorm(N,alpha.tr/(1-beta.tr),sqrt(tau.tr^2/(1-beta.tr^2)))
for(t in 1:nn){
    ws  = dnorm(y[t],x0,sig.tr)
    ind = sample(N,prob = ws,replace = TRUE)
    xs  = rnorm(N,alpha.tr+beta.tr*x0[ind],tau.tr)
    ws2 = dnorm(y[t],xs,sig.tr)/ws
    ind = sample(N,prob = ws2,replace = TRUE)
    x0  = xs[ind]
    xss2[t,] = c(q025(x0),mean(x0),median(x0),q975(x0))
}


lines(xss1[,2],col=4)
lines(xss2[,2],col=3)

c(sum((xss1[,2]-x.tr)^2),sum((xss1[,3]-x.tr)^2))/
c(sum((xss2[,2]-x.tr)^2),sum((xss2[,3]-x.tr)^2))


mean(xss1[,4]-xss1[,1])
mean(xss2[,4]-xss2[,1])



