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

lines(xss[,1],col='gray80')
lines(xss[,3],col='gray80')
lines(xss[,2],col=3,lwd=2)


###################
## parameter estimation
###################
runs  = 10
res   = list()

for(j in 1:runs)
{
    delta = runif(1,0.95,0.99) #0.95-0.99
    b2    = 1-((3*delta-1)/(2*delta))^2
    a     = sqrt(1-b2)

    th0   = cbind(rnorm(N,0,3),rnorm(N,0,3),rnorm(N,0,0.5))
    x0    = rnorm(N,th0[,1]/(1-bst2b(th0[,2])),sqrt(exp(th0[,3])^2/(1-bst2b(th0[,2])^2)))

    npar  = ncol(th0)
    xss   = matrix(NA,ncol=3,nrow=nn)
    thss  = array(NA,c(nn,npar,3))

    for(t in 1:nn){

        ws  = dnorm(y[t],th0[,1]+bst2b(th0[,2])*x0,sig.tr)
        ind = sample(N,prob = ws,replace = TRUE)

        L    = chol(b2*cov(th0)*(N-1)/N)
        m    = a*th0[ind,]+(1-a)*apply(th0[ind,],2,mean)[col(th0)]
        rmat = matrix(rnorm(N*ncol(th0)), nrow = N)
        ths  = rmat%*%L+m
        xs   = rnorm(N,ths[,1]+bst2b(ths[,2])*x0[ind],exp(ths[,3]))

        ws2  = dnorm(y[t],xs,sig.tr)/ws[ind]
        ind  = sample(N,prob = ws2,replace = TRUE)
        x0   = xs[ind]

        th0  = ths[ind,]
        xss[t,] = c(q025(x0),median(x0),q975(x0))
        thss[t,,1] = apply(th0,2,q025)
        thss[t,,2] = apply(th0,2,median)
        thss[t,,3] = apply(th0,2,q975)
    }

    res[[j]] = list(xs=xss,th=thss)
}


par(mfrow=c(2,2))
plot(x.tr,type='l',col=2,lwd=3)
for(j in 1:runs){
    lines(res[[j]]$xs[,2])
    lines(res[[j]]$xs[,1],col='gray80')
    lines(res[[j]]$xs[,3],col='gray80')
}

plot(res[[1]]$th[,1,2],type='l',ylim=c(-0.1,0.1))
for(j in 1:runs){
    lines(res[[j]]$th[,1,2])
    lines(res[[j]]$th[,1,1],col='gray80')
    lines(res[[j]]$th[,1,3],col='gray80')
}
abline(h=alpha.tr,col=2,lwd=3)


plot(res[[1]]$th[,2,2],type='l',ylim=c(0,3.5))
for(j in 1:runs){
    lines(res[[j]]$th[,2,2])
    lines(res[[j]]$th[,2,1],col='gray80')
    lines(res[[j]]$th[,2,3],col='gray80')
}
abline(h=b2bst(beta.tr),col=2,lwd=3)


plot(res[[1]]$th[,3,2],type='l',ylim=c(-1,0.5))
for(j in 1:runs){
    lines(res[[j]]$th[,3,2])
    lines(res[[j]]$th[,3,1],col='gray80')
    lines(res[[j]]$th[,3,3],col='gray80')
}
abline(h=log(tau.tr),col=2,lwd=3)


