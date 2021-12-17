rm(list=ls(all=TRUE))
library(moments)
#set.seed(1)
q025  = function(x) quantile(x,0.025)
q975  = function(x) quantile(x,0.975)
bst2b = function(x) 1/(1+exp(-x))
b2bst = function(x) log(x)-log(1-x)

data <- read.csv("data/IBM_weekly.csv")
Date = as.Date(data$Date[-1])
tmp  = (log(data$Adj.Close[-1])-log(data$Adj.Close[-length(data$Adj.Close)]))*100
lret = tmp-mean(tmp)

par(mfrow=c(1,1))
plot(Date,lret,type='l')
