rm(list=ls(all=TRUE))
#set.seed(1)
q025  = function(x) quantile(x,0.025)
q975  = function(x) quantile(x,0.975)
bst2b = function(x) 1/(1+exp(-x))
b2bst = function(x) log(x)-log(1-x)

