
#Survey Together Simulation Study
#3x3 Design         Nathaniel Dyrkton  Supervised by Paul Gustafson and Harlan Campbell
library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)
library(dclone)

#First step is to generate data under 3 conditions for phi (bias term)
#1: Phi is constant by  time, 2: Phi is linear in time, and 3: Phi follows a random walk 

inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/(1-x))
}



theta0 <- -2

sigmasq = rtruncnorm(1000,mean = 0.25, sd = sqrt(0.25),a = 0, b= Inf)

hist(sigmasq)

inv.logit(theta0)

quantile(inv.logit(rnorm(50000,mean = 0,sd = sqrt(0.511))),seq(0,1,0.1))-inv.logit(0)


result <- numeric(10000)

sigmasq = rtruncnorm(10000,mean = 0.5, sd = sqrt(0.25),a = 0, b= Inf)

for(i in 1:10000){
  
  result[i] <- inv.logit(rnorm(1,mean = 0,sd = sqrt(sigmasq[i])))-inv.logit(0) 
  
}

quantile(result,seq(0,1,0.05))

         