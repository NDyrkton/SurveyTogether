
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



#prior on sigmasq
sigmasq <-  rtruncnorm(10000,mean = 0, sd = sqrt(1),a = 0, b= Inf)

jump1 <- inv.logit(rtruncnorm(10000,0,sd = sqrt(sigmasq),a = 0, b = Inf))


quantile(jump1-0.5,seq(0,1,0.05))



x <- inv.logit(rnorm(10000,mean = 0,sd = sqrt(sigmasq)))
quantile(x-inv.logit(0),seq(0,1,0.05))         


sigmasq = rtruncnorm(100000,mean = 0, sd = sqrt(),a = 0, b= Inf)


pi <- rtruncnorm(10000,mean = 0,sd = sqrt(0.1),a = 0, b = Inf)

x2 <- exp(rnorm(100000,mean = 0,sd = sqrt(pi)))
quantile(x2-exp(0),seq(0,1,0.05))

x3 <- inv.logit(rnorm(10000,mean = 0, sd = 0.5))

quantile(x3,seq(0,1,0.05))
