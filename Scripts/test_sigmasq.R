

#prior simulations;
library(MCMCpack)
library(rjags)
library(truncnorm)

inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/(1-x))
}


#base sigmasq
#How much does sigmasq affect the positiverate?
sigmasq <-  rtruncnorm(1000000,mean = 0, sd = sqrt(1),a = 0, b= Inf)
jump.method <- inv.logit(rnorm(1000000,0,sd = sqrt(sigmasq)))
quantile(jump.method-0.5,seq(0,1,0.05))
####


#theta0 prior
theta0 <- rnorm(1000000,mean = 0,sd = sqrt(2))
quantile(inv.logit(theta0),seq(0,1,0.05))
inv.logit(quantile(theta0,seq(0,1,0.05)))

#check histogram
hist(inv.logit(theta0))


##########SECTION 3###################
#####################################
#theta0 prior
theta0 <- rnorm(1000000,mean = 0,sd = sqrt(1))
quantile(inv.logit(theta0),seq(0,1,0.05))
inv.logit(quantile(theta0,seq(0,1,0.05)))

#section 3 sigmasq
#How much does sigmasq affect the positiverate?
sigmasq <-  rtruncnorm(100000,mean = 0, sd = sqrt(0.1),a = 0, b= Inf)
jump.method <- inv.logit(rnorm(100000,0,sd = sqrt(sigmasq)))
quantile(jump.method-0.5,seq(0,1,0.05))
####


#test for linear model for phi -> section 3
gamma0 <- 1
gamma1 <- rnorm(1000000,mean = 0, sd = sqrt(0.01))

phi0 <- exp(gamma0)
phi1 <- exp(gamma0+gamma1)

quantile(phi1-phi0,seq(0,1,0.05))


#test for random walk model for phi -> section 3
pi <- rtruncnorm(100000,mean = 0,sd = sqrt(0.01),a = 0, b = Inf)

phi2 <- exp(rnorm(100000,mean = 0,sd = sqrt(pi)))
quantile(phi2-exp(0),seq(0,1,0.05))






