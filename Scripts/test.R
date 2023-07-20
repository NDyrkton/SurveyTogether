


library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)
library(dclone)

mod.linear.phi <- custommodel('
model{	
#likelihood

for (i in 1:T){
		phi[1,i] <- 1	
}


for (k in 2:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times[k,t])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.01)
positiverate[1]	<- ilogit(logitpositiverate[1])



for(t in 2:T){

	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],rho)T(logitpositiverate[t-1],)
	
	positiverate[t]	<- ilogit(logitpositiverate[t])
}


for (k in 1:K){
	for (t in 1:T){
	  
      	  
		Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
		#Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
	}
}

#priors
theta0 ~ dnorm(-2, 1);
rho ~ dgamma(5,12);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
	gamma1[k] ~ dnorm(0, 1/0.01);
}
}')



mod.linear.phi.2 <- custommodel('
model{	
#likelihood

for (i in 1:T){
		phi[1,i] <- 1	
}


for (k in 2:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times[k,t])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.01)
positiverate[1]	<- ilogit(logitpositiverate[1])



for(t in 2:T){


  #get cumulative positiverate
  
  for(i in 1:(t-1)){
  cposrate_t_1 <- cposrate_t_1 + ilogit(logitpositiverate[i])
  
  }
  #cannot jump higher than cumulative posrate
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],rho)T(,logit(1-cposrate_t_1))
	
	positiverate[t]	<- cposrate_t_1 + ilogit(logitpositiverate[t])
}

for(t in 1:T){
  
  posrate[t] <- dround(positiverate[t],3)
  
	P[t] ~ dbin(positiverate[t], N)
}


for (k in 1:K){
	for (t in 1:T){
	  
	  #Pt[t] <- dround(P[t],)
	  
		Y[k,t] ~ dbin(1-(1-(P[t]/N))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
rho ~ dgamma(5,12);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
	gamma1[k] ~ dnorm(0, 1/0.01);
}
}')


generate.dataset <- function(N= 10000, K =3, t = c(1:3), ns = rep(100,length(t)), phi = "constant"){
  Y <- matrix(NA,ncol = length(t),nrow = K)
  smalln <- t(matrix(rep(ns,K),ncol = K))
  theta_t <- numeric(length(t))
  logit_theta_t <- numeric(length(t))
  times <- t(matrix(rep(t,K),ncol = K))
  
  #priors on general parameters
  rho <- rgamma(1,shape = 3,rate = 1.5)
  theta0 <- rnorm(1,mean =-2, sd = 1)
  
  if(phi == "constant"){
    #generate param based on prior
    gamma0 <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1))) 
    #constant phi values
    phi <- exp(gamma0)
    
    
    logit_theta_t[1] <- rnorm(1,mean = theta0,sd = 1/10)
    theta_t[1] <- inv.logit(logit_theta_t[1])
    
    for(i in 2:length(t)){
      
      logit_theta_t[i] <- rtruncnorm(1,a = logit_theta_t[i-1], b = Inf, mean = logit_theta_t[i-1],sd = 1/rho^2)
      theta_t[i] <- inv.logit(logit_theta_t[i])
      
    }
    P_t <- rbinom(n = length(t),size = N, prob = theta_t)
    
    for(k in 1:K){
      
      for(i in 1:length(t)){
        Y_kt <- rnoncenhypergeom(n = 1, n1 = P_t[i],n2 = N-P_t[i], m1 = smalln[k,i], psi = phi[k])
        Y[k,i] <- Y_kt
        
      }
    }
    
    parameters <-  c(gamma0,rho,theta_t, theta0)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = "")
    theta.names <- paste(rep("theta",length(t)),1:length(t),sep = "")
    names(parameters) <-  c(gamma0.names,"rho",theta.names,"theta0")
    
    return(list(K=K, T=max(t), times=times, N=N, Y=Y, smalln=smalln, params = parameters))
    
  }else if(phi == "linear"){
    phi_kt <- matrix(NA,nrow =K,ncol = length(t))
    #priors
    gamma_0k <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1)))
    gamma_1k <- c(0,rnorm(K-1,mean = 0, sd = rep(0.1,K-1)))
    
    logit_theta_t[1] <- rnorm(1,mean = theta0,sd = 1/10)
    theta_t[1] <- inv.logit(logit_theta_t[1])
    
    for(i in 2:length(t)){
      logit_theta_t[i] <- rtruncnorm(1,a = logit_theta_t[i-1], b = Inf, mean = logit_theta_t[i-1],sd = 1/rho^2)
      theta_t[i] <- inv.logit(logit_theta_t[i])
      
    }
    
    P_t <- rbinom(n = length(t),size = N, prob = theta_t)
    
    for(k in 1:K){
      for(i in 1:length(t)){
        
        phi_kt[k,i] <- exp(gamma_0k[k] + gamma_1k[k]*t[i])
        
        Y_kt <- rnoncenhypergeom(n = 1, n1 = P_t[i],n2 = N-P_t[i], m1 = smalln[k,i], psi = phi_kt[k,i])
        Y[k,i] <- Y_kt
      }
      
    }  
    
    parameters <-  c(gamma_0k,gamma_1k,rho,theta_t, theta0)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = "")
    gamma1.names <- paste(rep("gamma1",K),1:K,sep = "")
    theta.names <- paste(rep("theta",length(t)),1:length(t),sep = "")
    names(parameters) <-  c(gamma0.names,gamma1.names,"rho",theta.names,"theta0")
    
    return(list(K=K, T=max(t), times=times, N=N, Y=Y, smalln=smalln,params = parameters))
    
  }else if(phi == "walk"){
    
    phi_kt <- matrix(NA,nrow =K,ncol = length(t))
    gamma_kt <- matrix(NA,nrow = K,ncol = length(t))
    
    gamma_0k <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1)))
    gamma_kt[1,] <- rep(0,length(t))
    #prior
    
    pi <- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = 1/10)
    
    
    
    logit_theta_t[1] <- rnorm(1,mean = theta0,sd = 1/10)
    theta_t[1] <- inv.logit(logit_theta_t[1])
    #first study is unbiased 
    
    phi_kt[1,] <- exp(gamma_kt[1,])
    gamma_kt[2:K,1] <- rnorm(K-1,mean = gamma_0k[2:K],sd = 1/10)
    phi_kt[2:K,1] <- exp(gamma_kt[2:K,1])
    
    for(i in 2:length(t)){
      
      gamma_kt[2:K,i] <- rnorm(K-1, mean = gamma_kt[2:K,i-1],sd = c(pi,pi))
      phi_kt[2:K,i] <- exp(gamma_kt[2:K,i])
      
      
      logit_theta_t[i] <- rtruncnorm(1,a = logit_theta_t[i-1], b = Inf, mean = logit_theta_t[i-1],sd = 1/rho^2)
      theta_t[i] <- inv.logit(logit_theta_t[i])
      
    }
    
    P_t <- rbinom(n = length(t),size = N, prob = theta_t)
    
    for(k in 1:K){
      
      for(i in 1:length(t)){
        
        
        Y_kt <- rnoncenhypergeom(n = 1, n1 = P_t[i],n2 = N-P_t[i], m1 = smalln[k,i], psi = phi_kt[k,i])
        Y[k,i] <- Y_kt
      }
    }
    parameters <- c(gamma_0k,as.numeric(gamma_kt),rho,pi,theta_t,theta0)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = "")
    gamma_kt.c <- expand.grid(1:K,1:length(t))
    gamma_kt.names <- paste(rep("gamma_",K*length(t)),as.character(gamma_kt.c$Var1),as.character(gamma_kt.c$Var2),sep = '')
    theta.names <- paste(rep("theta",length(t)),1:length(t),sep = "")
    names(parameters) <-  c(gamma0.names,gamma_kt.names,"rho","pi",theta.names,"theta0")
    
    
    return(list(K=K, T=max(t), times=times, N=N, Y=Y, smalln=smalln, params = parameters))
  }
  
}

inv.logit <- function(x){exp(x)/(1+exp(x))}

#N = 2500000, n = rep(25000,20)


my.dat <- generate.dataset(N = 255200373, n = rep(250,20),t = 1:20,phi = "linear")

y <- my.dat$params[grep("theta",names(my.dat$params))][-21]

plot(1:20,y)


my.dat2 <- my.dat

my.dat2$Y[1,2:5] <- NA

my.dat2$Y[2,c(3:6,15:18)] <- NA

my.dat2$Y[3,c(1:3,9:11)] <- NA


cl <- makePSOCKcluster(4)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")




line.linear <- jags.parfit(cl, my.dat, c("positiverate","rho"), mod.linear.phi,
                           n.chains=4,n.adapt = 100000,thin = 10, n.iter = 100000)

means.posrate <- summary(line.linear)$statistics[,1]
#check


gelman.diag(line.linear)
plot(1:20,y)
lines(means.posrate[-21])
lines(my.dat$Y[2,]/my.dat$smalln[2,], col = "blue")
lines(my.dat$Y[3,]/my.dat$smalln[3,], col = "red")





#test noncentral hypergeometric 
#binomial approximation

N <- 25520000
P <- N/3

phi <-  1.7

n <- 10000

param <- rnoncenhypergeom(10000,P,N-P,n,phi)
est <- rbinom(10000,n,1-(1-(P/N))^phi)

quantile(param)
quantile(est)


hist(param)
hist(est)

mean(param)
mean(est)

var(param)
var(est)




dnoncenhypergeom(min(param),P,N-P,n,phi)
dbinom(min(param),n,1-(1-(P/N))^phi)



#test small sample approximation

pos.rate <- 0.8732

P_t <- rbinom(1,size =25520000 ,prob = pos.rate)


surveys <- rbinom(1,size = 100,prob = 1-(1-(P_t/N))^1)

mean(surveys/100)
jags.mod <- custommodel("
                        
model{


  P ~ dbin(positiverate,N)


  Y ~ dbin(1-(1-(P/N)),smalln)
  
  #prior 
  
  positiverate ~ dunif(0.8,1)


}
                        
")


cl <- makePSOCKcluster(4)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


line.test <- jags.parfit(cl, list(Y=surveys,smalln = 100,N = N), c("positiverate"), jags.mod,
                           n.chains=4,n.adapt = 5000,thin = 5, n.iter =100000)


