library(MCMCpack)

set.seed(125)
posrate = 0.8
phi = 0.8
N <- 1000000
n = 1000
P <- rbinom(n = 1,size = N, prob = posrate)
Y <- rnoncenhypergeom(n = 1, n1 = P, N-P, m1 = n, psi = phi)


dbinom(Y,size = n, prob = (posrate*phi)/(1-posrate + (posrate*phi)))
dbinom(Y,size = n, prob =(1-(1-posrate)^phi))



dnoncenhypergeom(Y, n1 = P, N-P, m1 = n, psi = phi)
##ok looks better


#Survey Together Simulation Study
#3x3 Design         Nathaniel Dyrkton  Supervised by Paul Gustafson and Harlan Campbell
library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)
library(dclone)

#T = 10, 10 timpoints

#First step is to generate data under 3 conditions for phi (bias term)
#1: Phi is constant by  time, 2: Phi is linear in time, and 3: Phi follows a random walk 


inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/(1-x))
}

#now fixed to be consistent with notation in paper.


#How fixed to be consistent with notation in paper.
generate.dataset <- function(N= 10000000, K =3, t = c(1:5), phi = "constant"){
  Y <- matrix(NA,ncol = length(t),nrow = K)
  smalln <- matrix(0,ncol = length(t),nrow = K)
  #get default smalln
  
  for(k in 1:K){
    #k = 1 is n = 100
    smalln[k,] <- rep(100,length(t))
    
    if(k>1){
      #all other k is k = 1000
      smalln[k,] <- rep(1000,length(t))
    }
  }
  
  #initialize positive rate
  posrate_t <- numeric(length(t))
  theta_t <- numeric(length(t))
  times <- t(matrix(rep(t,K),ncol = K))
  
  #priors on general parameters
  sigmasq<- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = sqrt(0.1))
  theta0 <- rnorm(1,mean =0, sd = sqrt(1))
  
  if(phi == "constant"){
    #generate param based on prior
    gamma0 <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1))) 
    #constant phi values
    phi <- exp(gamma0)
    
    
    theta_t[1] <- rnorm(1,mean = theta0,sd = sqrt(sigmasq))
    posrate_t[1] <- inv.logit(theta_t[1])
    
    #random walk for positive rate
    for(i in 2:length(t)){
      
      theta_t[i] <- rnorm(1,mean = theta_t[i-1],sd = sqrt(sigmasq))
      posrate_t[i] <- inv.logit(theta_t[i])
      
    }
    
    #generate positive in population
    P_t <- rbinom(n = length(t),size = N, prob = posrate_t)
    
    for(k in 1:K){
      
      for(i in 1:length(t)){
        
        #generate survey sample
        
        Y_kt <- rnoncenhypergeom(n = 1, n1 = P_t[i],n2 = N-P_t[i], m1 = smalln[k,i], psi = phi[k])
        Y[k,i] <- Y_kt
        
      }
    }
    
    parameters <-  c(gamma0,sigmasq,posrate_t, theta0)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = "")
    posrate.names <- paste(rep("posrate",length(t)),1:length(t),sep = "")
    names(parameters) <-  c(gamma0.names,"sigmasq",posrate.names,"theta0")
    
    return(list(K=K, T=max(t), times=times, N=N, Y=Y, smalln=smalln, params = parameters))
    
  }else if(phi == "linear"){
    phi_kt <- matrix(NA,nrow =K,ncol = length(t))
    #priors
    gamma_0k <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1)))
    gamma_1k <- c(0,rnorm(K-1,mean = 0, sd = rep(sqrt(0.01),K-1)))
    
    theta_t[1] <- rnorm(1,mean = theta0,sd = sqrt(sigmasq))
    posrate_t[1] <- inv.logit(theta_t[1])
    
    for(i in 2:length(t)){
      theta_t[i] <- rnorm(1,mean = theta_t[i-1],sd = sqrt(sigmasq))
      posrate_t[i] <- inv.logit(theta_t[i])
      
    }
    
    P_t <- rbinom(n = length(t),size = N, prob = posrate_t)
    
    for(k in 1:K){
      for(i in 1:length(t)){
        
        phi_kt[k,i] <- exp(gamma_0k[k] + gamma_1k[k]*t[i])
        
        Y_kt <- rnoncenhypergeom(n = 1, n1 = P_t[i],n2 = N-P_t[i], m1 = smalln[k,i], psi = phi_kt[k,i])
        Y[k,i] <- Y_kt
      }
      
    }  
    
    parameters <-  c(gamma_0k,gamma_1k,sigmasq,posrate_t, theta0)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = "")
    gamma1.names <- paste(rep("gamma1",K),1:K,sep = "")
    posrate.names <- paste(rep("posrate",length(t)),1:length(t),sep = "")
    names(parameters) <-  c(gamma0.names,gamma1.names,"sigmasq",posrate.names,"theta0")
    
    return(list(K=K, T=max(t), times=times, N=N, Y=Y, smalln=smalln,params = parameters))
    
  }else if(phi == "walk"){
    
    phi_kt <- matrix(NA,nrow =K,ncol = length(t))
    gamma_kt <- matrix(NA,nrow = K,ncol = length(t))
    
    gamma_0k <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1)))
    gamma_kt[1,] <- rep(0,length(t))
    #prior
    
    pisq <- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = sqrt(1/100))
    
    
    
    theta_t[1] <- rnorm(1,mean = theta0,sd = sqrt(sigmasq))
    posrate_t[1] <- inv.logit(theta_t[1])
    #first study is unbiased 
    
    phi_kt[1,] <- exp(gamma_kt[1,])
    gamma_kt[2:K,1] <- rnorm(K-1,mean = gamma_0k[2:K],sd = sqrt(pisq))
    phi_kt[2:K,1] <- exp(gamma_kt[2:K,1])
    
    for(i in 2:length(t)){
      
      gamma_kt[2:K,i] <- rnorm(K-1, mean = gamma_kt[2:K,i-1],sd = sqrt(c(pisq,pisq)))
      phi_kt[2:K,i] <- exp(gamma_kt[2:K,i])
      
      
      theta_t[i] <- rnorm(1,mean = theta_t[i-1],sd = sqrt(sigmasq))
      posrate_t[i] <- inv.logit(theta_t[i])
      
    }
    
    P_t <- rbinom(n = length(t),size = N, prob = posrate_t)
    
    for(k in 1:K){
      
      for(i in 1:length(t)){
        
        Y_kt <- rnoncenhypergeom(n = 1, n1 = P_t[i],n2 = N-P_t[i], m1 = smalln[k,i], psi = phi_kt[k,i])
        Y[k,i] <- Y_kt
      }
    }
    parameters <- c(gamma_0k,as.numeric(gamma_kt),sigmasq,pisq,posrate_t,theta0)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = "")
    gamma_kt.c <- expand.grid(1:K,1:length(t))
    gamma_kt.names <- paste(rep("gamma_",K*length(t)),as.character(gamma_kt.c$Var1),as.character(gamma_kt.c$Var2),sep = '')
    posrate.names <- paste(rep("posrate",length(t)),1:length(t),sep = "")
    names(parameters) <-  c(gamma0.names,gamma_kt.names,"sigmasq","pisq",posrate.names,"theta0")
    
    
    return(list(K=K, T=max(t), times=times, N=N, Y=Y, smalln=smalln, params = parameters))
  }
  
}

mod.const.phi<- custommodel('
model{	
#likelihood

for (t in 1:T){
		phi[1,t] <- 1	}

for (k in 2:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1], 1/sigmasq)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}


for (k in 1:K){
	for (t in 1:T){
		
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1)
sigmasq ~ dnorm(0, 1/0.1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
}
}')


mod.const.phi.approx<- custommodel('
model{	
#likelihood

for (t in 1:T){
		phi[1,t] <- 1	}

for (k in 2:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1], 1/sigmasq)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}


for (k in 1:K){
	for (t in 1:T){
		
		Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1)
sigmasq ~ dnorm(0, 1/0.1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
}
}')

mod.const.phi.exact<- custommodel('
model{	
#likelihood

for (t in 1:T){
		phi[1,t] <- 1	}

for (k in 2:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1], 1/sigmasq)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}


for(t in 1:T){
	P[t] ~ dbin(positiverate[t], N)
}


for (k in 1:K){
	for (t in 1:T){
		
		Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
	}
}

#priors
theta0 ~ dnorm(0, 1)
sigmasq ~ dnorm(0, 1/0.1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
}
}')

#this function extracts the first survey only as "unbiased survey"
extract.unbiased <- function(datalist){
  K = 1
  T <- datalist$T
  new.list <- list(K = K, 
                   times = matrix(datalist$times[1,],ncol = T), N = datalist$N, T = T,
                   Y = matrix(datalist$Y[1,],ncol = T), smalln = matrix(datalist$smalln[1,],ncol = T))
  
  return(new.list)
}

#calculates the posterior mean at given timepoint (not used)
get.mean <- function(mcmc.obj,timepoint){
  return(summary(mcmc.obj)$statistics[,1][timepoint])
}

#calculates the posterior variance at given timepoint (not used)
get.var <- function(mcmc.obj,timepoint){
  return(summary(mcmc.obj)$statistics[,2][timepoint])
}

#gets 95% credible interval (not used)
get.CI <- function(line,var){
  
  quantiles <- summary(line)$quantile
  
  if(is.null(dim(quantiles))){
    
    lower.quantile <- quantiles[1]
    upper.quantile <- quantiles[5]
    
    return(list(Lower = lower.quantile,Upper = upper.quantile))
    
  }else{
    
    lower.quantile <- summary(line)$quantile[,1]
    upper.quantile <- summary(line)$quantile[,5]
    
    #extract variable of interest
    lower.quantile <- lower.quantile[grep(var,names(lower.quantile))] 
    upper.quantile <- upper.quantile[grep(var,names(upper.quantile))] 
    
    return(c(lower.quantile,upper.quantile))
    
  }
  
  
}

dcoptions("verbose"=F)#mute the output of dclone



#generate all datasets in advance
generate.data.replicates <- function(phi = "constant",t = 1:5,NN,K = 3){
  
  #function creates all data for simulation
  list.return <- list()
  for(i in 1:NN){
    gen.data <- generate.dataset(phi = phi,t = t)
    
    if(sum(gen.data$Y[2:K,]==1000 | gen.data$Y[2:K,]==0)>=1){
      #print if phi is unidentifiable
      print(paste("bad data on i =", i))
      
    }
    
    list.return[[i]] <- gen.data
  }
  
  return(list.return)
}

#check how many unidentifiable phi (not used)
check.bad.data <- function(data.list,t = 1:5,phi = 'constant',K = 3){
  
  for(i in 1:length(data.list)){
    
    if(sum(data.list[[i]]$Y[2:K,]==1000 | data.list[[i]]$Y[2:K,]==0)>=1){
      print(paste("bad data on i =", i))
      
      data.list[[i]] <- generate.dataset(phi = phi,t = t,K = K)
      
    }
  }
  
  return(data.list)
  
}

#how many replications
NN <- 100

set.seed(12345)
#generate the NN datasets
data.const <- generate.data.replicates(phi = "constant", NN = NN,t = 1:10)
data.linear  <- generate.data.replicates(phi = "linear", NN = NN,t = 1:10)
data.walk  <- generate.data.replicates(phi = "walk", NN = NN,t = 1:10)


#parallel 10 cores
cl <- makePSOCKcluster(3)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


chain1<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(1))
chain2<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(1+1))
chain3<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(1+2))

chains.init <- list(chain1,chain2,chain3)

data.const.phi <- data.const[[1]]

const.binom_limit <- jags.parfit(cl, data.const.phi, c("positiverate","gamma0"), mod.const.phi,
                           n.chains = 3,n.adapt = 25000,thin = 5, n.iter = 70000,inits = chains.init)

const.binom_approx <- jags.parfit(cl, data.const.phi, c("positiverate","gamma0"), mod.const.phi.approx,
                                  n.chains = 3,n.adapt = 25000,thin = 5, n.iter = 70000,inits = chains.init)

const.binom_exact <- jags.parfit(cl, data.const.phi, c("positiverate","gamma0"), mod.const.phi.exact,
                                 n.chains = 3,n.adapt = 25000,thin = 5, n.iter = 70000,inits = chains.init)


summary(const.binom_limit)
summary(const.binom_approx)
