#Survey Together Simulation Study
#3x3 Design  
#Simplified only using Normal Model
#Nathaniel Dyrkton  Supervised by Paul Gustafson and Harlan Campbell

#known sigma, unknown mu, unknown phi.

library(rjags)
library(MCMCpack)
library(ggplot2)
library(truncnorm)

#models

mod.const.phi<- '
model{	
#likelihood

for (i in 1:Ivec[1]){
		phi[1,i] <- 1	}

for (k in 2:K){
	for (i in 1:Ivec[k]){
		phi[k,i] <- gamma0[k]
	}
}
	
	
mu[1] ~ dnorm(mu0,1/10)
for(t in 2:T){
	mu[t]	~ dnorm(mu[t-1],rho)
}

for (k in 1:K){
	for (i in 1:Ivec[k]){
		
		Y[k,i] ~ dnorm(mu[times[k,i]] + phi[k,i], 10/sqrt(smalln[k,i]));
	}
}

#priors
mu0 ~ dnorm(75, 5);
rho ~ dnorm(0, 1/5)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 2);
}
}'

mod.linear.phi <- '
model{	
#likelihood

for (i in 1:Ivec[1]){
		phi[1,i] <- 1	}

for (k in 2:K){
	for (i in 1:Ivec[k]){
		phi[k,i] <- exp(gamma0[k] + gamma1[k]*times[k,i])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/10)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],rho)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

for(t in 1:T){
	P[t] ~ dbin(positiverate[t], N)
}

for (k in 1:K){
	for (i in 1:Ivec[k]){
		
		Y[k,i] ~ dhyper(P[times[k,i]], N-P[times[k,i]], smalln[k,i], phi[k,i]);
	}
}

#priors
theta0 ~ dnorm(0, 2);
rho ~ dnorm(0, 1/5)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 2);
	gamma1[k] ~ dnorm(0, 1/5);
}
}'

mod.walk.phi <- '
model{	
#likelihood

for (i in 1:Ivec[1]){
		phi[1,i] <- 1	
		gamma[1,i] <- 0
}


for (k in 2:K){

  gamma[k,1] ~ dnorm(gamma0[k],1/10)
  phi[k,1] <- exp(gamma[k,1])
  
	for (t in 2:T){
	  gamma[k,t] ~ dnorm(gamma[k,t-1], pi)
	  phi[k,t] <- exp(gamma[k,t])
	  
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/10)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],rho)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

for(t in 1:T){
	P[t] ~ dbin(positiverate[t], N)
}

for (k in 1:K){
	for (i in 1:Ivec[k]){
		
		Y[k,i] ~ dhyper(P[times[k,i]], N-P[times[k,i]], smalln[k,i], phi[k,i]);
	}
}

#priors
theta0 ~ dnorm(0, 2);
rho ~ dnorm(0, 1/5)T(0,);
pi ~ dnorm(0, 1/5)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 2);

}

}'








gen.norm.dat <- function(N= 10000, K =3, ts = c(1:5), ns = rep(100,length(ts)), phi = "constant"){
  sigma <- 5
  mu0 <- rnorm(1,mean = 75, sd = 5)
  len <- length(ts)
  mu <- numeric(len)
  rho <-  rinvgamma(1,shape = 2,scale = 2)
  Y <- matrix(NA,nrow = K,ncol = len)
  times <- t(matrix(rep(ts,K),ncol = K))
  smalln <- t(matrix(rep(ns,K),ncol = K))
  Ivec <- rep(3,K)
  
  
  if(phi == "constant"){
    phi <- c(0,rnorm(K-1,mean = 5,sd =1))
    
    #random walk for mu
    mu[1] <- rnorm(1,mean = mu0,sd = 5)
    
    for(t in 2:len){
      mu[t] <- rnorm(1,mean = mu[t-1], sd= rho)
      
    }
    # the biased surveys
    for(k in 1:K){
      Y[k,] <- rnorm(len,mean =mu + phi[k], sd = rep(10,len)/sqrt(ns)   )
    }
    
    parameters <- c(mu0,phi,rho,mu)
    phi.names <- paste(rep("phi",K),1:K,sep = '')
    mu.names <- paste(rep("mu",len),ts,sep = '')
    names(parameters) <- c("mu0",phi.names,"rho",mu.names)
    
    return(list(K=K, Ivec= Ivec,T=5, times=times, N=N, Y=Y, smalln=smalln, params = parameters))
    
  }else if(phi == "linear"){
    
    gamma0 <- c(0,rnorm(K-1,mean = 5,sd = c(5,5)))
    gamma1 <- c(0,rnorm(K-1,mean = 2,sd = c(1,1)))
    
    
    #random walk for mu
    mu[1] <- rnorm(1,mean = mu0,sd = 5)
    
    for(t in 2:len){
      mu[t] <- rnorm(1,mean = mu[t-1], sd= rho)
      
    }
    # the biased surveys
    for(k in 1:K){
      #vectorized      
      phi <- gamma0[k] + gamma1[k]*t

      Y[k,] <- rnorm(len,mean =mu + phi, sd = rep(10,len)/sqrt(ns)   )
    }
    
    parameters <- c(mu0,gamma0,gamma1,rho,mu)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = '')
    gamma1.names <- paste(rep("gamma1",K),1:K,sep = '')
    mu.names <- paste(rep("mu",len),ts,sep = '')
    names(parameters) <- c("mu0",gamma0.names,gamma1.names,"rho",mu.names)
    
    return(list(K=K, Ivec= Ivec,T=5, times=times, N=N, Y=Y, smalln=smalln, params = parameters))
    
  }else if(phi == "walk"){
    
    phi0 <- c(0, rnorm(K-1,mean = 5, sd = 3))
    phi_kt <- matrix(NA,nrow = K,ncol = len)
    phi_kt[1,] <- 0
    pi <- rinvgamma(1,shape = 1,scale = 1)
    
    
    #random walk for mu
    mu[1] <- rnorm(1,mean = mu0,sd = 5)
    
    #random walk for phi
    phi_kt[2:K,1] <- rnorm(K-1,mean = phi0[2:K],sd = c(2,2))
    

    for(t in 2:len){
      mu[t] <- rnorm(1,mean = mu[t-1], sd= rho)
      phi_kt[2:K,t] <- rnorm(K-1,mean = phi_kt[2:K,t-1],sd = c(pi,pi))
      
    }
    # the biased surveys
    for(k in 1:K){
      #vectorized
      Y[k,] <- rnorm(len,mean =mu + phi_kt[k,], sd = rep(10,len)/sqrt(ns)   )
    }
    ###
    parameters <- c(mu0,phi0,phi_kt,rho,mu)
    phi0.names <- paste(rep("phi0",K),1:K,sep = '')
    phi_kt.c <- expand.grid(1:K,1:len)
    phi_kt.names <- paste(rep("phi",K*len),as.character(phi_kt.c$Var1),as.character(phi_kt.c$Var2),sep = '')
    mu.names <- paste(rep("mu",len),ts,sep = '')
    
    names(parameters) <- c("mu0",phi0.names,phi_kt.names,"rho",mu.names)
    
    return(list(K=K, Ivec= Ivec,T=5, times=times, N=N, Y=Y, smalln=smalln, params = parameters))

  }else{
  
    print("phi type not recognized")
    return(NULL)
  }
}


generate.model.ests <- function(model.string, data.list, params ,n.chains =3, n.iter = 20000, thin = 10){
  
  jags.mod <- jags.model(textConnection(model.string), 
                         data = data.list, n.chains = n.chains, n.adapt = 10000,quiet = T)
  
  samps <- coda.samples(jags.mod, params, n.iter = n.iter, thin = 10, progress.bar = "none")
  #as of right now the estimate is a mean.
  bayes.est <- summary(samps)$statistics[,1]
  return(bayes.est)
}

extract.unbiased <- function(datalist){
  K = 1
  Ivec = 3
  T <- datalist$T
  new.list <- list(K = K, Ivec = Ivec, 
                   times = matrix(datalist$times[1,],ncol = T), N = datalist$N, T = T,
                   Y = matrix(datalist$Y[1,],ncol = T), smalln = matrix(datalist$smalln[1,],ncol = T))
  
  return(new.list)
}


#check
gen.norm.dat(phi = "constant")







