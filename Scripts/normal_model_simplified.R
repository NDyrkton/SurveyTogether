#Survey Together Simulation Study
#3x3 Design  
#Simplified only using Normal Model
#Nathaniel Dyrkton  Supervised by Paul Gustafson and Harlan Campbell

#known sigma, unkonwn mu, unkown phi.

library(rjags)
library(MCMCpack)
library(ggplot2)
library(truncnorm)

gen.norm.dat <- function(N= 10000, K =3, t = c(1:5), ns = rep(100,length(t)), phi = "constant"){
  sigma <- 5
  mu0 <- rnorm(1,mean = 75, sd = 5)
  len <- length(t)
  mu <- numeric(len)
  rho <-  rinvgamma(1,shape = 2,scale = 2)
  Y <- matrix(NA,nrow = K,ncol = len)
  times <- t(matrix(rep(t,K),ncol = K))
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
    mu.names <- paste(rep("mu",len),t,sep = '')
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
    mu.names <- paste(rep("mu",len),t,sep = '')
    names(parameters) <- c("mu0",gamma0.names,gamma1.names,"rho",mu.names)
    
    return(list(K=K, Ivec= Ivec,T=5, times=times, N=N, Y=Y, smalln=smalln, params = parameters))
    
  }else if(phi = "walk"){
    
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
    parameters <- c(mu0,phi0,phi,rho,mu)
    phi0.names <- paste(rep("phi0",K),1:K,sep = '')
    phi_kt.names <- paste(rep("phi",K),1:K,sep = '')
    mu.names <- paste(rep("mu",len),t,sep = '')
    names(parameters) <- c("mu0",gamma0.names,gamma1.names,"rho",mu.names)
    
    return(list(K=K, Ivec= Ivec,T=5, times=times, N=N, Y=Y, smalln=smalln, params = parameters))

  }else{
    
    
    print("phi type not recognized")
    return(NULL)
  }
  
  
  
  
}