
library(MCMCpack)
library(truncnorm)

inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/(1-x))
}

#now fixed to be consistent with notation in paper.
generate.dataset <- function(N= 10000, K =3, t = c(1:5), ns = rep(1000,length(t)), phi = "constant"){
  Y <- matrix(NA,ncol = length(t),nrow = K)
  smalln <- t(matrix(rep(ns,K),ncol = K))
  posrate_t <- numeric(length(t))
  theta_t <- numeric(length(t))
  times <- t(matrix(rep(t,K),ncol = K))
  
  #priors on general parameters
  sigmasq<- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = sqrt(0.5))
  theta0 <- rnorm(1,mean =0, sd = sqrt(0.5))
  
  if(phi == "constant"){
    #generate param based on prior
    gamma0 <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1))) 
    #constant phi values
    phi <- exp(gamma0)
    
    
    theta_t[1] <- rnorm(1,mean = theta0,sd = sqrt(sigmasq))
    posrate_t[1] <- inv.logit(theta_t[1])
    
    for(i in 2:length(t)){
      
      theta_t[i] <- rnorm(1,mean = theta_t[i-1],sd = sqrt(sigmasq))
      posrate_t[i] <- inv.logit(theta_t[i])
      
    }
    P_t <- rbinom(n = length(t),size = N, prob = posrate_t)
    
    for(k in 1:K){
      
      for(i in 1:length(t)){
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
    gamma_1k <- c(0,rnorm(K-1,mean = 0, sd = rep(0.1,K-1)))
    
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
    
    pisq <- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = sqrt(1/200))
    
    
    
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


generate.model.ests <- function(model.string, data.list, params ,n.chains =3, n.iter = 25000, thin = 10){
  
  jags.mod <- jags.model(textConnection(model.string), 
                         data = data.list, n.chains = n.chains, n.adapt = 10000,quiet = T)
  
  samps <- coda.samples(jags.mod, params, n.iter = n.iter, thin = 10, progress.bar = "none")
  #as of right now the estimate is a mean.
  bayes.est <- summary(samps)$statistics[,1]
  return(bayes.est)
}