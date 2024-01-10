
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
    smalln[k,] <- rep(100,length(t))
    
    if(k>1){
      smalln[k,] <- rep(1000,length(t))
    }
  }
  
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

#JAGS MODEL for constant phi
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
		
		Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1);
sigmasq ~ dnorm(0, 1/0.1)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);
}
}')

mod.linear.phi <- custommodel('
model{	
#likelihood

for (t in 1:T){
		phi[1,t] <- 1	}

for (k in 2:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times_centered[k,t])
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
theta0 ~ dnorm(0, 1);
sigmasq ~ dnorm(0, 1/0.1)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);
	gamma1[k] ~ dnorm(0, 1/0.01);
}
}')

mod.walk.phi <- custommodel('
model{	
#likelihood

for (i in 1:T){
		phi[1,i] <- 1	
}


for (k in 2:K){

  gamma[k,1] ~ dnorm(gamma0[k],1/pisq)
  phi[k,1] <- exp(gamma[k,1])
  
	for (t in 2:T){
	  gamma[k,t] ~ dnorm(gamma[k,t-1], 1/pisq)
	  phi[k,t] <- exp(gamma[k,t])
	  
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],1/sigmasq)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}


for (k in 1:K){
	for (t in 1:T){
		
		Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1);
sigmasq ~ dnorm(0, 1/0.1)T(0,);
pisq ~ dnorm(0, 1/0.01)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);

}

}')

#helper functions
#source("Scripts/helperfunctions.R")




extract.unbiased <- function(datalist){
  K = 1
  T <- datalist$T
  new.list <- list(K = K, 
                   times = matrix(datalist$times[1,],ncol = T), N = datalist$N, T = T,
                   Y = matrix(datalist$Y[1,],ncol = T), smalln = matrix(datalist$smalln[1,],ncol = T))
  
  return(new.list)
}

get.mean <- function(mcmc.obj,timepoint){
  return(summary(mcmc.obj)$statistics[,1][timepoint])
}

get.var <- function(mcmc.obj,timepoint){
  return(summary(mcmc.obj)$statistics[,2][timepoint])
}

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




generate.unidentifiable.data<- function(phi = "constant",t = 1:5,NN,K = 3){
  
  #function creates all data for simulation
  list.return <- list()
  j = 1
  for(i in 1:NN){

    gen.data <- generate.dataset(phi = phi,t = t)

    if(sum(gen.data$Y[2:K,]==1000 | gen.data$Y[2:K,]==0)>=1){
      print(paste("bad data on i =", i))

      list.return[[j]] <- gen.data
      j = j + 1
      
    }
    

  }
  
  return(list.return)
}

check.bad.data <- function(data.list,t = 1:5,phi = 'constant',K = 3){
  
  for(i in 1:length(data.list)){
    
    if(sum(data.list[[i]]$Y[2:K,]==1000 | data.list[[i]]$Y[2:K,]==0)>=1){
      print(paste("bad data on i =", i))
      
      data.list[[i]] <- generate.dataset(phi = phi,t = t,K = K)
      
    }
  }
  
  return(data.list)
  
}


NN <- 8000

set.seed(519523)
data.const <- generate.unidentifiable.data(phi = "constant", NN = NN,t = 1:10)
data.linear  <- generate.unidentifiable.data(phi = "linear", NN = NN,t = 1:10)
data.walk  <- generate.unidentifiable.data(phi = "walk", NN = NN,t = 1:10)

NN <- min(length(data.const), length(data.linear), length(data.walk))

cl <- makePSOCKcluster(10)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")

print("T = 10 Simulations")

run.simulation <- function(cl,data.const, data.linear, data.walk, NN = 500, ti = 5){
  
  error <- matrix(NA,nrow = NN, ncol = 9)
  colnames(error) <- c("const var","const bias","const 95 contain",
                       "linear var","linear bias", "linear 95 contain", 
                       "walk var","walk bias", "walk 95 contain")
  only.unbiased <-  matrix(NA, nrow = NN, ncol = 9)
  colnames(only.unbiased) <- c("const var","const bias","const 95 contain",
                               "linear var","linear bias", "linear 95 contain", 
                               "walk var","walk bias", "walk 95 contain")
  
  
  for(j in 1:NN){
    
    
    chain1<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(j))
    chain2<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(j+1))
    chain3<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(j+2))
    chain4<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(j+4))
    chain5<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(j+5))
    chain6<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(j+6))
    chain7<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(j+7))
    chain8<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(j+8))
    chain9<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(j+9))
    chain10<- list(.RNG.name = "base::Super-Duper", 
                   .RNG.seed = c(j+10))
    
    chains.init <- list(chain1,chain2,chain3,chain4,chain5,chain6,chain7,chain8,chain9,chain10)
    
    
    
    data.const.phi <- data.const[[j]]
    data.linear.phi <- data.linear[[j]]
    data.walk.phi <- data.walk[[j]]
    
    #center linear function
    data.const.phi$times_centered <- data.const.phi$times- median(data.const.phi$times)
    data.linear.phi$times_centered <- data.linear.phi$times - median(data.linear.phi$times)
    data.walk.phi$times_centered <- data.walk.phi$times- median(data.walk.phi$times)
    
    
    pos.rate.const <- data.const.phi$params[grep(paste0("posrate",ti),names(data.const.phi$params))]
    pos.rate.linear <- data.linear.phi$params[grep(paste0("posrate",ti),names(data.linear.phi$params))]
    pos.rate.walk <- data.walk.phi$params[grep(paste0("posrate",ti),names(data.walk.phi$params))]
    
    unbiased.const.phi <- extract.unbiased(data.const.phi)
    unbiased.linear.phi <- extract.unbiased(data.linear.phi)
    unbiased.walk.phi <- extract.unbiased(data.walk.phi)
    
    #data x model
    
    if(j %% 10 ==0) print(j)
    
    const.const <- jags.parfit(cl, data.const.phi, "positiverate", mod.const.phi,
                               n.chains = 10,n.adapt = 25000,thin = 5, n.iter = 70000,inits = chains.init)
    
    const.summary <- summary(const.const)
    const.lower <- const.summary$quantiles[,1][ti]
    const.upper <- const.summary$quantiles[,5][ti]
    
    error[j,"const var"] <- (const.summary$statistics[,2][ti])^2
    error[j,"const bias"] <- pos.rate.const - const.summary$statistics[,1][ti]
    error[j,"const 95 contain"] <- ifelse(pos.rate.const <= const.upper & pos.rate.const >= const.lower,TRUE,FALSE)
    
    
    linear.linear <- jags.parfit(cl, data.linear.phi, "positiverate", mod.linear.phi,
                                 n.chains = 10,n.adapt = 25000,thin = 5, n.iter = 70000,inits = chains.init)
    linear.summary <- summary(linear.linear)
    linear.lower <- linear.summary$quantiles[,1][ti]
    linear.upper <- linear.summary$quantiles[,5][ti]
    
    error[j,"linear var"] <- (linear.summary$statistics[,2][ti])^2
    error[j,"linear bias"] <- pos.rate.linear - linear.summary$statistics[,1][ti]
    error[j,"linear 95 contain"] <- ifelse(pos.rate.linear<= linear.upper & pos.rate.linear>= linear.lower,TRUE,FALSE)
    
    
    walk.walk <- jags.parfit(cl, data.walk.phi, "positiverate", mod.walk.phi,
                             n.chains = 10,n.adapt = 25000,thin = 5, n.iter = 70000,inits = chains.init)
    walk.summary <- summary(walk.walk)
    walk.lower <- walk.summary$quantiles[,1][ti]
    walk.upper <- walk.summary$quantiles[,5][ti]
    
    error[j,"walk var"] <- (walk.summary$statistics[,2][ti])^2
    error[j,"walk bias"] <- pos.rate.walk - walk.summary$statistics[,1][ti]
    error[j,"walk 95 contain"] <- ifelse(pos.rate.walk <= walk.upper & pos.rate.walk >= walk.lower,TRUE,FALSE)
    
    
    unb.const<- jags.parfit(cl, unbiased.const.phi, "positiverate", mod.const.phi,
                            n.chains = 10,n.adapt =20000,thin = 5, n.iter = 50000,inits = chains.init)
    
    unb.const.summary <- summary(unb.const)
    unb.const.lower <- unb.const.summary$quantiles[,1][ti]
    unb.const.upper <- unb.const.summary$quantiles[,5][ti]
    
    only.unbiased[j,"const var"] <- (unb.const.summary$statistics[,2][ti])^2
    only.unbiased[j,"const bias"] <- pos.rate.const - unb.const.summary$statistics[,1][ti]
    only.unbiased[j,"const 95 contain"] <- ifelse(pos.rate.const <= unb.const.upper & pos.rate.const >= unb.const.lower,TRUE,FALSE)
    
    
    unb.linear <- jags.parfit(cl, unbiased.linear.phi, "positiverate", mod.const.phi,
                              n.chains = 10,n.adapt = 20000,thin = 5, n.iter = 50000,inits = chains.init)
    unb.linear.summary <- summary(unb.linear)
    unb.linear.lower <- unb.linear.summary$quantiles[,1][ti]
    unb.linear.upper <- unb.linear.summary$quantiles[,5][ti]
    
    only.unbiased[j,"linear var"] <- (unb.linear.summary$statistics[,2][ti])^2
    only.unbiased[j,"linear bias"] <- pos.rate.linear - unb.linear.summary$statistics[,1][ti]
    only.unbiased[j,"linear 95 contain"] <- ifelse(pos.rate.linear <= unb.linear.upper & pos.rate.linear  >= unb.linear.lower,TRUE,FALSE)
    
    
    unb.walk <- jags.parfit(cl, unbiased.walk.phi, "positiverate", mod.const.phi,
                            n.chains = 10,n.adapt = 20000,thin = 5, n.iter = 50000,inits = chains.init)
    
    unb.walk.summary <- summary(unb.walk)
    unb.walk.lower <- unb.walk.summary$quantiles[,1][ti]
    unb.walk.upper <- unb.walk.summary$quantiles[,5][ti]
    
    only.unbiased[j,"walk var"] <- (unb.walk.summary$statistics[,2][ti])^2
    only.unbiased[j,"walk bias"] <- pos.rate.walk - unb.walk.summary$statistics[,1][ti]
    only.unbiased[j,"walk 95 contain"] <- ifelse(pos.rate.walk <= unb.walk.upper & pos.rate.walk  >= unb.walk.lower,TRUE,FALSE)
    
    
  }
  
  result.final <- rbind(error,only.unbiased)
  
  return(result.final)
  
  
}

#run the simulation
results.plot.final  <- run.simulation(cl,data.const,data.linear,data.walk,NN = 100,ti = 10)

#end parallel
stopCluster(cl)

write.csv(results.plot.final,"CheckCIT10N1000.csv",row.names = F)

models <- results.plot.final[1:100,]
unbiased <- results.plot.final[101:200,]
apply(models*100,2,mean)
apply(unbiased*100,2,mean)

