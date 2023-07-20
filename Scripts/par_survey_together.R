
#Survey Together Simulation Study
#3x3 Design         Nathaniel Dyrkton  Supervised by Paul Gustafson and Harlan Campbell
library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)
library(dclone)

#First step is to generate data under 3 conditions for phi (bias term)
#1: Phi is constant by  time, 2: Phi is linear in time, and 3: Phi follows a random walk 

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
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.01)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1], pow(rho,-2))
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
theta0 ~ dnorm(0, 1);
rho ~ dnorm(0, 1/0.01)T(0,);

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
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times[k,t])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.01)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1], pow(rho,-2))
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
theta0 ~ dnorm(0, 1);
rho ~ dnorm(0, 1/0.01)T(0,);

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

  gamma[k,1] ~ dnorm(gamma0[k],1/0.01)
  phi[k,1] <- exp(gamma[k,1])
  
	for (t in 2:T){
	  gamma[k,t] ~ dnorm(gamma[k,t-1], pow(pi,-2))
	  phi[k,t] <- exp(gamma[k,t])
	  
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.01)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],pow(rho,-2))
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
theta0 ~ dnorm(0, 1);
rho ~ dnorm(0, 1/0.01)T(0,);
pi ~ dnorm(0, 1/0.01)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);

}

}')

#helper functions

inv.logit <- function(x){
  exp(x)/(1+exp(x))
}


#try with t = 7
generate.dataset <- function(N= 10000, K =3, t = c(1:3), ns = rep(100,length(t)), phi = "constant"){
  Y <- matrix(NA,ncol = length(t),nrow = K)
  smalln <- t(matrix(rep(ns,K),ncol = K))
  theta_t <- numeric(length(t))
  logit_theta_t <- numeric(length(t))
  times <- t(matrix(rep(t,K),ncol = K))
  
  #priors on general parameters
  rho <- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = 1/10)
  theta0 <- rnorm(1,mean =0, sd = 1)
  
  if(phi == "constant"){
    #generate param based on prior
    gamma0 <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1))) 
    #constant phi values
    phi <- exp(gamma0)
    
    
    logit_theta_t[1] <- rnorm(1,mean = theta0,sd = 1/10)
    theta_t[1] <- inv.logit(logit_theta_t[1])
    
    for(i in 2:length(t)){
      
      logit_theta_t[i] <- rnorm(1,mean = logit_theta_t[i-1],sd = rho)
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
      logit_theta_t[i] <- rnorm(1,mean = logit_theta_t[i-1],sd = rho)
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
      
      
      logit_theta_t[i] <- rnorm(1,mean = logit_theta_t[i-1],sd = rho)
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

extract.unbiased <- function(datalist){
  K = 1
  T <- datalist$T
  new.list <- list(K = K, 
                   times = matrix(datalist$times[1,],ncol = T), N = datalist$N, T = T,
                   Y = matrix(datalist$Y[1,],ncol = T), smalln = matrix(datalist$smalln[1,],ncol = T))
  
  return(new.list)
}

get.mean <- function(mcmc.obj){
  return(summary(mcmc.obj)$statistics[,1][9])
}

dcoptions("verbose"=F)#mute the output


#start parallel
cl <- makePSOCKcluster(3)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


NN <- 250
set.seed(12345)

error <- matrix(NA,nrow = NN, ncol = 9)
colnames(error) <- c("const x const","const x linear", "const x walk",
                     "linear x const", "linear x linear","linear x walk",
                     "walk x const", "walk x linear", "walk x walk")

only.unbiased <-  matrix(NA, nrow = NN, ncol = 9)
colnames(only.unbiased) <- c("const x const","const x linear", "const x walk",
                             "linear x const", "linear x linear","linear x walk",
                             "walk x const", "walk x linear", "walk x walk")



for(j in 1:NN){
  
  data.const.phi <- generate.dataset(phi = 'constant',t = 1:9)
  data.linear.phi <- generate.dataset(phi = 'linear', t = 1:9)
  data.walk.phi <- generate.dataset(phi = 'walk', t= 1:9)
  
  pos.rate.const <- data.const.phi$params[grep("theta9",names(data.const.phi$params))]
  pos.rate.linear <- data.linear.phi$params[grep("theta9",names(data.linear.phi$params))]
  pos.rate.walk <- data.walk.phi$params[grep("theta9",names(data.walk.phi$params))]
  
  unbiased.const.phi <- extract.unbiased(data.const.phi)
  unbiased.linear.phi <- extract.unbiased(data.linear.phi)
  unbiased.walk.phi <- extract.unbiased(data.walk.phi)
  
  #data x model
  
  if(j %% 10 ==0) print(j)
  
  const.const <- jags.parfit(cl, data.const.phi[-8], "positiverate", mod.const.phi,
                             n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  error[j,"const x const"] <- pos.rate.const - get.mean(const.const)
  
  const.linear <- jags.parfit(cl, data.const.phi[-8], "positiverate", mod.linear.phi,
                              n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  error[j, "const x linear"] <- pos.rate.const - get.mean(const.linear)
  
  const.walk <- jags.parfit(cl, data.const.phi[-8], "positiverate", mod.walk.phi,
                            n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  error[j, "const x walk"] <- pos.rate.const - get.mean(const.walk)
  
  linear.const <- jags.parfit(cl, data.linear.phi[-8], "positiverate", mod.const.phi,
                              n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  error[j, "linear x const"] <- pos.rate.linear - get.mean(linear.const)
  
  linear.linear <- jags.parfit(cl, data.linear.phi[-8], "positiverate", mod.linear.phi,
                               n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  error[j, "linear x linear"] <- pos.rate.linear - get.mean(linear.linear)
  
  linear.walk <- jags.parfit(cl, data.linear.phi[-8], "positiverate", mod.walk.phi,
                             n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  error[j, "linear x walk"] <- pos.rate.linear - get.mean(linear.walk)
  
  walk.const <- jags.parfit(cl, data.walk.phi[-8], "positiverate", mod.const.phi,
                            n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  error[j, "walk x const"] <- pos.rate.walk - get.mean(walk.const)
  
  walk.linear <- jags.parfit(cl, data.walk.phi[-8], "positiverate", mod.linear.phi,
                             n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  error[j, "walk x linear"] <- pos.rate.walk - get.mean(walk.linear)
  
  walk.walk <- jags.parfit(cl, data.walk.phi[-8], "positiverate", mod.walk.phi,
                           n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  error[j, "walk x walk"] <- pos.rate.walk - get.mean(walk.walk)
  
  #one unbiased survey only.
  #data x model
  
  unb.const.const <- jags.parfit(cl, unbiased.const.phi, "positiverate", mod.const.phi,
                                 n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  only.unbiased[j, "const x const"] <- pos.rate.const - get.mean(unb.const.const)
  
  unb.const.linear <- jags.parfit(cl, unbiased.const.phi, "positiverate", mod.linear.phi,
                                  n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  only.unbiased[j, "const x linear"] <- pos.rate.const - get.mean(unb.const.linear)
  
  unb.const.walk <- jags.parfit(cl, unbiased.const.phi, "positiverate", mod.walk.phi,
                                n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  only.unbiased[j, "const x walk"] <- pos.rate.const - get.mean(unb.const.walk)
  
  unb.linear.const <- jags.parfit(cl, unbiased.linear.phi, "positiverate", mod.const.phi,
                                  n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  only.unbiased[j, "linear x const"] <- pos.rate.linear - get.mean(unb.linear.const)
  
  unb.linear.linear <- jags.parfit(cl, unbiased.linear.phi, "positiverate", mod.linear.phi,
                                   n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  only.unbiased[j, "linear x linear"] <- pos.rate.linear - get.mean(unb.linear.linear)
  
  unb.linear.walk <- jags.parfit(cl, unbiased.linear.phi, "positiverate", mod.walk.phi,
                                 n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  only.unbiased[j, "linear x walk"] <- pos.rate.linear - get.mean(unb.linear.walk)
  
  unb.walk.const <- jags.parfit(cl, unbiased.walk.phi, "positiverate", mod.const.phi,
                                n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  only.unbiased[j, "walk x const"] <- pos.rate.walk - get.mean(unb.walk.const)
  
  unb.walk.linear <- jags.parfit(cl, unbiased.walk.phi, "positiverate", mod.linear.phi,
                                 n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  only.unbiased[j, "walk x linear"] <- pos.rate.walk - get.mean(unb.walk.linear)
  
  unb.walk.walk <- jags.parfit(cl, unbiased.walk.phi, "positiverate", mod.walk.phi,
                               n.chains=3,n.adapt = 5000,thin = 10, n.iter = 25000)
  only.unbiased[j, "walk x walk"] <- pos.rate.walk - get.mean(unb.walk.walk)
  
  
}


result.RMSE <- apply(error,2,function(x){sqrt(mean((100*x)^2))})
result.RMSE.unb <- apply(only.unbiased,2,function(x){sqrt(mean((100*x)^2))})




results.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const","linear","walk"),3),RMSE = result.RMSE)
results.plot.unb <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const.1","linear.1","walk.1"),3),RMSE = result.RMSE.unb)

ggplot(data = results.plot, aes(x = data, y = RMSE,group = model,colour = model)) + geom_point() + geom_line() + theme_minimal() + 
  labs(x = "Data generation", y = "Root Mean Squared Error", title = paste("RMSE of NN = ",NN,"in 3x3 design, 9 timepoints")) + geom_point(data = results.plot.unb)+
  geom_line(data = results.plot.unb,linetype = "dashed",aes(colour = model,group=model)) + 
  scale_color_manual(values = c("const"="blue","const.1"="blue","linear" = "red","linear.1"="red","walk"="green","walk.1"="green"))
 


#end parallel
stopCluster(cl)


write.csv(results.plot,"Results/RMSE_with_bias.csv",row.names = T)
write.csv(results.plot.unb,"Results/RMSE_only_unbiased.csv",row.names = T)
write.csv(error,"Results/RMSE_total.csv",row.names = T)


