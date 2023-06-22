#Survey Together Simulation Study
#3x3 Design  
#Simplified only using Normal Model
#Nathaniel Dyrkton  Supervised by Paul Gustafson and Harlan Campbell

#known sigma, unknown mu, unknown phi.

library(rjags)
library(MCMCpack)
library(ggplot2)
library(truncnorm)
library(dclone)

#models

mod.const.phi<- custommodel('
model{	
#likelihood

for (i in 1:Ivec[1]){
		phi[1,i] <- 0	}

for (k in 2:K){
	for (i in 1:Ivec[k]){
		phi[k,i] <- gamma0[k]
	}
}
	
	
mu[1] ~ dnorm(mu0,pow(5,-2))
for(t in 2:T){
	mu[t]	~ dnorm(mu[t-1],rho)
}

for (k in 1:K){
	for (i in 1:Ivec[k]){
		
		Y[k,i] ~ dnorm(mu[times[k,i]] + phi[k,i], smalln[k,i]/100);
	}
}

#priors
mu0 ~ dnorm(75, pow(5,-2));
rho ~ dgamma(2,1/0.05)

gamma0[1]~ dnorm(0,1)
for (k in 2:K){
	gamma0[k] ~ dnorm(5, pow(3,-2));
}
}')

mod.linear.phi <- custommodel('
model{	
#likelihood

for (i in 1:Ivec[1]){
		phi[1,i] <- 0	}

for (k in 2:K){
	for (i in 1:Ivec[k]){
		phi[k,i] <- (gamma0[k] + gamma1[k]*times[k,i])
	}
}
	
	
mu[1] ~ dnorm(mu0,pow(5,-2))
for(t in 2:T){
	mu[t]	~ dnorm(mu[t-1],rho)
}

for (k in 1:K){
	for (i in 1:Ivec[k]){
		
		Y[k,i] ~ dnorm(mu[times[k,i]] + phi[k,i], smalln[k,i]/100);
	}
}

#priors
mu0 ~ dnorm(75, pow(5,-2));
rho ~ dgamma(2,1/0.05)

for (k in 2:K){
	gamma0[k] ~ dnorm(5, pow(3,-2));
	gamma1[k] ~ dnorm(2, 1)
}
}')

mod.walk.phi <- custommodel('
model{	
#likelihood

for (i in 1:Ivec[1]){
		phi[1,i] <- 0
		gamma[1,i] <- 0
}


for (k in 2:K){

  gamma[k,1] ~ dnorm(gamma0[k],1)
  phi[k,1] <- exp(gamma[k,1])
  
	for (t in 2:T){
	  gamma[k,t] ~ dnorm(gamma[k,t-1], pi)
	  phi[k,t] <- exp(gamma[k,t])
	  
	}
}
	
	
mu[1] ~ dnorm(mu0,pow(5,-2))
for(t in 2:T){
	mu[t]	~ dnorm(mu[t-1],rho)
}

for (k in 1:K){
	for (i in 1:Ivec[k]){
		
		Y[k,i] ~ dnorm(mu[times[k,i]] + phi[k,i], smalln[k,i]/100);
	}
}

#priors
mu0 ~ dnorm(75, pow(5,-2))
rho ~ dgamma(2,1/0.05)
pi ~ dgamma(2,1/0.5)


for (k in 2:K){
	gamma0[k] ~ dnorm(5, pow(3,-2));

}
}')




gen.norm.dat <- function(N= 10000, K =3, ts = c(1:5), ns = rep(100,length(ts)), phi = "constant"){
  sigma <- 5
  mu0 <- rnorm(1,mean = 75, sd = 5)
  len <- length(ts)
  mu <- numeric(len)
  rho <-  rgamma(1,shape = 2,scale = 0.05) 
  Y <- matrix(NA,nrow = K,ncol = len)
  times <- t(matrix(rep(ts,K),ncol = K))
  smalln <- t(matrix(rep(ns,K),ncol = K))
  Ivec <- rep(3,K)
  
  
  if(phi == "constant"){
    phi <- c(0,rnorm(K-1,mean = 5,sd =5))
    
    #random walk for mu
    mu[1] <- rnorm(1,mean = mu0,sd = 5)
    
    for(t in 2:len){
      mu[t] <- rnorm(1,mean = mu[t-1], sd= sqrt(1/rho))
      
    }
    # the biased surveys
    for(k in 1:K){
      Y[k,] <- rnorm(len,mean =mu + phi[k], sd = rep(10,len)/sqrt(ns)   )
    }
    
    parameters <- c(mu0,phi,rho,mu)
    phi.names <- paste(rep("gamma",K),1:K,sep = '')
    mu.names <- paste(rep("mu",len),ts,sep = '')
    names(parameters) <- c("mu0",phi.names,"rho",mu.names)
    
    return(list(K=K, Ivec= Ivec,T=5, times=times, N=N, Y=Y, smalln=smalln, params = parameters))
    
  }else if(phi == "linear"){
    
    gamma0 <- c(0,rnorm(K-1,mean = 5,sd = c(3,3)))
    gamma1 <- c(0,rnorm(K-1,mean = 2,sd = c(1,1)))
    
    
    #random walk for mu
    mu[1] <- rnorm(1,mean = mu0,sd = 5)
    
    for(t in 2:len){
      mu[t] <- rnorm(1,mean = mu[t-1], sd= sqrt(1/rho))
      
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
    
    phi0 <- c(0, rnorm(K-1,mean = 5, sd = 2))
    phi_kt <- matrix(NA,nrow = K,ncol = len)
    phi_kt[1,] <- 0
    pi <- rgamma(1,shape = 2,scale = 0.5)
    
    
    #random walk for mu
    mu[1] <- rnorm(1,mean = mu0,sd = 5)
    
    #random walk for phi
    phi_kt[2:K,1] <- rnorm(K-1,mean = phi0[2:K],sd = c(1,1))
    
    
    for(t in 2:len){
      mu[t] <- rnorm(1,mean = mu[t-1], sd= sqrt(1/rho))
      phi_kt[2:K,t] <- rnorm(K-1,mean = phi_kt[2:K,t-1],sd = sqrt(1/c(pi,pi)))
      
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

get.mean <- function(mcmc.obj){
  return(summary(mcmc.obj)$statistics[,1][5])
}



#run sim

#start parallel
cl <- makePSOCKcluster(3)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


NN <- 100
set.seed(11947194)

error <- matrix(NA,nrow = NN, ncol = 9)
colnames(error) <- c("const x const","const x linear", "const x walk",
                     "linear x const", "linear x linear","linear x walk",
                     "walk x const", "walk x linear", "walk x walk")

only.unbiased <-  matrix(NA, nrow = NN, ncol = 9)
colnames(only.unbiased) <- c("const x const","const x linear", "const x walk",
                             "linear x const", "linear x linear","linear x walk",
                             "walk x const", "walk x linear", "walk x walk")


dcoptions("verbose"=F)#mute the output

for(j in 1:NN){
  
  data.const.phi <- gen.norm.dat(phi = 'constant')
  data.linear.phi <- gen.norm.dat(phi = 'linear')
  data.walk.phi <- gen.norm.dat(phi = 'walk')
  
  pos.rate.const <- data.const.phi$params[grep("mu5",names(data.const.phi$params))]
  pos.rate.linear <- data.linear.phi$params[grep("mu5",names(data.linear.phi$params))]
  pos.rate.walk <- data.walk.phi$params[grep("mu5",names(data.walk.phi$params))]
  
  unbiased.const.phi <- extract.unbiased(data.const.phi)
  unbiased.linear.phi <- extract.unbiased(data.linear.phi)
  unbiased.walk.phi <- extract.unbiased(data.walk.phi)
  
  #data x model
  
  if(j %% 10 ==0) print(j)
  
  const.const <- jags.parfit(cl, data.const.phi[-8], "mu", mod.const.phi,
                             n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  error[j,"const x const"] <- pos.rate.const - get.mean(const.const)
  
  const.linear <- jags.parfit(cl, data.const.phi[-8], "mu", mod.linear.phi,
                              n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  error[j, "const x linear"] <- pos.rate.const - get.mean(const.linear)
  
  const.walk <- jags.parfit(cl, data.const.phi[-8], "mu", mod.walk.phi,
                            n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  error[j, "const x walk"] <- pos.rate.const - get.mean(const.walk)
  
  linear.const <- jags.parfit(cl, data.linear.phi[-8], "mu", mod.const.phi,
                              n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  error[j, "linear x const"] <- pos.rate.linear - get.mean(linear.const)
  
  linear.linear <- jags.parfit(cl, data.linear.phi[-8], "mu", mod.linear.phi,
                               n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  error[j, "linear x linear"] <- pos.rate.linear - get.mean(linear.linear)
  
  linear.walk <- jags.parfit(cl, data.linear.phi[-8], "mu", mod.walk.phi,
                             n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  error[j, "linear x walk"] <- pos.rate.linear - get.mean(linear.walk)
  
  walk.const <- jags.parfit(cl, data.walk.phi[-8], "mu", mod.const.phi,
                            n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  error[j, "walk x const"] <- pos.rate.walk - get.mean(walk.const)
  
  walk.linear <- jags.parfit(cl, data.walk.phi[-8], "mu", mod.linear.phi,
                             n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  error[j, "walk x linear"] <- pos.rate.walk - get.mean(walk.linear)
  
  walk.walk <- jags.parfit(cl, data.walk.phi[-8], "mu", mod.walk.phi,
                           n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  error[j, "walk x walk"] <- pos.rate.walk - get.mean(walk.walk)
  
  #one unbiased survey only.
  #data x model
  
  unb.const.const <- jags.parfit(cl, unbiased.const.phi, "mu", mod.const.phi,
                                 n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  only.unbiased[j, "const x const"] <- pos.rate.const - get.mean(unb.const.const)
  
  unb.const.linear <- jags.parfit(cl, unbiased.const.phi, "mu", mod.linear.phi,
                                  n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  only.unbiased[j, "const x linear"] <- pos.rate.const - get.mean(unb.const.linear)
  
  unb.const.walk <- jags.parfit(cl, unbiased.const.phi, "mu", mod.walk.phi,
                                n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  only.unbiased[j, "const x walk"] <- pos.rate.const - get.mean(unb.const.walk)
  
  unb.linear.const <- jags.parfit(cl, unbiased.linear.phi, "mu", mod.const.phi,
                                  n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  only.unbiased[j, "linear x const"] <- pos.rate.linear - get.mean(unb.linear.const)
  
  unb.linear.linear <- jags.parfit(cl, unbiased.linear.phi, "mu", mod.linear.phi,
                                   n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  only.unbiased[j, "linear x linear"] <- pos.rate.linear - get.mean(unb.linear.linear)
  
  unb.linear.walk <- jags.parfit(cl, unbiased.linear.phi, "mu", mod.walk.phi,
                                 n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  only.unbiased[j, "linear x walk"] <- pos.rate.linear - get.mean(unb.linear.walk)
  
  unb.walk.const <- jags.parfit(cl, unbiased.walk.phi, "mu", mod.const.phi,
                                n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  only.unbiased[j, "walk x const"] <- pos.rate.walk - get.mean(unb.walk.const)
  
  unb.walk.linear <- jags.parfit(cl, unbiased.walk.phi, "mu", mod.linear.phi,
                                 n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  only.unbiased[j, "walk x linear"] <- pos.rate.walk - get.mean(unb.walk.linear)
  
  unb.walk.walk <- jags.parfit(cl, unbiased.walk.phi, "mu", mod.walk.phi,
                               n.chains=3,n.adapt = 10000,thin = 10, n.iter = 25000)
  only.unbiased[j, "walk x walk"] <- pos.rate.walk - get.mean(unb.walk.walk)
  
  
}


result.RMSE <- apply(error,2,function(x){sqrt(mean((100*x)^2))})
result.RMSE.unb <- apply(only.unbiased,2,function(x){sqrt(mean((100*x)^2))})




results.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const","linear","walk"),3),RMSE = result.RMSE)
results.plot.unb <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const.1","linear.1","walk.1"),3),RMSE = result.RMSE.unb)

ggplot(data = results.plot, aes(x = data, y = RMSE,group = model,colour = model)) + geom_point() + geom_line() + theme_minimal() + 
  labs(x = "Data generation", y = "Root Mean Squared Error", title = "RMSE of NN = 500 in 3x3 design, 5 timepoints") + geom_point(data = results.plot.unb)+
  geom_line(data = results.plot.unb,linetype = "dashed",aes(colour = model,group=model)) + 
  scale_color_manual(values = c("const"="blue","const.1"="blue","linear" = "red","linear.1"="red","walk"="green","walk.1"="green"))



