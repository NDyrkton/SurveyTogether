
#Survey Together Simulation Study
#3x3 Design         Nathaniel Dyrkton  Supervised by Paul Gustafson and Harlan Campbell
library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)

#First step is to generate data under 3 conditions for phi (bias term)
#1: Phi is constant by  time, 2: Phi is linear in time, and 3: Phi follows a random walk 

#JAGS MODEL for constant phi
mod.constant.phi<- '
model{	
#likelihood

for (i in 1:Ivec[1]){
		phi[1,i] <- 1	}

for (k in 2:K){
	for (i in 1:Ivec[k]){
		phi[k,i] <- exp(gamma0[k])
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

#helper functions

inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

generate.dataset <- function(N= 10000, K =3, t = c(1,2,3), ns = c(100,100,100), phi = "constant"){
  Y <- matrix(NA,ncol = length(t),nrow = K)
  Ivec <- rep(3,K)
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

    
    parameters = c(gamma0,theta_t, theta0)
    names(parameters) = c("gamma01","gamma02","gamma03","theta_1","theta_2", "theta_3","theta0")
    
    return(list(K=K, Ivec= Ivec,T=3, times=times, N=N, Y=Y, smalln=smalln, params = parameters))
    
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

    parameters <- c(gamma_0k,gamma_1k,theta_t, theta0)
    names(parameters) <-  c("gamma01","gamma02","gamma03",
                            "gamma11","gamma12","gamma13" ,"theta_1","theta_2", "theta_3","theta0")
    return(list(K=K, Ivec= Ivec,T=3, times=times, N=N, Y=Y, smalln=smalln,params = parameters))
    
  }else if(phi == "walk"){
    
    phi_kt <- matrix(NA,nrow =K,ncol = length(t))
    gamma_kt <- matrix(NA,nrow = K,ncol = length(t))
    
    gamma_0k <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1)))
    gamma_kt[1,] <- c(0,0,0) 
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
    parameters <- c(gamma_0k,as.numeric(gamma_kt),theta_t,theta0)
    names(parameters) <- c("gamma01","gamma02","gamma03","gamma11",
                           "gamma21","gamma31","gamma12", "gamma22", "gamma32",
                           "gamma13","gamma23","gamma33", "theta_1","theta_2",
                           "theta_3","theta0")
    
    return(list(K=K, Ivec= Ivec,T=3, times=times, N=N, Y=Y, smalln=smalln, params = parameters))
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
  new.list <- list(K = K, Ivec = Ivec, 
                   times = matrix(datalist$times[1,],ncol = 3), N = datalist$N, T = 3,
                   Y = matrix(datalist$Y[1,],ncol = 3), smalln = matrix(datalist$smalln[1,],ncol = 3))
  
  return(new.list)
}



NN <- 300
set.seed(98301)

error <- matrix(NA,nrow = NN, ncol = 9)
colnames(error) <- c("const x const","const x linear", "const x walk",
                                           "linear x const", "linear x linear","linear x walk",
                                           "walk x const", "walk x linear", "walk x walk")

only.unbiased <-  matrix(NA, nrow = NN, ncol = 9)
colnames(only.unbiased) <- c("const x const","const x linear", "const x walk",
                     "linear x const", "linear x linear","linear x walk",
                     "walk x const", "walk x linear", "walk x walk")


for(j in 1:NN){
  
  data.const.phi <- generate.dataset(phi = 'constant')
  data.linear.phi <- generate.dataset(phi = 'linear')
  data.walk.phi <- generate.dataset(phi = 'walk')
  
  pos.rate.const <- data.const.phi$params["theta_3"]
  pos.rate.linear <- data.linear.phi$params["theta_3"]
  pos.rate.walk <- data.walk.phi$params["theta_3"]
  
  unbiased.const.phi <- extract.unbiased(data.const.phi)
  unbiased.linear.phi <- extract.unbiased(data.linear.phi)
  unbiased.walk.phi <- extract.unbiased(data.walk.phi)
  
  #data x model
  
  if(j %% 5 ==0) print(j)
  
  const.const <- generate.model.ests(model.string = mod.constant.phi, data.list = data.const.phi[-8], params = "positiverate")
  error[j,"const x const"] <- pos.rate.const - const.const[3]
  
  const.linear <- generate.model.ests(model.string = mod.linear.phi, data.list = data.const.phi[-8], params = "positiverate")
  error[j, "const x linear"] <- pos.rate.const - const.linear[3]
  
  const.walk <- generate.model.ests(model.string = mod.walk.phi, data.list = data.const.phi[-8], params = "positiverate")
  error[j, "const x walk"] <-  pos.rate.const - const.walk[3]
  
  linear.const <- generate.model.ests(model.string = mod.constant.phi, data.list = data.linear.phi[-8], params = "positiverate")
  error[j, "linear x const"] <- pos.rate.linear - linear.const[3]
  
  linear.linear <- generate.model.ests(model.string = mod.linear.phi, data.list = data.linear.phi[-8], params = "positiverate")
  error[j, "linear x linear"] <-pos.rate.linear - linear.linear[3]
  
  linear.walk <- generate.model.ests(model.string = mod.walk.phi, data.list = data.linear.phi[-8], params = "positiverate")
  error[j, "linear x walk"] <- pos.rate.linear - linear.walk[3]
  
  walk.const <- generate.model.ests(model.string = mod.constant.phi, data.list = data.walk.phi[-8], params = "positiverate")
  error[j, "walk x const"] <- pos.rate.walk - walk.const[3]
  
  walk.linear <- generate.model.ests(model.string = mod.linear.phi, data.list = data.walk.phi[-8], params = "positiverate")
  error[j, "walk x linear"] <- pos.rate.walk - walk.linear[3]
  
  walk.walk <- generate.model.ests(model.string = mod.walk.phi, data.list = data.walk.phi[-8], params = "positiverate")
  error[j, "walk x walk"] <- pos.rate.walk - walk.walk[3]
  
  #one unbiased survey only.
  #data x model
  
  unb.const.const <- generate.model.ests(model.string = mod.constant.phi, data.list = unbiased.const.phi, params = "positiverate")
  only.unbiased[j, "const x const"] <- pos.rate.const - unb.const.const[3]
  
  unb.const.linear <- generate.model.ests(model.string = mod.linear.phi, data.list = unbiased.const.phi, params = "positiverate")
  only.unbiased[j, "const x linear"] <- pos.rate.const - unb.const.linear[3]
  
  unb.const.walk <- generate.model.ests(model.string = mod.walk.phi, data.list = unbiased.const.phi, params = "positiverate")
  only.unbiased[j, "const x walk"] <- pos.rate.const - unb.const.walk[3]
  
  unb.linear.const <- generate.model.ests(model.string = mod.constant.phi, data.list = unbiased.linear.phi, params = "positiverate")
  only.unbiased[j, "linear x const"] <- pos.rate.linear - unb.linear.const[3]
  
  unb.linear.linear <- generate.model.ests(model.string = mod.linear.phi, data.list = unbiased.linear.phi, params = "positiverate")
  only.unbiased[j, "linear x linear"] <- pos.rate.linear - unb.linear.linear[3]
  
  unb.linear.walk <- generate.model.ests(model.string = mod.walk.phi, data.list = unbiased.linear.phi, params = "positiverate")
  only.unbiased[j, "linear x walk"] <- pos.rate.linear - unb.linear.walk[3]
  
  unb.walk.const <- generate.model.ests(model.string = mod.constant.phi, data.list = unbiased.walk.phi, params = "positiverate")
  only.unbiased[j, "walk x const"] <- pos.rate.walk - unb.walk.const[3]
  
  unb.walk.linear <- generate.model.ests(model.string = mod.linear.phi, data.list = unbiased.walk.phi, params = "positiverate")
  only.unbiased[j, "walk x linear"] <- pos.rate.walk - unb.walk.linear[3]
  
  unb.walk.walk <- generate.model.ests(model.string = mod.walk.phi, data.list = unbiased.walk.phi, params = "positiverate")
  only.unbiased[j, "walk x walk"] <- pos.rate.walk - unb.walk.walk[3]
  

}

result.RMSE <- apply(error,2,function(x){sqrt(mean((100*x)^2))})

result.RMSE.unb <- apply(only.unbiased,2,function(x){sqrt(mean((100*x)^2))})




results.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const","linear","walk"),3),RMSE = result.RMSE)
results.plot.unb <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const.1","linear.1","walk.1"),3),RMSE = result.RMSE.unb)

#results.final <- rbind(results)

ggplot(data = results.plot, aes(x = data, y = RMSE,group = model,colour = model)) + geom_point() + geom_line() + theme_minimal() + 
  labs(x = "Data generation", y = "Root Mean Squared Error", title = "RMSE of NN = 300 in 3x3 design") + geom_point(data = results.plot.unb)+
  geom_line(data = results.plot.unb,linetype = "dashed",aes(colour = model,group=model))




write.csv(results.plot,"C:/Users/ndyrk/OneDrive/Masters/Masters Project/RMSE_3x3.csv",row.names = T)

write.csv(results.plot.unb,"C:/Users/ndyrk/OneDrive/Masters/Masters Project/RMSE_3x3_unb.csv",row.names = T)

write.csv(error,"C:/Users/ndyrk/OneDrive/Masters/Masters Project/RMSE_total.csv",row.names = T)





###################################
###########
###   old code
### 
###########
##################















######OLD



#START SIMULATION
NN <- 20
set.seed(983012983)

positive.rate.squared.error <- matrix(NA,nrow = NN,ncol = 9)

colnames(positive.rate.squared.error) <- c("const x const","const x linear", "const x walk",
                                           "linear x const", "linear x linear","linear x walk",
                                           "walk x const", "walk x linear", "walk x walk")



for(j in 1:NN){
  data.const.phi <- generate.dataset(phi = 'constant')
  data.linear.phi <- generate.dataset(phi = 'linear')
  data.walk.phi <- generate.dataset(phi = 'walk')
  if(j %% 10 ==0) print(j)
  
  #run jags model
  
  #constant data, constant model
  jags.const.phi.const <- jags.model(textConnection(mod.constant.phi), 
                                     data = data.const.phi[-8],  
                                     n.chains = 3, n.adapt = 10000,quiet = T)
  
  params.const.phi <- c("theta0", "positiverate", "gamma0")	
  samps.const.phi.const <- coda.samples(jags.const.phi.const,
                                        params.const.phi, n.iter = 20000, thin = 10, progress.bar = "none")
  mean.const.phi.const <- summary(samps.const.phi.const)$statistics[,1]
  positive.rate.squared.error[j,"const x const"] <- (data.const.phi$params[6]-mean.const.phi.const[grep("positiverate[3]",names(mean.const.phi.const),fixed = T)])
  
  #linear data const model
  
  jags.linear.phi.const <- jags.model(textConnection(mod.constant.phi), 
                                      data = data.linear.phi[-8],  
                                      n.chains = 3, n.adapt = 10000,quiet = T)
  
  samps.linear.phi.const <- coda.samples(jags.linear.phi.const,
                                         params.const.phi, n.iter = 20000, thin = 10, progress.bar = "none")
  mean.linear.phi.const <- summary(samps.linear.phi.const)$statistics[,1]
  positive.rate.squared.error[j,"linear x const"] <- (data.linear.phi$params[9]-mean.linear.phi.const[grep("positiverate[3]",names(mean.linear.phi.const),fixed = T)])
  
  #walk data const model
  
  jags.walk.phi.const <- jags.model(textConnection(mod.constant.phi), 
                                    data = data.walk.phi[-8],  
                                    n.chains = 3, n.adapt = 10000,quiet = T)
  
  params.walk.phi <- c("theta0", "positiverate", "gamma0","phi")	
  samps.walk.phi.const <- coda.samples(jags.walk.phi.const,
                                       params.const.phi, n.iter = 20000, thin = 10, progress.bar = "none")
  
  mean.walk.phi.const <- summary(samps.walk.phi.const)$statistics[,1]
  positive.rate.squared.error[j,"walk x const"] <- (data.walk.phi$params[15]-mean.walk.phi.const[grep("positiverate[3]",names(mean.walk.phi.const),fixed = T)])
  
  
  #constant data, linear model
  
  jags.const.phi.linear <- jags.model(textConnection(mod.linear.phi), 
                                      data = data.const.phi[-8],  
                                      n.chains = 3, n.adapt = 10000,quiet = T)
  params.linear.phi <- c("theta0", "gamma0", "gamma1","positiverate")	
  samps.const.phi.linear <- coda.samples(jags.const.phi.linear, 
                                         params.linear.phi, n.iter = 20000, thin = 10, progress.bar = "none")
  mean.const.phi.linear <- summary(samps.const.phi.linear)$statistics[,1]
  
  positive.rate.squared.error[j,"const x linear"] <- (data.const.phi$params[6]-mean.const.phi.linear[grep("positiverate[3]",names(mean.const.phi.linear),fixed = T)])
  
  #linear data, linear model
  jags.linear.phi.linear <- jags.model(textConnection(mod.linear.phi), 
                                       data = data.linear.phi[-8],  
                                       n.chains = 3, n.adapt = 10000, quiet = T)
  samps.linear.phi.linear <- coda.samples(jags.linear.phi.linear,
                                          params.linear.phi, n.iter = 20000, thin = 10, progress.bar = "none")
  mean.linear.phi.linear <- summary(samps.linear.phi.linear)$statistics[,1]
  positive.rate.squared.error[j,"linear x linear"] <- (data.linear.phi$params[9]-mean.linear.phi.linear[grep("positiverate[3]",names(mean.linear.phi.linear),fixed = T)])
  
  #walk data, linear model
  jags.walk.phi.linear <- jags.model(textConnection(mod.linear.phi), 
                                     data = data.walk.phi[-8],  
                                     n.chains = 3, n.adapt = 10000, quiet = T)
  samps.walk.phi.linear <- coda.samples(jags.walk.phi.linear,
                                        params.linear.phi, n.iter = 20000, thin = 10, progress.bar = "none")
  mean.walk.phi.linear <- summary(samps.walk.phi.linear)$statistics[,1]
  positive.rate.squared.error[j,"walk x linear"] <- (data.walk.phi$params[15]-mean.walk.phi.linear[grep("positiverate[3]",names(mean.walk.phi.linear),fixed = T)])
  
  #constant data, walk model
  
  jags.const.phi.walk <- jags.model(textConnection(mod.walk.phi), 
                                    data = data.const.phi[-8],  
                                    n.chains = 3, n.adapt = 10000, quiet = T)
  samps.const.phi.walk <- coda.samples(jags.const.phi.walk,
                                       params.walk.phi, n.iter = 20000, thin = 10, progress.bar = "none")
  mean.const.phi.walk <- summary(samps.const.phi.walk)$statistics[,1]
  positive.rate.squared.error[j,"const x walk"] <- (data.const.phi$params[6]-mean.const.phi.walk[grep("positiverate[3]",names(mean.const.phi.walk),fixed = T)])
  
  #linear data, walk model
  
  jags.linear.phi.walk <- jags.model(textConnection(mod.walk.phi), 
                                     data = data.linear.phi[-8],  
                                     n.chains = 3, n.adapt = 10000, quiet = T)
  samps.linear.phi.walk <- coda.samples(jags.linear.phi.walk,
                                        params.walk.phi, n.iter = 20000, thin = 10, progress.bar = "none")
  mean.linear.phi.walk <- summary(samps.linear.phi.walk)$statistics[,1]
  positive.rate.squared.error[j,"linear x walk"] <- (data.linear.phi$params[9]-mean.linear.phi.walk[grep("positiverate[3]",names(mean.linear.phi.walk),fixed = T)])
  
  #walk data, walk model
  
  jags.walk.phi.walk <- jags.model(textConnection(mod.walk.phi), 
                                   data = data.walk.phi[-8],  
                                   n.chains = 3, n.adapt = 10000, quiet = T)
  samps.walk.phi.walk <- coda.samples(jags.walk.phi.walk,
                                      params.walk.phi, n.iter = 20000, thin = 10, progress.bar = "none")
  mean.walk.phi.walk <- summary(samps.walk.phi.walk)$statistics[,1]
  positive.rate.squared.error[j,"walk x walk"] <- (data.walk.phi$params[15]-mean.walk.phi.walk[grep("positiverate[3]",names(mean.walk.phi.walk),fixed = T)])
  
  
  
}

result.RMSE <- apply(positive.rate.squared.error,2,function(x){sqrt(mean((100*x)^2))})

library(ggplot2)
results.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const","linear","Walk"),3),RMSE = result.RMSE)

ggplot(data = results.plot, aes(x = data, y = RMSE,group = model,colour = model)) + geom_point() + geom_line() + theme_minimal() + 
  labs(x = "Data generation", y = "Root Mean Squared Error", title = "RMSE of NN = 100 in 3x3 design")

#write actual data set
write.csv(positive.rate.squared.error,"C:/Users/ndyrk/OneDrive/Masters/Masters Project/fullrun.csv",row.names = T)

write.csv(result.RMSE,"C:/Users/ndyrk/OneDrive/Masters/Masters Project/survey_together_simulation.csv",row.names = T)


