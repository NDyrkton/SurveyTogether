
#Survey Together Simulation Study
#3x3 Design         Nathaniel Dyrkton  Supervised by Paul Gustafson and Harlan Campbell
library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)
library(dclone)

#First step is to generate data under 3 conditions for phi (bias term)
#1: Phi is constant by  time, 2: Phi is linear in time, and 3: Phi follows a random walk 


inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/(1-x))
}


#How fixed to be consistent with notation in paper.
generate.dataset <- function(N= 10000, K =3, t = c(1:5), phi = "constant"){
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

for(t in 1:T){
	P[t] ~ dbin(positiverate[t], N)
}

for (k in 1:K){
	for (t in 1:T){
		
		Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
	}
}

#priors
theta0 ~ dnorm(0, 1/0.5);
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
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times[k,t])
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
theta0 ~ dnorm(0, 1/0.5);
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

for(t in 1:T){
	P[t] ~ dbin(positiverate[t], N)
}

for (k in 1:K){
	for (t in 1:T){
		
		Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
	}
}

#priors
theta0 ~ dnorm(0, 1/0.5);
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

dcoptions("verbose"=F)#mute the output of dlclone






generate.data.replicates <- function(phi = "constant",t = 1:5,NN,K = 3){
  
  #function creates all data for simulation
  list.return <- list()
  for(i in 1:NN){
    gen.data <- generate.dataset(phi = phi,t = t)
    
    if(sum(gen.data$Y[2:K,]==1000| gen.data$Y[2:K,]==0)>=1){
      print(paste("bad data on i =", i))
      gen.data <- generate.dataset(phi = phi,t = t,K = K)
    }
  
    list.return[[i]] <- gen.data
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


#narrow priors slightly


set.seed(1292374)
data.const <- generate.data.replicates(phi = "constant", NN = 500)
data.linear  <- generate.data.replicates(phi = "linear", NN = 500)
data.walk  <- generate.data.replicates(phi = "walk", NN = 500)


#check if bad data exists and replace.
data.const <- check.bad.data(data.const,phi = "constant")
data.linear <- check.bad.data(data.linear,phi = 'linear')
data.walk <- check.bad.data(data.walk,phi = 'walk')


#all checks pass



#start parallel
cl <- makePSOCKcluster(10)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")

NN <- 500

print("T = 5 simulations")

run.simulation <- function(cl,data.const, data.linear, data.walk, NN = 500, ti = 5){
  
  error <- matrix(NA,nrow = NN, ncol = 9)
  colnames(error) <- c("const x const","const x linear", "const x walk",
                       "linear x const", "linear x linear","linear x walk",
                       "walk x const", "walk x linear", "walk x walk")
  only.unbiased <-  matrix(NA, nrow = NN, ncol = 3)
  colnames(only.unbiased) <- c("const","linear", "walk")
  
  
  
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
    
    pos.rate.const <- data.const.phi$params[grep(paste0("posrate",ti),names(data.const.phi$params))]
    pos.rate.linear <- data.linear.phi$params[grep(paste0("posrate",ti),names(data.linear.phi$params))]
    pos.rate.walk <- data.walk.phi$params[grep(paste0("posrate",ti),names(data.walk.phi$params))]
    
    unbiased.const.phi <- extract.unbiased(data.const.phi)
    unbiased.linear.phi <- extract.unbiased(data.linear.phi)
    unbiased.walk.phi <- extract.unbiased(data.walk.phi)
    
    #data x model
    
    if(j %% 10 ==0) print(j)
    
    const.const <- jags.parfit(cl, data.const.phi[-8], "positiverate", mod.const.phi,
                               n.chains=4,n.adapt = 25000,thin = 5, n.iter = 60000,inits = chains.init)
    error[j,"const x const"] <- pos.rate.const - get.mean(const.const,ti)
    
    const.linear <- jags.parfit(cl, data.const.phi[-8], "positiverate", mod.linear.phi,
                                n.chains=4,n.adapt = 25000,thin = 5, n.iter = 60000,inits = chains.init)
    error[j, "const x linear"] <- pos.rate.const - get.mean(const.linear,ti)
    
    const.walk <- jags.parfit(cl, data.const.phi[-8], "positiverate", mod.walk.phi,
                              n.chains=4,n.adapt = 25000,thin = 5, n.iter = 60000,inits = chains.init)
    error[j, "const x walk"] <- pos.rate.const - get.mean(const.walk,ti)
    
    linear.const <- jags.parfit(cl, data.linear.phi[-8], "positiverate", mod.const.phi,
                                n.chains=4,n.adapt = 25000,thin = 5, n.iter = 60000,inits = chains.init)
    error[j, "linear x const"] <- pos.rate.linear - get.mean(linear.const,ti)
    
    linear.linear <- jags.parfit(cl, data.linear.phi[-8], "positiverate", mod.linear.phi,
                                 n.chains=4,n.adapt = 25000,thin = 5, n.iter = 60000,inits = chains.init)
    error[j, "linear x linear"] <- pos.rate.linear - get.mean(linear.linear,ti)
    
    linear.walk <- jags.parfit(cl, data.linear.phi[-8], "positiverate", mod.walk.phi,
                               n.chains=4,n.adapt = 25000,thin = 5, n.iter = 60000,inits = chains.init)
    error[j, "linear x walk"] <- pos.rate.linear - get.mean(linear.walk,ti)
    
    walk.const <- jags.parfit(cl, data.walk.phi[-8], "positiverate", mod.const.phi,
                              n.chains=4,n.adapt = 25000,thin = 5, n.iter = 60000,inits = chains.init)
    error[j, "walk x const"] <- pos.rate.walk - get.mean(walk.const,ti)
    
    walk.linear <- jags.parfit(cl, data.walk.phi[-8], "positiverate", mod.linear.phi,
                               n.chains=4,n.adapt = 25000,thin = 5, n.iter = 60000,inits = chains.init)
    error[j, "walk x linear"] <- pos.rate.walk - get.mean(walk.linear,ti)
    
    walk.walk <- jags.parfit(cl, data.walk.phi[-8], "positiverate", mod.walk.phi,
                             n.chains=4,n.adapt = 25000,thin = 5, n.iter = 60000,inits = chains.init)
    error[j, "walk x walk"] <- pos.rate.walk - get.mean(walk.walk,ti)
    
    #one unbiased survey only.
    #data x model
  
    
    unb.const<- jags.parfit(cl, unbiased.const.phi, "positiverate", mod.const.phi,
                            n.chains=4,n.adapt =20000,thin = 5, n.iter = 50000,inits = chains.init)
    only.unbiased[j, "const"] <- pos.rate.const - get.mean(unb.const,ti)
    
    unb.linear <- jags.parfit(cl, unbiased.linear.phi, "positiverate", mod.const.phi,
                              n.chains=4,n.adapt = 20000,thin = 5, n.iter = 50000,inits = chains.init)
    only.unbiased[j, "linear"] <- pos.rate.linear - get.mean(unb.linear,ti)
    
    unb.walk <- jags.parfit(cl, unbiased.walk.phi, "positiverate", mod.const.phi,
                            n.chains=4,n.adapt = 20000,thin = 5, n.iter = 50000,inits = chains.init)
    only.unbiased[j, "walk"] <- pos.rate.walk - get.mean(unb.walk,ti)
    
  }
  
  
  result.RMSE <- apply(error,2,function(x){sqrt(mean((100*x)^2))})
  result.RMSE.unb <- apply(only.unbiased,2,function(x){sqrt(mean((100*x)^2))})

  
  
  results.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3),c()), model = rep(c("const","linear","walk"),3),RMSE = result.RMSE)
  results.plot.unb <- data.frame(data = c("const","linear","walk"), model = rep("unbiased",3),RMSE = result.RMSE.unb)
  
  results.plot.final <- rbind(results.plot,results.plot.unb)
  
  return(results.plot.final)
  
  
}

#run the simulation

results.plot.final  <- run.simulation(cl,data.const,data.linear,data.walk,NN = 500)

ggplot(data = results.plot.final, aes(x = data, y = RMSE,group = model,colour = model)) + geom_point() + geom_line() + theme_minimal() + 
  labs(x = "Data generation", y = "Root Mean Squared Error", title = paste("RMSE of NN = ",NN,"in 3x3 design, 5 timepoints")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black"))
 


#end parallel
stopCluster(cl)


write.csv(results.plot.final,"RMSE_plot_t5_500.csv",row.names = F)


