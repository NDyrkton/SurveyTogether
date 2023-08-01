
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
sigmasq ~ dnorm(0, 1/0.5)T(0,);

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
sigmasq ~ dnorm(0, 1/0.5)T(0,);

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
sigmasq ~ dnorm(0, 1/0.5)T(0,);
pisq ~ dnorm(0, 1/0.01)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);

}

}')

#helper functions
source("Scripts/helperfunctions.R")




extract.unbiased <- function(datalist){
  K = 1
  T <- datalist$T
  new.list <- list(K = K, 
                   times = matrix(datalist$times[1,],ncol = T), N = datalist$N, T = T,
                   Y = matrix(datalist$Y[1,],ncol = T), smalln = matrix(datalist$smalln[1,],ncol = T))
  
  return(new.list)
}

get.mean <- function(mcmc.obj,timepoint){
  return(summary(mcmc.obj)$statistics[,1][5])
}

dcoptions("verbose"=F)#mute the output of dlclone



#start parallel
cl <- makePSOCKcluster(3)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


NN <- 10
set.seed(12345)

error <- matrix(NA,nrow = NN, ncol = 9)
colnames(error) <- c("const x const","const x linear", "const x walk",
                     "linear x const", "linear x linear","linear x walk",
                     "walk x const", "walk x linear", "walk x walk")

only.unbiased <-  matrix(NA, nrow = NN, ncol = 3)
colnames(only.unbiased) <- c("const","linear", "walk")


time = Sys.time()
for(j in 1:NN){
  
  
  chain1<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(j))
  chain2<- list(.RNG.name = "base::Super-Duper", 
                .RNG.seed = c(j+1))
  chain3<- list(.RNG.name = "base::Wichmann-Hill", 
                .RNG.seed = c(j+2))
  
  chains.init <- list(chain1,chain2,chain3)
  
  
  data.const.phi <- generate.dataset(phi = 'constant',t = 1:5,K = 3)
  data.linear.phi <- generate.dataset(phi = 'linear', t = 1:5,K = 3)
  data.walk.phi <- generate.dataset(phi = 'walk', t= 1:5,K = 3)
  
  pos.rate.const <- data.const.phi$params[grep("posrate5",names(data.const.phi$params))]
  pos.rate.linear <- data.linear.phi$params[grep("posrate5",names(data.linear.phi$params))]
  pos.rate.walk <- data.walk.phi$params[grep("posrate5",names(data.walk.phi$params))]
  
  unbiased.const.phi <- extract.unbiased(data.const.phi)
  unbiased.linear.phi <- extract.unbiased(data.linear.phi)
  unbiased.walk.phi <- extract.unbiased(data.walk.phi)
  
  #data x model
  
  if(j %% 10 ==0) print(j)
  
  const.const <- jags.parfit(cl, data.const.phi[-8], "positiverate", mod.const.phi,
                             n.chains=3,n.adapt = 10000,thin = 5, n.iter = 40000,inits = chains.init)
  error[j,"const x const"] <- pos.rate.const - get.mean(const.const)
  
  const.linear <- jags.parfit(cl, data.const.phi[-8], "positiverate", mod.linear.phi,
                              n.chains=3,n.adapt = 10000,thin = 5, n.iter = 40000,inits = chains.init)
  error[j, "const x linear"] <- pos.rate.const - get.mean(const.linear)
  
  const.walk <- jags.parfit(cl, data.const.phi[-8], "positiverate", mod.walk.phi,
                            n.chains=3,n.adapt = 10000,thin = 5, n.iter = 40000,inits = chains.init)
  error[j, "const x walk"] <- pos.rate.const - get.mean(const.walk)
  
  linear.const <- jags.parfit(cl, data.linear.phi[-8], "positiverate", mod.const.phi,
                              n.chains=3,n.adapt = 10000,thin = 5, n.iter = 40000,inits = chains.init)
  error[j, "linear x const"] <- pos.rate.linear - get.mean(linear.const)
  
  linear.linear <- jags.parfit(cl, data.linear.phi[-8], "positiverate", mod.linear.phi,
                               n.chains=3,n.adapt = 10000,thin = 5, n.iter = 40000,inits = chains.init)
  error[j, "linear x linear"] <- pos.rate.linear - get.mean(linear.linear)
  
  linear.walk <- jags.parfit(cl, data.linear.phi[-8], "positiverate", mod.walk.phi,
                             n.chains=3,n.adapt = 10000,thin = 5, n.iter = 40000,inits = chains.init)
  error[j, "linear x walk"] <- pos.rate.linear - get.mean(linear.walk)
  
  walk.const <- jags.parfit(cl, data.walk.phi[-8], "positiverate", mod.const.phi,
                            n.chains=3,n.adapt = 10000,thin = 5, n.iter = 40000,inits = chains.init)
  error[j, "walk x const"] <- pos.rate.walk - get.mean(walk.const)
  
  walk.linear <- jags.parfit(cl, data.walk.phi[-8], "positiverate", mod.linear.phi,
                             n.chains=3,n.adapt = 10000,thin = 5, n.iter = 40000,inits = chains.init)
  error[j, "walk x linear"] <- pos.rate.walk - get.mean(walk.linear)
  
  walk.walk <- jags.parfit(cl, data.walk.phi[-8], "positiverate", mod.walk.phi,
                           n.chains=3,n.adapt = 10000,thin = 5, n.iter = 40000,inits = chains.init)
  error[j, "walk x walk"] <- pos.rate.walk - get.mean(walk.walk)
  
  #one unbiased survey only.
  #data x model
  
  unb.const<- jags.parfit(cl, unbiased.const.phi, "positiverate", mod.const.phi,
                                 n.chains=3,n.adapt = 5000,thin = 5, n.iter = 30000,inits = chains.init)
  only.unbiased[j, "const"] <- pos.rate.const - get.mean(unb.const)
  
  unb.linear <- jags.parfit(cl, unbiased.linear.phi, "positiverate", mod.const.phi,
                                  n.chains=3,n.adapt = 5000,thin = 5, n.iter = 30000,inits = chains.init)
  only.unbiased[j, "linear"] <- pos.rate.linear - get.mean(unb.linear)
  
  unb.walk <- jags.parfit(cl, unbiased.walk.phi, "positiverate", mod.const.phi,
                                n.chains=3,n.adapt = 5000,thin = 5, n.iter = 30000,inits = chains.init)
  only.unbiased[j, "walk"] <- pos.rate.walk - get.mean(unb.walk)

}

Sys.time()-time

result.RMSE <- apply(error,2,function(x){sqrt(mean((100*x)^2))})
result.RMSE.unb <- apply(only.unbiased,2,function(x){sqrt(mean((100*x)^2))})

#1 hr for 50 iterations


#200 200 hours is 8 days


results.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3),c()), model = rep(c("const","linear","walk"),3),RMSE = result.RMSE)
results.plot.unb <- data.frame(data = c("const","linear","walk"), model = rep("unbiased",3),RMSE = result.RMSE.unb)

results.plot.final <- rbind(results.plot,results.plot.unb)

ggplot(data = results.plot.final, aes(x = data, y = RMSE,group = model,colour = model)) + geom_point() + geom_line() + theme_minimal() + 
  labs(x = "Data generation", y = "Root Mean Squared Error", title = paste("RMSE of NN = ",NN,"in 3x3 design, 5 timepoints")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black"))
 


#end parallel
stopCluster(cl)


write.csv(results.plot,"Results/RMSE_with_bias.csv",row.names = T)
write.csv(results.plot.unb,"Results/RMSE_only_unbiased.csv",row.names = T)
write.csv(error,"Results/RMSE_total.csv",row.names = T)


