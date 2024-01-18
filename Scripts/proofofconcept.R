
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
generate.dataset <- function(N= 10000, K =3, t = c(1:5), ns = rep(100,length(t)), phi = "constant"){
  Y <- matrix(NA,ncol = length(t),nrow = K)

  smalln <- if (is.matrix(ns)) ns  else t(matrix(rep(ns,K),ncol = K));
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
    gamma_1k <- c(0,rnorm(K-1,mean = 0, sd = rep(0.25,K-1)))
    
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
    
    pisq <- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = sqrt(0.1))
    
    
    
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
theta0 ~ dnorm(0, 1);
sigmasq ~ dnorm(0, 1/0.1)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);
	gamma1[k] ~ dnorm(0, 1/0.25);
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


cl <- makePSOCKcluster(4)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


set.seed(31525)

#generate dataset

new.n <- matrix(c(rep(100,10), rep(1000,10), rep(1000,10)), ncol = 10,nrow = 3,byrow = T)

data <- generate.dataset(N=100000,K = 3,t = c(1:10), phi = "linear",ns= new.n)


K_1 <- extract.surveys(data,1)

K_2 <- extract.surveys(data,2)

K_3 <- extract.surveys(data,3)




chain1<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(159))
chain2<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(260))
chain3<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(371))
chain4<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(482))

inits.chains <- list(chain1,chain2,chain3,chain4)

line.1 <- jags.parfit(cl, K_1, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=4,n.adapt = 50000,thin = 5, n.iter = 50000
                      ,inits = inits.chains)
line.2 <- jags.parfit(cl, K_2, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=4,n.adapt = 50000,thin = 5, n.iter = 50000
                      ,inits = inits.chains)
line.3 <- jags.parfit(cl, K_3, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=4,n.adapt = 50000,thin = 5, n.iter = 50000
                      ,inits = inits.chains)

line.full <- jags.parfit(cl, data, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=4,n.adapt = 200000,thin = 5, n.iter = 100000
                      ,inits = inits.chains)

posrate.1 <- get.point.est(line.1,'positiverate')
CI.1 <- get.CI(line.1, 'positiverate')


posrate.2 <- get.point.est(line.2,'positiverate')
CI.2 <- get.CI(line.2, 'positiverate')


posrate.3 <- get.point.est(line.3,'positiverate')
CI.3 <- get.CI(line.3, 'positiverate')


posrate.full <- get.point.est(line.full,'positiverate')
CI.full <- get.CI(line.full, 'positiverate')



frame <- data.frame(Survey = c(rep("(k=1)",10),rep("(k=2)",10),rep("(k=3)",10),rep("Method",10)),time = c(1:10,1:10,1:10,1:10), 
                               estimate = c(posrate.1,posrate.2,posrate.3,posrate.full),CI.L = c(CI.1$Lower,CI.2$Lower,CI.3$Lower,CI.full$Lower),
                               CI.U = c(CI.1$Upper,CI.2$Upper,CI.3$Upper,CI.full$Upper))
                    
frame.posrate <-  data.frame(time = c(1:10),posrate = data$params[grep("posrate",names(data$params))])

library(RColorBrewer)
ggplot(frame,aes(x = time, y = estimate, colour = Survey)) + geom_point() + geom_line() + 
  geom_ribbon(aes(ymin = CI.L,ymax = CI.U),alpha = 0.2) + theme_bw() + 
  labs(x = "Time", y = "Positive rate",title = "Example of the synthesis method on simulated data") + geom_point(data=frame.posrate,aes(x = time, y = posrate,colour = "True Positive rate")) +
  geom_line(data=frame.posrate,aes(x = time, y = posrate,colour = "True Positive rate"),linetype = "dashed") + scale_colour_brewer(palette = "Dark2")
  


mean((CI.1$Upper-CI.1$Lower)/(CI.full$Upper-CI.full$Lower))
