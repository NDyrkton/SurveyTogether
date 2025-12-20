
#Survey Together Proof of concept
# Nathaniel Dyrkton  Supervised by Paul Gustafson and Harlan Campbell
library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)
library(dclone)
library(dplyr)
library(RColorBrewer)


#T = 10, 10 timpoints

#First step is to generate data under 3 conditions for phi (bias term)
#1: Phi is constant by  time, 2: Phi is linear in time, and 3: Phi follows a random walk 


inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/(1-x))
}


get.point.est <- function(line,var,type = "median"){
  
  
  point.est <- summary(line)$statistics
  
  #breaks if number of time points is 1, checking null fix.
  
  if(is.null(dim(point.est))){
    
    return(point.est[1])
    
  }else{
    if(type == "mean"){
      means <- summary(line)$statistics[,1]
      return(means[grep(var,names(means))])
      
    }else if(type == "median"){
      
      medians <- summary(line)$quantiles[,3]
      return(medians[grep(var,names(medians))])
      
      
    }
    
  }
  
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
    
    return(list(Lower = lower.quantile,Upper = upper.quantile))
    
  }
  
  
}

extract.surveys <- function(datalist,row = 1){
  #function extracts given surveys from data list
  #same as extract.nona but does not shorten
  
  K <-  length(row)
  Y <- datalist$Y[row,]
  
  smalln <- datalist$smalln[row,]
  times <- datalist$times[row,]
  ti <- datalist$ti
  
  new.list <- list(K = K, 
                   times = matrix(times,ncol = ti), N = datalist$N, ti = ti,
                   Y = matrix(Y,ncol = ti), smalln = matrix(smalln,ncol = ti))
  return(new.list)
  
}


#now fixed to be consistent with notation in paper.
generate.dataset <- function(N= 50000, K =3, ti = c(1:5), phi = "constant"){
  Y <- matrix(NA,ncol = length(ti),nrow = K)
  smalln <- matrix(0,ncol = length(ti),nrow = K)
  #get default smalln
  
  for(k in 1:K){
    smalln[k,] <- rep(100,length(ti))
    
    if(k>1){
      smalln[k,] <- rep(1000,length(ti))
    }
  }
  
  posrate_t <- numeric(length(ti))
  theta_t <- numeric(length(ti))
  times <- t(matrix(rep(ti,K),ncol = K))
  
  #priors on general parameters
  sigmasq <- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = sqrt(1))
  theta0 <- rnorm(1,mean =0, sd = sqrt(2))
  
  if(phi == "constant"){
    #generate param based on prior
    gamma0 <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1))) 
    #constant phi values
    phi <- exp(gamma0)
    
    
    theta_t[1] <- rnorm(1,mean = theta0,sd = sqrt(sigmasq))
    posrate_t[1] <- inv.logit(theta_t[1])
    
    for(i in 2:length(ti)){
      
      theta_t[i] <- rnorm(1,mean = theta_t[i-1],sd = sqrt(sigmasq))
      posrate_t[i] <- inv.logit(theta_t[i])
      
    }
    P_t <- rbinom(n = length(ti),size = N, prob = posrate_t)
    
    for(k in 1:K){
      
      for(i in 1:length(ti)){
        Y_kt <- rnoncenhypergeom(n = 1, n1 = P_t[i],n2 = N-P_t[i], m1 = smalln[k,i], psi = phi[k])
        Y[k,i] <- Y_kt
        
      }
    }
    
    parameters <-  c(gamma0,sigmasq,posrate_t, theta0)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = "")
    posrate.names <- paste(rep("positiverate",length(ti)),1:length(ti),sep = "")
    names(parameters) <-  c(gamma0.names,"sigmasq",posrate.names,"theta0")
    
    return(list(K=K, ti =max(ti), times=times, N=N, Y=Y, smalln=smalln, params = parameters,times_centered = times-median(times)))
    
  }else if(phi == "linear"){
    phi_kt <- matrix(NA,nrow =K,ncol = length(ti))
    #priors
    gamma_0k <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1)))
    gamma_1k <- c(0,rnorm(K-1,mean = 0, sd = rep(sqrt(0.05),K-1)))
    
    theta_t[1] <- rnorm(1,mean = theta0,sd = sqrt(sigmasq))
    posrate_t[1] <- inv.logit(theta_t[1])
    
    for(i in 2:length(ti)){
      theta_t[i] <- rnorm(1,mean = theta_t[i-1],sd = sqrt(sigmasq))
      posrate_t[i] <- inv.logit(theta_t[i])
      
    }
    
    P_t <- rbinom(n = length(ti),size = N, prob = posrate_t)
    
    for(k in 1:K){
      for(i in 1:length(ti)){
        
        phi_kt[k,i] <- exp(gamma_0k[k] + gamma_1k[k]*ti[i])
        
        Y_kt <- rnoncenhypergeom(n = 1, n1 = P_t[i],n2 = N-P_t[i], m1 = smalln[k,i], psi = phi_kt[k,i])
        Y[k,i] <- Y_kt
      }
      
    }  
    
    parameters <-  c(gamma_0k,gamma_1k,sigmasq,posrate_t, theta0)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = "")
    gamma1.names <- paste(rep("gamma1",K),1:K,sep = "")
    posrate.names <- paste(rep("positiverate",length(ti)),1:length(ti),sep = "")
    names(parameters) <-  c(gamma0.names,gamma1.names,"sigmasq",posrate.names,"theta0")
    
    return(list(K=K, ti =max(ti), times=times, N=N, Y=Y, smalln=smalln,params = parameters, times_centered = times-median(times)))
    
  }else if(phi == "walk"){
    
    phi_kt <- matrix(NA,nrow =K,ncol = length(ti))
    gamma_kt <- matrix(NA,nrow = K,ncol = length(ti))
    
    gamma_0k <- c(0,rnorm(K-1,mean = 0, sd = rep(1,K-1)))
    gamma_kt[1,] <- rep(0,length(ti))
    #prior
    
    pisq <- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = sqrt(1))
    
    
    
    theta_t[1] <- rnorm(1,mean = theta0,sd = sqrt(sigmasq))
    posrate_t[1] <- inv.logit(theta_t[1])
    #first study is unbiased 
    
    phi_kt[1,] <- exp(gamma_kt[1,])
    gamma_kt[2:K,1] <- rnorm(K-1,mean = gamma_0k[2:K],sd = sqrt(pisq))
    phi_kt[2:K,1] <- exp(gamma_kt[2:K,1])
    
    for(i in 2:length(ti)){
      
      gamma_kt[2:K,i] <- rnorm(K-1, mean = gamma_kt[2:K,i-1],sd = sqrt(c(pisq,pisq)))
      phi_kt[2:K,i] <- exp(gamma_kt[2:K,i])
      
      
      theta_t[i] <- rnorm(1,mean = theta_t[i-1],sd = sqrt(sigmasq))
      posrate_t[i] <- inv.logit(theta_t[i])
      
    }
    
    P_t <- rbinom(n = length(ti),size = N, prob = posrate_t)
    
    for(k in 1:K){
      
      for(i in 1:length(ti)){
        
        Y_kt <- rnoncenhypergeom(n = 1, n1 = P_t[i],n2 = N-P_t[i], m1 = smalln[k,i], psi = phi_kt[k,i])
        Y[k,i] <- Y_kt
      }
    }
    parameters <- c(gamma_0k,as.numeric(gamma_kt),sigmasq,pisq,posrate_t,theta0)
    gamma0.names <- paste(rep("gamma0",K),1:K,sep = "")
    gamma_kt.c <- expand.grid(1:K,1:length(ti))
    gamma_kt.names <- paste(rep("gamma_",K*length(ti)),as.character(gamma_kt.c$Var1),as.character(gamma_kt.c$Var2),sep = '')
    posrate.names <- paste(rep("positiverate",length(ti)),1:length(ti),sep = "")
    names(parameters) <-  c(gamma0.names,gamma_kt.names,"sigmasq","pisq",posrate.names,"theta0")
    
    
    return(list(K=K, ti = max(ti), times=times, N=N, Y=Y, smalln=smalln, params = parameters, times_centered = times-median(times)))
  }
  
}

mod.base.walk  <- '
model{	
#likelihood

for (i in 1:ti){
		phi[1,i] <- 1	
}


for (k in 2:K){

  gamma[k,1] ~ dnorm(gamma0[k],1/pisq)
  phi[k,1] <- exp(gamma[k,1])
  
	for (t in 2:ti){
	  gamma[k,t] ~ dnorm(gamma[k,t-1], 1/pisq)
	  phi[k,t] <- exp(gamma[k,t])
	  
	}
}
	
	
theta[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(theta[1])
for(t in 2:ti){
	theta[t] ~ dnorm(theta[t-1],1/sigmasq)
	positiverate[t]	<- ilogit(theta[t])
}

#for(t in 1:ti){
	#P[t] ~ dbin(positiverate[t], N)
#}

for (k in 1:K){
	for (t in 1:ti){
		
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1/2);
sigmasq ~ dnorm(0, 1)T(0,);
pisq ~ dnorm(0, 1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1/10);

}

}'

mod.base.walk.wide.prior <- '
model{	
#likelihood

for (i in 1:ti){
		phi[1,i] <- 1	
}


for (k in 2:K){

  gamma[k,1] ~ dnorm(gamma0[k],1/pisq)
  phi[k,1] <- exp(gamma[k,1])
  
	for (t in 2:ti){
	  gamma[k,t] ~ dnorm(gamma[k,t-1], 1/pisq)
	  phi[k,t] <- exp(gamma[k,t])
	  
	}
}
	
	
theta[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(theta[1])
for(t in 2:ti){
	theta[t] ~ dnorm(theta[t-1],1/sigmasq);
	positiverate[t]	<- ilogit(theta[t])
}


for (k in 1:K){
	for (t in 1:ti){
		
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1/100);
sigmasq ~ dnorm(0, 1/100)T(0,);
pisq ~ dnorm(0, 1/100)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1/100);

}

}'

mod.linear.phi <-'
model{	
#likelihood

for (i in 1:ti){
		phi[1,i] <- 1	
}


for (k in 2:K){
	for (t in 1:ti){
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times[k,t])
	}
}
	
	
theta[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(theta[1])



for(t in 2:ti){

	theta[t] ~ dnorm(theta[t-1],1/sigmasq)
	#theta[t] ~ dunif(theta[t-1],theta[t-1]+rho)
	positiverate[t]	<- ilogit(theta[t])
}

#for(t in 1:ti){
          #normal approximation
# P[t] ~ dbin(positiverate[t], N)
#}


for (k in 1:K){
	for (t in 1:ti){

		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1/2);
sigmasq ~ dnorm(0, 1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1/10);
	gamma1[k] ~ dnorm(0, 1);
}
}'

mod.linear.phi.wide <-'
model{	
#likelihood

for (i in 1:ti){
		phi[1,i] <- 1	
}


for (k in 2:K){
	for (t in 1:ti){
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times[k,t])
	}
}
	
	
theta[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(theta[1])



for(t in 2:ti){

	theta[t] ~ dnorm(theta[t-1],1/sigmasq)
	positiverate[t]	<- ilogit(theta[t])
}



for (k in 1:K){
	for (t in 1:ti){
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1/100);
sigmasq ~ dnorm(0, 1/100)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1/100);
	gamma1[k] ~ dnorm(0, 1/100);
}
}'


mod.base.wide.prior <- '
model{	
#likelihood

for (i in 1:ti){
		phi[1,i] <- 1	
}


for (k in 2:K){

  gamma[k,1] ~ dnorm(gamma0[k],1/pisq)
  phi[k,1] <- exp(gamma[k,1])
  
	for (t in 2:ti){
	  gamma[k,t] ~ dnorm(gamma[k,t-1], 1/pisq)
	  phi[k,t] <- exp(gamma[k,t])
	  
	}
}
	
	
theta[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(theta[1])
for(t in 2:ti){
	theta[t] ~ dnorm(theta[t-1],1/sigmasq);
	positiverate[t]	<- ilogit(theta[t])
}


for (k in 1:K){
	for (t in 1:ti){
		
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1/100);
sigmasq ~ dnorm(0, 1/100)T(0,);
pisq ~ dnorm(0, 1/100)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1/100);

}

}'




n.chains <- 8
cl <- makePSOCKcluster(n.chains)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


set.seed(31525)

#generate dataset
data <- generate.dataset(N=100000,K = 3,t = c(1:10), phi = "linear")


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
chain5<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(149))
chain6<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(230))
chain7<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(321))
chain8<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(412))


inits.chains <- list(chain1,chain2,chain3,chain4,chain5,chain6,chain7,chain8)

line.1 <- jags.parfit(cl, K_1, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=n.chains,n.adapt = 5000, n.update = 10000, n.iter = 50000
                      ,inits = inits.chains)
line.2 <- jags.parfit(cl, K_2, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=n.chains,n.adapt = 5000, n.update = 10000, n.iter = 50000
                      ,inits = inits.chains)
line.3 <- jags.parfit(cl, K_3, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=n.chains,n.adapt = 5000, n.update = 10000, n.iter = 50000
                      ,inits = inits.chains)

line.linear.base <- jags.parfit(cl, data, c("positiverate","gamma0",'gamma1'), custommodel(mod.linear.phi),
                         n.chains=n.chains,n.adapt = 5000, n.update = 10000, n.iter = 500000
                         ,inits = inits.chains)

line.linear.wide <- jags.parfit(cl, data, c("positiverate","gamma0",'gamma1'), custommodel(mod.linear.phi.wide),
                                n.chains=n.chains,n.adapt = 5000, n.update = 10000, n.iter = 500000
                                ,inits = inits.chains)

line.walk.base <- jags.parfit(cl, data, c("positiverate",'gamma'), custommodel(mod.base.walk),
                                n.chains=n.chains,n.adapt = 5000, n.update = 10000, n.iter = 500000
                                ,inits = inits.chains)
line.walk.wide <- jags.parfit(cl, data, c("positiverate",'gamma'), custommodel(mod.base.walk.wide.prior),
                                n.chains=n.chains,n.adapt = 5000, n.update = 10000,thin = 5, n.iter = 500000
                                ,inits = inits.chains)

stopCluster(cl)


posrate.1 <- get.point.est(line.1,'positiverate')
CI.1 <- get.CI(line.1, 'positiverate')


posrate.2 <- get.point.est(line.2,'positiverate')
CI.2 <- get.CI(line.2, 'positiverate')


posrate.3 <- get.point.est(line.3,'positiverate')
CI.3 <- get.CI(line.3, 'positiverate')


posrate.linear.base <- get.point.est(line.linear.base,'positiverate')
CI.linear.base <- get.CI(line.linear.base, 'positiverate')

posrate.linear.wide <- get.point.est(line.linear.wide,'positiverate')
CI.linear.wide <- get.CI(line.linear.wide, 'positiverate')

posrate.walk.base <- get.point.est(line.walk.base,'positiverate')
CI.walk.base <- get.CI(line.walk.base, 'positiverate')


posrate.walk.wide <- get.point.est(line.walk.wide,'positiverate')
CI.walk.wide <- get.CI(line.walk.wide, 'positiverate')




estimates.frame <- data.frame(Survey = c(rep("Survey 1",10),rep("Survey 2",10),rep("Survey 3",10),rep("Linear Base",10), rep("Linear Wide",10),rep("Walk Base",10), rep("Walk Wide",10) ),
                              time = c(1:10,1:10,1:10,1:10,1:10,1:10,1:10), 
                    estimate = c(posrate.1,posrate.2,posrate.3,posrate.linear.base,posrate.linear.wide,posrate.walk.base,posrate.walk.wide),CI.L = c(CI.1$Lower,CI.2$Lower,CI.3$Lower,CI.linear.base$Lower ,
                                                                                                                                                     CI.linear.wide$Lower,CI.walk.base$Lower,CI.walk.wide$Lower),
                    CI.U = c(CI.1$Upper,CI.2$Upper,CI.3$Upper,CI.linear.base$Upper,CI.linear.wide$Upper,CI.walk.base$Upper,CI.walk.wide$Upper))

frame.posrate <-  data.frame(time = c(1:10),posrate = data$params[grep("positiverate",names(data$params))])


sensitivityPlots <- ggplot(estimates.frame %>% filter(!(Survey %in% c("Survey 2","Survey 3"))),aes(x = time, y = estimate, colour = Survey)) + geom_point() + geom_line(linewidth= 0.75) + 
  geom_ribbon(aes(ymin = CI.L,ymax = CI.U),alpha = 0.12) + theme_bw() + 
  labs(x = "Time", y = "Positive Rate",title = "Sensitivity of analysis to wider priors") + geom_point(data=frame.posrate,aes(x = time, y = posrate,colour = "True Positive Rate")) +
  geom_line(data=frame.posrate,aes(x = time, y = posrate,colour = "True Positive Rate"),linewidth = 1.20) +   scale_x_continuous(breaks = seq(1, 10, by = 1))  +
  scale_colour_manual(name = "Method Type",values = c("Survey 1" ="royalblue","Linear Base"="magenta" ,"Linear Wide"= "red3", "Walk Base"="orange", "Walk Wide"="limegreen" ,"True Positive Rate"="black"))   


ggsave("Figures/sensitivityPriors.png",plot = sensitivityPlots,width = 26, height = 14, units = 'cm')


#save the plot

#total width reduction
mean((CI.1$Upper-CI.1$Lower)/(CI.linear.base$Upper-CI.linear.base$Lower))
mean((CI.1$Upper-CI.1$Lower)/(CI.linear.wide$Upper-CI.linear.wide$Lower))
mean((CI.1$Upper-CI.1$Lower)/(CI.walk.base$Upper-CI.walk.base$Lower))
mean((CI.1$Upper-CI.1$Lower)/(CI.walk.wide$Upper-CI.walk.wide$Lower))







#generate dataset
set.seed(31525)


n.chains <- 8
cl <- makePSOCKcluster(n.chains)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


data <- generate.dataset(N=100000,K = 3,t = c(1:10), phi = "linear", n_unbiased = 10,n_biased = 100)


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
chain5<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(149))
chain6<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(230))
chain7<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(321))
chain8<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(412))




inits.chains <- list(chain1,chain2,chain3,chain4,chain5,chain6,chain7,chain8)

line.1 <- jags.parfit(cl, K_1, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=n.chains,n.adapt = 50000,thin = 5, n.iter = 50000
                      ,inits = inits.chains)
line.2 <- jags.parfit(cl, K_2, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=n.chains,n.adapt = 50000,thin = 5, n.iter = 50000
                      ,inits = inits.chains)
line.3 <- jags.parfit(cl, K_3, c("positiverate"), custommodel(mod.linear.phi),
                      n.chains=n.chains,n.adapt = 50000,thin = 5, n.iter = 50000
                      ,inits = inits.chains)

line.linear.base <- jags.parfit(cl, data, c("positiverate","gamma0",'gamma1'), custommodel(mod.linear.phi),
                                n.chains=n.chains,n.adapt = 200000,thin = 5, n.iter = 500000
                                ,inits = inits.chains)

line.linear.wide <- jags.parfit(cl, data, c("positiverate","gamma0",'gamma1'), custommodel(mod.linear.phi.wide),
                                n.chains=n.chains,n.adapt = 200000,thin = 5, n.iter = 500000
                                ,inits = inits.chains)

line.walk.base <- jags.parfit(cl, data, c("positiverate",'gamma'), custommodel(mod.base.walk),
                              n.chains=n.chains,n.adapt = 200000,thin = 5, n.iter = 500000
                              ,inits = inits.chains)
line.walk.wide <- jags.parfit(cl, data, c("positiverate",'gamma'), custommodel(mod.base.walk.wide.prior),
                              n.chains=n.chains,n.adapt = 200000,thin = 5, n.iter = 500000
                              ,inits = inits.chains)

stopCluster(cl)


posrate.1 <- get.point.est(line.1,'positiverate')
CI.1 <- get.CI(line.1, 'positiverate')


posrate.2 <- get.point.est(line.2,'positiverate')
CI.2 <- get.CI(line.2, 'positiverate')


posrate.3 <- get.point.est(line.3,'positiverate')
CI.3 <- get.CI(line.3, 'positiverate')


posrate.linear.base <- get.point.est(line.linear.base,'positiverate')
CI.linear.base <- get.CI(line.linear.base, 'positiverate')

posrate.linear.wide <- get.point.est(line.linear.wide,'positiverate')
CI.linear.wide <- get.CI(line.linear.wide, 'positiverate')

posrate.walk.base <- get.point.est(line.walk.base,'positiverate')
CI.walk.base <- get.CI(line.walk.base, 'positiverate')


posrate.walk.wide <- get.point.est(line.walk.wide,'positiverate')
CI.walk.wide <- get.CI(line.walk.wide, 'positiverate')




estimates.frame <- data.frame(Survey = c(rep("Survey 1",10),rep("Survey 2",10),rep("Survey 3",10),rep("Linear Base",10), rep("Linear Wide",10),rep("Walk Base",10), rep("Walk Wide",10) ),
                              time = c(1:10,1:10,1:10,1:10,1:10,1:10,1:10), 
                              estimate = c(posrate.1,posrate.2,posrate.3,posrate.linear.base,posrate.linear.wide,posrate.walk.base,posrate.walk.wide),CI.L = c(CI.1$Lower,CI.2$Lower,CI.3$Lower,CI.linear.base$Lower ,
                                                                                                                                                               CI.linear.wide$Lower,CI.walk.base$Lower,CI.walk.wide$Lower),
                              CI.U = c(CI.1$Upper,CI.2$Upper,CI.3$Upper,CI.linear.base$Upper,CI.linear.wide$Upper,CI.walk.base$Upper,CI.walk.wide$Upper))

frame.posrate <-  data.frame(time = c(1:10),posrate = data$params[grep("positiverate",names(data$params))])


sensitivityPlots <- ggplot(estimates.frame %>% filter(!(Survey %in% c("Survey 2","Survey 3"))),aes(x = time, y = estimate, colour = Survey)) + geom_point() + geom_line(linewidth= 0.75) + 
  geom_ribbon(aes(ymin = CI.L,ymax = CI.U),alpha = 0.12) + theme_bw() + 
  labs(x = "Time", y = "Positive Rate",title = "Sensitivity of analysis to wider priors") + geom_point(data=frame.posrate,aes(x = time, y = posrate,colour = "True Positive Rate")) +
  geom_line(data=frame.posrate,aes(x = time, y = posrate,colour = "True Positive Rate"),linewidth = 1.20) +   scale_x_continuous(breaks = seq(1, 10, by = 1))  +
  scale_colour_manual(name = "Method Type",values = c("Survey 1" ="royalblue","Linear Base"="magenta" ,"Linear Wide"= "red3", "Walk Base"="orange", "Walk Wide"="limegreen" ,"True Positive Rate"="black"))   


#total width reduction
mean((CI.1$Upper-CI.1$Lower)/(CI.linear.base$Upper-CI.linear.base$Lower))
mean((CI.1$Upper-CI.1$Lower)/(CI.linear.wide$Upper-CI.linear.wide$Lower))
mean((CI.1$Upper-CI.1$Lower)/(CI.walk.base$Upper-CI.walk.base$Lower))
mean((CI.1$Upper-CI.1$Lower)/(CI.walk.wide$Upper-CI.walk.wide$Lower))
#mean((CI.1$Upper-CI.1$Lower)/(CI.linear.base$Upper-CI.linear.base$Lower))






#now-cast
(CI.1$Upper[10]-CI.1$Lower[10])/(CI.full$Upper[10]-CI.full$Lower[10])
