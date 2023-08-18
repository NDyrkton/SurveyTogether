#Compare Confidence Intervals using unbiased
library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)
library(dclone)
library(dplyr)
source("Scripts/helperfunctions.R")
data.list <- readRDS("ddc_list.Rdata")


extract.unbiased.nona <- function(datalist,col = 1){
  K = 1
  Y <- datalist$Y[col,]

  smalln <- datalist$smalln[col,]
  smalln <- smalln[!is.na(Y)]
  
  Y <- Y[!is.na(Y)]
  T <- length(Y)
  times <- 1:T
  new.list <- list(K = K, 
                   times = matrix(times,ncol = T), N = datalist$N, T = T,
                   Y = matrix(Y,ncol = T), smalln = matrix(smalln,ncol = T))
  
  return(new.list)
}


get.point.est <- function(line,var){
  
  point.est <- summary(line)$statistics[,1]
  #variable of interest
  return(point.est[grep(var,names(point.est))])
}

get.CI <- function(line,var){
  
  lower.quantile <- summary(line)$quantile[,1]
  upper.quantile <- summary(line)$quantile[,5]
  
  #extract variable of interest
  lower.quantile <- lower.quantile[grep(var,names(lower.quantile))] 
  upper.quantile <- upper.quantile[grep(var,names(upper.quantile))] 
  
  return(list(Lower = lower.quantile,Upper = upper.quantile))
}


full.data.joined <- read.csv("Data/full_data_joined_ddc_vaccine.csv")
benchmark.data <- read.csv("Data/Benchmark_clean.csv")



mod.const.phi <- custommodel('
model{	
#likelihood

for (i in 1:T){
		phi[1,i] <- 1	
}


for (k in 2:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(logitpositiverate[1])



for(t in 2:T){

	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],1/sigmasq)T(logitpositiverate[t-1],)
	#logitpositiverate[t] ~ dunif(logitpositiverate[t-1],logitpositiverate[t-1]+rho)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

#for(t in 1:T){
          #normal approximation
# P[t] ~ dbin(positiverate[t], N)
#}


for (k in 1:K){
	for (t in 1:T){
    #		Y[k,t] ~ dbin(1-(1-(P[t]/N))^phi[k,t],smalln[k,t])
		Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
		#Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
	}
}

#priors
theta0 ~ dnorm(-2, 1);
sigmasq ~ dnorm(0, 1/5)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
}
}')


mod.linear.phi <- custommodel('
model{	
#likelihood

for (i in 1:T){
		phi[1,i] <- 1	
}


for (k in 2:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times[k,t])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(logitpositiverate[1])



for(t in 2:T){

	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],1/sigmasq)T(logitpositiverate[t-1],)
	#logitpositiverate[t] ~ dunif(logitpositiverate[t-1],logitpositiverate[t-1]+rho)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

#for(t in 1:T){
          #normal approximation
# P[t] ~ dbin(positiverate[t], N)
#}


for (k in 1:K){
	for (t in 1:T){
    #		Y[k,t] ~ dbin(1-(1-(P[t]/N))^phi[k,t],smalln[k,t])
		Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
		#Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
	}
}

#priors
theta0 ~ dnorm(-2, 1);
sigmasq ~ dnorm(0, 1/5)T(0,);

for (k in 2:K){
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
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],1/sigmasq)T(logitpositiverate[t-1],);
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

#for(t in 1:T){
	#P[t] ~ dbin(positiverate[t], N)
#}

for (k in 1:K){
	for (t in 1:T){
		
		Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
sigmasq ~ dnorm(0, 1/5)T(0,);
pisq ~ dnorm(0, 1/5)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);

}

}')



mod.linear.phi.2 <- custommodel('
model{	
#likelihood

for (i in 1:T){
		phi[1,i] <- 1	
}


for (k in 2:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k] + gamma1*times[k,t])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(logitpositiverate[1])



for(t in 2:T){

	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],1/sigmasq)T(logitpositiverate[t-1],)
	#logitpositiverate[t] ~ dunif(logitpositiverate[t-1],logitpositiverate[t-1]+rho)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

#for(t in 1:T){
          #normal approximation
# P[t] ~ dbin(positiverate[t], N)
#}


for (k in 1:K){
	for (t in 1:T){
    #		Y[k,t] ~ dbin(1-(1-(P[t]/N))^phi[k,t],smalln[k,t])
		Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
		#Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
	}
}

#priors
theta0 ~ dnorm(-2, 1);
sigmasq ~ dnorm(0, 1/5)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
}
	gamma1 ~ dnorm(0, 1/0.25);
}')




mod.walk.phi.2 <- custommodel('
model{	
#likelihood

for (i in 1:T){
		phi[1,i] <- 1	
		logodds[1,i] <- 1
}


gamma[1] <- 0
  
for(k in 2:K){
  logodds[k,1] <- gamma[1]+gamma0[k]
  phi[k,1] <- exp(gamma[1]+gamma0[k])
}
  
for(t in 2:T){

  gamma[t] ~ dnorm(gamma[t-1],1/pisq)
  
  for(k in 2:K){
    logodds[k,t] <- gamma[t] + logodds[k,t-1]
    phi[k,t] <- exp(logodds[k,t])
  
  }
  
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],1/sigmasq)T(logitpositiverate[t-1],);
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

#for(t in 1:T){
	#P[t] ~ dbin(positiverate[t], N)
#}

for (k in 1:K){
	for (t in 1:T){
		
		Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
sigmasq ~ dnorm(0, 1/5)T(0,);
pisq ~ dnorm(0, 1/5)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);

}



}')


#get unbiased data with no NAs
ipsos.dat <- extract.unbiased.nona(data.list,col = 1)
household.dat <- extract.unbiased.nona(data.list,col= 2)
facebook.dat <- extract.unbiased.nona(data.list,col = 3)


cl <- makePSOCKcluster(4)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")

chain1<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(159))
chain2<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(260))
chain3<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(371))
chain4<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(482))


#now run linear phi for each of the unbiased surveys


line.ipsos <- jags.parfit(cl, ipsos.dat, c("positiverate","gamma0","gamma1","sigmasq"), mod.linear.phi,
                           n.chains=4,n.adapt = 200000,thin = 5, n.iter = 50000,inits = list(chain1,chain2,chain3,chain4))

#check convergence
gelman.diag(line.ipsos)

point.ipsos <- get.point.est(line.ipsos,"positiverate")
CI.ipsos <- get.CI(line.ipsos,"positiverate")


line.household <- jags.parfit(cl, household.dat, c("positiverate","gamma0","gamma1","sigmasq"), mod.linear.phi,
                             n.chains=4,n.adapt = 200000,thin = 5, n.iter = 50000,inits = list(chain1,chain2,chain3,chain4))

gelman.diag(line.household)

point.household <- get.point.est(line.household,"positiverate")
CI.household <- get.CI(line.household,"positiverate")


line.facebook <- jags.parfit(cl, facebook.dat, c("positiverate","gamma0","gamma1","sigmasq"), mod.linear.phi,
                             n.chains=4,n.adapt = 200000,thin = 5, n.iter = 50000,inits = list(chain1,chain2,chain3,chain4))


gelman.diag(line.facebook)
point.facebook <- get.point.est(line.facebook,"positiverate")
CI.facebook <-  get.CI(line.facebook,"positiverate")



#run actual method (linear)

line.full <- jags.parfit(cl, data.list, c("positiverate","phi","sigmasq"), mod.linear.phi.2,
                         n.chains=4,n.adapt = 1000000,thin = 5, n.iter = 700000,inits = list(chain1,chain2,chain3,chain4))


#gelman.diag(line.full) 

point.full <- get.point.est(line.full,"positiverate")
CI.full <- get.CI(line.full,"positiverate")



#View(full.data.joined)
#in order facebook->household -> axios
full.data.joined$unb.lower.CI <- c(CI.facebook$Lower,CI.household$Lower,CI.ipsos$Lower)
full.data.joined$unb.upper.CI <- c(CI.facebook$Upper,CI.household$Upper,CI.ipsos$Upper)
full.data.joined$pointest <- c(point.facebook,point.household,point.ipsos)

facebook.end_date <- full.data.joined[full.data.joined$mode=="facebook",'end_date']


data.method <- data.frame(CI.Lower =CI.full$Lower,CI.Upper = CI.full$Upper,posrate = point.full,end_date = as.Date(facebook.end_date))





full.data.joined$end_date <- as.Date(full.data.joined$end_date)

#generate full plot
ggplot(data = full.data.joined,aes(x = end_date,y = posrate,col = mode,group = mode))+
  geom_point()+geom_line() +geom_errorbar(aes(ymin = CI.L,ymax = CI.U),col = 'black',lwd = 0.5)+
  geom_errorbar(aes(ymin=unb.lower.CI,ymax = unb.upper.CI)) + 
  geom_point(data = data.method,aes(x = end_date,y = posrate,group = "method",col = "method"))+
  geom_line(data = data.method,aes(x = end_date,y = posrate,group = "method",col = "method")) + 
  geom_errorbar(data = data.method,aes(ymin = CI.Lower,ymax = CI.Upper,group = "method",col = "method"),width = 1.75)+
  theme_minimal() + labs(x = "Date",y = "% Vaccinated",title = expression(paste("Plot of Method for linear ",phi))) +
  geom_ribbon(data = benchmark.data,aes(x = as.Date(date),y = posrate,ymin = posrate*0.98,ymax = posrate*1.02,group = "benchmark",col = "benchmark"),alpha = 0.2)





#try again with random walk for phi
line.full.walk <- jags.parfit(cl, data.list, c("positiverate","phi",'sigmasq'), mod.walk.phi.2,
                         n.chains=4,n.adapt = 700000,thin = 10, n.iter = 500000,inits = list(chain1,chain2,chain3,chain4))

#gelman.diag(line.full.walk)

point.walk <- get.point.est(line.full.walk,'positiverate')
CI.walk <- get.CI(line.full.walk,'positiverate')


#compare with constant phi

line.const <- jags.parfit(cl, data.list, c("positiverate","gamma0",'sigmasq'), mod.const.phi,
                          n.chains=4,n.adapt = 500000,thin = 10, n.iter = 400000,list(chain1,chain2,chain3,chain4))

gelman.diag(line.const)

point.const <- get.point.est(line.const,'positiverate')
CI.const <- get.CI(line.const,'positiverate')

#CI gain for constant phi vs walk phi
(CI.const$Upper-CI.const$Lower)/(CI.walk$Upper-CI.walk$Lower)


axios.ipsos <- full.data.joined %>% filter(mode == "ipsos_axios")

data.linear.walk <- data.frame(time = c(facebook.end_date,facebook.end_date,facebook.end_date,axios.ipsos$end_date),posrate = c(point.const,point.full,point.walk,axios.ipsos$posrate), CI.L=c(CI.const$Lower,CI.full$Lower,CI.walk$Lower,axios.ipsos$CI.L),
                               CI.U = c(CI.const$Upper,CI.full$Upper,CI.walk$Upper,axios.ipsos$CI.U),method = c(rep("Constant",20),rep("Linear",20),rep("Walk",20),rep("ipsos axios",11)))

data.linear.method <- data.linear.walk %>% filter(method != 'ipsos axios')

ggplot(data = data.linear.method,aes(x = as.Date(time), y = posrate,group = method,colour = method))+
  geom_point() + geom_line() + geom_errorbar(aes(ymin = CI.L,ymax = CI.U),width = 0.25) + theme_minimal()+
  labs(x= "Date",y ="% Vaccinated", title = expression(paste("Comparison of different models for ",phi)))

#plot gamma bias
#get gammas

const.gamma <- get.point.est(line.const,"gamma0")
const.CI <- get.CI(line.const,"gamma0")

point.linear.phi <- get.point.est(line.full,"phi")

linear.CI <- get.CI(line.full,"phi")
linear.CI.23.L<- linear.CI$Lower[-seq(1,60,by =3)]

linear.phi.CI.L.2 <-linear.CI.23.L[seq(1,40,by =2)] 
linear.phi.CI.L.3 <-linear.CI.23.L[-seq(1,40,by =2)] 

linear.CI.23.U <- linear.CI$Upper[-seq(1,60,by =3)]
linear.phi.CI.U.2 <- linear.CI.23.U[seq(1,40,by =2)]
linear.phi.CI.U.3 <- linear.CI.23.U[-seq(1,40,by =2)]

point.linear.23 <- point.linear.phi[-seq(1,60,by =3)]

point.linear.phi2 <- point.linear.23[seq(1,40,by =2)]
point.linear.phi3 <- point.linear.23[-seq(1,40,by =2)]





point.walk.phi <- get.point.est(line.full.walk,"phi")

walk.CI <- get.CI(line.full.walk,"phi")
walk.CI.23.L<- walk.CI$Lower[-seq(1,60,by =3)]

walk.phi.CI.L.2 <- walk.CI.23.L[seq(1,40,by =2)] 
walk.phi.CI.L.3 <- walk.CI.23.L[-seq(1,40,by =2)] 

walk.CI.23.U <- walk.CI$Upper[-seq(1,60,by =3)]

walk.phi.CI.U.2 <- walk.CI.23.U[seq(1,40,by =2)]
walk.phi.CI.U.3 <- walk.CI.23.U[-seq(1,40,by =2)]

point.walk.23 <- point.walk.phi[-seq(1,60,by =3)]

point.walk.phi2 <- point.walk.23[seq(1,40,by =2)]
point.walk.phi3 <- point.walk.23[-seq(1,40,by =2)]



#gamma2 <- gamma.walk[seq(1,40,by = 2)]
#gamma3 <- gamma.walk[-seq(1,40,by = 2)]

phi.data <- data.frame(method = c(rep("constant",20),rep("constant",20),rep('linear',20),rep('linear',20),rep('walk',20),rep('walk',20)),survey = c(rep("household",20),rep("facebook",20),rep("household",20),rep("facebook",20),rep("household",20),rep("facebook",20)),
                       phi = c(exp(rep(const.gamma[1],20)),exp(rep(const.gamma[2],20)),point.linear.phi2,point.linear.phi3, point.walk.phi2 ,point.walk.phi3),
                       time = c(1:20,1:20,1:20,1:20,1:20,1:20),
                       CI.L= c(exp(rep(const.CI$Lower[1],20)),exp(rep(const.CI$Lower[2],20)),linear.phi.CI.L.2,linear.phi.CI.L.3,walk.phi.CI.L.2,walk.phi.CI.L.3),
                       CI.U =c(exp(rep(const.CI$Upper[1],20)),exp(rep(const.CI$Upper[2],20)),linear.phi.CI.U.2,linear.phi.CI.U.3,walk.phi.CI.U.2,walk.phi.CI.U.3))

#phi.data <- phi.data %>% filter(method %in% c("walk",'constant'))

#phi.data.walk %>%
ggplot(data = phi.data,aes(x=time,y = phi,color = method,shape = survey))+geom_point() + geom_line()+ 
  geom_ribbon(aes(ymin = CI.L,ymax = CI.U),alpha = 0.2) +theme_minimal()+
    labs(x="timepoint",y = expression(phi), title = expression(paste("Plot of estimates of ",phi," by method and survey, parallel slopes",sep = " ")))


#exact effeciency gain
axios.na <- is.na(data.list$Y[1,])

gain <- 100*(axios.ipsos$unb.upper.CI-axios.ipsos$unb.lower.CI)/(CI.full$Upper[!axios.na]-CI.full$Lower[!axios.na])

mean(gain)


#walk vs linear efficiency

mean(100*(CI.walk$Upper-CI.walk$Lower)/(CI.full$Upper-CI.full$Lower))


