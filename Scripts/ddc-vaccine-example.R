
#Survey Together on real vaccine data, taken from https://github.com/vcbradley/ddc-vaccine-US

library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)
library(dclone)
library(dplyr)
library(forecast)

source("Scripts/helperfunctions.R")

fill.dates <- function(ref.date, vec.date, vec){
  
  return.vec <- numeric(length(ref.date))
  
  for(i in 1:length(vec.date)){
    
    tmp <- which( (vec.date[i]>= ref.date) & (vec.date[i] <= ref.date + 6))
    
    
    return.vec[tmp] <- vec[i] 
    
  }
  
  return.vec[which(return.vec==0)] <- NA
  return(return.vec)
  
}


#import data 

data.all <- read.csv("Data/full_data_ddc_vacine_us.csv")
data.benchmark <- read.csv("Data/benchmark-2021-05-26.csv")

data.all <- data.all %>% filter(end_date <= '2021-05-26' & pop == "US"& pct_error == 0)
data.all$Y <- round(data.all$n* data.all$pct_vaccinated)

data.benchmark <- data.benchmark %>% filter(state == "US") 

facebook <- data.all %>% filter(mode == "facebook")
household_pulse <- data.all %>% filter(mode == "household_pulse")
ipsos_axios <- data.all %>% filter(mode == "ipsos_axios")


facebook.Y <- facebook %>% group_by(end_date) %>% summarise(Y = max(Y), n = max(n))
household_pulse.Y <-  household_pulse %>% group_by(end_date) %>% summarise(Y = max(Y), n = max(n))
ipsos_axios.Y <- ipsos_axios %>% group_by(end_date) %>% summarise(Y = max(Y), n = max(n))

start.date <- min(data.all$end_date) ## 2021-01-09 #facebook only has date up to 2021-01-03
end.date <- max(data.benchmark$date) # 2021-05-25 #facebook ends at 2021-05-22

#ends on week of the 10th. so for each monday we start a new week (monday the 4th is the start of the week)
start.group <- as.Date("2021-01-09") #facebook has 20 weeks
end.group <- as.Date("2021-05-24")

data.benchmark <- data.benchmark %>% filter(date >= start.group & date <= end.group)
facebook.Y$end_date <- as.Date(facebook.Y$end_date)
household_pulse.Y$end_date <- as.Date(household_pulse.Y$end_date)
ipsos_axios.Y$end_date <- as.Date(ipsos_axios.Y$end_date)

Y <- matrix(NA,nrow =3, ncol = nrow(facebook.Y))


Y[3,] <-  facebook.Y$Y
Y[2,] <- fill.dates(ref.date = facebook.Y$end_date, vec.date = household_pulse.Y$end_date, vec = household_pulse.Y$Y)
Y[1,] <- fill.dates(ref.date = facebook.Y$end_date, vec.date = ipsos_axios.Y$end_date, vec = ipsos_axios.Y$Y)

n <- matrix(NA,nrow =3, ncol = nrow(facebook.Y))

n[3,] <- facebook.Y$n 
n[2,] <- fill.dates(ref.date = facebook.Y$end_date, vec.date = household_pulse.Y$end_date, vec = household_pulse.Y$n)
n[1,] <- fill.dates(ref.date = facebook.Y$end_date, vec.date = ipsos_axios.Y$end_date, vec = ipsos_axios.Y$n)


#need n

n[2,] <- round(na.interp(n[2,]))
n[1,] <- round(na.interp(n[1,]))

N <- 255200373
t <- 1:20
K = 3
times <- t(matrix(rep(t,K),ncol = K))


data.list <- list(K = 3, T = 20, N = N, times = times, Y = Y, smalln = n)


mod.linear.phi.2 <- custommodel('
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
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.01)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dunif(logitpositiverate[t-1],logitpositiverate[t-1]+rho)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

for(t in 1:T){
	P[t] ~ dbin(positiverate[t], N)
}


for (k in 1:K){
	for (t in 1:T){
	  
		Y[k,t] ~ dbin(1-(1-(P[t]/N))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
rho ~ dnorm(0, 1/5)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);
	gamma1[k] ~ dnorm(0, 1/0.001);
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
	  
		Y[k,t] ~ dbin(1-(1-(P[t]/N))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
rho ~ dnorm(0, 1/0.01)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);
	gamma1[k] ~ dnorm(0, 1/0.001);
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
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.001)
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
		
		Y[k,t] ~ dbin(1-(1-(P[t]/N))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
rho ~ dnorm(0, 1/2)T(0,);
pi ~ dnorm(0, 1/0.01)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);

}

}')

mod.walk.phi.2 <- custommodel('
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
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.001)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dunif(logitpositiverate[t-1],logitpositiverate[t-1]+rho)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

for(t in 1:T){
	P[t] ~ dbin(positiverate[t], N)
}

for (k in 1:K){
	for (t in 1:T){
		
		Y[k,t] ~ dbin(1-(1-(P[t]/N))^phi[k,t],smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
rho ~ dnorm(0, 1)T(0,);
pi ~ dnorm(0, 1/0.01)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);

}

}')





cl <- makePSOCKcluster(6)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")




line.linear <- jags.parfit(cl, data.list, c("positiverate","gamma0","gamma1","rho"), mod.linear.phi.2,
            n.chains=5,n.adapt = 200000,thin = 100, n.iter = 5000000)

means.posrate <- summary(line.linear)$statistics[,1]
#check
pos.rate <- means.posrate[grep("positiverate",names(means.posrate))]
CI.lower <- summary(line.linear)$quantile[,1][grep("positiverate",names(means.posrate))]
CI.upper <- summary(line.linear)$quantile[,5][grep("positiverate",names(means.posrate))]

facebook.pred <- facebook %>% group_by(end_date) %>% summarise(posrate = max(pct_vaccinated))
ipsos.pred <- ipsos_axios %>% group_by(end_date) %>% summarise(posrate = max(pct_vaccinated))
household_pulse.pred <- household_pulse %>% group_by(end_date) %>% summarise(posrate = max(pct_vaccinated))

my.data <- data.frame(dates= facebook.Y$end_date, posrate = pos.rate)

data.benchmark <- data.frame(dates = data.benchmark$date[which(data.benchmark$date %in% facebook.Y$end_date)], posrate = data.benchmark$pct_pop_vaccinated[which(data.benchmark$date %in% facebook.Y$end_date)])

data.bench.unique <- data.benchmark %>% group_by(dates) %>% summarise(posrate = max(posrate)) %>% transform(dates = as.Date(dates))

final.plot <- data.frame(dates = c(facebook.Y$end_date,data.bench.unique$dates,facebook.pred$end_date,ipsos.pred$end_date,household_pulse.pred$end_date), 
                         posrate = c(pos.rate,data.bench.unique$posrate,facebook.pred$posrate,ipsos.pred$posrate,household_pulse.pred$posrate),
                         source = c(rep("Method",20),rep("Benchmark",20),rep("facebook",20),rep("Ipsos-axios",length(ipsos.pred$posrate)),rep("Household-Pulse",length(household_pulse.pred$posrate))))

ggplot(data = final.plot,aes(x = as.Date(dates),y = posrate,colour = source)) + geom_point() + geom_line() + theme_minimal() + geom_crossbar(data = NULL,aes(ymin = CI.lower,ymax = CI.upper),inherit.aes = F)
  



line.walk <- jags.parfit(cl, data.list, c("positiverate","gamma0","rho"), mod.walk.phi.2,
                           n.chains=6,n.adapt = 150000,thin = 100, n.iter = 5000000)


means.posrate2 <- summary(line.walk)$statistics[,1]


CI.lower.2 <- summary(line.walk)$quantile[,1][grep("positiverate",names(means.posrate2))]
CI.upper.2 <- summary(line.walk)$quantile[,5][grep("positiverate",names(means.posrate2))]

pos.rate2 <- means.posrate2[grep("positiverate",names(means.posrate2))]

my.data2 <- data.frame(dates= facebook.Y$end_date, posrate = pos.rate2)


ggplot(data = my.data2,aes(x = as.Date(dates),y = posrate)) + geom_point(colour = 'blue') + geom_line(colour = "blue")+
  geom_point(data = data.bench.unique,aes(x = dates,y = posrate),colour = "grey") + geom_line(data = data.bench.unique,aes(x = dates,y = posrate),colour = "grey") +
  geom_ribbon(aes(ymin = CI.lower.2,ymax = CI.upper.2),alpha = 0.2) + theme_minimal() + labs(x = "Date",y = "Vaccine Rate") +
  









#provides errors



###fix

data.list2 <- generate.dataset(N =   255200373,t = 1:20, ns = rep(250000,20),phi= "linear")

data.list2$Y[2,c(3,5,10,12,18)] <- NA
data.list2$Y[3,c(1,2,8,11,20)] <- NA



line.linear2 <- jags.parfit(cl, data.list2[-8], c("positiverate","gamma0","gamma1","rho"), mod.linear.phi,
                           n.chains=4,n.adapt = 10000,thin = 10, n.iter = 50000)

summary(line.linear2)$statistics[,1]

data.list2$params

