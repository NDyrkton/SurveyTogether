
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

Y <- matrix(NA,nrow =4, ncol = nrow(facebook.Y))


Y[1,] <-  facebook.Y$Y
Y[2,] <- fill.dates(ref.date = facebook.Y$end_date, vec.date = household_pulse.Y$end_date, vec = household_pulse.Y$Y)
Y[3,] <- fill.dates(ref.date = facebook.Y$end_date, vec.date = ipsos_axios.Y$end_date, vec = ipsos_axios.Y$Y)

n <- matrix(NA,nrow =3, ncol = nrow(facebook.Y))

n[1,] <- facebook.Y$n 
n[2,] <- fill.dates(ref.date = facebook.Y$end_date, vec.date = household_pulse.Y$end_date, vec = household_pulse.Y$n)
n[3,] <- fill.dates(ref.date = facebook.Y$end_date, vec.date = ipsos_axios.Y$end_date, vec = ipsos_axios.Y$n)


#need n

n[2,] <- round(na.interp(n[2,]))
n[3,] <- round(na.interp(n[3,]))

N <- 255200373
t <- 1:20
K = 3
times <- t(matrix(rep(t,K),ncol = K))


data.list <- list(K = 3, T = 20, N = N, times = times, Y = Y, smalln = n)



mod.const.phi<- custommodel('
model{	
#likelihood


for (k in 1:K){
	for (t in 1:T){
		phi[k,t] <- gamma0[k]
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
theta0 ~ dnorm(0, 1/0.001);
rho ~ dnorm(0, 1/0.001)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1/0.05);
}
}')

mod.linear.phi <- custommodel('
model{	
#likelihood

for (k in 1:K){
	for (t in 1:T){
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times[k,t])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/10)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1], 1/rho)
	positiverate[t]	<- ilogit(logitpositiverate[t])
}

for(t in 1:T){
	#Pt[t] ~ dnorm(positiverate[t]*N, 1/ (positiverate[t]*(1-positiverate[t])*N))T(0,N);
	#P[t] <- dround(Pt[t],0)
	P[t] ~ dbin(positiverate[t], N)
	

}


for (k in 1:K){
	for (t in 1:T){
	  
		
		Y[k,t] ~ dhyper(ifelse(P[times[k,t]]< Y[k,t], Y[k,t],P[times[k,t]]), N-(ifelse(P[times[k,t]]< Y[k,t], Y[k,t],P[times[k,t]])), smalln[k,t], phi[k,t]);
		#we are getting Y> P[t] 
	}
}

#priors
theta0 ~ dnorm(10, 1/0.1);
rho ~ dnorm(0, 1/10)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(3, 1/0.001);
	gamma1[k] ~ dnorm(0.1, 1/0.001);
}
}')



cl <- makePSOCKcluster(4)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")




line.linear <- jags.parfit(cl, data.list, c("positiverate","gamma0","gamma1","rho"), mod.const.phi,
            n.chains=4,n.adapt = 10000,thin = 10, n.iter = 25000)

#provides errors



###fix

data.list2 <- generate.dataset(N =   255200373,t = 1:20, ns = rep(25000,20),phi= "linear")

data.list2$Y[2,c(3,5,10,12,18)] <- NA
data.list2$Y[3,c(1,2,8,11,20)] <- NA


#data.list.3 <- generate.dataset(phi= "linear")

line.linear2 <- jags.parfit(cl, data.list2, c("positiverate","gamma0","gamma1","rho"), mod.linear.phi,
                           n.chains=4,n.adapt = 1000,thin = 10, n.iter = 1000)




