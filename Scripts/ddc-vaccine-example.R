
#Survey Together on real vaccine data, taken from https://github.com/vcbradley/ddc-vaccine-US

library(MCMCpack)
library(rjags)
library(truncnorm)
library(ggplot2)
library(dclone)
library(dplyr)
library(forecast)
library(stringr)

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
data.all$end_date <- as.Date(data.all$end_date)
data.benchmark$date <- as.Date(data.benchmark$date)


#or.
full.data.joined <- data.all %>% group_by(mode,end_date) %>% summarise(Y = max(Y), n = max(n),CI.L = max(ci_2.5_samp),CI.U = max(ci_97.5_samp),posrate = max(pct_vaccinated))


facebook <- data.all %>% filter(mode == "facebook")
household_pulse <- data.all %>% filter(mode == "household_pulse")
ipsos_axios <- data.all %>% filter(mode == "ipsos_axios")


facebook.Y <- facebook %>% group_by(end_date) %>% summarise(Y = max(Y), n = max(n),CI.L = max(ci_2.5_samp),CI.U = max(ci_97.5_samp))
household_pulse.Y <-  household_pulse %>% group_by(end_date) %>% summarise(Y = max(Y), n = max(n),CI.L = max(ci_2.5_samp),CI.U = max(ci_97.5_samp))
ipsos_axios.Y <- ipsos_axios %>% group_by(end_date) %>% summarise(Y = max(Y), n = max(n),CI.L = max(ci_2.5_samp),CI.U = max(ci_97.5_samp))

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


#M is no more than 10x

data.list <- list(K = 3, T = 20, N = N, times = times, Y = Y, smalln = n)

#data.list$Y[1,] <- round(data.list$Y[1,]*100)
#data.list$smalln[1,] <- round(data.list$smalln[1,]*100)

#data.list$Y[2,] <- round(data.list$Y[2,]*10)
#data.list$smalln[2,] <- round(data.list$smalln[2,]*10)

#data.list$Y[3,] <- round(data.list$Y[3,]/100)
#data.list$smalln[3,] <- round(data.list$smalln[3,]/100)



#unbiased survey is present at all points

missing <- !is.na(Y[1,])
Y <-  Y[,missing]
t <- 1:11
K = 3
times <- t(matrix(rep(t,K),ncol = K))
n <- n[,missing]


data.list.full <- list(K = 3, T = 11, N = N, times = times, Y = Y, smalln = n)


#write.csv the full data to compare

write.csv(full.data.joined,"Data/full_data_joined_ddc_vaccine.csv",row.names = F)


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
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.01)
positiverate[1]	<- ilogit(logitpositiverate[1])



for(t in 2:T){

	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],rho)T(logitpositiverate[t-1],)
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
rho ~ dnorm(0, 1)T(0,);

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
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.01)
positiverate[1]	<- ilogit(logitpositiverate[1])



for(t in 2:T){

	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],rho)T(logitpositiverate[t-1],)
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
rho ~ dnorm(0, 1)T(0,);

for (k in 2:K){
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
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],rho)T(logitpositiverate[t-1],);
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
rho ~ dnorm(0, 1)T(0,);
pi ~ dnorm(0, 1/20)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);

}

}')



#create 4 parallel chains

cl <- makePSOCKcluster(4)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")

#run linear phi/ walk phi

line.linear <- jags.parfit(cl, data.list, c("positiverate","gamma"), mod.walk.phi,
            n.chains=4,n.adapt = 500000,thin = 10, n.iter = 100000)

means.posrate <- summary(line.linear)$statistics[,1]
pos.rate <- means.posrate[grep("positiverate",names(means.posrate))]



CI.lower <- summary(line.linear)$quantile[,1][grep("positiverate",names(means.posrate))]
CI.upper <- summary(line.linear)$quantile[,5][grep("positiverate",names(means.posrate))]



final.plot <- full.data.joined %>% select(mode,end_date,posrate,CI.L,CI.U)
preds <- data.frame(mode = rep("method",20),CI.L = CI.lower, CI.U = CI.upper,posrate = pos.rate,end_date = facebook.Y$end_date)
final.plot <- rbind(final.plot,preds)
#benchmark

data.bench.plot <- data.benchmark %>% filter(date %in% facebook.Y$end_date) %>% group_by(date) %>% summarise(posrate = max(pct_pop_vaccinated))
write.csv(data.bench.plot,"Data/Benchmark_clean.csv",row.names = F)

#only important one

final.plot <- final.plot %>% filter(mode %in% c('method',"ipsos_axios","facebook","household_pulse"))


ggplot(data = final.plot,aes(x = end_date,y = posrate,colour = mode)) + geom_point() + geom_line() + 
  theme_minimal() + geom_errorbar(aes(ymin = CI.L, ymax = CI.U)) +
  geom_point(data = data.bench.plot,aes(x = as.Date(date),y = posrate),colour = 'grey') + 
  geom_line(data = data.bench.plot,aes(x = as.Date(date),y = posrate),colour = 'grey') +
  labs(x = "Date", y = "Percentage Vaccinated",title = "Plot of survey estimates for random walk phi")
  

gelman.diag(line.linear) # still has not converged after 5 million burn-in



plot(line.linear)
gelman.plot(line.linear)


gammas <- means.posrate[1:40]
gamma2 <- gammas[seq(1,40,by = 2)]
gamma3 <- gammas[-seq(1,40,by = 2)]

gamma.data <- data.frame(time = rep(1:20,2), gamma = c(gamma2,gamma3),survey = c(rep("Household-Pulse",20),rep("Facebook",20)))


ggplot(data = gamma.data,aes(x = time, y = exp(gamma), col = survey)) + geom_point() + geom_line() +
  theme_minimal() + labs(x = "Time",y = expression("phi"),title = "Plot of Phi Over Time by Survey")



#export
saveRDS(data.list,file = "ddc_list.Rdata")
readRDS("ddc_list.Rdata")