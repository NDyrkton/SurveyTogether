#ddc-vaccine run for compute canada

library(dclone)
library(rjags)

data.list <- readRDS("ddc_list.Rdata")


#12 cores, 12 chains
cl <- makePSOCKcluster(12)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


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

#1 billion
adapt <- 1e9
iter <- 1e6 #1 million



line.linear <- jags.parfit(cl, data.list, c("positiverate","gamma0",'gamma1','rho','theta0'), mod.linear.phi.2,
                           n.chains=12,n.adapt = adapt, thin = 100, n.iter = iter)

means.posrate <- summary(line.linear)$statistics[,1]
#check

pos.rate <- means.posrate[grep("positiverate",names(means.posrate))]


CI.lower <- summary(line.linear)$quantile[,1][grep("positiverate",names(means.posrate))]
CI.upper <- summary(line.linear)$quantile[,5][grep("positiverate",names(means.posrate))]

print(gelman.diag(line.linear)) #see convergence

results <- data.frame(mean = means,posrate,CI.L = CI.lower, CI.U = CI.upper)

write.csv(results,file = "Run_ddc_example.csv")

stopCluster(cl)

