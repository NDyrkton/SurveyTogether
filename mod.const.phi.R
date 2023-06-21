
model{	
#likelihood

for (i in 1:Ivec[1]){
		phi[1,i] <- 1	}

for (k in 2:K){
	for (i in 1:Ivec[k]){
		phi[k,i] <- exp(gamma0[k])
	}
}
	
	
logitpositiverate[1] ~ dnorm(theta0,1/0.1)
positiverate[1]	<- ilogit(logitpositiverate[1])
for(t in 2:T){
	logitpositiverate[t] ~ dnorm(logitpositiverate[t-1],1/rho)
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
theta0 ~ dnorm(0, 1);
rho ~ dnorm(0, 1/10)T(0,);

for (k in 1:K){
	gamma0[k] ~ dnorm(0, 1);
}
}
