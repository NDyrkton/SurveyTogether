
mod.const.phi <- '
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
		#Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
		#Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
sigmasq ~ dnorm(0, 1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
}
}'


mod.linear.phi <-'
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
		#Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
		#Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
sigmasq ~ dnorm(0, 1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
	gamma1[k] ~ dnorm(0, 1/0.25);
}
}'


mod.walk.phi <- '
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
		
		#Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(-2, 1);
sigmasq ~ dnorm(0, 1)T(0,);
pisq ~ dnorm(0, 1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);

}

}'


mod.linear.phi.slope <-'
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
		#Y[k,t] ~ dbin(1-(1-(positiverate[t]))^phi[k,t],smalln[k,t])
		#Y[k,t] ~ dhyper(P[times[k,t]], N-P[times[k,t]], smalln[k,t], phi[k,t]);
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
		
	}
}

#priors
theta0 ~ dnorm(-2, 1);
sigmasq ~ dnorm(0, 1/5)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
}

gamma1 ~ dnorm(0, 1/0.25);


}'



mod.walk.phi.slope <- '
model{	
#likelihood

for (i in 1:T){
		phi[1,i] <- 0.98 + 0.001*i	
}


for(t in 1:T){

  gamma[t] ~ dnorm(gamma0[k],1/pisq)
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
pisq ~ dnorm(0, 1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);

}

}'



mod.walk.linear.phi <- '
model{	
#likelihood

for (i in 1:T){
		phi[1,i] <- 0.99 + 0.002*i
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
pisq ~ dnorm(0, 1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);

}

}'




