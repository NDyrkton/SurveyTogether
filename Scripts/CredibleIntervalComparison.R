#Compare Confidence Intervals using unbiased
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

ipsos.dat <- extract.unbiased.nona(data.list,col = 1)
household.dat <- extract.unbiased.nona(data.list,col=2)
facebook.dat <- extract.unbiased.nona(data.list,col = 3)


cl <- makePSOCKcluster(4)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")

#run linear phi/ walk phi

line.ipsos <- jags.parfit(cl, ipsos.dat, c("positiverate","gamma0","gamma1"), mod.linear.phi,
                           n.chains=4,n.adapt = 200000,thin = 5, n.iter = 50000)

#check convergence
gelman.diag(line.ipsos)

point.ipsos <- get.point.est(line.ipsos,"positiverate")
CI.ipsos <- get.CI(line.ipsos,"positiverate")


line.household <- jags.parfit(cl, household.dat, c("positiverate","gamma0","gamma1"), mod.linear.phi,
                             n.chains=4,n.adapt = 200000,thin = 5, n.iter = 50000)

gelman.diag(line.household)

point.household <- get.point.est(line.household,"positiverate")
CI.household <- get.CI(line.household,"positiverate")


line.facebook <- jags.parfit(cl, facebook.dat, c("positiverate","gamma0","gamma1"), mod.linear.phi,
                             n.chains=4,n.adapt = 200000,thin = 5, n.iter = 50000)


gelman.diag(line.facebook)
point.facebook <- get.point.est(line.facebook,"positiverate")
CI.facebook <-  get.CI(line.facebook,"positiverate")



#run actual method


line.full <- jags.parfit(cl, data.list, c("positiverate","gamma0","gamma1"), mod.linear.phi,
                         n.chains=4,n.adapt = 500000,thin = 10, n.iter = 800000)


gelman.diag(line.full) 

point.full <- get.point.est(line.full,"positiverate")
CI.full <- get.CI(line.full,"positiverate")



View(full.data.joined)
#in order facebook->household -> axios
full.data.joined$unb.lower.CI <- c(CI.facebook$Lower,CI.household$Lower,CI.ipsos$Lower)
full.data.joined$unb.upper.CI <- c(CI.facebook$Upper,CI.household$Upper,CI.ipsos$Upper)
full.data.joined$pointest <- c(point.facebook,point.household,point.ipsos)

facebook.end_date <- full.data.joined[full.data.joined$mode=="facebook",'end_date']


data.method <- data.frame(CI.Lower =CI.full$Lower,CI.Upper = CI.full$Upper,posrate = point.full,end_date = facebook.end_date)





full.data.joined$end_date <- as.Date(full.data.joined$end_date)
ggplot(data = full.data.joined,aes(x = end_date,y = posrate,col = mode,group = mode))+
  geom_point()+geom_line() +geom_errorbar(aes(ymin = CI.L,ymax = CI.U),col = 'black',lwd = 0.5)+
  geom_errorbar(aes(ymin=unb.lower.CI,ymax = unb.upper.CI)) + 
  geom_point(data = data.method,aes(x = end_date,y = posrate,group = "method",col = "method"))+
  geom_line(data = data.method,aes(x = end_date,y = posrate,group = "method",col = "method")) + 
  geom_errorbar(data = data.method,aes(ymin = CI.Lower,ymax = CI.Upper,group = "method",col = "method"),width = 2)+
  theme_minimal() + labs(x = "Date",y = "% Vaccinated",title = "Plot of loss of uncertainty linear Phi") +
  geom_ribbon(data = benchmark.data,aes(x = as.Date(date),y = posrate,ymin = posrate*0.98,ymax = posrate*1.02,group = "benchmark",col = "benchmark"),alpha = 0.2)


#get gammas

point.gamma0 <- get.point.est(line.full,"gamma0")
point.gamma1 <- get.point.est(line.full,"gamma1")

plot(1:20,)



