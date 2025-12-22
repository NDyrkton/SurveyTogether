library(foreach)
library(doParallel)
library(doRNG)
library(rjags)
library(coda)
library(MCMCpack)
library(truncnorm)


inv.logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/(1-x))
}

generate.dataset <- function(N= 10000000, K =3, ti = c(1:5), phi = "constant"){
  Y <- matrix(NA,ncol = length(ti),nrow = K)
  smalln <- matrix(0,ncol = length(ti),nrow = K)
  #get default smalln
  
  for(k in 1:K){
    smalln[k,] <- rep(100,length(ti))
    
    if(k>1){
      smalln[k,] <- rep(50000,length(ti))
    }
  }
  
  posrate_t <- numeric(length(ti))
  theta_t <- numeric(length(ti))
  times <- t(matrix(rep(ti,K),ncol = K))
  
  #priors on general parameters
  sigmasq <- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = sqrt(0.1))
  theta0 <- rnorm(1,mean =0, sd = sqrt(1))
  
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
    gamma_1k <- c(0,rnorm(K-1,mean = 0, sd = rep(sqrt(0.1),K-1)))
    
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
    
    pisq <- rtruncnorm(1,a = 0, b = Inf, mean = 0, sd = sqrt(0.1))
    
    
    
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

extract.unbiased <- function(datalist){
  K = 1
  ti <- datalist$ti
  new.list <- list(K = K, 
                   times = matrix(datalist$times[1,],ncol = ti), N = datalist$N, ti = ti,
                   Y = matrix(datalist$Y[1,],ncol = ti), smalln = matrix(datalist$smalln[1,],ncol = ti),
                   params = datalist$params)
  
  return(new.list)
}

check.bad.data <- function(data.list, ti, phi = 'constant', K = 3){
  
  max.n <- max(data.list$smalln)
  
  
  while(sum(data.list$Y[2:K,] == max.n | data.list$Y[2:K,] == 0)>=1){
  
    data.list <- generate.dataset(phi = phi,ti = ti,K = K)
    
  }
  
  return(data.list)
  
}


mod.const.phi <- ('
model{	
#likelihood

for (t in 1:ti){
		phi[1,t] <- 1	}

for (k in 2:K){
	for (t in 1:ti){
		phi[k,t] <- exp(gamma0[k])
	}
}
	
	
theta[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(theta[1])
for(t in 2:ti){
	theta[t] ~ dnorm(theta[t-1], 1/sigmasq)
	positiverate[t]	<- ilogit(theta[t])
}


for (k in 1:K){
	for (t in 1:ti){
		
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1);
sigmasq ~ dnorm(0, 1/0.1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
}
}')

mod.linear.phi <- ('
model{	
#likelihood

for (t in 1:ti){
		phi[1,t] <- 1	}

for (k in 2:K){
	for (t in 1:ti){
		phi[k,t] <- exp(gamma0[k] + gamma1[k]*times_centered[k,t])
	}
}
	
	
theta[1] ~ dnorm(theta0,1/sigmasq)
positiverate[1]	<- ilogit(theta[1])
for(t in 2:ti){
	theta[t] ~ dnorm(theta[t-1], 1/sigmasq)
	positiverate[t]	<- ilogit(theta[t])
}


for (k in 1:K){
	for (t in 1:ti){
		
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1);
sigmasq ~ dnorm(0, 1/0.1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);
	gamma1[k] ~ dnorm(0, 1/0.1);
}
}')

mod.walk.phi <- ('
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


for (k in 1:K){
	for (t in 1:ti){
		
		Y[k,t] ~ dbin(   (positiverate[t]*phi[k,t])/(1-positiverate[t] + (positiverate[t]*phi[k,t])),   smalln[k,t])
	}
}

#priors
theta0 ~ dnorm(0, 1);
sigmasq ~ dnorm(0, 1/0.1)T(0,);
pisq ~ dnorm(0, 1/0.1)T(0,);

for (k in 2:K){
	gamma0[k] ~ dnorm(0, 1);

}

}')


#cat(paste0("Running ", n_cores, " simultaneous 'MCMC Jags Fits' (Jobs).\n"))

#create scenarios
n_reps <- 1000
scenarios <- expand.grid(model = c("const","linear","walk"), dgm = c("const","linear","walk"), replicate_id = 1:n_reps)
scenarios$full.scenario <- paste(scenarios$dgm,"x",scenarios$model)
set.seed(42) 

models <- list(const = mod.const.phi, linear = mod.linear.phi, walk = mod.walk.phi)


#dat <- generate.dataset()

# run function

run.mod <- function(scenario, mod.name, model, dat, param, ti) {
  
  #truth
  posrate.t <- unname(dat$params[grep(paste0("positiverate",ti),names(dat$params))])
  
  
  model.txt <- textConnection(model)
  #fit the data
  model.fit <- suppressWarnings(rjags::jags.model(file =  model.txt,
    data = dat,
    n.chains = 4,      # 4 Chains running sequentially
    n.adapt = 5000,
    quiet = TRUE
  ))
  
  # burn in
  update(model.fit, n.iter = 1000, progress.bar = "none")
  
  # sample from posterior
  samps <- rjags::coda.samples(
    model.fit, 
    variable.names = param, 
    n.iter = 50000, 
    progress.bar = "none"
  )
  
  # calculate the median
  estimate <- unname(as.numeric((summary(samps)$quantiles[3])))
  #estimate <- unname(summary(samps)$statistics[1])
  close(model.txt)
  rm(samps)
  
  # Calculate R-hat for convergence check
  # rhat <- tryCatch({
  #   coda::gelman.diag(samps, multivariate = FALSE)$psrf[,1]
  # }, error = function(e) return(c(mu=NA, sigma=NA)))
  # 
  # return params
  return(data.frame(
    dgm = scenario,
    model = mod.name,
    est = estimate,
    posrate = posrate.t,
    bias = estimate-posrate.t
  ))
}

run.test <- run_iteration(scenario = scenarios[1,],model = mod.const.phi,dat = dat, param = "positiverate[5]",ti = 5)

run.test
# 5. Run Parallel Loop ----------------------------------------------------

# The %dorng% operator ensures reproducible unique data generation across all cores.


# setup parallel
#options(future.globals.maxSize = 2000 * 1024^2) # Increase to 2GB

cl <- makeCluster(19)
registerDoParallel(cl)

n_reps <- 100
scenarios <- expand.grid(model = c("constant","linear","walk","unbiased"), dgm = c("constant","linear","walk"), replicate_id = 1:n_reps)
scenarios$full.scenario <- paste(scenarios$dgm,"x",scenarios$model)


scenarios <- expand.grid(dgm = c("constant","linear","walk"), replicate_id = 1:n_reps)
#scenarios$full.scenario <- paste(scenarios$dgm,"x",scenarios$model)

models <- list(const = mod.const.phi, linear = mod.linear.phi, walk = mod.walk.phi)
#set.seed(42) 


time = 15
posrate.param <- paste("positiverate","[",time,"]")

time1 <- Sys.time()
results_matrix <- foreach(id = 1:nrow(scenarios), 
                          .packages = c("rjags", "coda","truncnorm","MCMCpack"), 
                          .combine = rbind,
                          .options.RNG = 20251220) %dorng% {
                            
                            #gc(verbose = FALSE)
                            #picked scenario
                            sc <- scenarios[id,]
                            
                            #data gen
                            dgm <- generate.dataset(ti = 1:time, phi = sc$dgm)
                            
                            dgm <- check.bad.data(dgm,ti = 1:time , K = 3, phi = sc$dgm)
                            
                            dgm.unbiased <- extract.unbiased(dgm)
                            
                            
    
                            #fit 4 models
                            fit.const <- run.mod(sc$dgm, mod.name = "const", model = mod.const.phi,dat = dgm, param = posrate.param , ti = time)
                            fit.linear <- run.mod(sc$dgm, mod.name = "linear" ,model = mod.linear.phi,dat = dgm, param = posrate.param , ti = time)                            
                            fit.walk <- run.mod(sc$dgm, mod.name = "walk",model = mod.walk.phi,dat = dgm, param = posrate.param , ti = time)                            
                            fit.unbiased <- run.mod(sc$dgm, mod.name = "unbiased", model = mod.const.phi, dat = dgm.unbiased, param = posrate.param , ti = time) 
                            
                            
                            
                            return.df <- do.call(rbind,list(fit.const, fit.linear, fit.walk, fit.unbiased))
                            
                            

                            # run the iteration
                            return(return.df)
                          }

stopCluster(cl)
time2 <- Sys.time()

print(time2-time1)

# 2/10
# 
# 0.2 mins per iterations
# 
# 0.2*1000

#200 mins for 1000 iterations

##
library(dplyr)
library(ggplot2)

#group the data

grouped.results <- results_matrix %>% group_by(dgm,model) %>% summarise(bias = mean(bias),mse = mean(bias^2))

ggplot(grouped.results,aes(x = dgm, y = bias, group = model,colour = model)) + geom_line() + geom_point()
ggplot(grouped.results,aes(x = dgm, y = mse, group = model,colour = model)) + geom_line() + geom_point()

results_df <- as.data.frame(results_matrix)
head(results_df)
