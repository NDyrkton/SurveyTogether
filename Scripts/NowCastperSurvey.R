#data with now casting at every point.
library(dclone)
library(MCMCpack)
library(dplyr)
library(rjags)

source("Scripts/JagsMods.R")
source("Scripts/helperfunctions.R")

#import data
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

extract.t <- function(list,t){
  return.list <- list
  
  return.list$T <- t
  
  return.list$Y <- matrix(return.list$Y[,1:t],ncol = t)
  return.list$smalln <- matrix(return.list$smalln[,1:t],ncol = t)
  return.list$times <- matrix(return.list$times[,1:t],ncol = t)
  
  return(return.list)
}




get.point.est <- function(line,var){
  
  
  point.est <- summary(line)$statistics
  
  #breaks if number of time points is 1, checking null fix.
  
  if(is.null(dim(point.est))){
    
    return(point.est[1])
    
  }else{
    means <- summary(line)$statistics[,1]
    
    return(means[grep(var,names(means))])
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




full.data.joined <- read.csv("Data/full_data_joined_ddc_vaccine.csv")
benchmark.data <- read.csv("Data/Benchmark_clean.csv")



ipsos.dat <- extract.unbiased.nona(data.list,col = 1)
household.dat <- extract.unbiased.nona(data.list,col= 2)
facebook.dat <- extract.unbiased.nona(data.list,col = 3)





cl <- makePSOCKcluster(4)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


ipsos.len <- length(ipsos.dat$Y)
fb.len <- length(household.dat$Y)
household.len <- length(household.dat$Y)


#now casts for each.
ipsos.posrates <- numeric(ipsos.len)
ipsos.CI <- list(CI.U=numeric(ipsos.len),CI.L = numeric(ipsos.len))

household.posrates <- numeric(household.len)
household.CI <- list(CI.U=numeric(household.len),CI.L = numeric(household.len))

facebook.posrates <- numeric(fb.len)
facebook.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))


full.posrates <- numeric(fb.len)
full.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))





ipsos.dat <- extract.unbiased.nona(data.list,col = 1)
household.dat <- extract.unbiased.nona(data.list,col= 2)
facebook.dat <- extract.unbiased.nona(data.list,col = 3)


for(t in 1:20){
  
  chain1<- list(.RNG.name = "base::Wichmann-Hill", 
                .RNG.seed = c(159+2*t))
  chain2<- list(.RNG.name = "base::Super-Duper", 
                .RNG.seed = c(260+3*t))
  chain3<- list(.RNG.name = "base::Wichmann-Hill", 
                .RNG.seed = c(371+4*t))
  chain4<- list(.RNG.name = "base::Super-Duper", 
                .RNG.seed = c(482+5*t))
  
  inits.chains <- list(chain1,chain2,chain3,chain4)
  
  
  data.t <- extract.t(data.list,t)
  

  ###ipsos run##
  if(t %in% c(1:ipsos.len)){
    
    ipsos.t <- extract.t(ipsos.dat,t)
    
    line.ipsos <- jags.parfit(cl,ipsos.t, c("positiverate"), custommodel(mod.linear.phi),
                              n.chains=4,n.adapt = 50000,thin = 5, n.iter = 50000,inits = inits.chains)
    
    ipsos.posrates[t] <- get.point.est(line.ipsos,"positiverate")[t]
    
    gelman.diag(line.ipsos)
    
    
    ipsos.CIs <- get.CI(line.ipsos,"positiverate")
    ipsos.CI$CI.U[t] <- ipsos.CIs$Upper[t]
    ipsos.CI$CI.L[t] <- ipsos.CIs$Lower[t]
    
    
    
  }
  
  if(t %in% c(1:household.len)){
    household.t <- extract.t(household.dat,t)
    
    line.household <- jags.parfit(cl,household.t, c("positiverate"), custommodel(mod.linear.phi),
                                  n.chains=4,n.adapt = 50000,thin = 5, n.iter = 50000,inits = inits.chains)
    
    household.posrates[t] <- get.point.est(line.household,"positiverate")[t]
    household.CIs <- get.CI(line.household,"positiverate")
    
    
    household.CI$CI.U[t] <- household.CIs$Upper[t]
    household.CI$CI.L[t] <- household.CIs$Lower[t]
    
  }
  
  
  ###facebook and full###
  
  facebook.t <- extract.t(facebook.dat,t)
  
  line.facebook <- jags.parfit(cl, facebook.t, c("positiverate"), custommodel(mod.linear.phi),
                               n.chains=4,n.adapt = 50000,thin = 5, n.iter = 50000,inits = inits.chains)
  
  facebook.posrates[t] <- get.point.est(line.facebook,"positiverate")[t]
  
  facebook.CIs <- get.CI(line.facebook,"positiverate")
  
  
  facebook.CI$CI.U[t] <- facebook.CIs$Upper[t]
  facebook.CI$CI.L[t] <- facebook.CIs$Lower[t]
  
  
  
  
  
  line.full <- jags.parfit(cl, data.t, c("positiverate","gamma0","gamma1"), custommodel(mod.linear.phi),
                           n.chains=4,n.adapt = 200000,thin = 5, n.iter = 550000,inits = inits.chains)
  
  if(any(gelman.diag(line.full)$psrf[,1] >= 1.1)){
    print("failed")
    
    line.full <- jags.parfit(cl, data.t, c("positiverate","gamma0"), custommodel(mod.linear.phi),
                             n.chains=4,n.adapt = 500000,thin = 5, n.iter = 700000,inits = inits.chains)
    
  }
  
  full.posrates[t] <- get.point.est(line.full,"positiverate")[t]
  
  full.CIs <- get.CI(line.full,"positiverate")
  
  
  full.CI$CI.U[t] <- full.CIs$Upper[t]
  full.CI$CI.L[t] <- full.CIs$Lower[t]
  
  
}

full.data.joined <- read.csv("Data/full_data_joined_ddc_vaccine.csv")
benchmark.data <- read.csv("Data/Benchmark_clean.csv")


full.data.joined$unb.lower.CI <- c(facebook.CI$Lower,household.CI$Lower,ipsos.CI$Lower)
full.data.joined$unb.upper.CI <- c(facebook.CI$Upper,household.CI$Upper,ipsos.CI$Upper)
full.data.joined$pointest <- c(facebook.posrates,household.posrates,ipsos.posrates)

method.dates <- full.data.joined$end_date[full.data.joined$mode=='facebook']

method.data <- data.frame(mode= rep("method",20),end_date = method.dates,pointest= full.posrates,CI.L=full.CI$CI.L,CI.U = full.CI$CI.U)

plot.data <- full.data.joined %>% select(mode,end_date,CI.L,CI.U,pointest) %>% rbind(method.data)


ggplot(plot.data,aes(x = as.Date(end_date),y = pointest, colour = mode)) + geom_point()+
  geom_line() + geom_errorbar(aes(ymin = CI.L,ymax = CI.U))+  theme_minimal() + labs(x = "Date",y = "% Vaccinated",title = expression(paste("Plot of Method for linear ",phi))) +
  geom_ribbon(data = benchmark.data,aes(x = as.Date(date),y = posrate,ymin = posrate*0.98,ymax = posrate*1.02,group = "benchmark",col = "benchmark"),alpha = 0.2)



