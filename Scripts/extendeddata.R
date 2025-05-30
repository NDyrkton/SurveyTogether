#redo on extended data.
library(dclone)
library(ggplot2)
library(dplyr)
library(lemon)
library(ggtext)
library(forecast)
library(rjags)

# fill in dates
fill.dates <- function(ref.date, vec.date, vec){
  
  return.vec <- numeric(length(ref.date))
  
  for(i in 1:length(vec.date)){
    
    tmp <- which( (vec.date[i]>= ref.date) & (vec.date[i] <= ref.date + 6))
    
    
    return.vec[tmp] <- vec[i] 
    
  }
  
  return.vec[which(return.vec==0)] <- NA
  return(return.vec)
  
}

#extract the unbiased survey from the data list
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


# point estimate to extract from jags
get.point.est <- function(line,var,type = "median"){
  
  
  point.est <- summary(line)$statistics
  
  #breaks if number of time points is 1, checking null fix.
  
  if(is.null(dim(point.est))){
    
    return(point.est[1])
    
  }else{
    if(type == "mean"){
      means <- summary(line)$statistics[,1]
      return(means[grep(var,names(means))])
      
    }else if(type == "median"){
      
      medians <- summary(line)$quantiles[,3]
      return(medians[grep(var,names(medians))])
      
      
    }

  }
  
}

#collect the 95% CIs from the mcmc lines
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
#extract the individual surveys from the data list
extract.surveys <- function(datalist,row = 1){
  #function extracts given surveys from data list
  #same as extract.nona but does not shorten
  
  K <-  length(row)
  Y <- datalist$Y[row,]
  
  smalln <- datalist$smalln[row,]
  times <- datalist$times[row,]
  T <- datalist$T
  
  new.list <- list(K = K, 
                   times = matrix(times,ncol = T), N = datalist$N, T = T,
                   Y = matrix(Y,ncol = T), smalln = matrix(smalln,ncol = T))
  return(new.list)
  
}




source("Scripts/JagsMods.R")
#source("Scripts/helperfunctions.R")


fb_df <- read.csv("Data/fb.csv")
cdc_df <- read.csv("Data/cdc.csv")
chp_df <- read.csv("Data/chp.csv")
ax_df <- read.csv("Data/ax.csv")
ax_n <- read.csv("Data/ax_n.csv")

fb_df$Y <- fb_df$vax*fb_df$n


#put in the n for axios ipsos
ax_df$n <- ax_n$n
#get Y
ax_df$Y <- ax_df$n*ax_df$vax

ax_df <- ax_df[nrow(ax_df):1,]


fb_df$ymd <- as.Date(fb_df$ymd)
cdc_df$ymd <- as.Date(cdc_df$ymd)
chp_df$ymd <- as.Date(chp_df$ymd)
ax_df$ymd <- as.Date(ax_df$ymd)



#chp  phase of interests are : phase 3, 3.1, and 3.2
 #weeks phase 3.2 weeks 39-34

chp_n_phase3.2 <- c(57064, 59833, 63536, 69114, 68799, 64562)
#weeks 33-28 ## min is april 14 2021
chp_n_phase3.1 <- c(66262, 68067, 70858, 78467, 68913)
#weeks  27-22march29-jan 18
chp_n_phase3 <- c(77104, 78306,77788 , 77122, 80567, 68348)

chp_n = c(chp_n_phase3.2,chp_n_phase3.2,chp_n_phase3)



#invert
chp_df$n <- chp_n[length(chp_n):1]
chp_df$Y <- chp_df$n*chp_df$vax
#taken from pdf

#start date is Jan-11-2021, max date is dec 4 2021

#round Y

fb_df$Y <- round(fb_df$Y)
chp_df$Y <- round(chp_df$Y)
ax_df$Y <- round(ax_df$Y)


#fb dates are the baseline
ref.dates <- fb_df$ymd


ax_y_full <- fill.dates(ref.dates, ax_df$ymd, ax_df$Y)
chp_y_full <- fill.dates(ref.dates, chp_df$ymd, chp_df$Y)

Y <- matrix(NA,ncol = length(ref.dates),nrow = 3)

Y[1,] <- ax_y_full
Y[2,] <- chp_y_full
Y[3,] <- fb_df$Y


n <-  matrix(NA,ncol = length(ref.dates),nrow = 3)

ax_n_full <- fill.dates(ref.dates, ax_df$ymd, ax_df$n)
chp_n_full <- fill.dates(ref.dates, chp_df$ymd, chp_df$n)

n[1,] <- round(na.interp(ax_n_full))
n[2,] <- round(na.interp(chp_n_full))
n[3,] <- fb_df$n

times <- matrix(rep(1:length(ref.dates),3),nrow = 3,byrow = T)
#data.list



data.list.extended <- list(K = 3, T = 48, N = 255200373,Y = Y, times =times, smalln = n)



#saveRDS(data.list.extended,"data.list.extended.Rdata")


#run starting here with pre-cals

#data.list.extended <- data.list


#START with full data

#get unbiased data with no NAs
ipsos.dat <- extract.unbiased.nona(data.list.extended,col = 1)
household.dat <- extract.unbiased.nona(data.list.extended,col= 2)
facebook.dat <- extract.unbiased.nona(data.list.extended,col = 3)


cl <- makePSOCKcluster(6)

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
chain5<- list(.RNG.name = "base::Wichmann-Hill", 
              .RNG.seed = c(334))
chain6<- list(.RNG.name = "base::Super-Duper", 
              .RNG.seed = c(569))

list.chains = list(chain1,chain2,chain3,chain4,chain5,chain6)



#now run linear phi for each of the unbiased surveys
line.ipsos <- jags.parfit(cl, ipsos.dat, c("positiverate","gamma0","gamma1","sigmasq"), custommodel(mod.linear.phi),
                          n.chains=6,n.adapt = 200000,thin = 5, n.iter = 500000,inits = list.chains)

#check convergence
gelman.diag(line.ipsos)

#get point estimates and CIs
point.ipsos <- get.point.est(line.ipsos,"positiverate",type = "median")
CI.ipsos <- get.CI(line.ipsos,"positiverate")


line.household <- jags.parfit(cl, household.dat, c("positiverate","gamma0","gamma1","sigmasq"), custommodel(mod.linear.phi),
                              n.chains=6,n.adapt = 200000,thin = 5, n.iter = 500000,inits = list.chains)

gelman.diag(line.household)

#point estimate and CIs
point.household <- get.point.est(line.household,"positiverate",type = "median")
CI.household <- get.CI(line.household,"positiverate")


line.facebook <- jags.parfit(cl, facebook.dat, c("positiverate","gamma0","gamma1","sigmasq"), custommodel(mod.linear.phi),
                             n.chains=6,n.adapt = 200000,thin = 5, n.iter = 500000,inits = list.chains)

#point estimates and CIs
gelman.diag(line.facebook)
point.facebook <- get.point.est(line.facebook,"positiverate",type = "median")
CI.facebook <-  get.CI(line.facebook,"positiverate")



#run method for three types of models for phi

line.const <- jags.parfit(cl, data.list.extended, c("positiverate","gamma0","sigmasq"), custommodel(mod.const.phi),
                          n.chains=6,n.adapt = 250000,thin = 5, n.iter = 500000,inits = list.chains)

line.linear <- jags.parfit(cl, data.list.extended, c("positiverate","sigmasq","phi"), custommodel(mod.linear.phi),
                         n.chains=6,n.adapt = 250000,thin = 5, n.iter = 500000,inits = list.chains)


line.walk <- jags.parfit(cl, data.list.extended, c("positiverate","gamma","sigmasq"), custommodel(mod.walk.phi),
                         n.chains=6,n.adapt = 250000,thin = 5, n.iter = 500000,inits = list.chains)


#constant model does not converge
#gelman.diag(line.const)

#gelman.diag(line.linear) #notrun if phi is included --- phi is 1 for k = 1, thus the funciton won't work
#gelman.diag(line.walk) 


#save all point estimates
point.const <-  get.point.est(line.const,"positiverate",type = "median")
CI.const <-  get.CI(line.const,"positiverate")

point.linear <- get.point.est(line.linear,"positiverate",type = "median")
CI.linear <- get.CI(line.linear,"positiverate")


point.walk <- get.point.est(line.walk,"positiverate",type = "median")
CI.walk <- get.CI(line.walk,"positiverate")



#include indicator for respective dfs
fb_df$mode <- rep("facebook",length(fb_df$ymd))
chp_df$mode <- rep("household-pulse",length(chp_df$ymd))
ax_df$mode <- rep("axios-ipsos",length(ax_df$ymd))

#save all estimates to dfs
fb_df$est <- point.facebook
fb_df$CI_L <- CI.facebook$Lower
fb_df$CI_U <- CI.facebook$Upper

chp_df$est <- point.household
chp_df$CI_L <- CI.household$Lower
chp_df$CI_U <- CI.household$Upper


ax_df$est <- point.ipsos
ax_df$CI_L <- CI.ipsos$Lower
ax_df$CI_U <- CI.ipsos$Upper

cdc_df$est <- cdc_df$vax

cdc_df$vax_lb2 <- cdc_df$vax*0.95
cdc_df$vax_ub2 <- cdc_df$vax*1.05


#data frame to compare positiverate estimates
compare.method <- data.frame(Method = c(rep("const",48),rep("linear",48),rep("walk",48)), ests = c(point.const,point.linear,point.walk),
                             CI.L=c(CI.const$Lower,CI.linear$Lower,CI.walk$Lower),CI.U=c(CI.const$Upper,CI.linear$Upper,CI.walk$Upper),dates = ref.dates)


#Figure.
ggplot(data = compare.method,aes(x=dates,y = ests,color = Method)) + geom_point() + 
  geom_line() + geom_errorbar(aes(ymin = CI.L,ymax = CI.U),width = 0) + theme_bw() + labs(x = NULL,y = "% of the population with at least one vaccine",title = expression(paste("Posterior estimates of % vaccinated by model for ",phi)))+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 0.90, by = 0.1), expand = expansion(mult = c(0, 0.05))) + scale_x_date(date_labels = "%b '%y", breaks = "1 month")


#prefered method (random walk model for phi) to be plotted against all single survey estimates
method_df <- data.frame(ymd = ref.dates, est = point.walk, CI_L = CI.walk$Lower, CI_U=CI.walk$Upper)



#Main plot to compare (inference plot)
fb_df %>% ggplot(aes(x = ymd, y = est)) + 
  geom_ribbon(data = cdc_df, aes(ymin = vax_lb2, ymax = vax_ub2), alpha = 0.3, color = "grey50") + 
  geom_pointline(aes(x = ymd, y = est),color = "#4891dc") +  geom_errorbar(aes(ymin = CI_L, ymax = CI_U), color = "#4891dc", width = 0)+ geom_pointline(data = ax_df,aes(x=ymd, y = est), color = "#cf7a30") + 
  geom_errorbar(data = ax_df, aes(ymin = CI_L, ymax = CI_U), color = "#cf7a30", width = 0) + 
  geom_pointline(data = chp_df,aes(x = ymd, y = est), color = "#69913b")  + geom_errorbar(data = chp_df, aes(ymin = CI_L, ymax = CI_U), color = "#69913b", width = 0) +
  geom_pointline(data = method_df,aes(x = ymd, y = est),color = "magenta") + geom_errorbar(data = method_df,aes(ymin = CI_L,ymax = CI_U),color = 'magenta',width = 0)+
  scale_x_date(date_labels = "%b '%y", breaks = "1 month") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 0.90, by = 0.1), expand = expansion(mult = c(0, 0.05))) + 
  theme_bw() + labs(x = NULL, y = "% of US Adults with 1+ dose Vaccination",title = "Synthesis method in comparison to single survey estimates") + 
  annotate("text", x = as.Date("2021-10-20"), y = 0.68, size = 3, label = "Axios-Ipsos", color = "#cf7a30") + 
  annotate("text", x = as.Date("2021-08-20"), y = 0.63, size = 3, label = "Method", color = "magenta")+
  annotate("text", x = as.Date("2021-08-01"), y = 0.87, size = 3, label = "Delphi-Facebook CTIS", color = "#4891dc") + 
  annotate("text", x = as.Date("2021-07-01"), y = 0.77, size = 3, label = "Census Household Pulse", color = "#69913b", angle = 10) + 
  annotate("label", x = as.Date("2021-05-01"), y = 0.53, size = 3, label = "CDC 18+\n(Retroactively updated)", angle = 5, color = "grey30", fill = "grey90", alpha = 0.6, label.size= 0, hjust = 0)


#calculate gain here:
gain <- (CI.ipsos$Upper-CI.ipsos$Lower)/(CI.walk$Upper[!is.na(data.list.extended$Y[1,])]-CI.walk$Lower[!is.na(data.list.extended$Y[1,])])

mean(gain)*100 #154.1535
median(gain)*100 #142.3487


gain.barplot <- data.frame(date = fb_df$ymd[!is.na(data.list.extended$Y[1,])], ratio = gain) 

ggplot(gain.barplot,aes(x = date,y =gain))+ geom_bar(stat = 'identity')+   scale_x_date(date_labels = "%b '%y", breaks = "1 month")+
  theme_bw() + labs(x = "Date", y = "Axios-Ipsos CI width/Synthesis CI width",title = "Width of 95% Credible Interval of Axios-Ipsos compared to the synthesis method") + geom_hline(yintercept = c(mean(gain),median(gain)),colour = c("blue",'red')) +
  annotate("text", x = as.Date("2021-10-5"), y = 1.65, size = 3, label = "Mean", color = "blue") + 
  annotate("text", x = as.Date("2021-10-25"), y = 1.35, size = 3, label = "Median", color = "red")


#collect estimates for phis by model

linear.phi <- get.point.est(line.linear,"phi")
CI.linear.phi <- get.CI(line.linear,"phi")


#Filter out cases where phi = 1, JAGS will automatically return all phi estimates (including for the unbiased survey)
linear.phi.2 <- linear.phi[seq(2,(48*3)-1,by = 3)]
linear.phi.3 <- linear.phi[seq(3,(48*3),by = 3)]
#.3 referes to k = 3, delphi-facebook, (third row of data.list.extended)

CI.linear.phi.2 <- lapply(CI.linear.phi,function(x){x[seq(2,(48*3)-1,by = 3)]})
CI.linear.phi.3 <- lapply(CI.linear.phi,function(x){x[seq(3,(48*3),by = 3)]})


point.const <- get.point.est(line.const,"gamma0")
CI.const <- get.CI(line.const,"gamma0")

point.gamma <- get.point.est(line.walk,"gamma")
gamma.CI <- get.CI(line.walk,"gamma")

point.gamma2 <- point.gamma[seq(1,96,by = 2)]
point.gamma3 <- point.gamma[-seq(1,96,by = 2)]


CI.L_gamma2 <- gamma.CI$Lower[seq(1,96,by = 2)]
CI.U_gamma2 <- gamma.CI$Upper[seq(1,96,by = 2)]

CI.L_gamma3 <- gamma.CI$Lower[-seq(1,96,by = 2)]
CI.U_gamma3 <- gamma.CI$Upper[-seq(1,96,by = 2)]


#concatenate into one large dataset
phi.dat <- data.frame(Survey = c(rep("Household-Pulse",48),rep("Delphi-Facebook",48),rep("Household-Pulse",48),rep("Delphi-Facebook",48),rep("Household-Pulse",48),rep("Delphi-Facebook",48)),
                      Method = c(rep("const",48),rep("const",48),rep("linear",48),rep("linear",48),rep("walk",48),rep("walk",48)),
                      phi = c(exp(rep(point.const[1],48)),exp(rep(point.const[2],48)),linear.phi.2,linear.phi.3,exp(point.gamma2),exp(point.gamma3)),
                      CI.L=c(exp(rep(CI.const$Lower[1],48)),exp(rep(CI.const$Lower[2],48)),CI.linear.phi.2$Lower,CI.linear.phi.3$Lower,exp(CI.L_gamma2),exp(CI.L_gamma3)),
                      CI.U =c(exp(rep(CI.const$Upper[1],48)),exp(rep(CI.const$Upper[2],48)),CI.linear.phi.2$Upper,CI.linear.phi.3$Upper,exp(CI.U_gamma2),exp(CI.U_gamma3)),
                      t= c(ref.dates,ref.dates,ref.dates,ref.dates,ref.dates,ref.dates))


#plot of comparison of phi estimates over time by method.
ggplot(data = phi.dat,aes(x = as.Date(t), y = phi, color = Method, shape = Survey)) + geom_point()+ 
  geom_line() + facet_grid(Survey~.)+ geom_ribbon(aes(ymin = CI.L,ymax = CI.U),alpha =0.25) + theme_bw() +labs(y = expression(phi[kt]), x = NULL,title = expression(paste(phi[kt]," by method and survey")))  +scale_x_date(date_labels = "%b '%y", breaks = "1 month") 


#summary(line.walk)
stopCluster(cl)


