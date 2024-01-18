#data with now casting at every point.
library(dclone)
library(MCMCpack)
library(dplyr)
library(rjags)#redo on extended data.
library(ggplot2)
library(dplyr)
library(lemon)
library(ggtext)
library(forecast)



source("Scripts/JagsMods.R")
#source("Scripts/helperfunctions.R")

#import data
data.list <- readRDS("data_extended.Rdata")

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




ipsos.dat <- extract.unbiased.nona(data.list,col = 1)
household.dat <- extract.unbiased.nona(data.list,col= 2)
facebook.dat <- extract.unbiased.nona(data.list,col = 3)


cl <- makePSOCKcluster(8)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")

dcoptions("verbose"=F)#mute the output of dclone

ipsos.len <- length(ipsos.dat$Y)
fb.len <- length(facebook.dat$Y)
household.len <- length(household.dat$Y)


#now casts for each parameter to be estimated.
ipsos.posrates <- numeric(ipsos.len)
ipsos.CI <- list(CI.U=numeric(ipsos.len),CI.L = numeric(ipsos.len))

household.posrates <- numeric(household.len)
household.CI <- list(CI.U=numeric(household.len),CI.L = numeric(household.len))

facebook.posrates <- numeric(fb.len)
facebook.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))

#all three surveys
full.posrates <- numeric(fb.len)
full.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))

#axios + fb
full.fb.posrates <- numeric(fb.len)
full.fb.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))

#axios + hp
full.hp.posrates <- numeric(fb.len)
full.hp.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))

sigmasq <- numeric(fb.len)
sigmasq.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))

phi.household <- numeric(fb.len)
phi.household.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))

phi.facebook <- numeric(fb.len)
phi.facebook.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))

pisq.est <- numeric(fb.len)
pisq.CI <- list(CI.U=numeric(fb.len),CI.L = numeric(fb.len))


data.list.fb <- extract.surveys(data.list,c(1,3))
data.list.hp <- extract.surveys(data.list,c(1,2))

for(t in 1:fb.len){
  print(t)
  
  chain1<- list(.RNG.name = "base::Wichmann-Hill", 
                .RNG.seed = c(159+2*t))
  chain2<- list(.RNG.name = "base::Super-Duper", 
                .RNG.seed = c(260+3*t))
  chain3<- list(.RNG.name = "base::Wichmann-Hill", 
                .RNG.seed = c(371+4*t))
  chain4<- list(.RNG.name = "base::Super-Duper", 
                .RNG.seed = c(482+5*t))
  chain5<- list(.RNG.name = "base::Wichmann-Hill", 
                .RNG.seed = c(159+6*t))
  chain6<- list(.RNG.name = "base::Super-Duper", 
                .RNG.seed = c(260+7*t))
  chain7<- list(.RNG.name = "base::Wichmann-Hill", 
                .RNG.seed = c(371+8*t))
  chain8<- list(.RNG.name = "base::Super-Duper", 
                .RNG.seed = c(482+9*t))
  
  inits.chains <- list(chain1,chain2,chain3,chain4,chain5,chain6,chain7,chain8)
  

  data.t <- extract.t(data.list,t)
  data.list.fb.t <- extract.t(data.list.fb,t)
  data.list.hp.t <- extract.t(data.list.hp,t)

  

  ###ipsos run##
  if(t %in% c(1:ipsos.len)){

    ipsos.t <- extract.t(ipsos.dat,t)

    line.ipsos <- jags.parfit(cl,ipsos.t, c("positiverate"), custommodel(mod.linear.phi),
                              n.chains=8,n.adapt = 50000,thin = 5, n.iter = 100000,inits = inits.chains)

    ipsos.posrates[t] <- get.point.est(line.ipsos,"positiverate")[t]

    gelman.diag(line.ipsos)


    ipsos.CIs <- get.CI(line.ipsos,"positiverate")
    ipsos.CI$CI.U[t] <- ipsos.CIs$Upper[t]
    ipsos.CI$CI.L[t] <- ipsos.CIs$Lower[t]



  }

  if(t %in% c(1:household.len)){
    household.t <- extract.t(household.dat,t)

    line.household <- jags.parfit(cl,household.t, c("positiverate"), custommodel(mod.linear.phi),
                                  n.chains=8,n.adapt = 50000,thin = 5, n.iter = 100000,inits = inits.chains)

    household.posrates[t] <- get.point.est(line.household,"positiverate")[t]
    household.CIs <- get.CI(line.household,"positiverate")


    household.CI$CI.U[t] <- household.CIs$Upper[t]
    household.CI$CI.L[t] <- household.CIs$Lower[t]

  }


###facebook and full###

  facebook.t <- extract.t(facebook.dat,t)

  line.facebook <- jags.parfit(cl, facebook.t, c("positiverate"), custommodel(mod.linear.phi),
                             n.chains=8,n.adapt = 50000,thin = 5, n.iter = 100000,inits = inits.chains)

  facebook.posrates[t] <- get.point.est(line.facebook,"positiverate")[t]

  facebook.CIs <- get.CI(line.facebook,"positiverate")


  facebook.CI$CI.U[t] <- facebook.CIs$Upper[t]
  facebook.CI$CI.L[t] <- facebook.CIs$Lower[t]

  
  
  

  line.full <- jags.parfit(cl, data.t, c("positiverate","sigmasq","gamma","pisq"), custommodel(mod.walk.phi),
                           n.chains=8,n.adapt = 350000,thin = 5, n.iter = 300000
                           ,inits = inits.chains)
  line.fb.axios <- jags.parfit(cl, data.list.fb.t, c("positiverate"), custommodel(mod.walk.phi),
                               n.chains=8,n.adapt = 300000,thin = 5, n.iter = 250000
                               ,inits = inits.chains)
  line.hp.axios <- jags.parfit(cl, data.list.hp.t, c("positiverate"), custommodel(mod.walk.phi),
                               n.chains=8,n.adapt = 300000,thin = 5, n.iter = 250000
                               ,inits = inits.chains)
  
  if(any(gelman.diag(line.full)$psrf[,1] >= 1.1)){
    print("failed")

    line.full <- jags.parfit(cl, data.t, c("positiverate","sigmasq","gamma","pisq"), custommodel(mod.walk.phi),
                             n.chains=8,n.adapt = 500000,thin = 5, n.iter = 600000,inits = inits.chains)

    print(gelman.diag(line.full))

  }

  full.posrates[t] <- get.point.est(line.full,"positiverate")[t]

  full.CIs <- get.CI(line.full,"positiverate")
  full.CI$CI.U[t] <- full.CIs$Upper[t]
  full.CI$CI.L[t] <- full.CIs$Lower[t]

  
  #for added fb, and added hp
  full.fb.posrates[t] <- get.point.est(line.fb.axios,"positiverate")[t]
  full.fb.CIs <- get.CI(line.fb.axios, "positiverate")
  full.fb.CI$CI.L[t] <- full.fb.CIs$Lower[t]
  full.fb.CI$CI.U[t] <- full.fb.CIs$Upper[t]
  
  
  #axios + hp
  full.hp.posrates[t] <- get.point.est(line.hp.axios,"positiverate")[t]
  full.hp.CIs <- get.CI(line.hp.axios,"positiverate")
  full.hp.CI$CI.L[t] <- full.hp.CIs$Lower[t]
  full.hp.CI$CI.U[t] <- full.hp.CIs$Upper[t]
  
  
  
 # get sigmasq ests
  sigmasq[t] <- get.point.est(line.full,'sigmasq')
  sigmasq.CIs <-  get.CI(line.full,'sigmasq')

  sigmasq.CI$CI.L[t] <- sigmasq.CIs$Lower
  sigmasq.CI$CI.U[t] <- sigmasq.CIs$Upper

  #phi ests
  phi.all <- tail(exp( get.point.est(line.full,'gamma')),2)
  phi.household[t] <- phi.all[1]
  phi.facebook[t] <- phi.all[2]



  phi.CIs <- lapply(get.CI(line.full,'gamma'),function(x){tail(exp(x),2)})

  phi.household.CI$CI.L[t] <- phi.CIs$Lower[1]
  phi.household.CI$CI.U[t] <- phi.CIs$Upper[1]

  phi.facebook.CI$CI.L[t] <- phi.CIs$Lower[2]
  phi.facebook.CI$CI.U[t] <- phi.CIs$Upper[2]
  #pi ests

  pisq.est[t] <- get.point.est(line.full,'pisq')
  pisq.CIs <- get.CI(line.full,'pisq')

  pisq.CI$CI.L[t] <- pisq.CIs$Lower
  pisq.CI$CI.U[t] <- pisq.CIs$Upper


}

#save results in data frame.
results.nowcast.full <- list(posrate = list(est = full.posrates,CI.L = full.CI$CI.L,CI.U = full.CI$CI.U),
                             sigmasq = list(est = sigmasq,CI.L = sigmasq.CI$CI.L,CI.U = sigmasq.CI$CI.U),
                             pisq = list(est = pisq.est,CI.L = pisq.CI$CI.L,CI.U = pisq.CI$CI.U),
                             phi = list(est.household = phi.household,est.facebook = phi.facebook, CI.household.L = phi.household.CI$CI.L, CI.household.U = phi.household.CI$CI.U, CI.facebook.L = phi.facebook.CI$CI.L, CI.facebook.U = phi.facebook.CI$CI.U))




fb_df$point.est <- facebook.posrates
fb_df$CI.L <- facebook.CI$CI.L
fb_df$CI.U <- facebook.CI$CI.U

chp_df$point.est <- household.posrates
chp_df$CI.L <- household.CI$CI.L
chp_df$CI.U <- household.CI$CI.U

ax_df$point.est <- ipsos.posrates
ax_df$CI.L <- ipsos.CI$CI.L
ax_df$CI.U <- ipsos.CI$CI.U

cdc_df$point.est <- cdc_df$vax



ref.dates <- fb_df$ymd




method_df <- data.frame(ymd = ref.dates, point.est = full.posrates, CI.L = full.CI$CI.L, CI.U=full.CI$CI.U)




fb_df %>% ggplot(aes(x = ymd, y = point.est)) + geom_errorbar(aes(ymin = CI.L,ymax = CI.U),width = 0)+
  geom_ribbon(data = cdc_df, aes(ymin = vax_lb, ymax = vax_ub), alpha = 0.3, color = "grey50") + 
  geom_pointline(color = "#4891dc") + geom_pointline(data = ax_df,aes(x = ymd, y = point.est), color = "#cf7a30") + 
  geom_errorbar(data = ax_df, aes(ymin = CI.L, ymax = CI.U), color = "#cf7a30", width = 0) + 
  geom_pointline(data = chp_df, color = "#69913b")  + geom_errorbar(data = chp_df,aes(ymin = CI.L, ymax = CI.U),width = 0)+
  geom_pointline(data = method_df,aes(x = ymd, y = point.est),color = "magenta") + geom_errorbar(data = method_df,aes(ymin = CI.L,ymax = CI.U),color = 'magenta',width = 0)+
  scale_x_date(date_labels = "%b '%y", breaks = "1 month") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 0.90, by = 0.1), expand = expansion(mult = c(0, 0.05))) + 
  theme_bw() + labs(x = NULL, y = "% of US Adults with 1+ dose Vaccination") + 
  annotate("text", x = as.Date("2021-10-20"), y = 0.68, size = 3, label = "Axios-Ipsos", color = "#cf7a30") + 
  annotate("text", x = as.Date("2021-08-20"), y = 0.63, size = 3, label = "Method", color = "magenta")+
  annotate("text", x = as.Date("2021-08-01"), y = 0.87, size = 3, label = "Delphi-Facebook CTIS", color = "#4891dc") + 
  annotate("text", x = as.Date("2021-07-01"), y = 0.77, size = 3, label = "Census Household Pulse", color = "#69913b", angle = 10) + 
  annotate("label", x = as.Date("2021-05-01"), y = 0.53, size = 3, label = "CDC 18+\n(Retroactively updated)", angle = 5, color = "grey30", fill = "grey90", alpha = 0.6, label.size= 0, hjust = 0)  + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(color = "black"), plot.caption = element_markdown(color = "grey40"))+ ggtitle("Now-cast performance of the synthesis method")

#Figure extends Bradley, Kurirwaki, Isakov, Sejdinovic, Meng, and Flaxman,<br> \"**Unrepresentative big surveys significantly overestimated US vaccine uptake**\" (_Nature_, Dec 2021, doi:10.1038/s41586-021-04198-4).<br> Article analyzed the period Jan-May 2021 with retroactively updated CDC numbers as of May 2021.<br> This figure extends the series up to December, with CDC's same series as of Dec 2021, with bands for potential +/- 2% error in CDC.<br> **Axios-Ipsos** (n = 1000 or so per point) shows +/- 3.4% 95 percent MOE, which is their modal value for the poll.<br> **Delphi-Facebook** (n = 250,000 per point) and **Census Household Pulse** (n = 75,000 per point) not shown.




sigmasq.df <- data.frame(ymd = ref.dates, point.est = sigmasq, CI.L = sigmasq.CI$CI.L,CI.U = sigmasq.CI$CI.U)
phi.df <- data.frame(ymd = c(ref.dates,ref.dates), survey = c(rep("household",48),rep("facebook",48)),
                     point.est = c(phi.household,phi.facebook),CI.L = c(phi.household.CI$CI.L,phi.facebook.CI$CI.L), CI.U = c(phi.household.CI$CI.U,phi.facebook.CI$CI.U))


ggplot(sigmasq.df, aes(x = ymd, y = point.est)) + geom_point() + geom_line() + geom_ribbon(aes(ymin = CI.L,ymax = CI.U),alpha = 0.2) + 
  theme_bw() + labs(x = "date", y = expression(paste("Estimates of ",sigma^2)),title = expression(paste("Estimates of ",sigma^2, "over time"))) 


ggplot(phi.df,aes(x = ymd, y= point.est, color = survey)) + geom_line() + geom_point() + geom_ribbon(aes(ymin =CI.L,ymax = CI.U),alpha = 0.2)+
  theme_bw() + ylim(0,2)

#mean gain
gain <- (ax_df$CI.U-ax_df$CI.L) / (method_df$CI.U[!is.na(data.list$Y[1,])]-method_df$CI.L[!is.na(data.list$Y[1,])])

mean(gain)
median(gain)

#79.9% method CI are 79.9% the size of axios alone


#get graph for pisq

pisq.df <- data.frame(ymd = ref.dates, point.est = pisq.est, CI.L = pisq.CI$CI.L,CI.U = pisq.CI$CI.U)
ggplot(pisq.df, aes(x = ymd, y = point.est)) + geom_point() + geom_line() + geom_ribbon(aes(ymin = CI.L,ymax = CI.U),alpha = 0.2) + 
  theme_bw() + labs(x = "date", y = expression(paste("Estimates of ",pi^2)),title = expression(paste("Estimates of ",pi^2, "over time"))) 




####################
write.csv(fb_df,"Data/fb_df_new.csv",row.names = F)
write.csv(chp_df,"Data/chp_df_new.csv",row.names = F)
write.csv(ax_df,"Data/ax_df_new.csv",row.names = F)
write.csv(method_df,"Data/nowcast_rw.csv",row.names = F)



#how much n gain from adding FB survey?

mean((full.fb.CI$CI.U[!is.na(data.list$Y[1,])]-full.fb.CI$CI.L[!is.na(data.list$Y[1,])])/(ax_df$CI.U-ax_df$CI.L))

mean((full.hp.CI$CI.U[!is.na(data.list$Y[1,])]-full.hp.CI$CI.L[!is.na(data.list$Y[1,])])/(ax_df$CI.U-ax_df$CI.L))

#about a 3.5 % gain on average
#CI now about 0.8554 the size of axios (FB)


#CI 0.9357 the size (hp)

#phat +_ Z_0.025*sqrt((phat(1-phat))/n)


#this means MOE is is also 96.78% the size.

reduction.fb <- (full.fb.CI$CI.U[!is.na(data.list$Y[1,])]-full.fb.CI$CI.L[!is.na(data.list$Y[1,])])/(ax_df$CI.U-ax_df$CI.L)
reduction.hp <- (full.hp.CI$CI.U[!is.na(data.list$Y[1,])]-full.hp.CI$CI.L[!is.na(data.list$Y[1,])])/(ax_df$CI.U-ax_df$CI.L)
reduction.full <- (full.CI$CI.U[!is.na(data.list$Y[1,])]-full.CI$CI.L[!is.na(data.list$Y[1,])])/(ax_df$CI.U-ax_df$CI.L)

#phat

MOE.ax <- (ax_df$CI.U-ax_df$CI.L)/2

phat <- ax_df$point.est
n.old <- (1.96^2 * phat*(1-phat))/MOE.ax^2
n.new <- (1.96^2 * phat*(1-phat))/(MOE.ax*reduction.fb)^2

mean(n.new-n.old) # gain of n = 609
median(n.new-n.old) 

n.new.hp <- (1.96^2 * phat*(1-phat))/(MOE.ax*reduction.hp)^2

mean(n.new.hp-n.old) # n = 245 gain
median(n.new.hp-n.old)
#combined

n.new.full <- (1.96^2 * phat*(1-phat))/(MOE.ax*reduction.full)^2


mean(n.new.full-n.old) # 924 gain

median(n.new.full-n.old)


n.df <- data.frame(n= c(n.new-n.old,n.new.hp-n.old,n.new.full-n.old), surveys = c(rep("Delphi-Facebook",23),rep("Household-Pulse",23),rep("Both Surveys",23)), date = c(ax_df$ymd,ax_df$ymd,ax_df$ymd))

ggplot(n.df,aes(x = date,y = n, fill = surveys)) + geom_bar(position='dodge',stat = "identity",col= 'black') + theme_bw() + 
  labs(x ="Date",y = "Number of iid samples gained compared to Ipsos-Axios",title = "Barplot of number of iid samples gained when including the biased surveys",fill = "Survey")+scale_x_date(date_labels = "%b '%y", breaks = "1 month") 

