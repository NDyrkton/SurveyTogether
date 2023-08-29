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
source("Scripts/helperfunctions.R")

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




ipsos.dat <- extract.unbiased.nona(data.list,col = 1)
household.dat <- extract.unbiased.nona(data.list,col= 2)
facebook.dat <- extract.unbiased.nona(data.list,col = 3)





cl <- makePSOCKcluster(4)

clusterEvalQ(cl, library(dclone))
load.module("lecuyer")
parLoadModule(cl,"lecuyer")


ipsos.len <- length(ipsos.dat$Y)
fb.len <- length(facebook.dat$Y)
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
dcoptions("verbose"=F)#mute the output of dlclone


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
  
  
  
  
  
  line.full <- jags.parfit(cl, data.t, c("positiverate"), custommodel(mod.walk.phi),
                           n.chains=4,n.adapt = 250000,thin = 5, n.iter = 200000
                           ,inits = inits.chains)
  
  if(any(gelman.diag(line.full)$psrf[,1] >= 1.1)){
    print("failed")
    
    line.full <- jags.parfit(cl, data.t, c("positiverate","gamma0"), custommodel(mod.walk.phi),
                             n.chains=4,n.adapt = 500000,thin = 5, n.iter = 500000,inits = inits.chains)
    
    gelman.diag(line.full)
    
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
  annotate("label", x = as.Date("2021-05-01"), y = 0.53, size = 3, label = "CDC 18+\n(Retroactively updated)", angle = 5, color = "grey30", fill = "grey90", alpha = 0.6, label.size= 0, hjust = 0) + 
  labs(caption = "<br>Figure extends Bradley, Kurirwaki, Isakov, Sejdinovic, Meng, and Flaxman,<br> \"**Unrepresentative big surveys significantly overestimated US vaccine uptake**\" (_Nature_, Dec 2021, doi:10.1038/s41586-021-04198-4).<br> Article analyzed the period Jan-May 2021 with retroactively updated CDC numbers as of May 2021.<br> This figure extends the series up to December, with CDC's same series as of Dec 2021, with bands for potential +/- 2% error in CDC.<br> **Axios-Ipsos** (n = 1000 or so per point) shows +/- 3.4% 95 percent MOE, which is their modal value for the poll.<br> **Delphi-Facebook** (n = 250,000 per point) and **Census Household Pulse** (n = 75,000 per point) not shown.") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(color = "black"), plot.caption = element_markdown(color = "grey40"))




#mean gain
mean(ax_df$vax_ub-ax_df$vax_lb)/mean(method_df$CI.U-method_df$CI.L)

# axios CI's are 60% larger



write.csv(fb_df,"Data/fb_df_new.csv",row.names = F)
write.csv(chp_df,"Data/chp_df_new.csv",row.names = F)
write.csv(ax_df,"Data/ax_df_new.csv",row.names = F)
write.csv(method_df,"Data/nowcast_rw.csv",row.names = F)







