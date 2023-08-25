#redo on extended data.
library(dclone)
library(ggplot2)
library(dplyr)
library(lemon)
library(ggtext)
library(forecast)
library(rjags)


fill.dates <- function(ref.date, vec.date, vec){
  
  return.vec <- numeric(length(ref.date))
  
  for(i in 1:length(vec.date)){
    
    tmp <- which( (vec.date[i]>= ref.date) & (vec.date[i] <= ref.date + 6))
    
    
    return.vec[tmp] <- vec[i] 
    
  }
  
  return.vec[which(return.vec==0)] <- NA
  return(return.vec)
  
}



source("Scripts/JagsMods.R")
source("Scripts/helperfunctions.R")


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


saveRDS(data.list.extended,"data.list.extended.Rdata")





#START with full data

#get unbiased data with no NAs
ipsos.dat <- extract.unbiased.nona(data.list.extended,col = 1)
household.dat <- extract.unbiased.nona(data.list.extended,col= 2)
facebook.dat <- extract.unbiased.nona(data.list.extended,col = 3)


cl <- makePSOCKcluster(4)

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


#now run linear phi for each of the unbiased surveys


line.ipsos <- jags.parfit(cl, ipsos.dat, c("positiverate","gamma0","gamma1","sigmasq"), custommodel(mod.linear.phi),
                          n.chains=4,n.adapt = 200000,thin = 5, n.iter = 50000,inits = list(chain1,chain2,chain3,chain4))

#check convergence
gelman.diag(line.ipsos)

point.ipsos <- get.point.est(line.ipsos,"positiverate")
CI.ipsos <- get.CI(line.ipsos,"positiverate")


line.household <- jags.parfit(cl, household.dat, c("positiverate","gamma0","gamma1","sigmasq"), custommodel(mod.linear.phi),
                              n.chains=4,n.adapt = 200000,thin = 5, n.iter = 50000,inits = list(chain1,chain2,chain3,chain4))

gelman.diag(line.household)

point.household <- get.point.est(line.household,"positiverate")
CI.household <- get.CI(line.household,"positiverate")


line.facebook <- jags.parfit(cl, facebook.dat, c("positiverate","gamma0","gamma1","sigmasq"), custommodel(mod.linear.phi),
                             n.chains=4,n.adapt = 200000,thin = 5, n.iter = 50000,inits = list(chain1,chain2,chain3,chain4))


gelman.diag(line.facebook)
point.facebook <- get.point.est(line.facebook,"positiverate")
CI.facebook <-  get.CI(line.facebook,"positiverate")



#run actual method (linear)

line.linear <- jags.parfit(cl, data.list.extended, c("positiverate","sigmasq","phi"), custommodel(mod.linear.phi),
                         n.chains=4,n.adapt = 1000000,thin = 5, n.iter = 700000,inits = list(chain1,chain2,chain3,chain4))


line.walk <- jags.parfit(cl, data.list.extended, c("positiverate","gamma","sigmasq"), custommodel(mod.walk.phi),
                         n.chains=4,n.adapt = 800000,thin = 5, n.iter = 700000,inits = list(chain1,chain2,chain3,chain4))


line.const <- jags.parfit(cl, data.list.extended, c("positiverate","gamma0","sigmasq"), custommodel(mod.const.phi),
                         n.chains=4,n.adapt = 1000000,thin = 5, n.iter = 700000,inits = list(chain1,chain2,chain3,chain4))


gelman.diag(line.linear)
gelman.diag(line.walk) 

#constant model seems like it is mis-specified, does not converge.
gelman.diag(line.const)

point.const <-  get.point.est(line.const,"positiverate")
CI.const <-  get.CI(line.const,"positiverate")

point.linear <- get.point.est(line.linear,"positiverate")
CI.linear <- get.CI(line.linear,"positiverate")


point.walk <- get.point.est(line.walk,"positiverate")
CI.walk <- get.CI(line.walk,"positiverate")


fb_df$mode <- rep("facebook",length(fb_df$ymd))
chp_df$mode <- rep("household-pulse",length(chp_df$ymd))
ax_df$mode <- rep("axios-ipsos",length(ax_df$ymd))

fb_df$CI_L <- CI.facebook$Lower
fb_df$CI_U <- CI.facebook$Upper

chp_df$CI_L <- CI.household$Lower
chp_df$CI_U <- CI.household$Upper

ax_df$CI_L <- CI.ipsos$Lower
ax_df$CI_U <- CI.ipsos$Upper


compare.method <- data.frame(method = c(rep("const",48),rep("linear",48),rep("walk",48)), ests = c(point.const,point.linear,point.walk),
                             CI.L=c(CI.const$Lower,CI.linear$Lower,CI.walk$Lower),CI.U=c(CI.const$Upper,CI.linear$Upper,CI.walk$Upper),dates = ref.dates)

ggplot(data = compare.method,aes(x=dates,y = ests,color = method)) + geom_point() + 
  geom_line() + geom_errorbar(aes(ymin = CI.L,ymax = CI.U),width = 1) + theme_minimal() + labs(x = 'Date',y = "% at least one vaccine",title = "Comparison of method by assumption on phi")


method_df <- data.frame(ymd = ref.dates, vax = point.walk, CI_L = CI.walk$Lower, CI_U=CI.walk$Upper)





fb_df %>% ggplot(aes(x = ymd, y = vax)) + 
  geom_ribbon(data = cdc_df, aes(ymin = vax_lb, ymax = vax_ub), alpha = 0.3, color = "grey50") + 
  geom_pointline(color = "#4891dc") +  geom_errorbar(aes(ymin = CI_L, ymax = CI_U), color = "#4891dc", width = 0)+ geom_pointline(data = ax_df, color = "#cf7a30") + 
  geom_errorbar(data = ax_df, aes(ymin = vax_lb, ymax = vax_ub), color = "#cf7a30", width = 0) + 
  geom_pointline(data = chp_df, color = "#69913b")  + geom_errorbar(data = chp_df, aes(ymin = CI_L, ymax = CI_U), color = "#69913b", width = 0) +
  geom_pointline(data = method_df,aes(x = ymd, y = vax),color = "magenta") + geom_errorbar(data = method_df,aes(ymin = CI_L,ymax = CI_U),color = 'magenta',width = 0)+
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


linear.phi <- get.point.est(line.linear,"phi")
CI.linear.phi <- get.CI(line.linear,"phi")

linear.phi.2 <- linear.phi[seq(2,(48*3)-1,by = 3)]
linear.phi.3 <- linear.phi[seq(3,(48*3),by = 3)]

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

phi.dat <- data.frame(survey = c(rep("household-pulse",48),rep("facebook",48),rep("household-pulse",48),rep("facebook",48),rep("household-pulse",48),rep("facebook",48)),
                      method = c(rep("const",48),rep("const",48),rep("linear",48),rep("linear",48),rep("walk",48),rep("walk",48)),
                      phi = c(exp(rep(point.const[1],48)),exp(rep(point.const[2],48)),linear.phi.2,linear.phi.3,exp(point.gamma2),exp(point.gamma3)),
                      CI.L=c(exp(rep(CI.const$Lower[1],48)),exp(rep(CI.const$Lower[2],48)),CI.linear.phi.2$Lower,CI.linear.phi.3$Lower,exp(CI.L_gamma2),exp(CI.L_gamma3)),
                      CI.U =c(exp(rep(CI.const$Upper[1],48)),exp(rep(CI.const$Upper[2],48)),CI.linear.phi.2$Upper,CI.linear.phi.3$Upper,exp(CI.U_gamma2),exp(CI.U_gamma3)),
                      t= c(1:48,1:48,1:48,1:48,1:48,1:48))

ggplot(data = phi.dat,aes(x = t, y = phi, color = method, shape = survey)) + geom_point()+ 
  geom_line() + facet_grid(survey~.)+ geom_ribbon(aes(ymin = CI.L,ymax = CI.U),alpha =0.25) + theme_minimal() +labs(y = expression(phi), x = "timepoints",title = "Plot of phi by method and survey, extended data")




ggplot(data = phi.dat,aes(x = t, y = phi, color = method, shape = survey)) + geom_point()+ 
  geom_line() + geom_ribbon(aes(ymin = CI.L,ymax = CI.U),alpha =0.25) + theme_minimal() +labs(y = expression(phi), x = "timepoints",title = "Plot of phi by method and survey, extended data")


#try with random walk ,same jumping?



