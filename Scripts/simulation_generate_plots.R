#create facet grid for 3 timepoint simulations
library(dplyr)
library(ggplot2)

calculate.MCSE <- function(MSE,dataset){
  
  MCSE <- numeric(ncol(dataset))
  names(MCSE) <- colnames(dataset)
  
  x <- dataset
  
  for(i in 1:ncol(x)){
    MCSE[i] <-   sqrt(    sum(((x[,i]^2)   -MSE[i])^2   )   /  (nrow(x)*(nrow(x)-1)))     
  }
  
  return(MCSE)
  
}



#########using N = 10.000,000

t5.error <- read.csv("Data/t5_sims_method.csv")
t5.unb <- read.csv("Data/t5_sims_only_unbiased.csv")

t10.error <- read.csv("Data/t10_sims_method.csv")
t10.unb <- read.csv("Data/t10_sims_only_unbiased.csv")

t15.error <- read.csv("Data/t15_sims_method.csv")
t15.unb<- read.csv("Data/t15_sims_only_unbiased.csv")

#the 3x3 MSE
t5.MSE <- apply(t5.error[1:1000,],2,function(x){mean(x^2)})
#unb MSE
t5.MSE.unb <- apply(t5.unb[1:1000,],2,function(x){mean(x^2)})

t5.MCSE <- calculate.MCSE(t5.MSE,t5.error[1:1000,])

t5.MCSE.unb <- calculate.MCSE(t5.MSE.unb,t5.unbiased[1:1000,])

t5.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const","linear","walk"),3),MSE = 100*t5.MSE,MCSE= 100*t5.MCSE)
t5.plot.unb <- data.frame(data = c("const","linear","walk"), model = rep("unbiased",3),MSE = 100*t5.MSE.unb,MCSE=100*t5.MCSE.unb)
t5.final <- rbind(t5.plot,t5.plot.unb)


#t10
t10.MSE <- apply(t10.error,2,function(x){mean(x^2)})
#unb MSE
t10.MSE.unb <- apply(t10.unb,2,function(x){mean(x^2)})
t10.MCSE <- calculate.MCSE(t10.MSE,t10.error)
t10.MCSE.unb <- calculate.MCSE(t10.MSE.unb,t10.unbiased)

t10.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const","linear","walk"),3),MSE = 100*t10.MSE,MCSE= 100*t10.MCSE)
t10.plot.unb <- data.frame(data = c("const","linear","walk"), model = rep("unbiased",3),MSE = 100*t10.MSE.unb,MCSE=100*t10.MCSE.unb)
t10.final <- rbind(t5.plot,t5.plot.unb)




#t10

t15.MSE <- apply(t15.error,2,function(x){mean(x^2)})
#unb MSE
t15.MSE.unb <- apply(t15.unb,2,function(x){mean(x^2)})
t15.MCSE <- calculate.MCSE(t15.MSE,t15.error)
t15.MCSE.unb <- calculate.MCSE(t15.MSE.unb,t15.unbiased)


t15.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const","linear","walk"),3),MSE = 100*t15.MSE, MCSE= 100*t15.MCSE)
t15.plot.unb <- data.frame(data = c("const","linear","walk"), model = rep("unbiased",3),MSE = 100*t15.MSE.unb,MCSE=100*t15.MCSE.unb)
t15.final <- rbind(t5.plot,t5.plot.unb)










####skip
t5 <- read.csv("Data/RMSE_plot_t5_1000.csv")
t10 <- read.csv("Data/RMSE_plot_t10_1000.csv")
t15 <- read.csv("Data/RMSE_plot_t15_1000.csv")

t5$time <- rep("5 timepoints",12)
t10$time <- rep("10 timepoints",12)
t15$time <- rep("15 timepoints",12)






all.time <- rbind(t5,t10,t15)

all.time$time <- factor(all.time$time,levels = c("5 timepoints","10 timepoints","15 timepoints"))

ggplot(data = all.time, aes(x = data, y = RMSE, group = model,colour = model)) + geom_point() + geom_line(linewidth = 0.75) + theme_bw() +
  facet_grid(.~time) + labs(x = "Data generation", y = "Root Mean Squared Error x 100", title = paste("RMSE of 1000 reptitions: 5, 10, and 15 timepoints")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black"))

#########using N = 100,000




CI.check <- read.csv("Data/CheckCIT10N1000.csv")

model <- CI.check[1:2000,]
unbiased <- CI.check[2001:4000,]

model %>% filter(const.95.contain==0) %>% View()

apply(100*model,2,mean)
apply(100*unbiased,2,mean)


##redone with N = 10,000,000, and new binomial approximation



t5 <- read.csv("Data/RMSE_1000_large_N_binom2.csv")
t10 <- read.csv("Data/RMSE_t10_1000_large_N_binom2.csv")
t15 <- read.csv("Data/RMSE_t15_1000_large_N_binom2.csv")

t5$time <- rep("5 time-points",12)
t10$time <- rep("10 time-points",12)
t15$time <- rep("15 time-points",12)



all.time <- rbind(t5,t10,t15)

all.time$time <- factor(all.time$time,levels = c("5 time-points","10 time-points","15 time-points"))

ggplot(data = all.time, aes(x = data, y = RMSE, group = model,colour = model)) + geom_point() + geom_line(linewidth = 0.75) + theme_bw() +
  facet_grid(.~time) + labs(x = "Data generation", y = "Root Mean Squared Error x 100", title = paste("RMSE of 1000 reptitions: 5, 10, and 15 time-points")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black"))

#########using N = 100,000


###including MCSE

t5.final$time <- rep("5 timepoints",12)
t10.final$time <- rep("10 timepoints",12)
t15.final$time <- rep("15 timepoints",12)




all.time <- rbind(t5.final,t10.final,t15.final)

all.time$time <- factor(all.time$time,levels = c("5 timepoints","10 timepoints","15 timepoints"))


ggplot(data = all.time, aes(x = data, y = MSE, group = model,colour = model)) + geom_point() + geom_line(linewidth = 0.75) + theme_bw() +
  facet_grid(.~time) + labs(x = "Data generation", y = "Mean Squared Error x 100", title = paste("MSE of 1000 reptitions: 5, 10, and 15 timepoints")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black")) + geom_ribbon(aes(ymin = MSE-1.96*MCSE,ymax = MSE+1.96*MCSE),alpha =0.2)


