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


###including MCSE error bars NN = 2,000 
t5 <- read.csv("Data/t5_sims_MCSE_2000.csv")
t10 <- read.csv("Data/t10_sims_MCSE_2000.csv")
t15 <- read.csv("Data/t15_sims_MCSE_2000.csv")

t5$time <- rep("5 time-points",12)
t10$time <- rep("10 time-points",12)
t15$time <- rep("15 time-points",12)


all.time <- rbind(t5,t10,t15)

all.time$time <- factor(all.time$time,levels = c("5 time-points","10 time-points","15 time-points"))

all.time$MSE <- all.time$MSE/100
all.time$MCSE <- all.time$MCSE/100


ggplot(data = all.time, aes(x = data, y = MSE, group = model,colour = model)) + geom_point() + geom_line(linewidth = 0.75) + theme_bw() +
  facet_grid(.~time) + labs(x = "Data generation", y = "Mean Squared Error", title = paste("MSE of 2000 reptitions: 5, 10, and 15 timepoints")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black")) + geom_ribbon(aes(ymin = MSE-1.96*MCSE,ymax = MSE+1.96*MCSE),alpha =0.2)





