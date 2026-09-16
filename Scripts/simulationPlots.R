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
t5 <- read.csv("Results/simulation_t5_summarised.csv")
t10 <- read.csv("Results/simulation_t10_summarised.csv")
t15 <- read.csv("Results/simulation_t15_summarised.csv")

t5$time <- rep("5 time-points",12)
t10$time <- rep("10 time-points",12)
t15$time <- rep("15 time-points",12)


all.time <- rbind(t5,t10,t15)

all.time$time <- factor(all.time$time,levels = c("5 time-points","10 time-points","15 time-points"))


bias.plot <- ggplot(data = all.time, aes(x = dgm, y = bias.total, group = model,colour = model)) + geom_point() + geom_line(linewidth = 0.75)  +
  facet_grid(.~time) + labs(x = "Data generation", y = "Bias", title = paste("Bias across 10,000 reptitions: 5, 10, and 15 time points")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black")) + geom_ribbon(aes(ymin = bias.total-1.96*mcse.bias,ymax = bias.total+1.96*mcse.bias),alpha =0.2)


mse.plot <- ggplot(data = all.time, aes(x = dgm, y = mse.total, group = model,colour = model)) + geom_point() + geom_line(linewidth = 0.75)  +
  facet_grid(.~time) + labs(x = "Data Generation", y = "MSE", title = paste("MSE across 10,000 reptitions: 5, 10, and 15 time points")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black")) + geom_ribbon(aes(ymin = mse.total-1.96*mcse.mse,ymax = mse.total+1.96*mcse.mse),alpha =0.2)

#save these plots
ggsave("Figures/simulation_bias_march2026.png",plot = bias.plot, width = 22, height = 14, unit = "cm")
ggsave("Figures/simulation_mse_march2026.png",plot = mse.plot, width = 22, height = 14, unit = "cm")


