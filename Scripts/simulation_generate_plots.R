#create facet grid for 3 timepoint simulations
library(dplyr)
library(ggplot2)


#########using N = 100,000


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


##redone with N = 10,000,000



t5 <- read.csv("Data/RMSE_1000_large_N.csv")
t10 <- read.csv("Data/RMSE_t10_1000_large_N.csv")
t15 <- read.csv("Data/RMSE_t15_1000_large_N.csv")

t5$time <- rep("5 timepoints",12)
t10$time <- rep("10 timepoints",12)
t15$time <- rep("15 timepoints",12)



all.time <- rbind(t5,t10,t15)

all.time$time <- factor(all.time$time,levels = c("5 timepoints","10 timepoints","15 timepoints"))

ggplot(data = all.time, aes(x = data, y = RMSE, group = model,colour = model)) + geom_point() + geom_line(linewidth = 0.75) + theme_bw() +
  facet_grid(.~time) + labs(x = "Data generation", y = "Root Mean Squared Error x 100", title = paste("RMSE of 1000 reptitions: 5, 10, and 15 timepoints")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black"))

#########using N = 100,000





