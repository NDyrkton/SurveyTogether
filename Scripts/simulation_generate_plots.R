#create facet grid for 3 timepoint simulations
library(dplyr)
library(ggplot2)

t5 <- read.csv("Data/RMSE_plot_t5_500.csv")
t10 <- read.csv("Data/RMSE_plot_t10_500.csv")
t15 <- read.csv("Data/RMSE_plot_t15_500.csv")

t5$time <- rep("5 timepoints",12)
t10$time <- rep("10 timepoints",12)
t15$time <- rep("15 timepoints",12)



all.time <- rbind(t5,t10,t15)

all.time$time <- factor(all.time$time,levels = c("5 timepoints","10 timepoints","15 timepoints"))

ggplot(data = all.time, aes(x = data, y = RMSE, group = model,colour = model)) + geom_point() + geom_line(linewidth = 0.75) + theme_minimal() +
  facet_grid(.~time) + labs(x = "Data generation", y = "Root Mean Squared Error (percentage)", title = paste("RMSE of 500 reptitions: 5, 10, and 15 timepoints")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black"))
