#View(results.plot.final)


MSE <- apply(results.plot.final[[1]],2,function(x){mean(x^2)})

calculate.MCSE <- function(MSE,dataset){
  
  MCSE <- numeric(ncol(dataset))
  names(MCSE) <- colnames(dataset)
  
  x <- dataset
  
  for(i in 1:ncol(x)){
    MCSE[i] <-   sqrt(sum((x^2-MSE[i])^2   /  (nrow(x)*(nrow(x)-1))))
  }
  
  return(MCSE)
  
}



apply(results.plot.final[[1]],2,length)



MCSE <- numeric(9)
names(MCSE) <- colnames(results.plot.final[[1]])

x <- results.plot.final[[1]]

for(i in 1:ncol(x)){
  
  MCSE[i] <-   sqrt(sum((x^2-MSE[i])^2   /  (nrow(x)*(nrow(x)-1))))
  
  
}

MSE1 <- 100*MSE
MCSE1 <- 100*MCSE


MSE.unb <- apply(results.plot.final[[2]],2,function(x){mean(x^2)})
MCSE.unb <- 100*calculate.MCSE(MSE.unb,results.plot.final[[2]])

MSE.unb <- MSE.unb*100





results.plot <- data.frame(data = c(rep(c("const"),3),rep(c("linear"),3), rep(c("walk"),3)), model = rep(c("const","linear","walk"),3),MSE = MSE1,MCSE= MCSE1)
results.plot.unb <- data.frame(data = c("const","linear","walk"), model = rep("unbiased",3),MSE = MSE.unb,MCSE=MCSE.unb)

final.plot <- rbind(results.plot,results.plot.unb)
View(final.plot)


ggplot(data = final.plot, aes(x = data, y = MSE, group = model,colour = model)) + geom_point() + geom_line(linewidth = 0.75) + theme_bw()  + labs(x = "Data generation", y = "Mean Squared Error x 100", title = paste("RMSE of 1000 reptitions: 5, 10, and 15 timepoints")) +
  scale_color_manual(values = c("const"="red","linear" = "green",walk="blue",unbiased = "black")) + geom_ribbon(aes(ymin = MSE-1.645*MCSE,ymax = MSE+1.645*MCSE,alpha = 0.1))

###
