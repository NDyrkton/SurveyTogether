#redo on extended data.
library(dclone)
library(ggplot2)
library(dplyr)
library(lemon)
library(ggtext)
library(forecast)


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

times <- matrix(1:length(ref.dates),nrow = 3)
data.list



data.list.extended <- list(K = 3, T = 48, N = 255200373, times =times, smalln = n)








fb_df %>% ggplot(aes(x = ymd, y = vax)) + 
  geom_ribbon(data = cdc_df, aes(ymin = vax_lb, ymax = vax_ub), alpha = 0.3, color = "grey50") + 
  geom_pointline(color = "#4891dc") + geom_pointline(data = ax_df, color = "#cf7a30") + 
  geom_errorbar(data = ax_df, aes(ymin = vax_lb, ymax = vax_ub), color = "#cf7a30", width = 0) + 
  geom_pointline(data = chp_df, color = "#69913b") + scale_x_date(date_labels = "%b '%y", breaks = "1 month") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 0.90, by = 0.1), expand = expansion(mult = c(0, 0.05))) + 
  theme_bw() + labs(x = NULL, y = "% of US Adults with 1+ dose Vaccination") + 
  annotate("text", x = as.Date("2021-10-20"), y = 0.68, size = 3, label = "Axios-Ipsos", color = "#cf7a30") + 
  annotate("text", x = as.Date("2021-08-01"), y = 0.87, size = 3, label = "Delphi-Facebook CTIS", color = "#4891dc") + 
  annotate("text", x = as.Date("2021-07-01"), y = 0.77, size = 3, label = "Census Household Pulse", color = "#69913b", angle = 10) + 
  annotate("label", x = as.Date("2021-05-01"), y = 0.53, size = 3, label = "CDC 18+\n(Retroactively updated)", angle = 5, color = "grey30", fill = "grey90", alpha = 0.6, label.size= 0, hjust = 0) + 
  labs(caption = "<br>Figure extends Bradley, Kurirwaki, Isakov, Sejdinovic, Meng, and Flaxman,<br> \"**Unrepresentative big surveys significantly overestimated US vaccine uptake**\" (_Nature_, Dec 2021, doi:10.1038/s41586-021-04198-4).<br> Article analyzed the period Jan-May 2021 with retroactively updated CDC numbers as of May 2021.<br> This figure extends the series up to December, with CDC's same series as of Dec 2021, with bands for potential +/- 2% error in CDC.<br> **Axios-Ipsos** (n = 1000 or so per point) shows +/- 3.4% 95 percent MOE, which is their modal value for the poll.<br> **Delphi-Facebook** (n = 250,000 per point) and **Census Household Pulse** (n = 75,000 per point) not shown.") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.text = element_text(color = "black"), plot.caption = element_markdown(color = "grey40"))








