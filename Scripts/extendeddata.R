#redo on extended data.
library(dclone)
library(ggplot2)
library(dplyr)
library(lemon)
library(ggtext)


source("Scripts/JagsMods.R")
source("Scripts/helperfunctions.R")


fb_df <- read.csv("Data/fb.csv")
cdc_df <- read.csv("Data/cdc.csv")
chp_df <- read.csv("Data/chp.csv")
ax_df <- read.csv("Data/ax.csv")

fb_df$ymd <- as.Date(fb_df$ymd)
cdc_df$ymd <- as.Date(cdc_df$ymd)
chp_df$ymd <- as.Date(chp_df$ymd)
ax_df$ymd <- as.Date(ax_df$ymd)



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


#fist get data into grouping by weeks


fill.dates <- function(ref.date, vec.date, vec){
  
  return.vec <- numeric(length(ref.date))
  
  for(i in 1:length(vec.date)){
    
    tmp <- which( (vec.date[i]>= ref.date) & (vec.date[i] <= ref.date + 6))
    
    
    return.vec[tmp] <- vec[i] 
    
  }
  
  return.vec[which(return.vec==0)] <- NA
  return(return.vec)
  
}

#facebook reference dates
ref.dates <- fb_df$ymd

chp.fill <- fill.dates(ref.dates,chp_df$ymd,vec = chp_df$vax)








