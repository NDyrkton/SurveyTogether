#original figure


#data with now casting at every point.
library(dplyr)
library(ggplot2)
library(dplyr)
library(lemon)
library(ggtext)
library(forecast)

data.list <- readRDS("data_extended.Rdata")

fb_df <- read.csv("Data/fb.csv")
cdc_df <- read.csv("Data/cdc.csv")
chp_df <- read.csv("Data/chp.csv")
ax_df <- read.csv("Data/ax.csv")
ax_n <- read.csv("Data/ax_n.csv")
ax_df <- ax_df[nrow(ax_df):1,]

cdc_df$vax_lb2 <- cdc_df$vax*0.95
cdc_df$vax_ub2 <- cdc_df$vax*1.05


fb_df$ymd <- as.Date(fb_df$ymd)
cdc_df$ymd <- as.Date(cdc_df$ymd)
chp_df$ymd <- as.Date(chp_df$ymd)
ax_df$ymd <- as.Date(ax_df$ymd)


#original graph, extends data from bradley et al. Includes 5% margin of error.

fb_df %>% ggplot(aes(x = ymd, y =vax)) +
  geom_ribbon(data = cdc_df, aes(ymin = vax_lb2, ymax = vax_ub2), alpha = 0.3, color = "grey50") + 
  geom_pointline(color = "#4891dc") + geom_pointline(data = ax_df,aes(x = ymd, y = vax), color = "#cf7a30") + 
  geom_errorbar(data = ax_df, aes(ymin = vax_lb, ymax = vax_ub), color = "#cf7a30", width = 0) + 
  geom_pointline(data = chp_df, color = "#69913b")  + geom_errorbar(data = chp_df,aes(ymin = vax-1.96*se, ymax = vax+1.96*se),color = "#69913b",width = 0)+
  scale_x_date(date_labels = "%b '%y", breaks = "1 month") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 0.90, by = 0.1), expand = expansion(mult = c(0, 0.05))) + 
  theme_bw() + labs(x = NULL, y = "% of US Adults with 1+ dose Vaccination") + 
  annotate("text", x = as.Date("2021-10-20"), y = 0.68, size = 3, label = "Axios-Ipsos", color = "#cf7a30") + 
  annotate("text", x = as.Date("2021-08-01"), y = 0.87, size = 3, label = "Delphi-Facebook CTIS", color = "#4891dc") + 
  annotate("text", x = as.Date("2021-07-01"), y = 0.77, size = 3, label = "Census Household Pulse", color = "#69913b", angle = 10) + 
  annotate("label", x = as.Date("2021-05-01"), y = 0.53, size = 3, label = "CDC 18+\n(Retroactively updated)", angle = 5, color = "grey30", fill = "grey90", alpha = 0.6, label.size= 0, hjust = 0)  + ggtitle("Large survey vaccine paradox (Bradley et al., 2021)")
