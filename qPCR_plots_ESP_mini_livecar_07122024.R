#set wd

#load libraries
library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)

#Edit raw qPCR output csv file
  #add column for sample ID
  #Change time stamps (UTC to PST)

#load edited csv file
qpcr_data <- read.csv("edited-2024-02-29 mini live car ESP 113119.csv")
head(qpcr_data)

#basic plot
ggplot(data=qpcr_data, aes(x = Sample, y = CT)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#boxplot
str(qpcr_data)

ggplot(data=qpcr_data, aes(x = Sample, y = CT)) +
  geom_boxplot(aes(color = Sample)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


#flip Y axis values, adjust Y axis range, group by technical replicate. 
qpcr_data$Sample <- factor(qpcr_data$Sample, levels = c("pos control", "11:07:19 AM", "12:07:06 PM", "1:07:06 PM", "3:07:06 PM", "7:07:05 PM", "1:07:07 AM", "1:07:05 PM", "NTC", "ext blank"))
str(qpcr_data)
ggplot(qpcr_data, aes(x = Sample, y = CT)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_y_reverse(limits=c(40, 15))

#add shaded area showing when the live car was deployed
qpcr_data$Sample <- factor(qpcr_data$Sample, levels = c("pos control", "11:07:19 AM", "12:07:06 PM", "1:07:06 PM", "3:07:06 PM", "7:07:05 PM", "1:07:07 AM", "1:07:05 PM", "NTC", "ext blank"))
str(qpcr_data)
ggplot(qpcr_data, aes(x = Sample, y = CT)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_y_reverse(limits=c(40, 15)) +
  annotate('rect', xmin=2.6, xmax=3.8, ymin=15, ymax=40, alpha=.2, fill='blue') +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))



