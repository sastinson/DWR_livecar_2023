#set wd

#load libraries
library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)

#Edit raw qPCR output csv file
#check amplification curves in QuantStudio, set CT threshold
#add column for sample ID
#Change time stamps (UTC to PST)

#load edited csv file
qpcr_data <- read.csv("redo ESP qPCR 081424.csv")
head(qpcr_data)

#dplyr to only plot live cars 1-9 and controls
qpcr_data2 <- qpcr_data %>% filter(livecar > 0) %>% filter(livecar!=10) %>% filter(biomass!="na")

#basic plot
ggplot(data=qpcr_data2, aes(x = biomass, y = Cq)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

#boxplot
ggplot(qpcr_data2, aes(x = biomass, y = Cq)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  scale_y_reverse(limits=c(41, 8))

#jitterplot
ggplot(qpcr_data2, aes(x = biomass, y = Cq)) +
  geom_jitter() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  scale_y_reverse(limits=c(41, 8))

#flip Y axis values, adjust Y axis range, group by technical replicate. 
qpcr_data2$biomass <- factor(qpcr_data2$biomass, levels = c("Larvae = 0", "Larva = 1", "Larvae = 50", "Larvae = 300", "NTC", "ext blank"))
str(qpcr_data2)

level_order <- c("NTC", "ext blank", "Larva = 1", "Larvae = 50", "Larvae = 300") 

#jitterplot, color by by livecar
ggplot(qpcr_data2, aes(x = factor(biomass, level = level_order), y = Cq)) +
  geom_jitter(width = 0.25, aes(colour = livecar)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + 
  scale_y_reverse(limits=c(42, 25))
