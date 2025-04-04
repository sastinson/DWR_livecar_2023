#load libraries
library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)


#load qPCR data from csv file
qpcr_data <- read.csv("2023-11-29_DWRLarvalSmelt_DryRun_qPCR.csv", stringsAsFactors = FALSE)
head(qpcr_data)

#subset
SR <- c("SR", "Control")
TN <- c("TN", "Control")
subset <- qpcr_data %>% filter(Sample.type %in% TN) 

#plot
ggplot(data=subset, aes(x = FilterID, y = Cq, label = FilterID, fill = factor(FilterID))) +
  geom_boxplot() +
  scale_y_continuous() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  