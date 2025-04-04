#ESP LE 2023
#set wd

#load libraries
library(ggplot2)
library(tidyr)
library(tidyverse)
library(dplyr)

#Edits to raw qPCR output csv file before loading into Rstudio:
#check amplification curves in QuantStudio, set CT threshold
#add column for sample ID in Excel
#Change "undetermined" Cq to 0
#Change time stamps (UTC to PST) in Excel

#load edited csv file
qpcr_data <- read.csv("redo ESP 41 qPCR 081424.csv")
head(qpcr_data)

#dplyr to only plot live cars 1-9 and controls
qpcr_data2 <- qpcr_data %>% filter(livecar > 0) %>% filter(livecar!=10) %>% filter(biomass!="na")

#std curve data - change Cq to relative conc

#Modified code from Avery/Katie at CFS below
#read in standard curve values, which were different for each plate in this study
std_curve <- readr::read_csv("ESP_stdcurves_082024.csv") %>% 
  rename(Alpha = "Log Intercept",
         Beta = "Log Slope") %>% dplyr::select(Plate, Alpha, Beta)

d_1 = qpcr_data2 %>% 
  left_join(std_curve, by = "Plate")
unique(d_1$Plate)

#calc conc rather than Cq
d_2 = d_1 %>% mutate(conc = exp((Cq-Alpha)/Beta))

d_2$biomass <- factor(d_2$biomass, levels = c("Larva = 1", "Larvae = 50", "Larvae = 300", "NTC", "ext blank"))

library(scales)

ggplot(d_2, aes(y = conc, x = biomass, colour = biomass)) +
  geom_jitter(alpha = 0.5, height = 0) +
  scale_color_manual(values = c("orange", "orange", "orange")) +
  labs(x = "Biomass", y = "Relative Concentration (ng/uL)") + 
  scale_y_continuous(labels = ~ sprintf(fmt = "%0.02e", .)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 9), 
        axis.title = element_text(size = 13), 
        legend.title = element_blank(), 
        legend.position = "none",
        panel.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.line = element_line(colour = "black"))

  