## DWR larval smelt 2023 live car analysis
# Avery Scherer
# 2024 Feb 08

library(artemis)
library(dplyr)

# because of issues with contamination in DWR's samples, we are only working with CFS samples here

# design:
# 3 biomasses (1, 50, 300)
# 2 distances (0-50 m, 51-100 m)
# 2 sampling methods (Sterivex, Smith Root)
# 3 reps per biomass
# 2 livecars per day so random block design across 4.5 days

######### set up ##################

# read in standard curve values, which were different for each plate in this study
std_curve <- readr::read_csv("1_data-raw/std_curve_2024-02-09.csv") %>% 
  rename(Alpha = "Log Intercept",
         Beta = "Log Slope") %>% 
  dplyr::select(Plate, Alpha, Beta)

d_1 = readRDS("1_data-raw/LarvalSmelt_subset.rds") %>% 
  # pull plate no from plate id
  mutate(Plate = stringr::str_sub(PlateID, start = 33, end = 34),
         Plate = as.numeric(ifelse(Plate == 10 | Plate == 11 | Plate == 12, Plate, stringr::str_sub(Plate, 1, 1)))) %>% 
  left_join(std_curve, by = "Plate")
unique(d_1$Plate)

# create version of dataset without the controls
d_2 = filter(d_1, Biomass != "0" & Biomass != "Control") %>% 
  # calc conc rather than Cq
  mutate(conc = exp((Cq-Alpha)/Beta),
         # code biomass as continuous and calc log version
         Biomass_c = as.numeric(Biomass),
         Biomass_t = log(Biomass_c),
         # create graphing labels for biomass and filter type
         Biomass_lables = case_when(Biomass == 1 ~ "larvae = 1",
                                    Biomass == 50 ~ "larvae = 50",
                                    Biomass == 300 ~ "larvae = 300"),
         Collector_labels = case_when(Collector == "CFS Smith Root 1.2um" ~ "Smith Root filters",
                                      Collector == "CFS Sterivex" ~ "Sterivex filters"))
# order factor levels
d_2$DistanceToLivecar <- factor(d_2$DistanceToLivecar, levels = c("50-0m", "100-50m"))
d_2$Biomass_lables <- factor(d_2$Biomass_lables, levels = c("larvae = 1", "larvae = 50", "larvae = 300"))
d_2$Collector_labels <- factor(d_2$Collector_labels, levels = c("Sterivex filters", "Smith Root filters"))

aspect_ratio <- 2.5

########## graph data ###################

# layer in the factors, starting with distance from livecar
(graph_distance <- ggplot(d_2, aes(y = conc, x = DistanceToLivecar, color = DistanceToLivecar)) +
    geom_jitter(alpha = 0.5, height = 0) +
    scale_color_manual(values = c("seagreen3", "deepskyblue2")) +
    labs(y = "Relative Concentration (ng/uL)") + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 9), 
          axis.title = element_text(size = 13), 
          legend.title = element_blank(), 
          legend.position = "none",
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")))
#ggsave(plot = graph_distance, height = 4, width = 2.5 * aspect_ratio, filename = "4_figs-and-tbls/AFS-calneva/graph-1_2024-04-01.jpg")

# add biomass
(graph_biomass <- ggplot(d_2, aes(y = conc, x = DistanceToLivecar, color = Biomass_lables)) +
    geom_jitter(alpha = 0.5, height = 0) +
    scale_color_manual(values = c("seagreen3", "deepskyblue2", "slateblue2")) +
    facet_grid(~Biomass_lables) + 
    labs(y = "Relative Concentration (ng/uL)") + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 9), 
          axis.title = element_text(size = 13), 
          legend.title = element_blank(), 
          legend.position = "none",
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")))
#ggsave(plot = graph_biomass, height = 3, width = 2.5 * aspect_ratio, filename = "4_figs-and-tbls/AFS-calneva/graph-2_2024-04-01.jpg")

# and finally filter type
(graph_filter <- ggplot(d_2, aes(y = conc, x = DistanceToLivecar, color = Collector)) +
    geom_jitter(alpha = 0.5, height = 0) +
    scale_color_manual(values = c("seagreen3", "deepskyblue3")) +
    facet_grid(Collector_labels~Biomass_lables) + 
    labs(y = "Relative Concentration (ng/uL)") + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 9), 
          axis.title = element_text(size = 13), 
          legend.title = element_blank(), 
          legend.position = "none",
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")))
#ggsave(plot = graph_filter, height = 4, width = 2.5 * aspect_ratio, filename = "4_figs-and-tbls/AFS-calneva/graph-3_2024-04-01.jpg")

########### model results #################

# final model selected in analysis script
m = eDNA_lmer(Cq ~ Collector + Biomass_t + DistanceToLivecar + 
                 (1|Day) + (1|FilterID),
               data = d_2,
               # for now use a general set of std curve values
               std_curve_alpha = d_2$Alpha,
               std_curve_beta = d_2$Beta,
               parallel_chains = 4L)

# same summary for use in plots
sum_m = as.data.frame(summary(m))
sum_m$Predictor = row.names(sum_m)
colnames(sum_m) = c("mean", "lwr", "median", "upr", "Predictor")
row2 = c("Intercept", "Filter type", "log(Biomass)", "Distance to Livecar", "sigma ln(eDNA)")
sum_m = cbind(sum_m,row2)

# set factor order - not sure why it plots backwards of this order
sum_m$row2 = factor(x = sum_m$row2, 
                      levels = c("Intercept", "sigma ln(eDNA)", "Filter type", "log(Biomass)", "Distance to Livecar"))

# stair step in the results, starting with sig effect of biomass
(graph_m_1 =
    ggplot(sum_m[sum_m$row2 != "Intercept", ], aes(x = median, y = row2, color = row2)) +
    geom_point(size = 2.5) +
    geom_segment(aes(xend = lwr, x = median, y = row2, yend = row2, color = row2)) +
    geom_segment(aes(xend = upr, x = median, y = row2, yend = row2, color = row2)) +
    geom_vline(xintercept = 0, color = "deeppink3", size = 0.75) +
    scale_color_manual(values = c("white", "white", "deepskyblue2", "white")) +
    labs(x = "Model Estimate (ln[eDNA])", y = NULL,
         subtitle = "Cq ~ Distance + log(Biomass) + Filter Type + (1|Day) + (1|FilterID)") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 10), 
          axis.title = element_text(size = 13),
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")))
#ggsave(plot = graph_m_1, height = 4, width = 2.5 * aspect_ratio, filename = "4_figs-and-tbls/AFS-calneva/graph-model-1_2024-04-01.jpg")

# add non sig effects of distance and filter type
(graph_m_2 =
    ggplot(sum_m[sum_m$row2 != "Intercept", ], aes(x = median, y = row2, color = row2)) +
    geom_point(size = 2.5) +
    geom_segment(aes(xend = lwr, x = median, y = row2, yend = row2, color = row2)) +
    geom_segment(aes(xend = upr, x = median, y = row2, yend = row2, color = row2)) +
    geom_vline(xintercept = 0, color = "deeppink3", size = 0.75) +
    scale_color_manual(values = c("white", "deepskyblue2", "ivory4", "deepskyblue2")) +
    labs(x = "Model Estimate (ln[eDNA])", y = NULL,
         subtitle = "Cq ~ Distance + log(Biomass) + Filter Type + (1|Day) + (1|FilterID)") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 10), 
          axis.title = element_text(size = 13),
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")))
#ggsave(plot = graph_m_2, height = 4, width = 2.5 * aspect_ratio, filename = "4_figs-and-tbls/AFS-calneva/graph-model-2_2024-04-01.jpg")

# and finally the large amount of variance
(graph_m_3 =
    ggplot(sum_m[sum_m$row2 != "Intercept", ], aes(x = median, y = row2, color = row2)) +
    geom_point(size = 2.5) +
    geom_segment(aes(xend = lwr, x = median, y = row2, yend = row2, color = row2)) +
    geom_segment(aes(xend = upr, x = median, y = row2, yend = row2, color = row2)) +
    geom_vline(xintercept = 0, color = "deeppink3", linewidth = 0.75) +
    scale_color_manual(values = c("deepskyblue2", "ivory4", "ivory4", "ivory4")) +
    labs(x = "Model Estimate (ln[eDNA])", y = NULL,
         subtitle = "Cq ~ Distance + log(Biomass) + Filter Type + (1|Day) + (1|FilterID)") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 10), 
          axis.title = element_text(size = 13),
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")))
#ggsave(plot = graph_m_3, height = 4, width = 2.5 * aspect_ratio, filename = "4_figs-and-tbls/AFS-calneva/graph-model-3_2024-04-01.jpg")

# same as the last plot, but all model results are grey (no highlight in blue) and biomass labeled (more accurately) abundance
# change biomass>abundance in summary
sum_m = as.data.frame(summary(m))
sum_m$Predictor = row.names(sum_m)
colnames(sum_m) = c("mean", "lwr", "median", "upr", "Predictor")
row2 = c("Intercept", "Filter type", "log(Abundance)", "Distance to Livecar", "sigma ln(eDNA)")
sum_m = cbind(sum_m,row2)
# reset factor order with abundance
sum_m$row2 = factor(x = sum_m$row2, 
                    levels = c("Intercept", "sigma ln(eDNA)", "Filter type", "log(Abundance)", "Distance to Livecar"))
(graph_m_4 =
    ggplot(sum_m[sum_m$row2 != "Intercept", ], aes(x = median, y = row2, color = row2)) +
    geom_point(size = 2.5) +
    geom_segment(aes(xend = lwr, x = median, y = row2, yend = row2, color = row2)) +
    geom_segment(aes(xend = upr, x = median, y = row2, yend = row2, color = row2)) +
    geom_vline(xintercept = 0, color = "deeppink3", linewidth = 0.75) +
    scale_color_manual(values = c("ivory4", "ivory4", "ivory4", "ivory4")) +
    labs(x = "Model Estimate (ln[eDNA])", y = NULL,
         subtitle = "Cq ~ Distance + log(Abundance) + Filter Type + (1|Day) + (1|FilterID)") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 10), 
          axis.title = element_text(size = 13),
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")))
#ggsave(plot = graph_m_4, height = 4, width = 2.5 * aspect_ratio, filename = "4_figs-and-tbls/graph-model_2024-05-13.jpg")

# within filter variation

# adjust concentration so Cq of 40 = conc of 0
d_3 <- d_2 %>% 
  mutate(conc_2 = ifelse(Cq == 40, 0, exp((Cq-Alpha)/Beta)))

(graph_filter <- ggplot(d_3, aes(y = conc_2, x = DistanceToLivecar, color = FilterID)) +
    geom_boxplot(alpha = 0.5) +
    facet_grid(Collector_labels~Biomass_lables) + 
    scale_color_grey() +
    labs(y = "Relative Concentration (ng/uL)") + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 13), 
          axis.text.y = element_text(size = 10), 
          axis.title = element_text(size = 13),
          plot.title = element_text(hjust = 0.5), 
          plot.subtitle = element_text(hjust = 0.5), 
          legend.position = "none",
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          axis.line = element_line(colour = "black")))
#ggsave(plot = graph_filter, height = 4, width = 2.5 * aspect_ratio, filename = "4_figs-and-tbls/AFS-calneva/graph-filter_2024-04-01.jpg")
