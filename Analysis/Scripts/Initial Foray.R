library(tidyverse)
library(vegan)
library(ggpubr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
flowercounts = read_csv("../Scratch/blacksand_flowercounts_2020.ar.data.csv")
geumexclusions = read_csv("../Scratch/blacksand_geumexclusions_2020.ar.data.csv")
geumindividualvisits = read_csv("../Scratch/blacksand_geumindividualvisits_2020.ar.data.csv")
topography = read_csv("../Scratch/blacksand_topography_2020.ar.data.csv")
visitors = read_csv("../Scratch/blacksand_visitors_2020.ar.data.csv")
head(dataframe)


# Create a contingency table
pollinator_data_summary <- visitors %>%
  group_by(date, site, subplot, treatment, visitor_nomenclature) %>%
  summarise(visit_count = sum(visitornumber)) %>%
  spread(key = visitor_nomenclature, value = visit_count, fill = 0)

# View the contingency table (optional)
head(pollinator_data_summary)

pollinator_data_summary$richness <- specnumber(pollinator_data_summary[,-c(1:3)])  # Excluding date, site, and subplot columns

# Calculate Shannon Diversity Index
pollinator_data_summary$simpson_diversity <- diversity(pollinator_data_summary[,-c(1:4)], index = "simpson")

# Create a data frame for plotting
diversity_df <- data.frame(date = pollinator_data_summary$date, Simpson = pollinator_data_summary$simpson_diversity)

# Plot Shannon Diversity Index over time
library(ggplot2)
ggplot(diversity_df, aes(x = date, y = Simpson)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Simpson Diversity Index over Time", x = "Date", y = "Simpson Diversity Index")


# Add statistical comparisons (e.g., t-test)
ggplot(pollinator_data_summary, aes(x = treatment, y = simpson_diversity, fill = treatment)) +
  geom_boxplot() +
  facet_wrap(~site) +
  stat_compare_means(method = "t.test") +  # t-test between treatments
  theme_minimal() +
  labs(title = "Simpson Diversity Index by Treatment and Site", x = "Treatment", y = "Shannon Diversity Index") +
  scale_fill_manual(values = c("Control" = "skyblue", "Early" = "orange")) +
  ylim(0, 1)