library(tidyverse)
library(vegan)
library(ggpubr)
library(broom)
library(AER)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
flowercounts = read_csv("../Scratch/blacksand_flowercounts_2020.ar.data.csv")
geumexclusions = read_csv("../Scratch/blacksand_geumexclusions_2020.ar.data.csv")
geumindividualvisits = read_csv("../Scratch/blacksand_geumindividualvisits_2020.ar.data.csv")
topography = read_csv("../Scratch/blacksand_topography_2020.ar.data.csv")
visitors = read_csv("../Scratch/blacksand_visitors_2020.ar.data.csv")
snowmelt = read_csv("../Scratch/snowmelt.csv")

# Clean snowmelt data to match pollinator data

names(snowmelt) = c("LTER", "site", "year", "treatment", "subplot", "melt_doy")
snowmelt$site = case_match(snowmelt$site, "East_Knoll" ~ "EastKnoll", .default = snowmelt$site)
snowmelt$subplot = case_match(snowmelt$subplot, "C"~"E", "B"~"C", .default = snowmelt$subplot)

# ILLEGAL simulation of snowmelt day for the missing B and D subplots that were present
# in Rose-Person's data but missing from the ITEX snowmelt data. Need to find a better
# and responsible way to account for this. 
simulated_snowmelt <- snowmelt %>%
  # Group by site, year, and treatment for averaging
  group_by(site, year, treatment) %>%
  # Create rows for blocks B and D
  summarise(
    melt_doy_B = mean(melt_doy[subplot %in% c("A", "C")], na.rm = TRUE),
    melt_doy_D = mean(melt_doy[subplot %in% c("C", "E")], na.rm = TRUE)
  ) %>%
  # Reshape the data into the right format (replicate site, year, and treatment for B and D)
  pivot_longer(cols = starts_with("melt_doy"), names_to = "subplot", values_to = "melt_doy") %>%
  mutate(subplot = case_when(
    subplot == "melt_doy_B" ~ "B",
    subplot == "melt_doy_D" ~ "D",
    TRUE ~ subplot
  ))

# Combine with the original data
snowmelt <- bind_rows(snowmelt, simulated_snowmelt)

# Collapse species observations into counts at each subplot grouped by date, site, treatment.
pollinator_data_summary <- visitors %>%
  group_by(date, site, subplot, treatment, visitor_nomenclature) %>%
  summarise(visit_count = sum(visitornumber)) %>%
  spread(key = visitor_nomenclature, value = visit_count, fill = 0)

head(pollinator_data_summary)

# Calculate alpha diversity
pollinator_data_summary$richness <- specnumber(pollinator_data_summary[,-c(1:4)])  # Excluding date, site, and subplot columns

# Calculate Simpson diversity
pollinator_data_summary$simpson_diversity <- diversity(pollinator_data_summary[,-c(1:4)], index = "simpson")

# Add year variable for pollinator data
pollinator_data_summary$year = format(as.Date(pollinator_data_summary$date, format="%Y-%m-%d"),"%Y")

# Bind snowmelt and pollinator data
combined = merge(pollinator_data_summary, snowmelt, by = c("year", "subplot", "site", "treatment"), all.y = TRUE)
combined = select(combined, LTER, year, date, site, subplot, treatment, richness, simpson_diversity, melt_doy, everything())


# Justify treatment as an instrumental variable by showing its correlation with snowmelt time/snowmelt advancement. 
#############################

# Default lm regressing outcome on IV and other relevant factors. 
justify = lm(melt_doy ~ treatment + site + subplot, data = snowmelt)

# F test to see if the instrument explains enough of the explanatory variable
justify_ftest <- waldtest(justify, .~.-treatment)
print(justify_ftest)

# P-value is significant but F statistic is not super high so mid justification of IV


# IVreg function
diversity_iv = combined %>%
  ivreg(simpson_diversity ~ melt_doy + subplot + site | treatment + 
          subplot + site, data = .)

diversity_iv
  
