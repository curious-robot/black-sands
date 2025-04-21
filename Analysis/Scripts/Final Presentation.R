# Setup ----

library(tidyverse)
library(vegan)
library(ggpubr)
library(broom)
library(AER)
library(lme4)
library(MatchIt)
library(fixest)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
flowercounts = read_csv("../Scratch/blacksand_flowercounts_2020.ar.data.csv")
geumexclusions = read_csv("../Scratch/blacksand_geumexclusions_2020.ar.data.csv")
geumindividualvisits = read_csv("../Scratch/blacksand_geumindividualvisits_2020.ar.data.csv")
topography = read_csv("../Scratch/blacksand_topography_2020.ar.data.csv")
visitors = read_csv("../Scratch/blacksand_visitors_2020.ar.data.csv")
snowmelt = read_csv("../Scratch/snowmelt.csv")

# Clean snowmelt data to match pollinator data -----

names(snowmelt) = c("LTER", "site", "year", "treatment", "subplot", "melt_doy")
snowmelt$site = case_match(snowmelt$site, "East_Knoll" ~ "EastKnoll", .default = snowmelt$site)
snowmelt$subplot = case_match(snowmelt$subplot, "C"~"E", "B"~"C", .default = snowmelt$subplot)

# Simulate Snowmelt ------
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


# Get Species Observations -----------
# Collapse species observations into counts at each subplot grouped by date, site, treatment.
pollinator_data_summary <- visitors %>%
  mutate(year = as.character(format(as.Date(date, format="%Y-%m-%d"), "%Y"))) %>%
  group_by(year, site, subplot, treatment, visitor_nomenclature) %>%
  summarise(visit_count = sum(visitornumber)) %>%
  spread(key = visitor_nomenclature, value = visit_count, fill = 0)

# Calculate alpha diversity
pollinator_data_summary$richness <- specnumber(pollinator_data_summary[,-c(1:4)])  # Excluding date, site, and subplot columns

# Calculate Simpson diversity
pollinator_data_summary$simpson_diversity <- diversity(pollinator_data_summary[,-c(1:4)], index = "simpson")

# Combine Everything
# Bind snowmelt and pollinator data
combined <- pollinator_data_summary %>%
  merge(snowmelt, by = c("year", "subplot", "site", "treatment"), all.y = TRUE) %>%
  merge(topography, by = c("subplot", "site", "treatment"), all.y = TRUE) %>%
  mutate(
    simpson_diversity = ifelse(is.na(simpson_diversity), 0, simpson_diversity),
    richness = ifelse(is.na(richness), 0, richness)
  ) %>%
  filter(!is.na(year), !is.na(LTER)) %>%
  select(LTER, year, site, subplot, treatment, treated, richness, simpson_diversity, melt_doy, slope, elevation, tpi3, tpi11, tpi31, everything())


# Justify Instrument -----
# Justify treatment as an instrumental variable by showing its correlation with snowmelt time/snowmelt advancement. 

# Default lm regressing treatment variable on IV and other relevant factors. 
justify = lmer(melt_doy ~ treatment + (1 | subplot), data = combined)
justify = lm(melt_doy ~ treatment + slope + elevation, data = combined)
justify = feols(melt_doy ~ treatment | site, data = combined)

# F test to see if the instrument explains enough of the explanatory variable
wald(justify, ~ treatment)

snowmelt_effect = ggplot(combined, aes(x = treatment, y = melt_doy, fill = treatment)) +
  geom_boxplot() +
  labs(title = "Snowmelt Day Comparison Between Control and Treatment",
       x = "Treatment",
       y = "Snowmelt Day of Year") +
  scale_fill_manual(values = c("Control" = "cyan", "Early" = "gray")) +  # Optional: choose your colors
  theme_minimal() +
  facet_wrap(~ site) # Create a separate subplot for each site
snowmelt_effect
# P-value is significant but F statistic is not super high so mid justification of IV

# Implement Instrument --------------

# IVreg function
diversity_iv = ivreg(simpson_diversity ~ melt_doy + site + subplot | treatment + 
          subplot + site, data = combined)

diversity_iv_est = coef(summary(diversity_iv))["melt_doy", "Estimate"]


# OLS Method -----------

diversity_ols = combined %>% 
  lm(simpson_diversity ~ melt_doy + slope + elevation + tpi11, data = .)

diversity_ols_est = coef(summary(diversity_ols))["melt_doy", "Estimate"]


# Comparison of Methods -----------
both_methods = data.frame(c("OLS", "IV"), c(diversity_ols_est, diversity_iv_est))
colnames(both_methods) = c("method", "estimate")

### Reshape the data-------
### Visualize distribution of estimates---------
ggplot(both_methods, aes(x = method, y = estimate)) + 
  geom_boxplot() +
  theme_classic() +
  ylim(0, 0.002) +
  labs(y = "Estimated effect of snow melt doy",
       x = "Estimation method") 

combined$LTER = NULL


# Plots ------------
ggplot(combined[combined$year == "2020",], aes(x = site, y = simpson_diversity, fill = treatment)) +
  geom_boxplot() +
  labs(
    title = "Simpson Diversity by Site in 2020",
    x = "Site",
    y = "Simpson Diversity Index"
  ) +
  theme_minimal()


# Remove subplot from regression since they're nested and reduce power
# Remove East Knoll
# 2-Way Fixed effect on site and data collection day
# Include error in slideshow

