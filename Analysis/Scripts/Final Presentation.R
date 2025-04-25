# Setup ----

library(tidyverse)
library(vegan)
library(ggpubr)
library(broom)
library(AER)
library(lme4)
library(MatchIt)
library(fixest)
library(mgcv)
library(gratia)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
flowercounts = read_csv("../Scratch/blacksand_flowercounts_2020.ar.data.csv")
topography = read_csv("../Scratch/blacksand_topography_2020.ar.data.csv")
visitors = read_csv("../Scratch/blacksand_visitors_2020.ar.data.csv")
#snowmelt = read_csv("../Scratch/snowmelt.csv")
snowdepth_raw = read_csv("../Scratch/black_sand_snow_depth.cf.data.csv")

# Clean snowdepth data ------
snowdepth = snowdepth_raw %>%
  setNames(c("site", "treatment", "subplot", "distance_from_top", "date", "depth", "comments")) %>%
  filter(treatment != "", !is.na(treatment), depth != "NaN") %>%
  mutate(
    # rename subplot levels
    subplot = case_when(
      subplot == "Top" ~ "A",
      subplot == "Topmiddle" ~ "B",
      subplot == "Middle" ~ "C",
      subplot == "Bottommiddle" ~ "D",
      subplot == "Bottom" ~ "E",
      .default = NA
    ),
    # rename East Knoll
    site = ifelse(site == "East_Knoll", "EastKnoll", site),
    # get day-of-year for each date
    doy = yday(date),
    # convert snow depth
    depth = as.numeric(depth),
    # extract year from date
    year = format(as.Date(date, format="%Y/%m/%d"),"%Y")
  ) %>%
  # bind snowdepth data with topography data
  left_join(topography, by = c("site", "subplot", "treatment")) %>%
  mutate(treatment = as.factor(treatment))

snowdepth_soddie = filter(snowdepth, site != "Soddie")
snowdepth_lefty = filter(snowdepth, site != "Lefty")
snowdepth_three = filter(snowdepth, site != "Lefty" & site != "Soddie")


snowmelt = snowdepth %>%
  mutate(
    plot = paste(site, subplot)
  ) %>%
  group_by(year, site, subplot, treatment) %>%
  summarise(
    melt_doy = if (any(depth == 0, na.rm = TRUE)) {
      min(doy[depth == 0], na.rm = TRUE)
    } else {
      NA
    },
    .groups = "drop",
    tpi11 = first(tpi11),
    slope = first(slope),
    elevation = first(elevation),
    plot = first(plot)
  ) %>%
  filter(!is.na(melt_doy)) %>%
  mutate(year = as.numeric(year))

# Clean flower specie data ----------

# Assess treatment as an IV Using GAM --------------
# non linear gam model to assess the effect of treatment on the relationship between doy and depth.
gam_model = gam(
  # treatment is included as an "additive factor" to account for varried intercepts (overall deeper snow in treatment v control)
  depth ~ s(doy, by = treatment) + site + treatment, 
  data = snowdepth, 
  method = "REML"
)

summary(gam_model)

# non linear gam model to assess the relationship between doy and depth under "null" conditions (treatment not included in model)

gam_null <- gam(
  depth ~ s(doy) + treatment + site,
  data = snowdepth,
  method = "REML"
)

summary(gam_null)

# anova demonstrates that the model where treatment is allowed to affect depth ~ doy spline is better than null model
anova(gam_model, gam_null, test = "Chisq")

# data frame for fitting values
newdata_all <- expand.grid(
  doy = seq(min(snowdepth$doy), max(snowdepth$doy), length.out = 100),
  treatment = levels(factor(snowdepth$treatment)),
  site = levels(factor(snowdepth$site))
)

# predict fitted values
newdata_all$fit <- predict(gam_model, newdata = newdata_all, type = "response")

### plot predicted splines and snowdepth ~ doy data ------------
ggplot() +
  geom_point(data = snowdepth, aes(x = doy, y = depth, color = treatment), alpha = 0.2) +
  geom_line(data = newdata_all, aes(x = doy, y = fit, color = treatment), size = 1.1) +
  facet_wrap(~site) +
  labs(title = "Observed vs Predicted Snow Depth",
       subtitle = "GAM Predictions with Raw Data",
       x = "Day of Year", y = "Snow Depth (cm)") +
  theme_minimal()

### plot predicted values only.-----------
ggplot(newdata_all, aes(x = doy, y = fit, color = treatment)) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~site) +
  labs(title = "Predicted Snow Depth over Day of Year",
       subtitle = "Split by Site and Treatment",
       x = "Day of Year", y = "Predicted Snow Depth (cm)") +
  theme_minimal()



# Get Species Observations -----------
# Collapse species observations into counts at each subplot grouped by date, site, treatment.
pollinators <- visitors %>%
  group_by(date, site, subplot, treatment, visitor_nomenclature) %>%
  summarise(visit_count = sum(visitornumber), .groups = "drop") %>%
  pivot_wider(names_from = visitor_nomenclature, values_from = visit_count, values_fill = 0) %>%
  mutate(
    richness = replace_na(specnumber(across(where(is.numeric), .names = "richness")), 0),
    simpson = replace_na(diversity(across(where(is.numeric)), index = "simpson"), 0 ),
    obs_doy = yday(date),
    treatment = factor(treatment),
    year = year(date)
  ) %>%
  left_join(snowmelt, by = c("year", "site", "subplot", "treatment"))

# Justify Instrument Using linear model -----
# Justify treatment as an instrumental variable by showing its correlation with snowmelt time/snowmelt advancement. 
# Default lm regressing treatment variable on IV and other relevant factors. 
#justify = lmer(melt_doy ~ treatment + (1 | subplot), data = snowmelt)
justify = lm(melt_doy ~ treatment + site + subplot + year, data = snowmelt)
#justify = feols(simpson_diversity ~ 1 | melt_doy ~ treatment, data = combined)

# F test to see if the instrument explains enough of the explanatory variable
summary(justify)

snowmelt_effect = ggplot(snowmelt, aes(x = factor(year), y = melt_doy, fill = treatment)) +
  geom_boxplot() +
  labs(title = "Snowmelt Day Comparison Between Control and Treatment",
       x = "Treatment",
       y = "Snowmelt Day of Year") +
  scale_fill_manual(values = c("Control" = "cyan", "Early" = "gray")) +  # Optional: choose your colors
  theme_minimal() +
  facet_wrap(~ site) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
snowmelt_effect
# P-value is significant but F statistic is not super high so mid justification of IV
# Implement Instrument Using GAM ------------

gam_pollinator <- gam(
  # treatment is included as an "additive factor" to account for varried intercepts (overall deeper snow in treatment v control)
  simpson ~ s(doy, by = treatment) + site + treatment, 
  data = pollinators, 
  method = "REML"
)

summary(gam_pollinator)

# non linear gam model to assess the relationship between doy and depth under "null" conditions (treatment not included in model)

gam_pollinator_null <- gam(
  simpson ~ s(doy) + site + year + treatment, 
  data = pollinators, 
  method = "REML"
)

summary(gam_pollinator_null)

# anova demonstrates that the model where treatment is allowed to affect depth ~ doy spline is better than null model
anova(gam_pollinator, gam_pollinator_null, test = "Chisq")

# data frame for figure
plotdata <- expand.grid(
  doy = seq(min(pollinators$doy), max(pollinators$doy), length.out = 100),
  treatment = levels(factor(pollinators$treatment)),
  site = levels(factor(pollinators$site))
)

# predict fitted values
plotdata$fit = predict(gam_pollinator, newdata = plotdata, type = "response")

# plot values.
ggplot(plotdata, aes(x = doy, y = fit, color = treatment)) +
  geom_line(size = 1.1) +
  geom_point(data = pollinators, aes(x = doy, y = simpson, color = treatment), alpha = 0.2) +
  facet_wrap(~site) +
  labs(title = "Predicted Pollinator Simpson Diversity",
       subtitle = "Split by Site and Treatment",
       x = "Day of Year", y = "Predicted Simpson Diversity") +
  theme_minimal()

# IV with LSDV TWFE ----------

twfe_lsdv = lm(richness ~ treatment + site + subplot, data = pollinators)
summary(twfe_lsdv)

tpi_snowmelt = lm(melt_doy ~ tpi11 + site + subplot, data = snowmelt)
summary(tpi_snowmelt)

ggplot() +
  geom_point(data = snowmelt, aes(x = tpi11, y = melt_doy, color = treatment), alpha = 0.2) +
  facet_wrap(~site) +
  labs(title = "Predicted Pollinator Simpson Diversity",
       subtitle = "Split by Site and Treatment",
       x = "Day of Year", y = "Predicted Simpson Diversity") +
  theme_minimal()

ggplot(filter(snowmelt, year == 2020), aes(x = tpi11, fill = factor(treatment))) +
  geom_density(alpha = 0.5) +
  labs(x = "Subplot (A-E) TPI11", y = "Density", fill = "Treatment") +
  facet_wrap(~subplot) +
  theme_minimal()

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

