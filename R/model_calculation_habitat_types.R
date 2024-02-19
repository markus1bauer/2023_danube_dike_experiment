# Grassland experiment on dikes
# Classification of plots to habitat types ####
# Just descriptive analysis

# Markus Bauer
# 2024-02-19



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)

### Start ###
rm(list = ls())

### Load data ###
sites <- read_csv(
  here("data", "processed", "data_processed_sites_spatial.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types =
    cols(
      .default = "?",
      plot = "f",
      site = "f",
      sand_ratio = "f",
      substrate_depth = col_factor(levels = c("30", "15")),
      target_type = col_factor(levels = c("hay_meadow", "dry_grassland")),
      seed_density = "f",
      exposition = col_factor(levels = c("north", "south")),
      survey_year = "c"
    )
) %>%
  ### Exclude data of seed mixtures
  filter(survey_year != "seeded") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    botanist_year = str_c(survey_year, botanist, exposition, sep = " "),
    botanist_year = factor(botanist_year),
    n = esy,
    id = factor(id)
  ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, survey_year, botanist_year, n
  )



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics #################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Data exploration ##########################################################


### a Graphs of raw data -------------------------------------------------------
library(broom)
data <- sites %>%
  group_by(survey_year_fct, n) %>%
  count() %>%
  mutate(ratio = round(nn/288, digits = 2)) %>%
  print(n = 30)

data %>%
  ggplot(aes(y = ratio, x = survey_year_fct)) +
  facet_grid(~ n) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 0.4))
