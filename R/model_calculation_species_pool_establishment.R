# Grassland experiment on dikes
# Calculate the species pool and the overall establishment of species ####

# Markus Bauer
# 2023-02-18



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)

### Start ###
rm(list = ls())

### Load data ###
traits <- read_csv(
  here("data", "processed", "data_processed_traits.csv"),
  col_names = TRUE,
  na = c("", "NA", "na"),
  col_types = cols(.default = "?")
  )

base::load(file = here("outputs", "models", "model_nmds.Rdata"))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Calculate #################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Check size of species pool for seeding ###
traits %>%
  select(seeded_hay_meadow, seeded_dry_grassland) %>%
  mutate(across(where(is.character), ~as.numeric(.x))) %>%
  summarise(across(where(is.numeric), ~sum(.x)))

### Check species number used for NMDS ###
ordi$species %>%
  as.data.frame() %>%
  rownames_to_column(var = "name") %>%
  as_tibble()

### How many species of the species pool established? ###
traits %>%
  select(
    name, seeded_hay_meadow, seeded_dry_grassland,
    starts_with("total_established")
    ) %>%
  pivot_longer(
    cols = starts_with("total_estab"), names_to = "year", values_to = "value"
    ) %>%
  pivot_longer(
    cols = starts_with("seeded_"), names_to = "type", values_to = "seeded"
    ) %>%
  filter(seeded == 1 & value > 0) %>%
  group_by(type, year) %>%
  count()

### Average establishment rate of all species ###
traits %>%
  select(name, seeded, starts_with("rate")) %>%
  filter(seeded == 1) %>%
  pivot_longer(
    cols = starts_with("rate_"), names_to = "year", values_to = "value"
    ) %>%
  group_by(year) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE)
    )
