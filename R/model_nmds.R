# Model for ordination 84-93-18 ####
# Markus Bauer
# Citation: Markus Bauer & Harald Albrecht (2020) Basic and Applied Ecology 42, 15-26
# https://doi.org/10.1016/j.baae.2019.11.003



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Load data ###
sites <- read_csv("data_processed_sites.csv",
                  col_names = TRUE,
                  na = c("na", "NA", ""),
                  col_types =
                    cols(
                      .default = "?",
                      id = "f",
                      plot = "f",
                      site = "f",
                      exposition = col_factor(levels = c("north", "south")),
                      sand_ratio = "f",
                      substrate_depth = "f",
                      target_type = "c",
                      seed_density = "f"
                    )) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year
  ) %>%
  mutate(survey_year_fct = factor(survey_year))

species <- read_csv("data_processed_species.csv",
                  col_names = TRUE,
                  na = c("na", "NA", ""),
                  col_types = cols(
                    .default = "d",
                    name = "f",
                  )) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(sites, by = "id") %>%
  column_to_rownames("id") %>%
  mutate(across(where(is.numeric), ~replace(., 0, NA)))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### 1 NMDS ####################################################################

#### a ordination -------------------------------------------------------------
set.seed(123)
(ordi <- metaMDS(species,
                 dist = "bray", binary = FALSE, autotransform = TRUE,
                 try = 99, previous.best = TRUE, na.rm = TRUE))
save(ordi, file = here("outputs", "models", "model_fcs_simple.Rdata"))
stressplot(ordi) # stress: 0.207; Non-metric fit R² =.957

#### b environmental factors --------------------------------------------------
dist <- dist(ordi, method = "euclidean")
(ef <- envfit(ordi ~ surveyYear_fac + exposition + sandRatio + substrateDepth +
                targetType + seedDensity,
              data = sites,
              permu = 999,
              na.rm = TRUE))
plot(ordi, type = "n")
plot(ef, add = TRUE, p. = .05)
text(ordi, dis = "sites", cex = .7)
ordiellipse(ordi, sites$surveyYear_fac, kind = "sd", draw = "lines", label = TRUE)


