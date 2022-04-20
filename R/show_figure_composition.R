# Dike grassland experiment
# Show_figure species composition ####
# Markus Bauer
# 2022-04-05


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation #########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(vegan)
library(adespatial)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Load data ###
sites <- read_csv("data_processed_sites.csv",
                  col_names = TRUE,
                  na = c("na", "NA", ""), col_types =
                    cols(
                      .default = "?",
                      id = "f",
                      plot = "f",
                      block = "f",
                      exposition = "f",
                      sandRatio = "d",
                      substrateDepth = "d",
                      targetType = "f",
                      seedDensity = "d"
                    )) %>%
  filter(surveyYear == 2021 & !str_detect(id, "C")) %>%
  select(
    id, plot, block, exposition, sandRatio, substrateDepth, targetType,
    seedDensity, surveyYear, vegetationCov
  ) %>%
  mutate(
    exposition_numeric = as.double(exposition),
    targetType_numeric = as.double(targetType),
    block_numeric = as.double(block)
  )

species <- read_csv("data_processed_species.csv",
                    col_names = TRUE,
                    na = c("na", "NA", ""), col_types =
                      cols(
                        .default = "d",
                        name = "f"
                      )
) %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = name, values_from = "value") %>%
  semi_join(sites, by = "id") %>%
  arrange(id) %>%
  column_to_rownames("id")

sites <- sites %>%
  column_to_rownames("id")

sites_soil <- sites %>%
  select(sandRatio, substrateDepth)
sites_seedmix <- sites %>%
  select(targetType_numeric, seedDensity)
sites_space <- sites %>%
  select(exposition_numeric)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Beta diversity #####################################################

### * check collinearity ####
data <- sites %>%
  select(where(is.numeric), -ends_with("numeric"), -surveyYear)
GGally::ggpairs(data, lower = list(continuous = "smooth_loess"))
#--> no relationship has r > 0.7 (Dormann et al. 2013 Ecography)

### * Calculate: Baselga presence-absence ####
beta <- beta.div.comp(species, coef = "BS", quant = FALSE)
beta$Note
beta$part # total = 0.325, substitution = 0.279, subsets = 0.046
beta_total <- beta$D %>% # sörensen dissimilarity
  as.matrix() %>%
  as.data.frame()
beta_substitution <- beta$repl %>% # replacement / simpson dissimilarity
  as.matrix()
beta_subsets <- beta$rich %>% # nestedness
  as.matrix()


## 3 db-RDA: replacemnet ######################################################

### a full model --------------------------------------------------------------

m1 <- dbrda(
  beta_substitution ~ sandRatio + substrateDepth + exposition +
    seedDensity + targetType,
  data = sites
)
anova(m1, permutations = how(nperm = 9999)) # P = 1e-04
(r2adj <- RsquareAdj(m1)$adj.r.squared) # R2adj = .398


## 3 db-RDA: nestedness #######################################################

### a full model --------------------------------------------------------------

m1 <- dbrda(
  beta_subsets ~ sandRatio + substrateDepth +
    seedDensity + targetType + exposition,
  data = sites
)
anova(m1, permutations = how(nperm = 9999)) # P = 1e-04

