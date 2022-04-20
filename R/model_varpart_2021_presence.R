# Dike grassland experiment
# Variation partitioning of 2021 (presence-absence data) ####
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


## 2 Variation partitioning: overall ##########################################

m1_total_varpart <- varpart(beta_total, sites_soil, sites_seedmix, sites_space)
plot(
  m1_total_varpart,
  Xnames = c("Substrate", "Seedmix", "Exposition"),
  cutoff = 0.01, digits = 1, bg = NA, id.size = 1
)


## 3 db-RDA: replacemnet ######################################################

### a full model --------------------------------------------------------------

m1 <- dbrda(
  beta_substitution ~ sandRatio + substrateDepth + exposition +
    seedDensity + targetType,
  data = sites
)
anova(m1, permutations = how(nperm = 9999)) # P = 1e-04
(r2adj <- RsquareAdj(m1)$adj.r.squared) # R2adj = .398

### b forward selection -------------------------------------------------------

#### * Soil ####
m1 <- dbrda(
  beta_substitution ~ sandRatio + substrateDepth, data = sites
)
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(beta_substitution,
                   sites_soil,
                   adjR2thresh = r2adj,
                   nperm = 9999)
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_soil))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_soil_selected <- sites %>%
  select(sandRatio)

#### * Seedmix ####
m1 <- dbrda(beta_substitution ~ targetType + seedDensity, data = sites)
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(beta_substitution,
                   sites_seedmix,
                   adjR2thresh = r2adj,
                   nperm = 9999)
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_seedmix))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_seedmix_selected <- sites %>%
  select(targetType_numeric)

#### * Exposition ####
m1 <- dbrda(beta_substitution ~ exposition, data = sites)
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(beta_substitution,
                   sites_space,
                   #adjR2thresh = r2adj,
                   nperm = 9999)
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_space))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_space_selected <- sites %>%
  select(exposition_numeric)

### c Variation partitioning -------------------------------------------------

m1_substitution_varpart <- varpart(
  beta_substitution,
  sites_soil_selected, sites_seedmix_selected, sites_space_selected
)
tiff(
  here("outputs", "figures", "figure_2021_replacement_800dpi_12x12cm.tiff"),
  res = 72, width = 12, height = 12, units = "cm", compression = "none"
)
plot(
  m1_substitution_varpart,
  Xnames = c("Site", "Seedmix", "Exposition"),
  cutoff = 0.01, digits = 2, bg = NA
)
dev.off()

### d partial db-RDA --------------------------------------------------------

#### * Seedmix / targetType ####
m1_substitution <- dbrda(
  beta_substitution ~ targetType +
    Condition(exposition + sandRatio),
  data = sites
)
anova(m1_substitution, permutations = how(nperm = 9999)) # p = 1e-04***
RsquareAdj(m1_substitution) # R2adj = .030

#### * Exposition / exposition ####
m1_substitution <- dbrda(
  beta_substitution ~ exposition +
    Condition(sandRatio + targetType),
  data = sites
)
anova(m1_substitution, permutations = how(nperm = 9999)) # p = 1e-04 ***
RsquareAdj(m1_substitution) # R2adj = .224


## 3 db-RDA: nestedness #######################################################

### a full model --------------------------------------------------------------

m1 <- dbrda(
  beta_subsets ~ sandRatio + substrateDepth +
    seedDensity + targetType + exposition,
  data = sites
)
anova(m1, permutations = how(nperm = 9999)) # P = 1e-04

### b forward selection -------------------------------------------------------

#### * Soil ####
m1 <- dbrda(
  beta_subsets ~ sandRatio + substrateDepth, data = sites
)
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(beta_subsets,
                   sites_soil,
                   #adjR2thresh = r2adj,
                   nperm = 9999)
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_soil))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_soil_selected <- sites %>%
  select()

#### * Seedmix ####
m1 <- dbrda(beta_subsets ~ targetType + seedDensity, data = sites)
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(beta_subsets,
                   sites_seedmix,
                   #adjR2thresh = r2adj,
                   nperm = 9999)
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_seedmix))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_seedmix_selected <- sites %>%
  select(seedDensity) #Dummy

#### * exposition ####
m1 <- dbrda(beta_subsets ~ exposition, data = sites)
r2adj <- RsquareAdj(m1)$adj.r.squared
sel <- forward.sel(beta_subsets,
                   sites_space,
                   adjR2thresh = r2adj,
                   nperm = 9999)
sel$p_adj <- p.adjust(sel$pvalue, method = "holm", n = ncol(sites_space))
sel # https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
sites_space_selected <- sites %>%
  select(exposition_numeric)

### c Variation partitioning -------------------------------------------------

m1_subsets_varpart <- varpart(
  beta_subsets, sites_seedmix_selected, sites_space_selected
)
tiff(
  here("outputs", "figures", "figure_2021_nestedness_800dpi_12x12cm.tiff"),
  res = 72, width = 12, height = 12, units = "cm", compression = "none"
)
plot(
  m1_subsets_varpart,
  Xnames = c("Seedmix", "Exposition"),
  cutoff = 0.01, digits = 2, bg = NA
)
dev.off()

### d partial db-RDA --------------------------------------------------------

#### * Exposition / exposition ####
m1_subsets <- dbrda(
  beta_subsets ~ exposition +
    Condition(sandRatio),
  data = sites
)
anova(m1_subsets, permutations = how(nperm = 9999)) # p = 0.013 *
RsquareAdj(m1_subsets) # R2adj = .025

