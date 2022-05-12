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
                      block = "f",
                      exposition = col_factor(levels = c("north", "south")),
                      sandRatio = "f",
                      substrateDepth = "f",
                      targetType = "c",
                      seedDensity = "f"
                    )) %>%
  select(
    id, plot, block, exposition, sandRatio, substrateDepth, targetType,
    seedDensity, surveyYear
  ) %>%
  filter(targetType != "non_ffh") %>%
  mutate(
    targetType = if_else(targetType == "mixed", "hay_meadow", targetType),
    surveyYear_fac = as.character(surveyYear),
    surveyYear_fac = if_else(block == "C", "reference", surveyYear_fac),
    surveyYear_fac = factor(surveyYear_fac)
  )

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


### 2 PERMANOVA ###############################################################

dist <- dist(ordi, method = "bray")
(disp.year <- betadisper(dist, sites$surveyYear_fac * sites$exposition * sites$targetType))
permutest(disp.year) # similar dispersion -> PERMANOVA possible
(permanova.year <- adonis(species ~ year,
                          data = sites, 
                          strata = sites$plot,
                          permutations = 999,
                          method = "bray"))
densityplot(permustats(permanova.year)) # significant


contrasts(edataT$year) <- cbind(c(1, 0, -1), c(0, 1, -1))
Year.ortho <- model.matrix(~ edataT$year)[, -1]
year8418 <- Year.ortho[, 1]
year9318 <- Year.ortho[, 2]
m3 <- adonis(vdataT ~ year8418+year9318,
             data = edataT, strata = edataT$plot,
             permutations = 999, method = "bray")
m3 <- adonis(vdataT ~ year,
             data = edataT, strata = edataT$plot,
             permutations = 999, method = "bray")
m3 #beide sig.; R2 = 0.15 und 0.09

