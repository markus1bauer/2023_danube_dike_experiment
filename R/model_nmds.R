# Dike grassland field experiment
# NMDS ordination ####
# Model building

# Markus Bauer
# 2022-11-15



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d", "ordi")))
setwd(here("data", "processed"))

#### * Load sites data ####

sites_experiment <- read_csv("data_processed_sites.csv",
                             col_names = TRUE, na = c("na", "NA", ""),
                             col_types = cols(.default = "?")) %>%
  mutate(reference = survey_year)
sites_splot <- read_csv("data_processed_sites_splot.csv",
                        col_names = TRUE, na = c("na", "NA", ""),
                        col_types = cols(
                          .default = "?",
                          survey_year = "c"
                        )) %>%
  mutate(
    reference = if_else(
      esy == "E12a", "+Reference", if_else(
        esy == "E22", "+Reference", "other"
      )
    ),
    target_type = if_else(
      esy == "E12a", "dry_grassland", if_else(
        esy == "E22", "hay_meadow", "other"
      )
    ),
    exposition = "other"
  )
sites_bauer <- read_csv("data_processed_sites_bauer.csv",
                        col_names = TRUE, na = c("na", "NA", ""),
                        col_types = cols(
                          .default = "?",
                          survey_year = "c"
                        )) %>%
  filter(exposition == "south" | exposition == "north") %>%
  mutate(
    reference = if_else(
      esy == "R1A", "+Reference", if_else(
        esy == "R22", "+Reference", if_else(
          esy == "R", "Grassland", if_else(
            esy == "?", "no", if_else(
              esy == "+", "no", if_else(
                esy == "R21", "Grassland", if_else(
                  esy == "V38", "-Reference", "other"
                )
              )
            )
          )
        )
      )
    ),
    target_type = if_else(
      esy == "R1A", "dry_grassland", if_else(
        esy == "R22", "hay_meadow", "other"
      )
    )
  )
sites <- sites_experiment %>%
  bind_rows(sites_splot, sites_bauer) %>%
  select(
    id, esy, reference,
    exposition, sand_ratio, substrate_depth, target_type, seed_density,
    survey_year, longitude, latitude, elevation, plot_size
  ) %>%
  arrange(id)


#### * Load species data ####

species_experiment <- read_csv("data_processed_species.csv",
                               col_names = TRUE, na = c("na", "NA", ""),
                               col_types = cols(.default = "?"))
species_splot <- read_csv("data_processed_species_splot.csv",
                          col_names = TRUE, na = c("na", "NA", ""),
                          col_types = cols(.default = "?"))
species_bauer <- read_csv("data_processed_species_bauer.csv",
                          col_names = TRUE, na = c("na", "NA", ""),
                          col_types = cols(.default = "?"))

### Exclude rare species (< 0.5% accumulated cover in all plots)
data <- species_experiment %>%
  full_join(species_splot, by = "name") %>%
  full_join(species_bauer, by = "name") %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  group_by(name) %>%
  summarise(total_cover_species = sum(value)) %>%
  filter(total_cover_species < 0.5)

species <- species_experiment %>%
  full_join(species_splot, by = "name") %>%
  full_join(species_bauer, by = "name") %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  filter(!(name %in% data$name)) %>% # use 'data' to filter
  pivot_wider(names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(sites, by = "id")

sites <- sites %>%
  semi_join(species, by = "id")

species <- species %>%
  column_to_rownames("id")

rm(list = setdiff(ls(), c(
  "sites", "species", "theme_mb", "vegan_cov_ellipse", "ordi"
)))

#### * Model ####

set.seed(12)
(ordi <- metaMDS(species, binary = TRUE,
                 try = 50, previous.best = TRUE, na.rm = TRUE))
#Wisonsin sqrt transformation, stress type 1
stressplot(ordi) # stress: 0.205; Non-metric fit R² =.958



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
save(ordi, file = here("outputs", "models", "model_nmds.Rdata"))
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


