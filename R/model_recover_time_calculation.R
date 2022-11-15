# Dike grassland field experiment
# Recover time (presence-absence data) ####
# Markus Bauer
# 2022-05-24



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(brms)
library(DHARMa)
library(bayesplot)
library(loo)
library(tidybayes)
library(emmeans)

### Start ###
rm(list = ls())

### Load data ###
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

load(file = here("outputs", "models", "model_nmds.Rdata"))

sites2 <- sites %>%
  mutate(nmds1 = ordi$points[, 1], nmds2 = ordi$points[, 2]) %>%
  group_by(exposition, target_type)
  mutate(data = sites %>% filter(reference == "+Reference"), mean = mean(nmds1))
  
rm(list = setdiff(ls(), c(
  "sites", "theme_mb", "vegan_cov_ellipse", "ordi"
)))




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics #################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Data exploration ##########################################################


### a Graphs of raw data -------------------------------------------------------

plot1 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = sand_ratio)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Sand ratio [vol%] (Data of 2021)")
plot2 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = substrate_depth)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Substrate depth [cm] (Data of 2021)")
plot3 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = target_type)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Target type (Data of 2021)")
plot4 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = seed_density)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Seed density [g/m²] (Data of 2021)")
(plot1 + plot2) / (plot3 + plot4)
plot1 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = exposition)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Exposition (Data of 2021)")
plot2 <- ggplot(sites, aes(y = n, x = survey_year_fct)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Survey year")
plot3 <- ggplot(sites %>% filter(survey_year == 2021), aes(y = n, x = site)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Blocks (Data of 2021)")
plot4 <- ggplot(sites, aes(y = n, x = botanist_year)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Botanists and survey year") +
  theme(axis.text.x = element_text(angle = 90))
(plot1 + plot2) / (plot3 + plot4)
ggplot(data = sites,
       aes(x = sand_ratio, y = n, color = target_type)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_grid(exposition ~ survey_year_fct) +
  labs(title = "Target type vs. sand ratio [vol%]")
ggplot(data = sites,
       aes(x = sand_ratio, y = n, color = substrate_depth)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_grid(exposition ~ survey_year_fct) +
  labs(title = "Substrate depth [cm] vs. sand ratio [vol%]")
ggplot(data = sites,
       aes(x = seed_density, y = n, color = target_type)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_grid(exposition ~ survey_year_fct) +
  labs(title = "Target type vs. Seed density [g/m²]")


### b Outliers, zero-inflation, transformations? ----------------------------

sites %>% group_by(exposition) %>% count(site)
boxplot(sites$n)
ggplot(sites, aes(x = exposition, y = n)) + geom_quasirandom()
ggplot(sites, aes(x = n)) + geom_histogram(binwidth = 0.03)
ggplot(sites, aes(x = n)) + geom_density()
ggplot(sites, aes(sqrt(n))) + geom_density()



## 2 Model building ###########################################################


### a Preparation -------------------------------------------------------------

### Posssible priors ###
get_prior(n ~ target_type + exposition + sand_ratio + survey_year_fct +
            seed_density + substrate_depth +
            (1 | site/plot) + (1 | botanist_year),
          data = sites)
### Example of normal distribution
ggplot(data = data.frame(x = c(-2, 2)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0.1, sd = 1))
### Example of cauchy distribution
ggplot(data = data.frame(x = c(-2, 2)), aes(x = x)) +
  stat_function(fun = dcauchy, n = 101, args = list(location = 0, scale = 1))
### Example of a student t distribution
ggplot(data.frame(x = c(-2, 2)), aes(x = x)) +
  stat_function(fun = dstudent_t, args = list(df = 3, mu = 0, sigma = 2.5))

### a Models ------------------------------------------------------------------

### Specifications ###
iter = 10000
chains = 4
thin = 2
priors <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(0.1, 1)", class = "b", coef = "sand_ratio25"),
  set_prior("normal(0.2, 1)", class = "b", coef = "sand_ratio50"),
  set_prior("normal(0.1, 1)", class = "b", coef = "expositionsouth"),
  set_prior("normal(0.1, 1)", class = "b", coef = "survey_year_fct2019"),
  set_prior("normal(0.2, 1)", class = "b", coef = "survey_year_fct2020"),
  set_prior("normal(0.3, 1)", class = "b",coef = "survey_year_fct2021"),
  set_prior("cauchy(0, 1)", class = "sigma")
)
### Models ###
m_simple <- brm(n ~ target_type + exposition + sand_ratio + survey_year_fct +
                  seed_density + substrate_depth +
                  (1 | site/plot) + (1 | botanist_year),
                data = sites, 
                family = gaussian("identity"),
                prior = priors,
                chains = chains,
                iter = iter,
                thin = thin,
                warmup = floor(iter / 2),
                save_pars = save_pars(all = TRUE),
                cores = parallel::detectCores(),
                seed = 121)

m_full <- brm(n ~ (target_type + exposition + sand_ratio + survey_year_fct +
                     seed_density + substrate_depth)^3 +
                target_type:sand_ratio:exposition:survey_year_fct +
                substrate_depth:sand_ratio:exposition:survey_year_fct +
                (1 | site/plot) + (1 | botanist_year),
              data = sites, 
              family = gaussian("identity"),
              prior = priors,
              chains = chains,
              iter = iter,
              thin = thin,
              warmup = floor(iter / 2),
              save_pars = save_pars(all = TRUE),
              cores = parallel::detectCores(),
              seed = 123)

m1 <- brm(n ~ (target_type + exposition + sand_ratio + survey_year_fct)^4 +
            substrate_depth + seed_density +
            (1 | site/plot) + (1 | botanist_year),
          data = sites, 
          family = gaussian("identity"),
          prior = priors,
          chains = chains,
          iter = iter,
          thin = thin,
          warmup = floor(iter / 2),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)

m2 <- brm(n ~ (target_type + exposition + sand_ratio + survey_year_fct)^4 +
            substrate_depth + seed_density +
            sand_ratio:substrate_depth +
            substrate_depth:exposition +
            seed_density:exposition +
            (1 | site/plot) + (1 | botanist_year),
          data = sites, 
          family = gaussian("identity"),
          prior = priors,
          chains = chains,
          iter = iter,
          thin = thin,
          warmup = floor(iter / 2),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)

m3 <- brm(n ~ (target_type + exposition + sand_ratio + survey_year_fct)^2 +
            substrate_depth + seed_density +
            substrate_depth:sand_ratio +
            seed_density:exposition +
            target_type:exposition:survey_year_fct +
            sand_ratio:exposition:survey_year_fct +
            seed_density:exposition:survey_year_fct +
            (1 | site/plot) + (1 | botanist_year),
          data = sites,
          family = gaussian("identity"),
          prior = priors,
          chains = chains,
          iter = iter,
          thin = thin,
          warmup = floor(iter / 2),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)


m3_flat <- brm(n ~ (target_type + exposition + sand_ratio + survey_year_fct)^2 +
                 substrate_depth + seed_density +
                 substrate_depth:sand_ratio +
                 seed_density:exposition +
                 target_type:exposition:survey_year_fct +
                 sand_ratio:exposition:survey_year_fct +
                 seed_density:exposition:survey_year_fct +
                 (1 | site/plot) + (1 | botanist_year),
               data = sites,
               family = gaussian("identity"),
               prior = c(
                 set_prior("cauchy(0, 1)", class = "sigma")
               ),
               chains = chains,
               iter = iter,
               thin = thin,
               warmup = floor(iter / 2),
               save_pars = save_pars(all = TRUE),
               cores = parallel::detectCores(),
               seed = 123)

### Save ###
save(m_simple, file = here("outputs", "models", "model_recover_simple.Rdata"))
save(m_full, file = here("outputs", "models", "model_recover_full.Rdata"))
save(m1, file = here("outputs", "models", "model_recover_1.Rdata"))
save(m2, file = here("outputs", "models", "model_recover_2.Rdata"))
save(m3, file = here("outputs", "models", "model_recover_3.Rdata"))
save(m3_flat, file = here("outputs", "models", "model_recover_3_flat.Rdata"))
