# Dike grassland field experiment
# Persistence ####
# Model building

# Markus Bauer
# 2022-11-15



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(brms)

### Start ###
rm(list = ls())

### Load data ###
sites <- read_csv(
  here("data", "processed", "data_processed_sites_temporal.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types =
    cols(
      .default = "?",
      plot = "f",
      site = "f",
      sand_ratio = "f",
      substrate_depth = col_factor(levels = c("30", "15")),
      target_type = col_factor(levels = c(
        "hay_meadow", "dry_grassland"
      )),
      seed_density = "f",
      exposition = col_factor(levels = c(
        "north", "south"
      )),
      survey_year = "c"
    )
  ) %>%
  ### Exclude data of seed mixtures
  filter(presabu == "presence") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    id = factor(id),
    n = persistence
  ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, survey_year, n
  )



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics #################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Data exploration ##########################################################


### a Graphs of raw data -------------------------------------------------------

plot1 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = sand_ratio)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  facet_grid(~ survey_year_fct) +
  labs(title = "Sand ratio [vol%]")
plot2 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = substrate_depth)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  facet_grid(~ survey_year_fct) +
  labs(title = "Substrate depth [cm]")
plot3 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = target_type)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  facet_grid(~ survey_year_fct) +
  labs(title = "Target type")
plot4 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = seed_density)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  facet_grid(~ survey_year_fct) +
  labs(title = "Seed density [g/m²]")
(plot1 + plot2) / (plot3 + plot4)
plot1 <- ggplot(sites %>% filter(survey_year == 2021),
                aes(y = n, x = exposition)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  facet_grid(~ survey_year_fct) +
  labs(title = "Exposition")
plot2 <- ggplot(sites, aes(y = n, x = survey_year_fct)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  labs(title = "Survey year")
plot3 <- ggplot(sites %>% filter(survey_year == 2021), aes(y = n, x = site)) +
  geom_quasirandom(color = "grey") + geom_boxplot(fill = "transparent") +
  facet_grid(~ survey_year_fct) +
  labs(title = "Blocks")
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


### b Outliers, zero-inflation, transformations? ------------------------------

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
            (1 | site/plot),
          data = sites)
### Example of normal distribution
ggplot(data = data.frame(x = c(-40, 40)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 10)) +
  expand_limits(y = 0)
### Example of cauchy distribution
ggplot(data = data.frame(x = c(-40, 40)), aes(x = x)) +
  stat_function(fun = dcauchy, n = 101, args = list(location = 0, scale = 10)) +
  expand_limits(y = 0)
### Example of a student t distribution
ggplot(data.frame(x = c(-40, 40)), aes(x = x)) +
  stat_function(fun = dstudent_t, args = list(df = 3, mu = 0, sigma = 10)) +
  expand_limits(y = 0)


### b Models ------------------------------------------------------------------

### Specifications ###
iter = 10000
chains = 4
thin = 2
priors <- c(
  set_prior("normal(0, 20)", class = "b"),
  set_prior("normal(-2.5, 20)", class = "b", coef = "sand_ratio25"),
  set_prior("normal(-5, 20)", class = "b", coef = "sand_ratio50"),
  set_prior("normal(5, 20)", class = "b", coef = "expositionsouth"),
  set_prior("normal(-2.5, 20)", class = "b", coef = "survey_year_fct2019"),
  set_prior("normal(-5, 20)", class = "b", coef = "survey_year_fct2020"),
  set_prior("normal(-7.5, 20)", class = "b", coef = "survey_year_fct2021"),
  set_prior("cauchy(0, 10)", class = "sigma")
)
### Models ###
m_simple <- brm(n ~ target_type + exposition + sand_ratio + survey_year_fct +
                  seed_density + substrate_depth +
                  (1 | site/plot),
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
                (1 | site/plot),
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
            (1 | site/plot),
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
            (1 | site/plot),
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
            (1 | site/plot),
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


m1_flat <- brm(n ~ (target_type + exposition + sand_ratio + survey_year_fct)^4 +
                 substrate_depth + seed_density +
                 (1 | site/plot),
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


### c Save ---------------------------------------------------------------------

save(m_simple, file = here(
  "outputs", "models", "model_persistence_simple.Rdata"
  ))
save(m_full, file = here("outputs", "models", "model_persistence_full.Rdata"))
save(m1, file = here("outputs", "models", "model_persistence_1.Rdata"))
save(m2, file = here("outputs", "models", "model_persistence_2.Rdata"))
save(m3, file = here("outputs", "models", "model_persistence_3.Rdata"))
save(m1_flat, file = here(
  "outputs", "models", "model_persistence_1_flat.Rdata"
  ))
