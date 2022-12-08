# Dike grassland field experiment
# Favourable conservation status ####
# Model building

# Markus Bauer
# 2022-12-06



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
  here("data", "processed", "data_processed_sites.csv"),
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
    n = fcs_target,
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


### b Outliers, zero-inflation, transformations? ----------------------------

sites %>% group_by(exposition) %>% count(site)
boxplot(sites$n)
ggplot(sites, aes(x = exposition, y = n)) + geom_quasirandom()
ggplot(sites, aes(x = n)) + geom_histogram(binwidth = 0.03)
ggplot(sites, aes(x = n)) + geom_density()



## 2 Model building ###########################################################


### a Possible priors ----------------------------------------------------------

get_prior(n ~ target_type + exposition + sand_ratio + survey_year_fct +
            seed_density + substrate_depth +
            (1 | site/plot) + (1 | botanist_year),
          data = sites)

ggplot(data = data.frame(x = c(-5, 5)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 2)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for Intercept")
ggplot(data = data.frame(x = c(-5, 5)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0.3, sd = 2)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for treatments")
ggplot(data = data.frame(x = c(-5, 5)), aes(x = x)) +
  stat_function(fun = dcauchy, n = 101, args = list(location = 0, scale = 1)) +
  expand_limits(y = 0) +
  ggtitle("Cauchy distribution")
ggplot(data.frame(x = c(-5, 5)), aes(x = x)) +
  stat_function(fun = dstudent_t, args = list(df = 3, mu = 0, sigma = 2.5)) +
  expand_limits(y = 0) +
  ggtitle(expression(Student~italic(t)*"-distribution"))


### b Model specifications -----------------------------------------------------

# NUTS sampler used
iter = 10000
chains = 4
thin = 2
seed = 123
warmup = floor(iter / 2)
priors <- c(
  set_prior("normal(0, 2)", class = "Intercept"),
  set_prior("normal(0, 2)", class = "b"),
  set_prior("normal(0.1, 2)", class = "b", coef = "sand_ratio25"),
  set_prior("normal(0.2, 2)", class = "b", coef = "sand_ratio50"),
  set_prior("normal(0.1, 2)", class = "b", coef = "survey_year_fct2019"),
  set_prior("normal(0.2, 2)", class = "b", coef = "survey_year_fct2020"),
  set_prior("normal(0.3, 2)", class = "b",coef = "survey_year_fct2021"),
  set_prior("cauchy(0, 1)", class = "sigma")
)


### c Models ------------------------------------------------------------------

m_simple <- brm(
  n ~ sand_ratio + target_type + exposition + survey_year_fct +
    substrate_depth + seed_density +
    botanist_year + (1 | site/plot),
  data = sites, 
  family = gaussian("identity"),
  prior = priors,
  chains = chains,
  iter = iter,
  thin = thin,
  control = list(max_treedepth = 13),
  warmup = warmup,
  save_pars = save_pars(all = TRUE),
  cores = parallel::detectCores(),
  seed = seed
)

m_full <- brm(
  n ~ sand_ratio * target_type * exposition * survey_year_fct +
    substrate_depth * seed_density +
    substrate_depth:exposition +
    seed_density:exposition +
    substrate_depth:survey_year_fct +
    seed_density:survey_year_fct +
    substrate_depth:exposition:survey_year_fct +
    seed_density:exposition:survey_year_fct +
    botanist_year + (1 | site/plot),
  data = sites, 
  family = gaussian("identity"),
  prior = priors,
  chains = chains,
  iter = iter,
  thin = thin,
  control = list(max_treedepth = 13),
  warmup = warmup,
  save_pars = save_pars(all = TRUE),
  cores = parallel::detectCores(),
  seed = seed
)

m1 <- brm(
  n ~ sand_ratio * substrate_depth * exposition * survey_year_fct +
    target_type + seed_density +
    botanist_year + (1 | site/plot),
  data = sites, 
  family = gaussian("identity"),
  prior = priors,
  chains = chains,
  iter = iter,
  thin = thin,
  control = list(max_treedepth = 13),
  warmup = floor(iter / 2),
  save_pars = save_pars(all = TRUE),
  cores = parallel::detectCores(),
  seed = seed
)

m2 <- brm(
  n ~ sand_ratio * target_type * exposition * survey_year_fct +
    substrate_depth + seed_density +
    substrate_depth:exposition +
    seed_density:exposition +
    substrate_depth:survey_year_fct +
    seed_density:survey_year_fct +
    botanist_year + (1 | site/plot),
  data = sites, 
  family = gaussian("identity"),
  prior = priors,
  chains = chains,
  iter = iter,
  thin = thin,
  control = list(max_treedepth = 13),
  warmup = warmup,
  save_pars = save_pars(all = TRUE),
  cores = parallel::detectCores(),
  seed = seed
)

m3 <- brm(
  n ~ (sand_ratio + target_type + seed_density + substrate_depth) *
    exposition * survey_year_fct +
    botanist_year + (1 | site/plot),
  data = sites,
  family = gaussian("identity"),
  prior = priors,
  chains = chains,
  iter = iter,
  thin = thin,
  control = list(max_treedepth = 13),
  warmup = warmup,
  save_pars = save_pars(all = TRUE),
  cores = parallel::detectCores(),
  seed = seed
)

m2_flat <- brm(
  n ~ sand_ratio * target_type * exposition * survey_year_fct +
    substrate_depth + seed_density +
    substrate_depth:exposition +
    seed_density:exposition +
    substrate_depth:survey_year_fct +
    seed_density:survey_year_fct +
    botanist_year + (1 | site/plot),
  data = sites, 
  family = gaussian("identity"),
  prior = c(
    set_prior("normal(0, 4)", class = "Intercept"),
    set_prior("normal(0, 4)", class = "b"),
    set_prior("cauchy(0, 1)", class = "sigma")
  ),
  chains = chains,
  iter = iter,
  thin = thin,
  control = list(max_treedepth = 13),
  warmup = warmup,
  save_pars = save_pars(all = TRUE),
  cores = parallel::detectCores(),
  seed = seed
)

m2_prior <- brm(
  n ~ sand_ratio * target_type * exposition * survey_year_fct +
    substrate_depth + seed_density +
    substrate_depth:exposition +
    seed_density:exposition +
    substrate_depth:survey_year_fct +
    seed_density:survey_year_fct +
    botanist_year + (1 | site/plot),
  data = sites, 
  family = gaussian("identity"),
  prior = priors,
  sample_prior = "only",
  chains = chains,
  iter = iter,
  thin = thin,
  control = list(max_treedepth = 13),
  warmup = warmup,
  save_pars = save_pars(all = TRUE),
  cores = parallel::detectCores(),
  seed = seed
)


### d Save ---------------------------------------------------------------------

save(m_simple, file = here("outputs", "models", "model_fcs_simple.Rdata"))
save(m_full, file = here("outputs", "models", "model_fcs_full.Rdata"))
save(m1, file = here("outputs", "models", "model_fcs_1.Rdata"))
save(m2, file = here("outputs", "models", "model_fcs_2.Rdata"))
save(m3, file = here("outputs", "models", "model_fcs_3.Rdata"))
save(m2_flat, file = here("outputs", "models", "model_fcs_2_flat.Rdata"))
save(m2_prior, file = here("outputs", "models", "model_fcs_2_prior.Rdata"))
