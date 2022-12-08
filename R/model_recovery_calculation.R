# Dike grassland field experiment
# Recovery completeness ####
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
  here("data", "processed", "data_processed_sites_nmds.csv"),
  col_names = TRUE, na = c("na", "NA", ""),
  col_types = cols(
    .default = "?",
    id = "f",
    plot = "f",
    site = "f",
    exposition = col_factor(levels = c("north", "south", "other")),
    sand_ratio = "f",
    substrate_depth = col_factor(levels = c("30", "15")),
    target_type = col_factor(levels = c(
      "hay_meadow", "dry_grassland", "other"
      )),
    seed_density = "f"
    )
  ) %>%
  filter(reference == "2018" | reference == "2019" | reference == "2020" |
           reference == "2021") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    survey_year = as.numeric(survey_year),
    botanist_year = str_c(survey_year, botanist, exposition, sep = " "),
    botanist_year = factor(botanist_year),
    n = recovery_time
    )

rm(list = setdiff(ls(), c("sites")))



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
            (1 | site/plot) + (1|botanist_year),
          data = sites)
ggplot(data = data.frame(x = c(-2, 0)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = -1, sd = 1)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for Intecept")
ggplot(data = data.frame(x = c(-2, 2)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = .8)) +
  expand_limits(y = 0) +
  ggtitle("Normal distribution for treatments")
ggplot(data = data.frame(x = c(-2, 2)), aes(x = x)) +
  stat_function(fun = dcauchy, n = 101, args = list(location = 0, scale = 1)) +
  expand_limits(y = 0) +
  ggtitle("Cauchy distribution") # See Lemoine 2019 https://doi.org/10.1111/oik.05985
ggplot(data.frame(x = c(-2, 2)), aes(x = x)) +
  stat_function(fun = dstudent_t, args = list(df = 3, mu = 0, sigma = 2.5)) +
  expand_limits(y = 0) +
  ggtitle(expression(Student~italic(t)*"-distribution")) # Software standard


### b Model specifications -----------------------------------------------------

# NUTS sampler used
iter = 20000
chains = 4
thin = 2
seed = 123
warmup = floor(iter / 2)
priors <- c(
  set_prior("normal(-1, 1)", class = "Intercept"),
  set_prior("normal(0, .8)", class = "b"),
  set_prior("normal(0.1, .8)", class = "b", coef = "sand_ratio25"),
  set_prior("normal(0.2, .8)", class = "b", coef = "sand_ratio50"),
  set_prior("normal(-0.1, .8)", class = "b", coef = "expositionsouth"),
  set_prior("normal(0.1, .8)", class = "b", coef = "survey_year_fct2019"),
  set_prior("normal(0.2, .8)", class = "b", coef = "survey_year_fct2020"),
  set_prior("normal(0.3, .8)", class = "b",coef = "survey_year_fct2021"),
  set_prior("normal(0.3, .8)", class = "sd"),
  set_prior("cauchy(0, 1)", class = "sigma")
)


### c Models ------------------------------------------------------------------

m_simple <- brm(
  n ~ sand_ratio + target_type + exposition + survey_year_fct +
    substrate_depth + seed_density +
    (1 | site/plot) + (1 | botanist_year),
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
    botanist_year +
    (1 | site/plot),
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
    botanist_year +
    (1 | site/plot),
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
    botanist_year +
    (1 | site/plot),
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
    botanist_year +
    (1 | site/plot),
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
    substrate_depth * seed_density +
    substrate_depth:exposition +
    seed_density:exposition +
    substrate_depth:survey_year_fct +
    seed_density:survey_year_fct +
    substrate_depth:exposition:survey_year_fct +
    seed_density:exposition:survey_year_fct +
    botanist_year +
    (1 | site/plot),
  data = sites, 
  family = gaussian("identity"),
  prior = c(
    set_prior("normal(0, 3)", class = "b"),
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
    substrate_depth * seed_density +
    substrate_depth:exposition +
    seed_density:exposition +
    substrate_depth:survey_year_fct +
    seed_density:survey_year_fct +
    substrate_depth:exposition:survey_year_fct +
    seed_density:exposition:survey_year_fct +
    botanist_year +
    (1 | site/plot),
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

save(m_simple, file = here("outputs", "models", "model_recovery_simple.Rdata"))
save(m_full, file = here("outputs", "models", "model_recovery_full.Rdata"))
save(m1, file = here("outputs", "models", "model_recovery_1.Rdata"))
save(m2, file = here("outputs", "models", "model_recovery_2.Rdata"))
save(m3, file = here("outputs", "models", "model_recovery_3.Rdata"))
save(m2_flat, file = here("outputs", "models", "model_recovery_2_flat.Rdata"))
save(m2_prior, file = here("outputs", "models", "model_recovery_2_prior.Rdata"))
