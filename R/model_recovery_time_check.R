# Dike grassland field experiment
# Recover time ####
# Model check

# Markus Bauer
# 2022-11-15



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
sites <- read_csv(here("data", "processed", "data_processed_sites_nmds.csv"),
                  col_names = TRUE, na = c("na", "NA", ""),
                  col_types = cols(
                    .default = "?",
                    id = "f",
                    plot = "f",
                    site = "f",
                    exposition = col_factor(levels = c("north", "south", "other")),
                    sand_ratio = "f",
                    substrate_depth = col_factor(levels = c("30", "15")),
                    target_type = col_factor(levels = c("dry_grassland",
                                                        "hay_meadow", "other")),
                    seed_density = "f"
                  )) %>%
  filter(reference == "2018" | reference == "2019" | reference == "2020" |
           reference == "2021") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    survey_year = as.numeric(survey_year),
    botanist_year = str_c(survey_year, botanist, sep = " "),
    n = recovery_time
  )

### Load models ###
base::load(file = here("outputs", "models", "model_recovery_simple.Rdata"))
base::load(file = here("outputs", "models", "model_recovery_full.Rdata"))
base::load(file = here("outputs", "models", "model_recovery_1.Rdata"))
base::load(file = here("outputs", "models", "model_recovery_2.Rdata"))
base::load(file = here("outputs", "models", "model_recovery_3.Rdata"))
base::load(file = here("outputs", "models", "model_recovery_3_flat.Rdata"))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics #################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Model check ###############################################################


### a Model comparison ---------------------------------------------------------

m_1 <- m1
m_2 <- m2
rm(list = setdiff(ls(), c("sites", "m_1", "m_2")))
m_1$formula
m_2$formula
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot)) 
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot)) 
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)


### b Model check -----------------------------------------------------------

#### * DHARMa ####
DHARMa.helpers::dh_check_brms(m_1, integer = TRUE)
DHARMa.helpers::dh_check_brms(m_2, integer = TRUE)

#### * Preparation ####
posterior1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(
    variable = c(
      "b_sand_ratio25",
      "b_sand_ratio50",
      "b_substrate_depth15",
      "b_target_typedry_grassland",
      "b_seed_density8",
      "b_expositionsouth",
      "b_survey_year_fct2019",
      "b_survey_year_fct2020",
      "b_survey_year_fct2021",
      "sd_site__Intercept",
      "sd_site:plot__Intercept",
      "sd_botanist_year__Intercept",
      "sigma"
    )
  )
posterior2 <- m_2 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(
    variable = c(
      "b_sand_ratio25",
      "b_sand_ratio50",
      "b_substrate_depth15",
      "b_target_typedry_grassland",
      "b_seed_density8",
      "b_expositionsouth",
      "b_survey_year_fct2019",
      "b_survey_year_fct2020",
      "b_survey_year_fct2021",
      "sd_site__Intercept",
      "sd_site:plot__Intercept",
      "sd_botanist_year__Intercept",
      "sigma"
    )
  )
hmc_diagnostics1 <- nuts_params(m_1)
hmc_diagnostics2 <- nuts_params(m_2)
y <- sites$n
yrep1 <- posterior_predict(m_1, draws = 500)
yrep2 <- posterior_predict(m_2, draws = 500)
loo1 <- loo(m_1, save_psis = TRUE, moment_match = FALSE)
loo2 <- loo(m_2, save_psis = TRUE, moment_match = FALSE)
draws1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::summarize_draws() %>%
  filter(str_starts(variable, "b_"))
draws2 <- m_2 %>%
  posterior::as_draws() %>%
  posterior::summarize_draws() %>%
  filter(str_starts(variable, "b_"))

#### * Samling efficency/effectiveness (Rhat and EFF) ####
range(draws1$rhat)
range(draws2$rhat)
range(draws1$ess_bulk)
range(draws2$ess_bulk)
range(draws1$ess_tail)
range(draws2$ess_tail)

#### * MCMC diagnostics ####
mcmc_trace(posterior1, np = hmc_diagnostics1)
mcmc_trace(posterior2, np = hmc_diagnostics2)
mcmc_pairs(posterior1, off_diag_args = list(size = 1.2))
mcmc_pairs(posterior2, off_diag_args = list(size = 1.2))
mcmc_scatter(m_1, np = hmc_diagnostics1, size = 1,
             pars = c("b_survey_year_fct2020", "b_survey_year_fct2019"))
mcmc_scatter(m_2, np = hmc_diagnostics2, size = 1,
             pars = c("b_survey_year_fct2020", "b_survey_year_fct2019"))
mcmc_parcoord(posterior1, np = hmc_diagnostics1)
mcmc_parcoord(posterior2, np = hmc_diagnostics2)

#### * Posterior predictive check ####
#### Kernel density
p1 <- ppc_dens_overlay(y, yrep1[1:50, ])
p2 <- ppc_dens_overlay(y, yrep2[1:50, ])
p1 / p2
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$site)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$site)
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$exposition)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$exposition)
p1 / p2
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$survey_year_fct)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$survey_year_fct)
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$target_type)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$target_type)
p1 / p2
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$seed_density)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$seed_density)
p1 / p2
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$sand_ratio)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$sand_ratio)
p1 / p2
p1 <- ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$substrate_depth)
p2 <- ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$substrate_depth)
p1 / p2
#### Histograms of statistics skew
p1 <- ppc_stat(y, yrep1, binwidth = 0.001)
p2 <- ppc_stat(y, yrep2, binwidth = 0.001)
p1 / p2
ppc_stat_grouped(y, yrep1, group = sites$site, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$site, binwidth = 0.001)
p1 <- ppc_stat_grouped(y, yrep1, group = sites$exposition, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$exposition, binwidth = 0.001)
p1 / p2
ppc_stat_grouped(y, yrep1, group = sites$survey_year_fct, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$survey_year_fct, binwidth = 0.001)
p1 <- ppc_stat_grouped(y, yrep1, group = sites$target_type, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$target_type, binwidth = 0.001)
p1 / p2
p1 <- ppc_stat_grouped(y, yrep1, group = sites$seed_density, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$seed_density, binwidth = 0.001)
p1 / p2
p1 <- ppc_stat_grouped(y, yrep1, group = sites$sand_ratio, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$sand_ratio, binwidth = 0.001)
p1 / p2
p1 <- ppc_stat_grouped(y, yrep1, group = sites$substrate_depth, binwidth = 0.001)
p2 <- ppc_stat_grouped(y, yrep2, group = sites$substrate_depth, binwidth = 0.001)
p1 / p2
#### LOO-PIT plots
p1 <- ppc_loo_pit_overlay(y, yrep1, lw = weights(loo1$psis_object))
p2 <- ppc_loo_pit_overlay(y, yrep2, lw = weights(loo2$psis_object))
p1 / p2
plot(loo1)
plot(loo2)


#### * Autocorrelation ####
mcmc_acf(posterior1, lags = 10)
mcmc_acf(posterior2, lags = 10)



## 2 Output of chosen model ####################################################


### a Model output ------------------------------------------------------------

prior_summary(m_1, all = FALSE)
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot) ) 
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
draws1
mcmc_intervals(
  posterior1,
  prob = 0.66,
  prob_outer = 0.95,
  point_est = "mean"
)
mcmc_intervals(
  posterior2,
  prob = 0.66,
  prob_outer = 0.95,
  point_est = "mean"
)


### b Effect sizes ------------------------------------------------------------

(emm <- emmeans(m_1, revpairwise ~ target_type + sand_ratio |
                  exposition | survey_year_fct, type = "response"))

write.csv(draws1, here("outputs", "statistics", "table_recovery_time.csv"))
