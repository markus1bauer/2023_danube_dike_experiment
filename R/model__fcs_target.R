# Dike grassland field experiment
# Favourable conservation status ####
# Markus Bauer
# 2022-10-18



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
setwd(here("data", "processed"))
rm(list = setdiff(ls(), c(
  "m_simple", "m_full", "m3", "m3_flat", "m_simple_flat", "m3_iter"
  )))

### Load data ###
sites <- read_csv("data_processed_sites.csv",
                  col_names = TRUE, na = c("na", "NA", ""), col_types =
                    cols(
                      .default = "?",
                      plot = "f",
                      site = "f",
                      survey_year = "d",
                      exposition = "f",
                      sand_ratio = "f",
                      substrate_depth = "f",
                      seed_density = "f",
                      target_type = "f"
                    )) %>%
  filter(survey_year != "seeded") %>%
  mutate(
    n = fcs_target,
    survey_year_fct = factor(survey_year),
    botanist_year = str_c(botanist, survey_year, sep = " "),
    botanist_year = factor(botanist_year),
    id = factor(id)
  ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, survey_year, botanist_year, n
  )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Data exploration ########################################################


### a Graphs -----------------------------------------------------------------

#simple effects:
ggplot(sites %>% filter(survey_year == 2021), aes(y = n, x = exposition)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites %>% filter(survey_year == 2020), aes(y = n, x = exposition)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites %>% filter(survey_year == 2021), aes(y = n, x = target_type)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(survey_year == 2021), aes(y = n, x = seed_density)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(survey_year == 2021), aes(y = n, x = substrate_depth)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(survey_year == 2021), aes(y = n, x = sand_ratio)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites %>% filter(survey_year == 2021), aes(y = n, x = site)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = botanist_year)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites, aes(y = n, x = survey_year_fct)) +
  geom_boxplot() + geom_quasirandom()
#2way
ggplot(sites %>% filter(survey_year == 2021), aes(x = exposition, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ target_type)
ggplot(sites %>% filter(survey_year == 2021), aes(x = sand_ratio, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ target_type)
ggplot(sites %>% filter(survey_year == 2021), aes(x = substrate_depth, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ target_type)
ggplot(sites %>% filter(survey_year == 2021), aes(x = sand_ratio, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ exposition)
ggplot(sites, aes(x = survey_year_fct, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ site)
#3way
ggplot(sites, aes(x = exposition, y = n, color = sand_ratio)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ target_type)
ggplot(sites %>% filter(survey_year == 2021),
       aes(x = exposition, y = n, color = sand_ratio)) +
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ target_type)
ggplot(sites %>% filter(survey_year == 2021),
       aes(x = exposition, y = n, color = substrate_depth)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ target_type)
ggplot(sites %>% filter(survey_year == 2021),
       aes(x = substrate_depth, y = n, color = sand_ratio)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ target_type)
ggplot(sites %>% filter(survey_year == 2021),
       aes(x = substrate_depth, y = n, color = sand_ratio)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ exposition)
ggplot(sites, aes(x = factor(survey_year), y = n, color = target_type)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ exposition) # Plot for figure
#4way
ggplot(sites %>% filter(survey_year == 2021),
       aes(x = exposition, y = n, color = sand_ratio)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_grid(substrate_depth ~ target_type)
ggplot(data = sites,
       aes(x = survey_year_fct, y = n, color = target_type)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_grid(exposition ~ sand_ratio)


#### b Outliers, zero-inflation, transformations? ----------------------------

dotchart((sites$n), groups = factor(sites$exposition),
         main = "Cleveland dotplot")
sites %>% group_by(exposition) %>% count(site)
boxplot(sites$n)
plot(table((sites$n)), type = "h",
     xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(x = n)) + geom_density()
ggplot(sites, aes(sqrt(n))) + geom_density()



## 2 Model building ###########################################################


### a models -----------------------------------------------------------------

iter = 10000
chains = 4
thin = 2
priors <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(0.1, 1)", class = "b", coef = "sand_ratio25"),
  set_prior("normal(0.2, 1)", class = "b", coef = "sand_ratio50"),
  set_prior("normal(0.1, 1)", class = "b", coef = "expositionnorth"),
  set_prior("normal(0.1, 1)", class = "b", coef = "survey_year_fct2019"),
  set_prior("normal(0.2, 1)", class = "b", coef = "survey_year_fct2020"),
  set_prior("normal(0.3, 1)", class = "b",coef = "survey_year_fct2021"),
  set_prior("cauchy(0, 1)", class = "sigma")
)

### Posssible priors ###
get_prior(n ~ target_type + exposition + sand_ratio + survey_year_fct +
            seed_density + substrate_depth +
            (1 | site/plot) + (1 | botanist_year),
          data = sites)
### Example of normal distribution ###
ggplot(data = data.frame(x = c(-2, 2)), aes(x = x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0.1, sd = 1))
### Example of cauchy distribution ###
ggplot(data = data.frame(x = c(-2, 2)), aes(x = x)) +
  stat_function(fun = dcauchy, n = 101, args = list(location = 0, scale = 1))
### Example of a student t distribution ###
ggplot(data.frame(x = c(-2, 2)), aes(x = x)) +
  stat_function(fun = dstudent_t, args = list(df = 3, mu = 0, sigma = 2.5))


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

m_full_flat <- brm(n ~ (target_type + exposition + sand_ratio + survey_year_fct +
                      seed_density + substrate_depth)^3 +
                 target_type:sand_ratio:exposition:survey_year_fct +
                 substrate_depth:sand_ratio:exposition:survey_year_fct +
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


### b comparison ------------------------------------------------------------

m_1 <- m1
m_2 <- m3
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot) + (1 | botanist_year)) 
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot) + (1 | botanist_year)) 
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)


### c model check -----------------------------------------------------------

#### * DHARMa ####
createDHARMa(
  simulatedResponse = t(posterior_predict(m_1)),
  observedResponse = sites$n,
  fittedPredictedResponse = apply(t(posterior_epred(m_1)), 1, mean),
  integerResponse = TRUE
  ) %>%
  plot()
createDHARMa(
  simulatedResponse = t(posterior_predict(m_2)),
  observedResponse = sites$n,
  fittedPredictedResponse = apply(t(posterior_epred(m_2)), 1, mean),
  integerResponse = TRUE
  ) %>%
  plot()

#### * Preparation ####
posterior1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(
    variable = c(
      "b_sand_ratio25",
      "b_sand_ratio50",
      "b_substrate_depth30",
      "b_target_typedry_grassland",
      "b_seed_density8",
      "b_expositionnorth",
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
      "b_substrate_depth30",
      "b_target_typedry_grassland",
      "b_seed_density8",
      "b_expositionnorth",
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
mcmc_scatter(m_1,
             pars = c("b_survey_year_fct2020", "b_survey_year_fct2019"),
             np = hmc_diagnostics1,
             size = 1)
mcmc_scatter(m_2,
             pars = c("b_survey_year_fct2020", "b_survey_year_fct2019"),
             np = hmc_diagnostics2,
             size = 1)
mcmc_parcoord(posterior1, np = hmc_diagnostics1)
mcmc_parcoord(posterior2, np = hmc_diagnostics2)

#### * Posterior predictive check ####
#### Kernel density
ppc_dens_overlay(y, yrep1[1:50, ])
ppc_dens_overlay(y, yrep2[1:50, ])
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$site)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$site)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$exposition)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$exposition)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$survey_year_fct)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$survey_year_fct)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$target_type)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$target_type)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$seed_density)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$seed_density)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$sand_ratio)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$sand_ratio)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$substrate_depth)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$substrate_depth)
#### Histograms of statistics skew
ppc_stat(y, yrep1, binwidth = 0.001)
ppc_stat(y, yrep2, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$site, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$site, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$exposition, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$exposition, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$survey_year_fct, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$survey_year_fct, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$target_type, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$target_type, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$seed_density, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$seed_density, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$sand_ratio, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$sand_ratio, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$substrate_depth, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$substrate_depth, binwidth = 0.001)
#### LOO-PIT plots
ppc_loo_pit_overlay(y, yrep1, lw = weights(loo1$psis_object))
ppc_loo_pit_overlay(y, yrep2, lw = weights(loo2$psis_object))
plot(loo1)
plot(loo2)

#### * Autocorrelation ####
mcmc_acf(posterior1, lags = 10)
mcmc_acf(posterior2, lags = 10)



## 3 Chosen model output #####################################################


### a Model output ------------------------------------------------------------

prior_summary(m_1)
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot) + (1 | botanist_year)) 
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | site/plot) + (1 | botanist_year)) 
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
draws1
mcmc_intervals(
  posterior1,
  prob = 0.5,
  prob_outer = 0.95,
  point_est = "mean"
)
sjPlot::plot_model(m_1, m_2, type = "est")
sjPlot::plot_model(m_2, type = "est")
sjPlot::plot_model(
  m_1, m_2, type = "est", ppd = TRUE,
  terms = c(
    "sand_ratio", "substrate_depth", "target_type",
    "seed_density",  "exposition", "survey_year_fct"
    )
  )


### b Effect sizes ------------------------------------------------------------

(emm <- emmeans(m_1, revpairwise ~ target_type + sand_ratio |
                  exposition | survey_year_fct, type = "response"))
(emm <- emmeans(m_1, revpairwise ~ substrate_depth + sand_ratio |
                  exposition | survey_year_fct, type = "response"))
(emm <- emmeans(m_1, revpairwise ~ seed_density |
                  exposition | survey_year_fct, type = "response"))

### Save ###
save(m3, file = here("data", "processed", "model_fcs_3.Rdata"))
save(m3_flat, file = here("data", "processed", "model_fcs_3_flat.Rdata"))
