# Dike grassland experiment
# Persistence Presence-Absence ####
# Markus Bauer
# 2022-05-13



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(ggbeeswarm)
library(brms)
library(DHARMa)
library(bayesplot)
library(tidybayes)
library(emmeans)

### Start ###
#rm(list = ls())
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
                      exposition = col_factor(levels = c("north", "south")),
                      sandRatio = "f",
                      substrateDepth = "f",
                      targetType = "f",
                      seedDensity = "f"
                    )) %>%
  filter(
    !str_detect(id, "C") & presabu == "presence" & surveyYear != "seeded"
  ) %>%
  select(
    id, plot, block, exposition, sandRatio, substrateDepth, targetType,
    seedDensity, surveyYear, B, C, D
  ) %>%
  pivot_longer(cols = c(B, C, D), names_to = "index", values_to = "n",
               names_transform = as.factor) %>%
  mutate(
    index = fct_rev(index),
    surveyYear = as.double(surveyYear),
    surveyYear_fac = as.factor(surveyYear)
  )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Data exploration ########################################################

### a Graphs -----------------------------------------------------------------
#simple effects:
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(y = n, x = exposition)) +
  geom_boxplot() + geom_quasirandom() # 
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(y = n, x = targetType)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(y = n, x = seedDensity)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(y = n, x = substrateDepth)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(y = n, x = sandRatio)) +
  geom_boxplot() + geom_quasirandom() # 
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(y = n, x = block)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites %>% filter(index == "B"), aes(y = n, x = surveyYear_fac)) +
  geom_boxplot() + geom_quasirandom() #
#2way
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(x = targetType, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ exposition)
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(x = sandRatio, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ exposition)
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(x = substrateDepth, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ exposition)
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(x = seedDensity, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ exposition) # 
#3way
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(x = sandRatio, y = n, color = substrateDepth)) + 
  geom_boxplot() + geom_quasirandom(dodge.width = .8, alpha = .5) +
  facet_wrap(~ exposition)
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(x = seedDensity, y = n, color = targetType)) + 
  geom_boxplot() + geom_quasirandom(dodge.width = .8, alpha = .5) +
  facet_wrap(~ exposition)
ggplot(sites %>% filter(surveyYear == 2021, index == "B"),
       aes(x = sandRatio, y = n, color = targetType)) + 
  geom_boxplot() + geom_quasirandom(dodge.width = .8, alpha = .5) +
  facet_wrap(~ exposition)
ggplot(sites %>% filter(index == "B"),
       aes(x = factor(surveyYear), y = n, color = sandRatio)) + 
  geom_boxplot() + geom_quasirandom(dodge.width = .8, alpha = .5) +
  facet_wrap(~ exposition) # Plot for figure
#4way
ggplot(sites %>% filter(index == "B"),
       aes(x = sandRatio, y = n, color = targetType)) + 
  geom_boxplot() + geom_quasirandom(dodge.width = .8, alpha = .5) +
  facet_grid(exposition ~ factor(surveyYear))

##### b Outliers, zero-inflation, transformations? ----------------------------
dotchart((sites$n), groups = factor(sites$index),
         main = "Cleveland dotplot")
sites %>% count(block)
boxplot(sites$n)
plot(table((sites$n)), type = "h",
     xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(x = n)) + geom_density()
ggplot(sites, aes(x = sqrt(n))) + geom_density()


## 2 Model building ###########################################################

### a models -----------------------------------------------------------------

iter = 10000
chains = 4
thin = 2

m_simple <- brm(n ~ index + exposition + surveyYear_fac +
                  targetType + sandRatio + seedDensity + substrateDepth +
                  (1 | block/plot),
                data = sites, 
                family = gaussian("identity"),
                prior = c(
                  set_prior("normal(0, 0.6)", class = "b"),
                  set_prior("cauchy(0, 1)", class = "sigma")
                ),
                chains = chains,
                iter = iter,
                thin = thin,
                warmup = floor(iter / 2),
                save_pars = save_pars(all = TRUE),
                cores = parallel::detectCores(),
                seed = 123)

m_full <- brm(n ~ (index + exposition + surveyYear_fac +
                     targetType + sandRatio + seedDensity + substrateDepth)^3 +
                targetType:exposition:surveyYear_fac:index +
                sandRatio:exposition:surveyYear_fac:index +
                (1 | block/plot),
              data = sites, 
              family = gaussian("identity"),
              prior = c(
                set_prior("normal(0, 0.6)", class = "b"),
                set_prior("cauchy(0, 1)", class = "sigma")
              ),
              chains = chains,
              iter = iter,
              thin = thin,
              warmup = floor(iter / 2),
              save_pars = save_pars(all = TRUE),
              cores = parallel::detectCores(),
              seed = 123)

m31 <- brm(n ~ (index + exposition + surveyYear_fac +
                 targetType + sandRatio)^3 + seedDensity + substrateDepth +
            targetType:exposition:surveyYear_fac:index +
            sandRatio:exposition:surveyYear_fac:index +
            (1 | block/plot),
          data = sites, 
          family = gaussian("identity"),
          prior = c(
            set_prior("normal(0, 0.6)", class = "b"),
            set_prior("cauchy(0, 1)", class = "sigma")
          ),
          chains = chains,
          iter = iter,
          thin = thin,
          control = list(adapt_delta = 0.99, max_treedepth = 12),
          warmup = floor(iter / 2),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)

m_full_flat <- brm(n ~ (index + exposition + surveyYear_fac +
                          targetType + sandRatio + seedDensity + substrateDepth)^3 +
                     targetType:exposition:surveyYear_fac:index +
                     sandRatio:exposition:surveyYear_fac:index +
                     (1 | block/plot),
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

m3_flat <- brm(n ~ (index + exposition + surveyYear_fac +
                      targetType + sandRatio)^3 + seedDensity + substrateDepth +
                 targetType:exposition:surveyYear_fac:index +
                 sandRatio:exposition:surveyYear_fac:index +
                 (1 | block/plot),
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

m_1 <- m3_flat
m_2 <- m3
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | block/plot) + (1 | botanist_year)) 
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | block/plot) + (1 | botanist_year)) 
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)

### c model check -----------------------------------------------------------

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

### * Preparation for evaluation ####
posterior1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(
    variable = c(
      "b_substrateDepth30",
      "b_seedDensity8",
      "b_targetTypedry_grassland",
      "b_sandRatio25",
      "b_sandRatio50",
      "b_expositionsouth",
      "sd_block__Intercept",
      "sd_block:plot__Intercept",
      "sigma"
    )
  )
posterior2 <- m_2 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(
    variable = c(
      "b_substrateDepth30",
      "b_seedDensity8",
      "b_targetTypedry_grassland",
      "b_sandRatio25",
      "b_sandRatio50",
      "b_expositionsouth",
      "sd_block__Intercept",
      "sd_block:plot__Intercept",
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

### * Samling efficency/effectiveness (Rhat and EFF) ####
range(draws1$rhat)
range(draws2$rhat)
range(draws1$ess_bulk)
range(draws2$ess_bulk)
range(draws1$ess_tail)
range(draws2$ess_tail)

### * MCMC diagnostics ####
mcmc_trace(posterior1, np = hmc_diagnostics1)
mcmc_trace(posterior2, np = hmc_diagnostics2)
mcmc_pairs(posterior1, off_diag_args = list(size = 1.2))
mcmc_pairs(posterior2, off_diag_args = list(size = 1.2))
mcmc_scatter(m_1,
             pars = c("b_surveyYear_fac2020", "b_surveyYear_fac2019"),
             np = hmc_diagnostics1,
             size = 1)
mcmc_scatter(m_2,
             pars = c("b_surveyYear_fac2020", "b_surveyYear_fac2019"),
             np = hmc_diagnostics2,
             size = 1)
mcmc_parcoord(posterior1, np = hmc_diagnostics1)
mcmc_parcoord(posterior2, np = hmc_diagnostics2)

### * Posterior predictive check ####
#### Kernel density
ppc_dens_overlay(y, yrep1[1:50, ])
ppc_dens_overlay(y, yrep2[1:50, ])
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$block)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$block)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$exposition)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$exposition)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$surveyYear_fac)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$surveyYear_fac)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$targetType)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$targetType)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$seedDensity)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$seedDensity)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$sandRatio)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$sandRatio)
ppc_dens_overlay_grouped(y, yrep1[1:50, ], group = sites$substrateDepth)
ppc_dens_overlay_grouped(y, yrep2[1:50, ], group = sites$substrateDepth)
#### Histograms of statistics skew
ppc_stat(y, yrep1, binwidth = 0.001)
ppc_stat(y, yrep2, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$block, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$block, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$exposition, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$exposition, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$surveyYear_fac, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$surveyYear_fac, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$targetType, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$targetType, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$seedDensity, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$seedDensity, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$sandRatio, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$sandRatio, binwidth = 0.001)
ppc_stat_grouped(y, yrep1, group = sites$substrateDepth, binwidth = 0.001)
ppc_stat_grouped(y, yrep2, group = sites$substrateDepth, binwidth = 0.001)
#### LOO-PIT plots
ppc_loo_pit_overlay(y, yrep1, lw = weights(loo1$psis_object))
ppc_loo_pit_overlay(y, yrep2, lw = weights(loo2$psis_object))
plot(loo1)
plot(loo2)

### * Autocorrelation ####
mcmc_acf(posterior1, lags = 10)
mcmc_acf(posterior2, lags = 10)


## 3 Chosen model output #####################################################

### a Model output ------------------------------------------------------------

prior_summary(m3_flat)
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | block/plot) + (1 | botanist_year)) 
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 | block/plot) + (1 | botanist_year)) 
bayes_R2(m_1, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
bayes_R2(m_2, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
draws1
mcmc_intervals(
  posterior2,
  prob = 0.5,
  prob_outer = 0.95,
  point_est = "mean"
)
sjPlot::plot_model(m3, type = "est")
sjPlot::plot_model(m3_flat, type = "est")
sjPlot::plot_model(m_1, type = "pred", ppd = TRUE, terms = c(
  "sandRatio", "index", "surveyYear_fac", "exposition"
))
sjPlot::plot_model(m_1, type = "pred", ppd = TRUE, terms = c(
  "targetType", "index", "surveyYear_fac", "exposition"
))


### b Effect sizes ------------------------------------------------------------

(emm <- emmeans(m_1, revpairwise ~ targetType + sandRatio |
                  exposition | surveyYear_fac | index, type = "response"))
(emm <- emmeans(m_1, revpairwise ~ substrateDepth + sandRatio |
                  exposition | surveyYear_fac < index, type = "response"))
(emm <- emmeans(m_1, revpairwise ~ seedDensity |
                  exposition | surveyYear_fac | index, type = "response"))

### Save ###
write.csv(draws2, here("outputs", "statistics", "table_persistence_mX.csv"))
