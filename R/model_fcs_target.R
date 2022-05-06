# Dike grassland experiment
# Community weighted mean of specific leaf area ####
# Markus Bauer
# 2022-04-11



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
#rm(list = ls())
setwd(here("data", "processed"))

### Load data ###
sites <- read_csv("data_processed_sites.csv",
                  col_names = TRUE, na = c("na", "NA"), col_types =
                    cols(
                      .default = "?",
                      block = "f",
                      exposition = "f",
                      sandRatio = "f",
                      substrateDepth = "f",
                      seedDensity = "f",
                      targetType = "f"
                    )) %>%
  filter(
    !str_detect(id, "C") & presabu == "presence" & surveyYear != "seeded"
  ) %>%
  mutate(
    n = fcs_target,
    surveyYear_fac = factor(surveyYear),
    targetType = factor(targetType),
    botanist_year = str_c(botanist, surveyYear, sep = " "),
    botanist_year = factor(botanist_year),
    surveyYear = as.double(surveyYear)
  ) %>%
  select(
    id, plot, block, exposition, sandRatio, substrateDepth, targetType,
    seedDensity, surveyYear_fac, surveyYear, botanist_year, n
  )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Data exploration ########################################################

### a Graphs -----------------------------------------------------------------
#simple effects:
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = exposition)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites %>% filter(surveyYear == 2020), aes(y = n, x = exposition)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = targetType)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = seedDensity)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = substrateDepth)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = sandRatio)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = block)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = botanist_year)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites, aes(y = n, x = surveyYear_fac)) +
  geom_boxplot() + geom_quasirandom()
#2way
ggplot(sites %>% filter(surveyYear == 2021), aes(x = exposition, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021), aes(x = sandRatio, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021), aes(x = substrateDepth, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021), aes(x = sandRatio, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ exposition)
ggplot(sites, aes(x = surveyYear_fac, y = n)) + 
  geom_boxplot() + geom_quasirandom() +
  facet_wrap(~ block)
#3way
ggplot(sites, aes(x = exposition, y = n, color = sandRatio)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = exposition, y = n, color = sandRatio)) +
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = exposition, y = n, color = substrateDepth)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = substrateDepth, y = n, color = sandRatio)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = substrateDepth, y = n, color = sandRatio)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ exposition)
ggplot(sites, aes(x = factor(surveyYear), y = n, color = targetType)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_wrap(~ exposition) # Plot for figure
#4way
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = exposition, y = n, color = sandRatio)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_grid(substrateDepth ~ targetType)
ggplot(data = sites,
       aes(x = surveyYear_fac, y = n, color = targetType)) + 
  geom_quasirandom(alpha = 0.5, dodge.width = 0.8) +
  geom_boxplot(fill = "transparent") +
  facet_grid(exposition ~ sandRatio)


##### b Outliers, zero-inflation, transformations? ----------------------------
dotchart((sites$n), groups = factor(sites$exposition),
         main = "Cleveland dotplot")
sites %>% group_by(exposition) %>% count(block)
boxplot(sites$n)
plot(table((sites$n)), type = "h",
     xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(x = n)) + geom_density()
ggplot(sites, aes(sqrt(n))) + geom_density()


## 2 Model building ###########################################################

### a models -----------------------------------------------------------------

iter = 1000#10000
chains = 3#4
thin = 2

m_simple1 <- brm(n ~ targetType + exposition + sandRatio + surveyYear_fac +
                   seedDensity + substrateDepth +
                   (1 | block/plot) + (1 | botanist_year),
                 data = sites, 
                 family = gaussian("identity"),
                 prior = c(
                   set_prior("normal(0, 3)", class = "b"),
                   set_prior("normal(0, 3)", class = "b", coef = "expositionnorth"),
                   set_prior("cauchy(0, 1)", class = "sigma")
                 ),
                 chains = chains,
                 iter = iter,
                 thin = thin,
                 warmup = floor(iter / 2),
                 save_pars = save_pars(all = TRUE),
                 cores = parallel::detectCores(),
                 seed = 123)

m_full1 <- brm(n ~ (targetType + exposition + sandRatio + surveyYear_fac +
                      seedDensity + substrateDepth)^3 +
                 targetType:sandRatio:exposition:surveyYear_fac +
                 substrateDepth:sandRatio:exposition:surveyYear_fac +
                 (1 | block/plot) + (1 | botanist_year),
               data = sites, 
               family = gaussian("identity"),
               prior = c(
                 set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 3)", class = "b", coef = "expositionnorth"),
                 set_prior("cauchy(0, 1)", class = "sigma")
               ),
               chains = chains,
               iter = iter,
               thin = thin,
               warmup = floor(iter / 2),
               save_pars = save_pars(all = TRUE),
               cores = parallel::detectCores(),
               seed = 123)

m31 <- brm(n ~ (targetType + exposition + sandRatio + surveyYear_fac)^2 +
             substrateDepth:sandRatio +
             seedDensity:exposition +
             targetType:exposition:surveyYear_fac +
             sandRatio:exposition:surveyYear_fac +
             seedDensity:exposition:surveyYear_fac +
             (1 | block/plot) + (1 | botanist_year),
           data = sites, 
           family = gaussian("identity"),
           prior = c(
             set_prior("normal(0, 3)", class = "b"),
             set_prior("normal(0, 3)", class = "b", coef = "expositionnorth"),
             set_prior("cauchy(0, 1)", class = "sigma")
           ),
           chains = chains,
           iter = iter,
           thin = thin,
           warmup = floor(iter / 2),
           save_pars = save_pars(all = TRUE),
           cores = parallel::detectCores(),
           seed = 123)

m_full1_flat <- brm(n ~ (targetType + exposition + sandRatio + surveyYear_fac +
                      seedDensity + substrateDepth)^3 +
                 targetType:sandRatio:exposition:surveyYear_fac +
                 substrateDepth:sandRatio:exposition:surveyYear_fac +
                 (1 | block/plot) + (1 | botanist_year),
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

m31_flat <- brm(n ~ (targetType + exposition + sandRatio + surveyYear_fac)^2 +
             substrateDepth:sandRatio +
             seedDensity:exposition +
             targetType:exposition:surveyYear_fac +
             sandRatio:exposition:surveyYear_fac +
             seedDensity:exposition:surveyYear_fac +
             (1 | block/plot) + (1 | botanist_year),
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

m4 <- lmerTest::lmer(n ~ (targetType + sandRatio + exposition +
                       surveyYear_fac)^2 +
                       seedDensity + substrateDepth +
                       seedDensity:targetType +
                       substrateDepth:sandRatio +
                       seedDensity:targetType:exposition +
                       substrateDepth:sandRatio:exposition +
                       sandRatio:surveyYear_fac:exposition +
                       sandRatio:exposition:targetType +
                       surveyYear_fac:exposition:targetType +
                       (1 | block/plot),
                     data = sites,
                     REML = FALSE);lme4::isSingular(m4);simulateResiduals(m4, plot = TRUE);car::Anova(m4, type = 2)




### b comparison ------------------------------------------------------------

m_1 <- m31
m_2 <- m21

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

posterior1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::subset_draws(
    variable = c(
      "b_substrateDepth30",
      "b_seedDensity8",
      "b_targetTypedry_grassland",
      "b_sandRatio25",
      "b_sandRatio50",
      "b_expositionnorth",
      "sd_block__Intercept",
      "sd_block:plot__Intercept",
      "sd_botanist_year__Intercept",
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
      "b_expositionnorth",
      "sd_block__Intercept",
      "sd_block:plot__Intercept",
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

### Samling efficency/effectiveness (Rhat and EFF) ###
draws1 <- m_1 %>%
  posterior::as_draws() %>%
  posterior::summarize_draws() %>%
  filter(str_starts(variable, "b_"))
draws2 <- m_2 %>%
  posterior::as_draws() %>%
  posterior::summarize_draws() %>%
  filter(str_starts(variable, "b_"))
range(draws1$rhat)
range(draws2$rhat)
range(draws1$ess_bulk)
range(draws2$ess_bulk)
range(draws1$ess_tail)
range(draws2$ess_tail)

### MCMC diagnostics ###
mcmc_trace(posterior1, np = hmc_diagnostics)
mcmc_trace(posterior2, np = hmc_diagnostics)
mcmc_pairs(posterior1, off_diag_args = list(size = 1.2))
mcmc_pairs(posterior2, off_diag_args = list(size = 1.2))
mcmc_scatter(m1,
             pars = c("b_surveyYear_fac2020", "b_surveyYear_fac2019"),
             np = hmc_diagnostics,
             size = 1)
mcmc_scatter(m2,
             pars = c("b_surveyYear_fac2020", "b_surveyYear_fac2019"),
             np = hmc_diagnostics,
             size = 1)
mcmc_parcoord(posterior1, np = hmc_diagnostics)
mcmc_parcoord(posterior2, np = hmc_diagnostics)

### Posterior predictive check ###
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
### Autocorrelation ###
mcmc_acf(posterior1, lags = 10)
mcmc_acf(posterior2, lags = 10)


## 3 Chosen model output #####################################################

### a Model output ------------------------------------------------------------

prior_summary(m)
bayes_R2(m, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 + exposition | block/plot) + (1 | botanist_year)) 
bayes_R2(m, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
draws
mcmc_intervals(
  posterior,
  prob = 0.5,
  prob_outer = 0.95,
  point_est = "mean"
)
sjPlot::plot_model(m, type = "est")
sjPlot::plot_model(m, type = "pred", terms = c(
  "sandRatio", "substrateDepth", "surveyYear_fac", "exposition"
))
sjPlot::plot_model(m, type = "pred", terms = c(
  "targetType", "sandRatio", "surveyYear_fac", "exposition"
  ))
sjPlot::plot_model(m, type = "pred", terms = c(
  "sandRatio", "seedDensity", "surveyYear_fac", "exposition"
))

### b Effect sizes ------------------------------------------------------------

(emm <- emmeans(m, revpairwise ~ targetType + sandRatio |
                  exposition | surveyYear_fac, type = "response"))
(emm <- emmeans(m, revpairwise ~ substrateDepth + sandRatio |
                  exposition | surveyYear_fac, type = "response"))
(emm <- emmeans(m, revpairwise ~ seedDensity |
                  exposition | surveyYear_fac, type = "response"))
### Save ###
write.csv(draws, here("outputs", "statistics", "table_fcs_target.csv"))
areas
ggsave("figure_fcs_target_300dpi_14x14cm.tiff", 
       dpi = 300, width = 14, height = 14, units = "cm",
       path = here("outputs", "figures"))

m2_flat %>%
  posterior::as_draws() %>%
  posterior::subset_draws(
    variable = c(
      "b_expositionnorth",
      "b_sandRatio25",
      "b_sandRatio50",
      "b_targetTypedry_grassland",
      "b_surveyYear_fac2019",
      "b_surveyYear_fac2020",
      "b_surveyYear_fac2021",
      "b_substrateDepth30",
      "b_seedDensity8",
      "sd_block__Intercept",
      "sigma"
    )
  ) %>%
  bayesplot::mcmc_areas(posterior,
                        prob = .5,
                        prob_outer = .95,
                        point_est = "mean")

ggsave("figure_fcs_target_flat_300dpi_14x14cm.tiff", 
       dpi = 300, width = 14, height = 14, units = "cm",
       path = here("outputs", "figures"))
