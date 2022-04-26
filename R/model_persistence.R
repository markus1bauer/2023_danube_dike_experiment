# Dike grassland experiment
# Persistence Presence-Absence ####
# Markus Bauer
# 2022-04-25



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
  filter(!str_detect(id, "C")) %>%
  mutate(n = D,
         surveyYear_fac = factor(surveyYear),
         botanist_year = str_c(botanist, surveyYear, sep = " "),
         botanist_year = factor(botanist_year)) %>%
  select(id, plot, block, surveyYear, surveyYear_fac, botanist, botanist_year,
         exposition, sandRatio, substrateDepth, seedDensity, targetType, n)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Data exploration ########################################################

### a Graphs -----------------------------------------------------------------
#simple effects:
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = exposition)) +
  geom_boxplot() + geom_quasirandom() # 
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = targetType)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = seedDensity)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = substrateDepth)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = sandRatio)) +
  geom_boxplot() + geom_quasirandom() # 
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = block)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites, aes(y = n, x = botanist_year)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites, aes(y = n, x = surveyYear_fac)) +
  geom_boxplot() + geom_quasirandom() #
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
  facet_wrap(~ exposition) # sig
#3way
ggplot(sites, aes(x = exposition, y = n, color = sandRatio)) + 
  geom_boxplot() +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = exposition, y = n, color = sandRatio)) + 
  geom_boxplot() +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = exposition, y = n, color = substrateDepth)) + 
  geom_boxplot() +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = substrateDepth, y = n, color = sandRatio)) + 
  geom_boxplot() +
  facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = substrateDepth, y = n, color = sandRatio)) + 
  geom_boxplot() +
  facet_wrap(~ exposition)
ggplot(sites, aes(x = factor(surveyYear), y = n, color = targetType)) + 
  geom_boxplot() +
  facet_wrap(~ exposition) # Plot for figure
#4way
ggplot(sites,
       aes(x = exposition, y = n, color = sandRatio, fill = targetType)) + 
  geom_boxplot() +
  facet_wrap(~ factor(surveyYear))

##### b Outliers, zero-inflation, transformations? ----------------------------
dotchart((sites$n), groups = factor(sites$exposition),
         main = "Cleveland dotplot")
sites %>% count(block)
boxplot(sites$n)
plot(table((sites$n)), type = "h",
     xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(x = n)) + geom_density()
ggplot(sites, aes(x = log(2 - n))) + geom_density()


## 2 Model building ###########################################################

### a models -----------------------------------------------------------------

iter = 1000
chains = 4
thin = 2

m1 <- brm(n ~ targetType + exposition + sandRatio + surveyYear_fac +
            seedDensity + substrateDepth +
            (1 + exposition | block/plot) + (1 | botanist_year),
          data = sites, 
          family = gaussian("identity"),
          #prior = c(
          #  set_prior("normal(0, 1)", class = "b"),
          #  set_prior("normal(0.5, 2)", class = "b", coef = "expositionnorth"),
          #  set_prior("cauchy(0, 1)", class = "sigma")
          #),
          chains = chains,
          iter = iter,
          thin = thin,
          warmup = floor(iter / 2),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)

m2 <- brm(n ~ targetType + exposition * sandRatio * surveyYear_fac +
            seedDensity + substrateDepth +
            (1 + exposition | block/plot) + (1 | botanist_year),
          data = sites, 
          family = gaussian("identity"),
          prior = c(
            set_prior("normal(0, 1)", class = "b"),
            set_prior("normal(0.5, 2)", class = "b", coef = "expositionnorth"),
            set_prior("cauchy(0, 1)", class = "sigma")
          ),
          chains = chains,
          iter = iter,
          thin = thin,
          warmup = floor(iter / 2),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)

m2_flat <- brm(n ~ targetType + exposition * sandRatio * surveyYear_fac +
                 seedDensity + substrateDepth +
                 (1 + exposition | block/plot) + (1 | botanist_year),
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

m3 <- brm(n ~ (surveyYear_fac + exposition + sandRatio + targetType) +
            seedDensity + substrateDepth +
            (1 + exposition | block/plot) +
            (1 | botanist:surveyYear_fac) +
            (1 | surveyYear_fac:exposition:sandRatio) +
            (1 | surveyYear_fac:exposition:targetType),
          data = sites, 
          family = gaussian("identity"),
          prior = c(
            set_prior("normal(0, 1)", class = "b"),
            set_prior("normal(0.5, 2)", class = "b", coef = "expositionnorth"),
            set_prior("cauchy(0, 1)", class = "sigma")
          ),
          chains = chains,
          iter = iter,
          thin = thin,
          warmup = floor(iter / 2),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)

m4 <- lmerTest::lmer(n ~ targetType + sandRatio * exposition * surveyYear_fac +
                       seedDensity + substrateDepth +
                       (1 | block/plot),
                     data = sites,
                     REML = FALSE);lme4::isSingular(m4);simulateResiduals(m4, plot = TRUE)




### b comparison ------------------------------------------------------------

bayes_factor(m1, m2, log = TRUE)
bayes_factor(m2, m2_flat, log = TRUE)
m <- m1

### c model check -----------------------------------------------------------

createDHARMa(
  simulatedResponse = t(posterior_predict(m)),
  observedResponse = sites$n,
  fittedPredictedResponse = apply(t(posterior_epred(m)), 1, mean),
  integerResponse = TRUE
) %>%
  plot()

posterior <- m %>%
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
  )

### Samling efficency/effectiveness (Rhat and EFF) ###
(draws <- m %>%
    posterior::as_draws() %>%
    posterior::summarize_draws() %>%
    filter(str_starts(variable, "b_")))

### Trace plots ###
bayesplot::mcmc_trace(posterior,
                      facet_args = list(ncol = 2, labeller = label_parsed))

### Posterior distributions ###
(areas <- bayesplot::mcmc_areas(posterior,
                                prob = .5,
                                prob_outer = .95,
                                point_est = "mean"))

### Autocorrelation ###
rstan::stan_ac(m$fit)

### Posterior predictive check ###
pp_check(m, ndraws = 50, type = "dens_overlay")
pp_check(m, ndraws = 50, type = "dens_overlay_grouped", group = "block")
pp_check(m, ndraws = 50, type = "dens_overlay_grouped", group = "exposition")
pp_check(m, ndraws = 50, type = "dens_overlay_grouped", group = "sandRatio")
pp_check(m, ndraws = 50, type = "dens_overlay_grouped", group = "targetType")
pp_check(m, ndraws = 50, type = "dens_overlay_grouped",
         group = "surveyYear_fac")



## 3 Chosen model output #####################################################

### a Model output ------------------------------------------------------------

prior_summary(m)
bayes_R2(m, probs = c(0.05, 0.5, 0.95),
         re_formula =  ~ (1 + exposition | block/plot) + (1 | botanist_year)) 
bayes_R2(m, probs = c(0.05, 0.5, 0.95),
         re_formula = 1 ~ 1)
draws
sjPlot::plot_model(m, type = "est")
sjPlot::plot_model(m, type = "int")
sjPlot::plot_model(m, type = "pred", terms = "exposition")
sjPlot::plot_model(m, type = "pred", terms = "targetType")
sjPlot::plot_model(m, type = "pred", terms = c("exposition", "sandRatio"))

### b Effect sizes ------------------------------------------------------------

(emm <- emmeans(m, revpairwise ~ exposition, type = "response"))
(emm <- emmeans(m, revpairwise ~ targetType, type = "response"))

### Save ###
write.csv(draws, here("outputs", "statistics", "table_cwm_abu_sla.csv"))
areas
ggsave("figure_cwm_abu_sla_300dpi_14x14cm.tiff", 
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

ggsave("figure_cwm_abu_sla_flat_300dpi_14x14cm.tiff", 
       dpi = 300, width = 14, height = 14, units = "cm",
       path = here("outputs", "figures"))
