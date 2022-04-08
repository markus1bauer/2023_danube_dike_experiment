# Model for community weighted mean of specific leaf area ####
# Markus Bauer



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(brms)
library(DHARMa)
library(MCMCvis)
library(emmeans)

### Start ###
rm(list = ls())
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
                    )
) %>%
  filter(!str_detect(id, "C")) %>%
  mutate(n = cwmAbuSla,
         surveyYear_fac = factor(surveyYear),
         botanist_year = str_c(botanist, surveyYear, sep = " "),
         botanist_year = factor(botanist_year)) %>%
  select(id, plot, block, surveyYear, surveyYear_fac, botanist, botanist_year,
         exposition, sandRatio, substrateDepth, seedDensity, targetType, n)




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Data exploration ########################################################

#### a Graphs -----------------------------------------------------------------
#simple effects:
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = exposition)) +
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
  facet_wrap(~ exposition)
#4way
ggplot(sites %>% filter(surveyYear == 2021), aes(x = exposition, y = n,
                  color = sandRatio, fill = substrateDepth)) + 
  geom_boxplot() +
  facet_wrap(~ targetType)

##### b Outliers, zero-inflation, transformations? ----------------------------
dotchart((sites$n), groups = factor(sites$exposition),
         main = "Cleveland dotplot")
sites %>% count(block)
boxplot(sites$n)
plot(table((sites$n)), type = "h",
     xlab = "Observed values", ylab = "Frequency")
ggplot(sites, aes(n)) + geom_density()
ggplot(sites, aes(log(n))) + geom_density()


## 2 Model building ###########################################################

#### a models -----------------------------------------------------------------
iter = 10000
chains = 4

m1 <- brm(n ~ 1 + 
           (1 | block/plot) + (1 | botanist_year),
         data = sites, 
         warmup = 500,
         chains = chains,
         iter = iter,
         prior = set_prior("cauchy(0,1)", class = "sigma"),
         control = list(adapt_delta = 0.99),
         save_pars = save_pars(all = TRUE),
         cores = parallel::detectCores(),
         seed = 123)
summary(m1)
m2 <- brm(n ~ (exposition + substrateDepth + sandRatio + targetType +
                 seedDensity) + 
            (1 | block/plot) + (1 | botanist_year),
          data = sites, 
          warmup = 500,
          chains = chains,
          iter = iter,
          prior = set_prior("cauchy(0,1)", class = "sigma"),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)
summary(m2) # check convergence by ESS and PSRF
m3 <- brm(n ~ (exposition + substrateDepth + sandRatio + targetType +
                 seedDensity) + 
            (1 + targetType | block/plot) + (1 + targetType | botanist_year),
          data = sites, 
          warmup = 500,
          chains = chains,
          iter = iter,
          prior = set_prior("cauchy(0,1)", class = "sigma"),
          control = list(adapt_delta = 0.99),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)
summary(m3) # check convergence by ESS and PSRF
m4 <- brm(n ~ (exposition + substrateDepth + sandRatio + seedDensity) *
                 targetType + 
            (1 + targetType | block/plot) + (1 + targetType | botanist_year),
          data = sites, 
          warmup = 500,
          chains = chains,
          iter = iter,
          prior = set_prior("cauchy(0,1)", class = "sigma"),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)
summary(m4) # check convergence by ESS and PSRF
m5 <- brm(n ~ exposition * substrateDepth * sandRatio + seedDensity +
            targetType + 
            (1 + targetType | block/plot) + (1 + targetType | botanist_year),
          data = sites, 
          warmup = 500,
          chains = chains,
          iter = iter,
          prior = set_prior("cauchy(0,1)", class = "sigma"),
          control = list(adapt_delta = 0.99),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)
summary(m5) # check convergence by ESS and PSRF
m6 <- brm(n ~ exposition * substrateDepth * sandRatio + seedDensity +
            targetType + botanist_year +
            (1 + targetType | block/plot),
          data = sites, 
          warmup = 500,
          chains = chains,
          iter = iter,
          prior = c(
            set_prior("cauchy(0,1)", class = "sigma"),
            set_prior("normal(0.1,2)", class = "b", coef = "targetType")
            ),
          control = list(adapt_delta = 1.1),
          save_pars = save_pars(all = TRUE),
          cores = parallel::detectCores(),
          seed = 123)
summary(m6) # check convergence by ESS and PSRF
prior_summary(m6)
plot(hypothesis(m6, "targetTypedry_grassland > 0"))
# A general guideline suggests that values less than 1.05 are good, between 1.05 and 1.10 are ok, and above 1.10 have not converged well.
coda::effectiveSize(m3)
plot(m)
gelman.diag(Samples)
correlationPlot(Samples)
MCMCtrace(m6, 
          #params = c('beta[1]', 'beta[2]', 'beta[3]'), 
          ISB = FALSE, 
          exact = TRUE,
          Rhat = TRUE,
          n.eff = TRUE,
          pdf = TRUE)
MCMCplot(m6)#, 
         #params = "beta", 
         #ci = c(50, 80), 
         #HPD = FALSE, 
         #ref_ovl = TRUE)

#### b comparison ------------------------------------------------------------
brms::bayes_factor(m4, m5, log = FALSE)
rm(m3)

#### c model check -----------------------------------------------------------
pp_check(m2, nsamples = 50)
pp_check(m2, nsamples = 50, type = "scatter_avg_grouped",
         group = "block")
pp_check(m2, type = "loo_intervals")
pp_check(m2, type = "loo_pit")

model_check <- createDHARMa(
  simulatedResponse = t(posterior_predict(m6)),
  observedResponse = sites$n,
  fittedPredictedResponse = apply(t(posterior_epred(m6)), 1, mean),
  integerResponse = TRUE
  ) %>%
  plot()


## 3 Chosen model output #####################################################

### Model output -------------------------------------------------------------
bayes_R2(m2, re_formula =  ~ (1 | block/plot), probs = c(0.05, 0.5, 0.95) ) 
bayes_R2(m2, re_formula = 1 ~ 1, probs = c(0.05, 0.5, 0.95) ) 
# conditional Bayes R2 (Gelman et al. 2019)
summary(m2)
marginalPlot(m2$fit)
MCMCvis::MCMCplot(m2)

### Effect sizes --------------------------------------------------------------
(emm <- emmeans(m2, revpairwise ~ exposition, type = "response"))

### Save ###
table <- tid(m2$fit, conf.int = TRUE, conf.method = "HPDinterval", rhat = TRUE, ess = TRUE)
write.csv(table, here("outputs", "statistics", "table_cwmAbuSla.csv"))
          