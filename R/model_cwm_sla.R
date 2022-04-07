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
library(BayesianTools)
library(DHARMa)
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
         surveyYear_fac = factor(surveyYear)) %>%
  select(id, plot, block, surveyYear, surveyYear_fac, exposition,
         sandRatio, substrateDepth, seedDensity, targetType, n)
    



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Data exploration ########################################################

#### a Graphs -----------------------------------------------------------------
#simple effects:
ggplot(sites, aes(y = n, x = exposition)) +
  geom_boxplot() + geom_quasirandom()
ggplot(sites %>% filter(surveyYear == 2021), aes(y = n, x = targetType)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites, aes(y = n, x = seedDensity)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites, aes(y = n, x = substrateDepth)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites, aes(y = n, x = sandRatio)) +
  geom_boxplot() + geom_quasirandom() 
ggplot(sites, aes(y = n, x = block)) +
  geom_boxplot() + geom_quasirandom() 
#2way
ggplot(sites %>% filter(surveyYear == 2021), aes(x = exposition, y = n)) + 
  geom_boxplot() + geom_quasirandom() + facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021), aes(x = sandRatio, y = n)) + 
  geom_boxplot() + geom_quasirandom() + facet_wrap(~ targetType)
ggplot(sites, aes(x = substrateDepth, y = n)) + 
  geom_boxplot() + geom_quasirandom() + facet_wrap(~ targetType)
ggplot(sites, aes(x = sandRatio, y = n)) + 
  geom_boxplot() + geom_quasirandom() + facet_wrap(~ exposition)
#3way
ggplot(sites, aes(x = exposition, y = n, color = sandRatio)) + 
  geom_boxplot() + facet_wrap(~ targetType)
ggplot(sites %>% filter(surveyYear == 2021),
       aes(x = exposition, y = n, color = sandRatio)) + 
  geom_boxplot() + facet_wrap(~ targetType)
ggplot(sites, aes(x = exposition, y = n, color = substrateDepth)) + 
  geom_boxplot() + facet_wrap(~ targetType)
ggplot(sites, aes(x = substrateDepth, y = n, color = sandRatio)) + 
  geom_boxplot() + facet_wrap(~ targetType)
ggplot(sites, aes(x = factor(surveyYear), y = n, color = targetType)) + 
  geom_boxplot() + facet_wrap(~ exposition)
#4way
ggplot(sites, aes(x = exposition, y = n, color = sandRatio, fill = substrateDepth)) + 
  geom_boxplot() + facet_wrap(~ targetType)

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
#random structure
m1a <- brm(n ~ 1 + (surveyYear | block/plot), data = sites)
VarCorr(m1a) # convergence problems
m1b <- brm(n ~ 1 + (surveyYear | block), data = sites)
VarCorr(m1b) # convergence problems
m1c <- brm(n ~ 1 + (1 | block/plot), data = sites)
VarCorr(m1c)
m1d <- brm(n ~ 1 + (1 | block/plot), data = sites,)
VarCorr(m1d)
#fixed effects
m2 <- brm((n) ~ (exposition + substrateDepth + sandRatio + targetType +
                    seedDensity) +
             exposition:sandRatio + exposition:targetType +
             (1 | block/plot) + (1 | surveyYear_fac),
          data = sites, 
          family = gaussian,
          cores = parallel::detectCores(), 
          chains = 4,
          iter = 10000,
          control = list(adapt_delta = 0.99, max_treedepth = 12), 
          seed = 123, 
          prior = set_prior("normal(0,2)", class="b"))
#fixed and site and year effects
m3 <- brm(sqrt(n) ~ (exposition + substrateDepth + sandRatio + targetType +
                    seedDensity) + surveyYear_fac +
             exposition:sandRatio + exposition:targetType +
             (1 | block/plot), data = sites)
m4 <- brm(sqrt(n) ~ (exposition + substrateDepth + sandRatio + targetType +
                        seedDensity + surveyYear_fac)^3 +
             (1 | block/plot), data = sites)


#### b comparison ------------------------------------------------------------
anova(m2, m3, m4)
rm(m1a, m1b, m1c, m1d, m4)

#### c model check -----------------------------------------------------------
Samples <- coda.samples(jagsModel, variable.names = para.names, n.iter = 5000)

x = getSample(Samples)
# note - yesterday, we calcualted the predictions from the parameters
# here we observe them direct - this is the normal way to calcualte the 
# posterior predictive distribution
posteriorPredDistr = x[,5:(4+599)]
posteriorPredSim = x[,(5+599):(4+2*599)]


sim <- createDHARMa(simulatedResponse = t(posteriorPredSim), 
                   observedResponse = Owls$SiblingNegotiation, 
                   fittedPredictedResponse = apply(posteriorPredDistr, 2, median), 
                   integerResponse = TRUE)
plot(sim)
testDispersion(res)
testZeroInflation(res)
plotResiduals(Dat$Veg,res$scaledResiduals)

testSpatialAutocorrelation(res)

## 3 Chosen model output #####################################################

### Model output -------------------------------------------------------------
summary(m)
plot(m)
correlationPlot(m)
marginalPlot(m)

### Effect sizes --------------------------------------------------------------
(emm <- emmeans(m3, revpairwise ~ exposition, type = "response"))
plot(emm, comparison = T)
(emm <- emmeans(m3, revpairwise ~ side, type = "response"))

### Save ###
table <- broom::tidy(car::Anova(m3, type = 2))
setwd(here("data/tables"))
write.csv(table, "table_anova_cwmAbuSla.csv")
