# Model for community weighted mean of specific leaf area ####
# Markus Bauer



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(lmerTest)
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
ggplot(sites, aes(y = n, x = targetType)) +
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
ggplot(sites, aes(x = exposition, y = n)) + 
  geom_boxplot() + geom_quasirandom() + facet_wrap(~ targetType)
ggplot(sites, aes(x = sandRatio, y = n)) + 
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
m1a <- lmer(n ~ 1 + (surveyYear | block/plot), sites, REML = FALSE)
VarCorr(m1a) # convergence problems
m1b <- lmer(n ~ 1 + (surveyYear | block), sites, REML = FALSE)
VarCorr(m1b) # convergence problems
m1c <- lmer(n ~ 1 + (1 | block/plot), sites, REML = FALSE)
VarCorr(m1c)
m1d <- lmer(n ~ 1 + (1 | block/plot), sites, REML = FALSE)
VarCorr(m1d)
#fixed effects
m2 <- lmer((n) ~ (exposition + substrateDepth + sandRatio + targetType +
                    seedDensity) +
             exposition:sandRatio + exposition:targetType +
             (1 | block/plot) + (1 | surveyYear_fac), sites, REML = FALSE)
simulateResiduals(m2, plot = TRUE)
#fixed and site and year effects
m3 <- lmer(sqrt(n) ~ (exposition + substrateDepth + sandRatio + targetType +
                    seedDensity) + surveyYear_fac +
             exposition:sandRatio + exposition:targetType +
             (1 | block/plot), sites, REML = FALSE)
simulateResiduals(m3, plot = TRUE)
m4 <- lmer(sqrt(n) ~ (exposition + substrateDepth + sandRatio + targetType +
                        seedDensity + surveyYear_fac)^3 +
             (1 | block/plot), sites, REML = FALSE)
simulateResiduals(m4, plot = TRUE)


#### b comparison ------------------------------------------------------------
anova(m2, m3, m4)
rm(m1a, m1b, m1c, m1d, m4)

#### c model check -----------------------------------------------------------
simulationOutput <- simulateResiduals(m3, plot = FALSE)
plotResiduals(simulationOutput$scaledResiduals, sites$surveyYearF)
plotResiduals(simulationOutput$scaledResiduals, sites$locationAbb)
plotResiduals(simulationOutput$scaledResiduals, sites$block)
plotResiduals(simulationOutput$scaledResiduals, sites$plot)
plotResiduals(simulationOutput$scaledResiduals, sites$side)
plotResiduals(simulationOutput$scaledResiduals, sites$exposition)
plotResiduals(simulationOutput$scaledResiduals, sites$PC1)
plotResiduals(simulationOutput$scaledResiduals, sites$PC2)
plotResiduals(simulationOutput$scaledResiduals, sites$PC3)


## 3 Chosen model output #####################################################

### Model output -------------------------------------------------------------
MuMIn::r.squaredGLMM(m4)
0.411 /  0.425
VarCorr(m3)
sjPlot::plot_model(m3, type = "re", show.values = TRUE)
car::Anova(m2, type = 3)

### Effect sizes --------------------------------------------------------------
(emm <- emmeans(m3, revpairwise ~ exposition, type = "response"))
plot(emm, comparison = T)
(emm <- emmeans(m3, revpairwise ~ side, type = "response"))

### Save ###
table <- broom::tidy(car::Anova(m3, type = 2))
setwd(here("data/tables"))
write.csv(table, "table_anova_cwmAbuSla.csv")
