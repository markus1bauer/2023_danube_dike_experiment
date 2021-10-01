# Opinion paper ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(lmerTest)
library(DHARMa)
library(emmeans)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv("data_processed_sites.csv", col_names = T, na = c("na", "NA"), col_types = 
                    cols(
                      .default = "?",
                      block = "f",
                      surveyYear = "f"
                    )) %>%
  select(id, block, surveyYear, speciesRichness, targetRichness, vegetationCov, bioMass, accumulatedCov, graminoidCov)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Data exploration #####################################################################################

#### a Graphs ---------------------------------------------------------------------------------------------
#simple effects:
ggplot(sites, aes(y = targetRichness, x = log(bioMass + 1), color = block)) + 
  geom_point() + 
  geom_smooth(method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.1)))
ggsave("opinion_paper_exp_biomass.tiff", 
       dpi = 300, width = 8, height = 7, units = "cm",
       path = here("outputs/figures"))
ggplot(sites, aes(y = log(targetRichness), x = vegetationCov, color = surveyYear)) + 
  geom_point() + 
  geom_smooth(method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.1)))
ggsave("opinion_paper_exp_vegetationCov.tiff", 
       dpi = 300, width = 8, height = 7, units = "cm",
       path = here("outputs/figures"))
ggplot(sites, aes(y = log(targetRichness), x = accumulatedCov, color = surveyYear)) + 
  geom_point() + 
  geom_smooth(method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.1)))
ggplot(sites, aes(y = targetRichness, x = log(graminoidCov + 1), color = surveyYear)) + 
  geom_point() + 
  geom_smooth(method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.1)))
ggsave("opinion_paper_exp_graminoidCov.tiff", 
       dpi = 300, width = 8, height = 7, units = "cm",
       path = here("outputs/figures"))


## 2 Model building ################################################################################

#### a models ----------------------------------------------------------------------------------------


#### b comparison -----------------------------------------------------------------------------------------
anova(m2, m3, m4)
rm(m1a, m1b, m1c, m1d, m4)

#### c model check -----------------------------------------------------------------------------------------
simulationOutput <- simulateResiduals(m3, plot = F)
plotResiduals(simulationOutput$scaledResiduals, sites$surveyYearF)
plotResiduals(simulationOutput$scaledResiduals, sites$locationAbb)
plotResiduals(simulationOutput$scaledResiduals, sites$block)
plotResiduals(simulationOutput$scaledResiduals, sites$plot)
plotResiduals(simulationOutput$scaledResiduals, sites$side)
plotResiduals(simulationOutput$scaledResiduals, sites$exposition)
plotResiduals(simulationOutput$scaledResiduals, sites$PC1)
plotResiduals(simulationOutput$scaledResiduals, sites$PC2)
plotResiduals(simulationOutput$scaledResiduals, sites$PC3)


## 3 Chosen model output ################################################################################

### Model output ---------------------------------------------------------------------------------------------
MuMIn::r.squaredGLMM(m5)
0.411 /  0.425
MuMIn::r.squaredGLMM(m2) # to see which R2m is due to 'location' and 'surveyYearF'
0.015 / 0.241
rm(m2)
VarCorr(m3)
sjPlot::plot_model(m3, type = "re", show.values = T)
car::Anova(m5, type = 2)

### Effect sizes -----------------------------------------------------------------------------------------
(emm <- emmeans(m3, revpairwise ~ exposition, type = "response"))
plot(emm, comparison = T)
(emm <- emmeans(m3, revpairwise ~ side, type = "response"))

### Save ###
table <- tidy(car::Anova(m2, type = 2))
setwd(here("data/tables"))
write.csv2(table, "table_anova_cwmAbuHeight.csv")
