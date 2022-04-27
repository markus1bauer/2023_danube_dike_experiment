# Dike grassland experiment
# Show_figure species composition ####
# Markus Bauer
# 2022-04-27


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
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
                      exposition = "f",
                      sandRatio = "d",
                      substrateDepth = "d",
                      targetType = "f",
                      seedDensity = "d"
                    )) %>%
  filter(!str_detect(id, "C")) %>%
  select(
    id, plot, block, exposition, sandRatio, substrateDepth, targetType,
    seedDensity, surveyYear, cwmAbuSla
  ) %>%
  mutate(
    n = cwmAbuSla,
    surveyYear_fac = factor(surveyYear)
  )

data <- sites

### * Model ####


### * Functions ####
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_text(angle = 0, hjust = 0.5, size = 9,
                             color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9,
                              color = "black"),
    axis.line = element_line(),
    legend.key = element_rect(fill = "white"),
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(graph_a <- ggplot() +
   geom_beeswarm(
     aes(y = n, x = surveyYear_fac, color = targetType),
     data = data,
     alpha = 0.5,
     dodge.width = 0.8
   ) +
   geom_boxplot(
     aes(y = n, x = surveyYear_fac, fill = targetType),
     data = data,
     alpha = 0.5
   ) +
   #facet_grid(~targetType) +
    #geom_hline(
     # yintercept = c(
      #  mean(sites$y),
       # mean(sites$y) + 0.5 * sd(sites$y),
        #mean(sites$y) - 0.5 * sd(sites$y)
      #),
      #linetype = c(1, 2, 2),
      #color = "grey70"
    #) +
    #scale_y_continuous(limits = c(0, .92), breaks = seq(-100, 400, .1)) +
    #scale_shape_manual(values = c("circle", "circle open")) +
    #labs(x = "", y = expression(Temporal ~ "beta" ~ diversity ~ "[" * italic("D")[sor] * "]")) +
    theme_mb())

### Save ###
ggsave(here("outputs", "figures", "sla_800dpi_8x8cm.tiff"),
       dpi = 800, width = 8, height = 8, units = "cm")


