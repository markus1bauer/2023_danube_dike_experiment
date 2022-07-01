# Dike grassland experiment
# Show_figure species composition ####
# Markus Bauer
# 2022-05-025


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
                      id = "c",
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
    seedDensity, surveyYear, accumulatedCov
  ) %>%
  mutate(
    n = accumulatedCov,
    surveyYear_fac = factor(surveyYear),
    targetType = factor(targetType)
  )
table(sites$block, sites$surveyYear)
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
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Adapted seedmix #################################################################

(graph_a <- ggplot() +
   geom_quasirandom(
     aes(y = n, x = sandRatio, color = targetType),
     data = sites,
     alpha = 0.5,
     dodge.width = 0.8,
     cex = .5
   ) +
   geom_hline(
     yintercept = c(75, 100),
     linetype = c("dashed", "solid", "dashed", "solid", "dashed", "solid",
                  "dashed", "solid", "dashed", "solid", "dashed", "solid",
                  "dashed", "solid", "dashed", "solid"),
     size = .3,
     color = "grey70"
   ) +
   geom_boxplot(
     aes(y = n, x = sandRatio, fill = targetType),
     data = sites,
     alpha = 0.5
   ) +
   facet_grid(
     exposition ~ surveyYear_fac,
     labeller = as_labeller(
       c(south = "South", north = "North",
         "2018" = "2018", "2019" = "2019", "2020" = "2020", "2021" = "2021")
     )
   ) +
   scale_y_continuous(limits = c(0, 125), breaks = seq(-100, 400, 25)) +
   scale_color_manual(labels = c("Hay meadow", "Dry grassland"),
                      values = c("#00BFC4", "#F8766D")) +
   scale_fill_manual(labels = c("Hay meadow", "Dry grassland"),
                     values = c("#00BFC4", "#F8766D")) +
   labs(
     x = "Sand ratio [%]", fill = "", color = "",
     y = expression(Accumulated ~ cover ~ "[%]")
   ) +
   theme_mb())

### Save ###
ggsave(here("outputs", "figures", "figure_box_accumulatedCov_800dpi_27x9cm.tiff"),
       dpi = 800, width = 27, height = 9, units = "cm")
ggsave(here("outputs", "figures", "figure_box_accumulatedCov_800dpi_16.5x14cm.tiff"),
       dpi = 800, width = 16.5, height = 14, units = "cm")


## 2 Seed density #################################################################

(graph_b <- ggplot() +
    geom_quasirandom(
      aes(y = n, x = seedDensity, color = targetType),
      data = sites,
      alpha = 0.5,
      dodge.width = 0.8,
      cex = .5
    ) +
    geom_hline(
      yintercept = c(75, 100),
      linetype = c("dashed", "solid", "dashed", "solid", "dashed", "solid",
                   "dashed", "solid", "dashed", "solid", "dashed", "solid",
                   "dashed", "solid", "dashed", "solid"),
      size = .3,
      color = "grey70"
    ) +
    geom_boxplot(
      aes(y = n, x = seedDensity, fill = targetType),
      data = sites,
      alpha = 0.5
    ) +
    facet_grid(
      exposition ~ surveyYear_fac,
      labeller = as_labeller(
        c(south = "South", north = "North",
          "2018" = "2018", "2019" = "2019", "2020" = "2020", "2021" = "2021")
      )
    ) +
    scale_y_continuous(limits = c(0, 125), breaks = seq(-100, 400, 25)) +
    scale_color_manual(labels = c("Hay meadow", "Dry grassland"),
                       values = c("#00BFC4", "#F8766D")) +
    scale_fill_manual(labels = c("Hay meadow", "Dry grassland"),
                      values = c("#00BFC4", "#F8766D")) +
    labs(
      x = "Seed density [g/m²]", fill = "", color = "",
      y = expression(Accumulated ~ cover ~ "[%]")
    ) +
    theme_mb())

### Save ###
ggsave(here("outputs", "figures", "figure_box_accumulatedCov_seedDensity_800dpi_27x9cm.tiff"),
       dpi = 800, width = 27, height = 9, units = "cm")
ggsave(here("outputs", "figures", "figure_box_accumulatedCov_seedDensity_800dpi_16.5x14cm.tiff"),
       dpi = 800, width = 16.5, height = 14, units = "cm")
