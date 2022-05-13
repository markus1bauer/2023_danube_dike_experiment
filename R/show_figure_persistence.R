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
#rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
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

### * Model ####
m <- m3

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


## 1 Boxplots #################################################################

(graph_a <- ggplot() +
    geom_quasirandom(
      aes(y = n, x = surveyYear_fac, color = sandRatio),
      data = sites,
      method = "pseudorandom",
      alpha = 0.2,
      dodge.width = 0.8,
      cex = .5
    ) +
    geom_boxplot(
      aes(y = n, x = surveyYear_fac, fill = sandRatio),
      data = sites,
      alpha = 0.5
    ) +
    facet_grid(
      index ~ exposition,
      labeller = as_labeller(
        c(south = "South", north = "North",
          B = "Losses", C = "Gains", D = "Total")
        )
      ) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(-100, 400, .25)) +
    scale_color_manual(values = c("#990000", "#CC6600", "#FFFF00")) +
    scale_fill_manual(values = c("#990000", "#CC6600", "#FFFF00")) +
    labs(
      x = "", fill = "Sand ratio [%]", color = "Sand ratio [%]",
      y = expression(Persistence)
      ) +
    theme_mb())

### Save ###
ggsave(here("outputs", "figures", "figure_box_persistence_800dpi_25x9cm.tiff"),
       dpi = 800, width = 25, height = 9, units = "cm")


## 2 Marginal effects ##########################################################

(graph_b <- ggeffects::ggemmeans(
  m,
  terms = c(
    "surveyYear_fac", "targetType", "exposition", "index"
  ),
  ci.lvl = 0.95,
  type = "fixed",
  typical = "mean",
  back.transform = TRUE,
  ppd = TRUE
) %>%
  mutate(facet = fct_relevel(facet, "north", "south")) %>%
  ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    size = .3,
    color = "grey70"
  ) +
  geom_point(
    aes(
      x = x, y = predicted, color = group
    ),
    position = position_dodge(.5)
  ) +
  geom_pointrange(
    aes(
      x = x, y = predicted, color = group, ymin = conf.low, ymax = conf.high
    ),
    position = position_dodge(.5)
  ) +
  facet_grid(
    index ~ exposition,
    labeller = as_labeller(
      c(south = "South", north = "North",
        B = "Losses", C = "Gains", D = "Total")
    )
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(-100, 400, .25)) +
  scale_color_manual(labels = c("Hay meadow", "Dry grassland"),
                     values = c("#00BFC4", "#F8766D")) +
  scale_fill_manual(labels = c("Hay meadow", "Dry grassland"),
                    values = c("#00BFC4", "#F8766D")) +
  labs(
    x = "", color = "",
    y = expression(Persistence)
  ) +
  theme_mb())

### Save ###
ggsave(here("outputs", "figures",
            "figure_stat_persistence_800dpi_16.5x14cm.tiff"),
       dpi = 800, width = 16.5, height = 14, units = "cm")
