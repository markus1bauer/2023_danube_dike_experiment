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
  mutate(
    n = fcs_target,
    surveyYear_fac = factor(surveyYear),
    targetType = factor(targetType)
  ) %>%
  select(
    id, plot, block, exposition, sandRatio, substrateDepth, targetType,
    seedDensity, surveyYear_fac, n
  )

### * Model ####
m <- m31

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
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(graph_a <- ggplot() +
    geom_quasirandom(
      aes(y = n, x = surveyYear_fac, color = targetType),
      data = sites,
      alpha = 0.5,
      dodge.width = 0.8,
      cex = .5
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      size = .3,
      color = "grey70"
    ) +
    geom_boxplot(
      aes(y = n, x = surveyYear_fac, fill = targetType),
      data = sites,
      alpha = 0.5
    ) +
    facet_grid(
      exposition ~ sandRatio,
      labeller = as_labeller(
        c(south = "South", north = "North",
        "0" = "0% Sand", "25" = "25% Sand", "50" = "50% Sand")
      )
    ) +
    scale_y_continuous(limits = c(-2.8, 1.9), breaks = seq(-100, 400, 1)) +
    scale_color_manual(labels = c("Hay meadow", "Dry grassland"),
                       values = c("#00BFC4", "#F8766D")) +
    scale_fill_manual(labels = c("Hay meadow", "Dry grassland"),
                      values = c("#00BFC4", "#F8766D")) +
    labs(
      x = "", fill = "", color = "",
      y = expression(
        Favourable ~ Conservation ~ Status ~ "(FCS)"
        )
      ) +
    theme_mb())

(graph_b <- ggeffects::ggemmeans(
      m,
      terms = c(
        "surveyYear_fac", "substrateDepth", "exposition", "sandRatio"
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
      facet ~ panel,
      labeller = as_labeller(
        c(south = "South", north = "North",
          "0" = "0% Sand", "25" = "25% Sand", "50" = "50% Sand")
      )
    ) +
    scale_y_continuous(limits = c(-1.9, 1.1), breaks = seq(-100, 400, 1)) +
    scale_color_manual(labels = c("Hay meadow", "Dry grassland"),
                       values = c("#00BFC4", "#F8766D")) +
    scale_fill_manual(labels = c("Hay meadow", "Dry grassland"),
                      values = c("#00BFC4", "#F8766D")) +
    labs(
      x = "", color = "",
      y = expression(
        Favourable ~ Conservation ~ Status ~ "(FCS)"
      )
    ) +
    theme_mb())

### Save ###
ggsave(here("outputs", "figures",
            "figure_stat_fcs_target_800dpi_16.5x14cm.tiff"),
       dpi = 800, width = 16.5, height = 14, units = "cm")
