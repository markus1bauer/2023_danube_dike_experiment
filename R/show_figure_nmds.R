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
                      exposition = col_factor(levels = c("north", "south")),
                      sandRatio = "f",
                      substrateDepth = "f",
                      targetType = "c",
                      seedDensity = "f"
                    )) %>%
  select(
    id, plot, block, exposition, sandRatio, substrateDepth, targetType,
    seedDensity, surveyYear, NMDS1, NMDS2
  ) %>%
  filter(targetType != "non_ffh") %>%
  mutate(
    targetType = if_else(targetType == "mixed", "hay_meadow", targetType),
    surveyYear_fac = as.character(surveyYear),
    surveyYear_fac = if_else(block == "C", "reference", surveyYear_fac)
  )

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
    legend.position = "bottom",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


(graph_b <- ggplot() +
   geom_point(
     aes(y = NMDS2, x = NMDS1, color = surveyYear_fac, shape = surveyYear_fac),
     data = sites,
     alpha = 1,
     cex = 2
   ) +
   facet_grid(
     exposition ~ targetType,
     labeller = as_labeller(
       c(south = "South", north = "North",
         "dry_grassland" = "Dry grassland", "hay_meadow" = "Hay meadow")
     )
   ) +
   coord_fixed() +
   scale_color_manual(values = c(
     "yellow1", "orange2", "firebrick3", "deeppink4", "royalblue", "black"
   )) +
   scale_shape_manual(values = c(
     "circle", "circle", "circle", "circle", "square", "square open"
   )) +
   labs(
     x = "NMDS1", fill = "", color = "", shape = "",
     y = "NMDS2"
   ) +
   theme_mb())


### Save ###
ggsave(here("outputs", "figures", "figure_nmds_800dpi_16.5x16cm.tiff"),
       dpi = 800, width = 16.5, height = 16, units = "cm")
