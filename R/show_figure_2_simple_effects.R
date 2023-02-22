# Dike grassland field experiment
# Simple effects ####
# Show figure 2

# Markus Bauer
# 2022-12-12



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(patchwork)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot #######################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



graph_a + theme(legend.position = "none") +
   graph_b + theme(legend.position = "none") +
  (graph_c + theme(legend.position = c(.19, .88))) +
  plot_layout(guides = "keep") +
  plot_annotation(tag_levels = "A", tag_prefix = "", tag_suffix = "") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

### Save ###

ggsave(
  here("outputs", "figures", "figure_2_800dpi_16.5x5cm.tiff"),
  dpi = 800, width = 16.5, height = 5, units = "cm"
)
