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



graph_a +
  annotate("text", x = -.29, y = 5.71, label = "A", fontface = 2) +
  theme(legend.position = "none") +
  graph_b +
  annotate("text", x = -.07, y = 5.71, label = "B", fontface = 2) +
  theme(legend.position = "none") +
  (graph_c +
     annotate("text", x = -.5, y = 5.71, label = "C", fontface = 2) +
     theme(legend.position = c(.19, .88))) +
  plot_layout(guides = "keep")

### Save ###

ggsave(
  here("outputs", "figures", "figure_3_300dpi_16.5x5cm.tiff"),
  dpi = 300, width = 16.5, height = 5, units = "cm"
)
