# Show map of the Danube dike experiment ####
# Markus Bauer
# Citation: Markus Bauer, Jakob Huber, Katharina Beck, Johannes Kollmann (xxxx)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(tidyverse)
library(sf)
library(ggmap)
library(ggrepel)
library(RColorBrewer)
library(patchwork)

### Start ###
rm(list = ls())
setwd("Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/processed/shp_files")

### Load data ###
ger <- st_read("germany.shp")
load("background_google_experiment.rda")
load("background_toner_experiment.rda")
load("background_terrain_experiment.rda")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
themeMB <- function(){
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid = element_line(colour = NA),
    text  = element_text(size = 10, color = "black"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(), 
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(.5, 0, 0, 0, "cm")
  )
}


## 1 Preparation ##############################################################################

### a Map of project site -----------------------------------------------------------------------
sitesGraph <- ggmap(background_google_experiment) + 
   coord_sf(crs = st_crs(4326)) +
   ggspatial::annotation_scale(location = "br", pad_y = unit(0.6, "cm"), pad_x = unit(0.7, "cm"),
                               width_hint = 0.4, height = unit(0.2, "cm")) +
   ggspatial::annotation_north_arrow(location = "br", pad_y = unit(1.1, "cm"), pad_x = unit(0.6, "cm"), 
                                     which_north = "true", style = ggspatial::north_arrow_fancy_orienteering(), height = unit(1, "cm"), width = unit(1, "cm")) +
   themeMB()

### b Germany -----------------------------------------------------------------------
gerGraph <- ggplot() +
  geom_sf(data = ger, fill = "transparent", colour = "black") +
  geom_point(aes(x = 12.886, y = 48.839), size = 1, colour = "white") +
  themeMB() +
  theme(
    plot.background = element_blank()
  )

### c Inset -----------------------------------------------------------------------
sitesGraph + inset_element(gerGraph, 
                           left = .01, 
                           bottom = .1, 
                           right = .3, 
                           top = .4, 
                           on_top = T)


# 2 Save ##############################################################################

ggsave("figure_1_map_(300dpi_8x11cm).tiff", 
       dpi = 300, width = 8, height = 11, units = "cm",
       path = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/outputs/figures")
