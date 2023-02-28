# Dike grassland field experiment
# Show map of the Danube dike experiment ####
# Show figure 1

# Markus Bauer
# 2023-02-09



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(sf)
library(raster)
library(ggmap)
library(tmap)
library(RColorBrewer)
library(ggrepel)
library(patchwork)
library(grid)

### Start ###
rm(list = ls())
setwd(here("data", "processed", "spatial"))

### Load data ###
germany <- st_read("germany.shp")
danube <- st_read("danube.shp")
dikes <- st_read("dikes.shp")
blocks <- st_read("blocks.shp")
blocks_table <- read_csv(
  "blocks_table.csv",
  col_names = TRUE,
  col_types = cols(block = col_factor())
)
load("bg_stamen_terrain.rda")

### Functions ###
theme_mb <- function() {
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



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Map with ggmap ###########################################################


### a Map of project site -----------------------------------------------------

(sites_graph <- ggmap(bg_stamen_terrain) +
   geom_point(
     aes(x = x, y = y), data = blocks_table, size = 2, color = "black"
     ) +
   geom_text_repel(
     aes(label = block, x = x, y = y),
     data = blocks_table,
     min.segment.length = 1,
     color = "black",
     direction = "y",
     nudge_y = .001
   ) +
   coord_sf(crs = st_crs(4326)) +
   ggspatial::annotation_scale(
     location = "br",
     pad_y = unit(0.6, "cm"),
     pad_x = unit(0.7, "cm"),
     width_hint = 0.4,
     height = unit(0.15, "cm"),
     text_col = "black"
   ) +
   ggspatial::annotation_north_arrow(
     location = "br",
     pad_y = unit(1.1, "cm"),
     pad_x = unit(0.6, "cm"),
     which_north = "true",
     style = ggspatial::north_arrow_fancy_orienteering(
       text_col = "black"
     ),
     height = unit(.8, "cm"),
     width = unit(.8, "cm")
   ) +
   theme_mb())


### b Germany -----------------------------------------------------------------

germany_graph <- ggplot() +
  geom_sf(data = germany, fill = "transparent", colour = "black") +
  geom_point(aes(x = 12.886, y = 48.839), size = 1, colour = "black") +
  theme_mb() +
  theme(
    plot.background = element_blank()
  )


### c Inset --------------------------------------------------------------------

sites_graph + inset_element(
  germany_graph,
  left = .01,
  bottom = .1,
  right = .3,
  top = .4,
  on_top = TRUE
)



# 2 Map with tmap ##############################################################


### a Map of project site ------------------------------------------------------

location_experiment <- st_point(c(12.886, 48.839)) %>%
  st_sfc(crs = 4326)
tmap_mode("plot")
bbox <- st_bbox(c(ymin = 48.835, ymax = 48.845, xmin = 12.87, xmax = 12.895))
(tmap <- tm_shape(dikes, bbox = bbox) +
    tm_lines(lwd = 3) +
    tm_shape(danube) +
    tm_fill(col = "grey") +
    tm_text("river") +
    tm_shape(blocks) +
    tm_fill("red", ) +
    tm_text("block", size = 1, col = "black", ymod = 0.7) +
    tm_compass(position = c("right", "bottom"), size = 2) +
    tm_scale_bar(
      position = c("right", "bottom", with = 0.4), breaks = c(.25, .5)
      ) +
    tm_layout(frame = FALSE))
(tmap_germany <- tm_shape(germany) +
    tm_borders(col = "black") +
    tm_shape(location_experiment) +
    tm_dots(size = .075) +
    tm_layout(frame = FALSE, bg.color = "transparent"))
(tmap_germany_white <- tm_shape(germany) +
    tm_borders(col = "white") +
    tm_shape(location_experiment) +
    tm_dots(size = .075, col = "white") +
    tm_layout(frame = FALSE, bg.color = "transparent"))


### b Save ---------------------------------------------------------------------

tmap_save(
  tmap,
  insets_tm = tmap_germany,
  insets_vp = viewport(
    x = unit(1.1, "cm"), y = unit(1.1, "cm"),
    width = unit(1, "cm"), height = unit(1.5, "cm")
  ),
  dpi = 300, width = 8, height = 5, units = "cm",
  filename = paste0(
    here("outputs", "figures"), "/", "figure_1_map_tmap_300dpi_8x5cm.tiff"
    )
)

tmap_save(
  tmap_germany_white,
  dpi = 300, width = 2, height = 3.5, units = "cm", bg = "transparent",
  device = tiff,
  type = "cairo",
  filename = paste0(
    here("outputs", "figures"), "/", "figure_1_germany_300dpi_2x3.5cm.tiff"
  )
)
