# Show map of the Danube dike experiment ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
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
setwd(here("data","processed", "spatial"))

### Load data ###
germany <- st_read("germany.shp")
danube <- st_read("danube.shp")
dikes <- st_read("dikes.shp")
blocks <- st_read("blocks.shp")
blocks2 <- read_csv2("blocks2.csv", col_names = T, col_types = 
                      cols(block = col_factor())
)
bg_google_satellite2 <- raster("bg_google_satellite2.tif")
load("bg_google_satellite.rda")
load("bg_stamen_terrain.rda")



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


## 1 Map with ggmap ##############################################################################

### a Map of project site -----------------------------------------------------------------------
(sitesGraph <- ggmap(bg_google_satellite) + 
   geom_point(data = blocks2, aes(x = x, y = y), size = 2, color = "white") +
   geom_text_repel(data = blocks2, aes(label = block, x = x, y = y),
                   min.segment.length = 1, 
                   color = "white", 
                   direction = "y", 
                   nudge_y = .001) +
   coord_sf(crs = st_crs(4326)) +
   ggspatial::annotation_scale(location = "br", 
                               pad_y = unit(0.6, "cm"), 
                               pad_x = unit(0.7, "cm"),
                               width_hint = 0.4, 
                               height = unit(0.15, "cm"),
                               text_col = "white") +
   ggspatial::annotation_north_arrow(location = "br", 
                                     pad_y = unit(1.1, "cm"), 
                                     pad_x = unit(0.6, "cm"), 
                                     which_north = "true", 
                                     style = ggspatial::north_arrow_fancy_orienteering(
                                       text_col = "white"
                                       ), 
                                     height = unit(.8, "cm"), 
                                     width = unit(.8, "cm")) +
   themeMB())

### b Germany -----------------------------------------------------------------------
gerGraph <- ggplot() +
  geom_sf(data = germany, fill = "transparent", colour = "white") +
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

### d Save -----------------------------------------------------------------------
ggsave("figure_1_map_300dpi_8x11cm.tiff", 
       dpi = 300, width = 8, height = 11, units = "cm",
       path = here("outputs", "figures")
)


# 2 Map with tmap ##############################################################################

### a Map of project site -----------------------------------------------------------------------
tmap_mode("plot")
bbox <- st_bbox(c(ymin = 48.835, ymax = 48.845, xmin = 12.87, xmax = 12.895))
#tm_shape(bg_google_satellite2, bbox = bbox) +
  #tm_raster()
(tmap <- tm_shape(dikes, bbox = bbox) +
  tm_lines() +
  tm_shape(danube) +
  tm_fill(col = "grey") +
  tm_text("river") +
  tm_shape(blocks) +
  tm_fill("red") +
  tm_text("block", size = 1, col = "black", ymod = 0.7) +
  tm_compass(position = c("right", "bottom"), size = 2) +
  tm_scale_bar(position = c("right", "bottom", with = 0.4)) +
  tm_layout(frame = F))
tmap_ger <- tm_shape(germany) +
  tm_borders(col = "black") +
  tm_layout(frame = F)

### b Save -----------------------------------------------------------------------
tmap_save(tmap,
          insets_tm = tmap_ger,
          insets_vp = viewport(x = unit(3.1, "cm"),
                               y = unit(3.3, "cm"),
                               width = unit(3, "cm"),
                               height = unit(4, "cm")), 
          filename = paste0(here("outputs/figures"), "/", "figure_1_map_tmap_(300dpi_8x11cm).tiff"),
          dpi = 300)
