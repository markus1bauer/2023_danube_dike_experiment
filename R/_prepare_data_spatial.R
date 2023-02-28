# Grassland experiment on dikes
# Prepare spatial data ####

# Markus Bauer
# 2023-02-28



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(raster)
library(ggmap)
library(sf)
library(mapview)
library(mapedit)

### Start ###
rm(list = ls())



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Load shp files #############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Base maps #################################################################


germany <- raster::getData(
  'GADM', country = 'DEU', level = 0, download =  TRUE
) %>% 
  st_as_sf() %>%
  st_set_crs(4326)

bg_stamen_terrain <- get_map(
  source = "stamen",
  maptype = "terrain",
  location= c(left = 12.87, bottom = 48.835, right = 12.895, top = 48.845),
  zoom = 14,
  scale = 1
)
ggmap(bg_stamen_terrain)




## 2 Load and transform shp files ##############################################


dikes <- st_read(here("data", "raw", "spatial", "Deich.shp")) %>%
  st_transform(crs = 4326) %>%
  st_crop(ymin = 48.835, ymax = 48.845, xmin = 12.87, xmax = 12.895)

blocks <- st_read(here("data", "raw", "spatial", "blocks.kml")) %>%
  group_by(Name) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_convex_hull() %>%
  rename(block = Name) %>%
  mutate(
    block = as_factor(block),
    block = fct_recode(
      block, "6" = "1", "5" = "2", "4" = "3", "3" = "4", "2" = "5", "1" = "6"
      )
    )

blocks_table <- blocks %>%
  st_centroid()
blocks_table <- bind_cols(
  as_tibble(as_factor(st_drop_geometry(blocks_table))),
  as_tibble(st_coordinates(blocks_table))
  ) %>%
  rename(x = X, y = Y)



## 3 Digitize shp files ########################################################


### Digitized only once: ###
#data <- mapview() %>% editMap()
#mapview(data$finished)
#danube <- data$finished %>% st_as_sf()
#plot(st_geometry(danube))


  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



write_csv(
  blocks_table,
  file = here("data", "processed", "spatial", "blocks_table.csv")
)
st_write(
  germany,
  layer = "germany.shp",
  driver = "ESRI Shapefile",
  delete_layer = TRUE,
  dsn = here("data", "processed", "spatial")
)
### Digitized shp saved only once: ###
#st_write(
#  danube,
#  layer = "danube.shp",
#  driver = "ESRI Shapefile",
#  delete_layer = TRUE,
#  dsn = here("data", "processed", "spatial")
#)
st_write(
  dikes,
  layer = "dikes.shp",
  driver = "ESRI Shapefile",
  delete_layer = TRUE,
  dsn = here("data", "processed", "spatial")
)
st_write(
  blocks,
  layer = "blocks.shp",
  driver = "ESRI Shapefile",
  delete_layer = TRUE,
  dsn = here("data", "processed", "spatial")
)
save(
  bg_stamen_terrain,
  file = paste0(
    here("data", "processed", "spatial"), "/", "bg_stamen_terrain.rda"
    )
)
