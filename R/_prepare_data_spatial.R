# Grassland experiment on dikes
# Prepare spatial data ####

# Markus Bauer
# 2023-02-09



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
setwd(here("data", "raw", "spatial"))
register_google(key = "AIzaSyB5nQU_dgB_kPsQkk-_cq7pA0g1-2qka4E")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Load shp files #############################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Base maps #################################################################


germany <- raster::getData(
  'GADM', country = 'DEU', level = 0, download =  TRUE
) %>% 
  st_as_sf() %>%
  st_set_crs(4326)

bg_google_satellite <- get_map(
  location = c(lon = 12.886, lat = 48.839), 
  zoom = 15,
  scale = 1,
  maptype = "satellite",
  color = "color",
  source = "google"
  )
ggmap(bg_google_satellite)

bg_stamen_terrain <- get_map(
  location= c(left = 12.87, bottom = 48.835, right = 12.895, top = 48.845),
  zoom = 14, 
  scale = 1,
  maptype = "terrain-background",
  source = "stamen"
)
ggmap(bg_stamen_terrain)

mgmap <- as.matrix(bg_google_satellite)
vgmap <- as.vector(mgmap)
vgmaprgb <- col2rgb(vgmap)
gmapr <- matrix(vgmaprgb[1, ], ncol = ncol(mgmap), nrow = nrow(mgmap))
gmapg <- matrix(vgmaprgb[2, ], ncol = ncol(mgmap), nrow = nrow(mgmap))
gmapb <- matrix(vgmaprgb[3, ], ncol = ncol(mgmap), nrow = nrow(mgmap))
raster <- brick(raster(gmapr), raster(gmapg), raster(gmapb))
rm(gmapr, gmapg, gmapb)
rgmaprgbGM <- raster
projection(rgmaprgbGM) <- CRS("+init=epsg:3857")
rprobextSpDF <- as(
  extent(unlist(attr(bg_google_satellite, which = "bb"))[c(2, 4, 1, 3)]),
  "SpatialPolygons"
  )
projection(rprobextSpDF) <- CRS("+init=epsg:4326")
rprobextGM <- spTransform(rprobextSpDF, CRS("+init=epsg:3857"))
extent(rgmaprgbGM) <- c(rprobextGM@bbox[1, ], rprobextGM@bbox[2, ])
bg_google_satellite2 <- rgmaprgbGM
rm(mgmap, raster, rgmaprgbGM, rprobextGM, rprobextSpDF, vgmaprgb, vgmap)



## 2 Transform shp files #######################################################


dikes <- st_read("Deich.shp") %>%
  st_transform(crs = 4326) %>%
  st_crop(ymin = 48.835, ymax = 48.845, xmin = 12.87, xmax = 12.895)
### change existing shp-file ###
#data <- mapview(st_geometry(blocks)) %>%
  #editMap(targetLayerId = "blocks")
#mapview(data$drawn)

blocks <- st_read("blocks.kml") %>%
  group_by(Name) %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_convex_hull %>%
  rename(block = Name) %>%
  mutate(block = as_factor(block)) %>%
  mutate(block = fct_recode(
    block,
    "6" = "1", "5" = "2", "4" = "3", "3" = "4", "2" = "5", "1" = "6"
    ))

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
#danube <- data$finished %>%
  #st_as_sf() %>%
  #rename(river = X_leaflet_id) %>%
  #mutate(river = str_replace(as.character(river), ".", "Danube")) %>%
  #mutate(river = str_extract(river, "Danube")) %>%
  #st_crop(ymin = 48.835, ymax = 48.845, xmin = 12.87, xmax = 12.895)
#plot(st_geometry(danube))


  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



write_csv2(blocks_table, file = "blocks2.csv")
writeRaster(
  bg_google_satellite2,
  file = paste0(
    here("data", "processed", "spatial"), "/", "bg_google_satellite2.grd"
  ),
  overwrite = TRUE, 
  datatype = "raster"
)
st_write(
  germany,
  layer = "germany.shp",
  driver = "ESRI Shapefile",
  delete_layer = TRUE,
  dsn = here("data", "processed", "spatial")
)
### Digitized only once: ###
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
  bg_google_satellite,
  file = paste0(
    here("data", "processed", "spatial"), "/", "bg_google_satellite.rda"
    )
)
save(
  bg_stamen_terrain,
  file = paste0(
    here("data", "processed", "spatial"), "/", "bg_stamen_terrain.rda"
    )
)
