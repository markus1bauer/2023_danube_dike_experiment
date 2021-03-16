# Prepare spatial data ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(tidyverse)
library(sf)
library(ggmap)

### Start ###
rm(list = ls())
setwd("Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/raw")
register_google(key = "AIzaSyB5nQU_dgB_kPsQkk-_cq7pA0g1-2qka4E")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Load shp files ##########################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Base map #################################################################################################

data <- raster::getData('GADM', country = 'DEU', level = 0, download =  F)
data <- st_as_sf(data)
ger <- st_set_crs(data, 4326)

background_google_experiment <- get_map(
  location = c(lon = 12.886, lat = 48.839), 
  zoom = 15,
  scale = 1,
  maptype = "satellite",
  color = "color",
  source = "google"
  )
ggmap(background_google_experiment)

background_toner_experiment <- get_map(
  location = c(left = 12.87, bottom = 48.835, right = 12.895, top = 48.845),
  zoom = 14, 
  scale = 1,
  maptype = "toner",
  source = "stamen"
)
ggmap(background_toner_experiment)

background_terrain_experiment <- get_map(
  location= c(left = 12.87, bottom = 48.835, right = 12.895, top = 48.845),
  zoom = 14, 
  scale = 1,
  maptype = "terrain",
  source = "stamen"
)
ggmap(background_terrain_experiment)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


st_write(germany, layer = "germany.shp", driver = "ESRI Shapefile",
         dsn = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/processed/shp_files")
save(background_google_experiment, file = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/processed/shp_files/background_google_experiment.rda")
save(background_toner_experiment, file = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/processed/shp_files/background_toner_experiment.rda")
save(background_terrain_experiment, file = "Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/processed/shp_files/background_terrain_experiment.rda")
