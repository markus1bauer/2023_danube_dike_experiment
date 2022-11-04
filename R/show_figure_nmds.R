# Dike grassland experiment
# Show_figure species composition ####
# Markus Bauer
# 2022-10-28



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data", "processed"))

### Functions ###
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

vegan_cov_ellipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#### * Load sites data ####

sites_experiment <- read_csv("data_processed_sites.csv",
                             col_names = TRUE, na = c("na", "NA", ""),
                             col_types = cols(.default = "?")) %>%
  mutate(
    reference = survey_year,
    target_type = if_else(
      target_type == "dry_grassland", "Dry grassland", if_else(
        target_type == "hay_meadow", "Hay meadow", "other"
        )
      )
    )
sites_splot <- read_csv("data_processed_sites_splot.csv",
                             col_names = TRUE, na = c("na", "NA", ""),
                             col_types = cols(
                               .default = "?",
                               survey_year = "c"
                               )) %>%
  mutate(
    reference = if_else(
      esy == "E12a", "Dry grassland", if_else(
        esy == "E22", "Hay meadow", "other"
        )
      ),
    target_type = if_else(
      esy == "E12a", "Dry grassland", if_else(
        esy == "E22", "Hay meadow", "other"
      )
    ),
    exposition = "other"
    )

sites_bauer <- read_csv("data_processed_sites_bauer.csv",
                             col_names = TRUE, na = c("na", "NA", ""),
                             col_types = cols(
                               .default = "?",
                               survey_year = "c"
                               )) %>%
  filter(exposition == "south" | exposition == "north") %>%
  mutate(
    reference = if_else(
      esy == "R12", "Dry grassland", if_else(
        esy == "R22", "Hay meadow", if_else(
          esy == "R", "Grassland", if_else(
            esy == "?", "no", "Fail"
          )
        )
      )
    ),
    target_type = if_else(
      esy == "R12", "Dry grassland", if_else(
        esy == "R22", "Hay meadow", "other"
        )
      )
    )
sites <- sites_experiment %>%
  bind_rows(sites_splot, sites_bauer) %>%
  select(
    id, esy, reference,
    exposition, sand_ratio, substrate_depth, target_type, seed_density,
    survey_year, longitude, latitude, elevation, plot_size
    )

#### * Load species data ####

species_experiment <- read_csv("data_processed_species.csv",
                          col_names = TRUE, na = c("na", "NA", ""),
                          col_types = cols(.default = "?"))
species_splot <- read_csv("data_processed_species_splot.csv",
                          col_names = TRUE, na = c("na", "NA", ""),
                          col_types = cols(.default = "?"))
species_bauer <- read_csv("data_processed_species_bauer.csv",
                        col_names = TRUE, na = c("na", "NA", ""),
                        col_types = cols(.default = "?"))

### Exclude rare species (< 0.5% accumulated cover in all plots)
data <- species_experiment %>%
  full_join(species_splot, by = "name") %>%
  full_join(species_bauer, by = "name") %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  group_by(name) %>%
  summarise(total_cover_species = sum(value)) %>%
  filter(total_cover_species < 0.5)

species <- species_experiment %>%
  full_join(species_splot, by = "name") %>%
  full_join(species_bauer, by = "name") %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  filter(!(name %in% data$name)) %>% # use 'data' to filter
  pivot_wider(names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(sites, by = "id")

sites <- sites %>%
  semi_join(species, by = "id")

species <- species %>%
  column_to_rownames("id")

rm(list = setdiff(ls(), c(
  "sites", "species", "theme_mb", "vegan_cov_ellipse"
  )))

#### * Model ####

set.seed(12)
(ordi <- metaMDS(species, binary = TRUE,
                 try = 50, previous.best = TRUE, na.rm = TRUE))
#Wisonsin sqrt transformation, stress type 1
stressplot(ordi) # stress: 0.205; Non-metric fit R² =.958
data_nmds <- read_csv("data_processed_nmds.csv",
                               col_names = TRUE, na = c("na", "NA", ""),
                               col_types = cols(.default = "?"))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### * Preparation ####

ellipses <- tibble()

data_nmds <-  sites %>%
  select(id, reference, exposition, target_type) %>% # modify group
  mutate(group_type = as_factor(reference), # modify group
         NMDS1 = ordi$points[, 1],
         NMDS2 = ordi$points[, 2]) %>%
  group_by(group_type) %>%
  mutate(mean1 = mean(NMDS1),
         mean2 = mean(NMDS2))
#write_csv(data_nmds, here("data", "processed", "data_processed_nmds.csv"))

for (group in levels(data_nmds$group_type)) {
  
  ellipses_calc <- data_nmds %>%
    filter(group_type == group) %>%
    with(
      cov.wt(
        cbind(NMDS1, NMDS2),
        wt = rep(1 / length(NMDS1), length(NMDS1))
      )
    )
  
  ellipses <-
    vegan_cov_ellipse(
      cov = ellipses_calc$cov, center = ellipses_calc$center
    ) %>%
    as_tibble() %>%
    bind_cols(group_type = group) %>%
    bind_rows(ellipses)
  
  data_ellipses <- ellipses
  
}

#### * Plot ####

(graph_a <- ggplot() +
   geom_point(
     aes(y = NMDS2, x = NMDS1, color = group_type, shape = group_type),
     data = data_nmds,
     cex = 2
   ) +
   facet_grid(
     exposition ~ target_type,
     labeller = as_labeller(
       c(south = "South", north = "North",
         "dry_grassland" = "Dry grassland", "hay_meadow" = "Hay meadow")
     )
   ) +
   coord_fixed() +
   theme_mb())

   geom_path(
     aes(x = NMDS1, y = NMDS2, linetype = group_type, color = group_type),
     data = data_ellipses %>% filter(group_type != "no"),
     size = 1,
     show.legend = FALSE
   ) +
   scale_color_manual(
     values = c(
       "orange1", "firebrick2", "deeppink3", "mediumpurple4",
       "royalblue", "black"
     )
   ) +
   scale_shape_manual(
     values = c(
       "circle open", "circle open", "circle open", "circle open",
       "square", "square open"
     )
   ) +
   scale_alpha_manual(values = c(.3, .3, .3, .3, .7, .6)) +
   labs(
     x = "NMDS1", y = "NMDS2", fill = "", color = "", shape = "", alpha = ""
   ) +
   theme_mb())


### Save ###
ggsave(here("outputs", "figures", "figure_nmds_800dpi_16.5x16cm.tiff"),
       dpi = 800, width = 16.5, height = 16, units = "cm")
