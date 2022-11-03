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

veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#### * Load sites data ####

sites_experiment <- read_csv("data_processed_sites.csv",
                             col_names = TRUE, na = c("na", "NA", ""),
                             col_types = cols(.default = "?"))
sites_splot <- read_csv("data_processed_sites_splot.csv",
                             col_names = TRUE, na = c("na", "NA", ""),
                             col_types = cols(.default = "?"))
sites_bauer <- read_csv("data_processed_sites_bauer.csv",
                             col_names = TRUE, na = c("na", "NA", ""),
                             col_types = cols(.default = "?"))
sites <- sites_experiment %>%
  bind_rows(sites_splot, sites_bauer)

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
  


#### * Model ####

set.seed(123)
(ordi <- metaMDS(species,
                 dist = "bray", binary = FALSE, autotransform = TRUE,
                 try = 99, previous.best = TRUE, na.rm = TRUE))
stressplot(ordi) # stress: 0.207; Non-metric fit R² =.957




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


ellipses <- tibble()

data1 <-  sites %>%
  mutate(group_type = str_c(surveyYear_fac, exposition, targetType,
                            sep = "."),
         group_type = factor(group_type))

for(group in levels(data1$group_type)) {
  
  data2 <- data1 %>%
    filter(group_type == group) %>%
      with(
        cov.wt(
        cbind(NMDS1, NMDS2),
        wt = rep(1 / length(NMDS1), length(NMDS1))
        )
      )
  
  ellipses <- 
    veganCovEllipse(cov = data2$cov, center = data2$center) %>%
    as_tibble() %>%
    bind_cols(group = group) %>%
    bind_rows(ellipses)
  
  data <- ellipses %>%
    separate(
      group,
      sep = "\\.",
      c("surveyYear_fac", "exposition", "targetType")
      )
  
}


(graph_a <- ggplot() +
    geom_point(
      aes(y = NMDS2, x = NMDS1, color = surveyYear_fac,
          shape = surveyYear_fac, alpha = surveyYear_fac),
      data = sites,
      cex = 2
    ) +
    geom_path(
      aes(x = NMDS1, y = NMDS2, color = surveyYear_fac),
      data = data,
      size = 1
    ) +
    facet_grid(
      exposition ~ targetType,
      labeller = as_labeller(
        c(south = "South", north = "North",
          "dry_grassland" = "Dry grassland", "hay_meadow" = "Hay meadow")
      )
    ) +
    coord_fixed() +
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
