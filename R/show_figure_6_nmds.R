# Dike grassland field experiment
# NMDS ordination ####
# Show figure 6

# Markus Bauer
# 2022-11-15



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d", "ordi")))
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
    legend.text = element_text(size = 9),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}

vegan_cov_ellipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi / npoints
  circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(circle %*% chol(cov)))
}

#### Load sites data and model ###

sites_nmds <- read_csv(
  "data_processed_sites_nmds.csv",
  col_names = TRUE,
  na = c("na", "NA", ""),
  col_types = cols(.default = "?")
)

sites_nmds %>%
  group_by(source, esy) %>%
  filter(esy == "R1A" | esy == "R22" | esy == "V38") %>%
  count()

base::load(file = here("outputs", "models", "model_nmds.Rdata"))
ordi



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### 1 Preparation #############################################################


ellipses <- tibble()

### Copy references to both expositions ###
data <- sites_nmds %>%
  filter(exposition == "other") %>%
  mutate(exposition = "north")
data <- data %>%
  mutate(exposition = "south") %>%
  bind_rows(data)
sites_nmds <- sites_nmds %>%
  filter(!(exposition == "other")) %>%
  bind_rows(data)

### Copy negative references to both target_types ###
data <- sites_nmds %>%
  filter(target_type == "other") %>%
  mutate(target_type = "hay_meadow")
data <- data %>%
  mutate(target_type = "dry_grassland") %>%
  bind_rows(data)
sites_nmds <- sites_nmds %>%
  filter(!(target_type == "other")) %>%
  bind_rows(data) %>%
  filter(!(reference == "negative_reference" & exposition == "north"))

sites_nmds <- sites_nmds %>%
  ### Select variables to plot ###
  select(
    id, NMDS1, NMDS2, reference, exposition, target_type # modify group
    ) %>%
  mutate(
    reference = if_else(
      reference == "positive_reference", "+Reference", if_else(
        reference == "negative_reference", "-Reference", reference
      )
    ),
    group_type = str_c(
      reference, exposition, target_type, sep = "." # modify group
      )
    ) %>%
  group_by(group_type) %>%
  mutate(
    mean1 = mean(NMDS1),
    mean2 = mean(NMDS2),
    group_type = factor(group_type),
    reference = factor(reference)
  )

### Calculate ellipses for plot ###
for (group in levels(sites_nmds$group_type)) {
  
  ellipses_calc <- sites_nmds %>%
    filter(group_type == group) %>%
    with(
      cov.wt(
        cbind(NMDS1, NMDS2),
        wt = rep(1 / length(NMDS1), length(NMDS1))
      )
    )
  
  ellipses <- vegan_cov_ellipse(
    cov = ellipses_calc$cov, center = ellipses_calc$center
    ) %>%
    as_tibble() %>%
    bind_cols(group_type = group) %>%
    bind_rows(ellipses)
  
  data_ellipses <- ellipses %>%
    separate(
      group_type,
      sep = "\\.",
      c("reference", "exposition", "target_type")
    )
}



#### 2 Plot ###################################################################


(graph_a <- ggplot() +
   geom_point(
     aes(y = NMDS2, x = NMDS1,
         color = reference, shape = reference, alpha = reference),
     data = sites_nmds,
     cex = 2
   ) +
   facet_grid(
     exposition ~ target_type,
     labeller = as_labeller(
       c(south = "South", north = "North",
         "dry_grassland" = "Dry grassland", "hay_meadow" = "Hay meadow")
     )
   ) +
   geom_vline(xintercept = 0, linetype = "dashed") +
   geom_hline(yintercept = 0, linetype = "dashed") +
   geom_path(
     aes(x = NMDS1, y = NMDS2, color = reference),
     data = data_ellipses,
     size = 1,
     show.legend = FALSE
   ) +
   coord_fixed() +
   scale_shape_manual(
     values = c(
       "triangle", "square",
       "circle open", "circle open", "circle open", "circle open",
       "square open"
     )
   ) +
   scale_color_manual(
     values = c(
       "grey", "royalblue",
       "orange1", "firebrick2", "deeppink3", "mediumpurple4",
       "black"
     )
   ) +
   scale_alpha_manual(values = c(.7, .7, .3, .3, .3, .3, .6)) +
   labs(
     x = "NMDS1", y = "NMDS2", fill = "", color = "", shape = "", alpha = ""
   ) +
   theme_mb())


### Save ###

ggsave(here("outputs", "figures", "figure_6_nmds_800dpi_16.5x16cm.tiff"),
       dpi = 800, width = 16.5, height = 16, units = "cm")
graph_a + theme(legend.position = "right")
ggsave(here("outputs", "figures", "figure_6_nmds_800dpi_18x13cm.tiff"),
       dpi = 800, width = 18, height = 13, units = "cm")
