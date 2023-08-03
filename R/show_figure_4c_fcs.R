# Dike grassland field experiment
# Favourable conservation status ####
# Show figure 4

# Markus Bauer
# 2022-12-12



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(tidybayes))
library(ggbeeswarm)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data", "processed"))

### Load data ###
sites <- read_csv(
  here("data", "processed", "data_processed_sites_spatial.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types =
    cols(
      .default = "?",
      plot = "f",
      site = "f",
      sand_ratio = "f",
      substrate_depth = col_factor(levels = c("30", "15")),
      target_type = col_factor(levels = c("hay_meadow", "dry_grassland")),
      seed_density = "f",
      exposition = col_factor(levels = c("north", "south")),
      survey_year = "c"
    )
) %>%
  ### Exclude data of seed mixtures
  filter(survey_year != "seeded") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    botanist_year = str_c(survey_year, botanist, exposition, sep = " "),
    botanist_year = factor(botanist_year),
    n = fcs_target,
    id = factor(id)
  ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, survey_year, botanist_year, n
  )

### * Model ####
base::load(file = here("outputs", "models", "model_fcs_2.Rdata"))

model <- sites %>%
  add_epred_draws(m2, allow_new_levels = TRUE)

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
    legend.text = element_text(size = 9),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## 1 Marginal effects #########################################################


(p1 <- ggplot(data = sites) +
   geom_quasirandom(
     aes(y = n, x = sand_ratio, color = target_type),
     alpha = 0.2,
     dodge.width = 0.9,
     cex = .5
   ) +
   geom_hline(
     yintercept = 0,
     linetype = "dashed",
     linewidth = .3,
     color = "black"
   ) +
   tidybayes::stat_pointinterval(
     aes(y = .epred, x = sand_ratio, color = target_type),
     data = model,
     point_interval = median_qi,
     .width = c(0.66, 0.95),
     point_size = 2,
     position = "dodge"
   ) +
   facet_grid(
     exposition ~ survey_year_fct,
     labeller = as_labeller(
       c(south = "South", north = "North",
         "2018" = "2018", "2019" = "2019", "2020" = "2020", "2021" = "2021")
     )
   ) +
   scale_y_continuous(limits = c(-2.8, 2.0), breaks = seq(-100, 400, 1)) +
   scale_color_manual(
     breaks = c("hay_meadow", "dry_grassland"),
     labels = c("Hay meadow", "Dry grassland"),
     values = c("#00BFC4", "#F8766D")
     ) +
   labs(
     x = "Sand ratio [%]",
     y = expression(Favourable ~ Conservation ~ Status ~ "(FCS)"),
     color = ""
   ) +
   theme_mb())

### Save ###

ggsave(
  here("outputs", "figures", "figure_4c_fcs_600dpi_24x8cm.tiff"),
  dpi = 600, width = 24, height = 8, units = "cm"
  )

p1 + theme(legend.position = "bottom")
ggsave(
  here("outputs", "figures", "figure_4c_fcs_600dpi_16.5x8cm.tiff"),
  dpi = 600, width = 16.5, height = 8, units = "cm"
  )

(graph_c <- p1 +
    theme(
      strip.background.x = element_blank(),
      strip.text.x = element_blank()
    ))
