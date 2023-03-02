# Dike grassland field experiment
# Favourable Conservation Status (FCS) ####
# Show figure 3

# Markus Bauer
# 2022-12-12



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(ggrepel)
library(tidybayes)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))

### Load data ###
sites <- read_csv(
  here("data", "processed", "data_processed_sites_temporal.csv"),
  col_names = TRUE,
  na = c("na", "NA", ""),
  col_types = cols(
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
  filter(presabu == "presence") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    botanist_year = str_c(survey_year, botanist, exposition, sep = " "),
    botanist_year = factor(botanist_year),
    id = factor(id),
    n = persistence / 100
  ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, survey_year, botanist_year, n
  )

### * Model ####
load(file = here("outputs", "models", "model_persistence_2.Rdata"))

model <- sites %>%
  tidybayes::add_epred_draws(m2, allow_new_levels = TRUE)

### Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "white"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_text(
      angle = 0, hjust = 0.5, size = 9, color = "black"
    ),
    axis.title = element_text(
      angle = 0, hjust = 0.5, size = 9, color = "black"
    ),
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


## 1 Coefficients #############################################################


get_variables(m2)

(p1 <- m2 %>%
    tidybayes::gather_draws(
      b_sand_ratio25, b_sand_ratio50, b_substrate_depth15,
      b_target_typedry_grassland, b_seed_density8
    ) %>%
    mutate(
      .variable = as.factor(.variable),
      .variable = fct_relevel(
        .variable, "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth15",
        "b_target_typedry_grassland", "b_seed_density8"
      ),
      .variable = fct_relevel(.variable, rev),
      .variable = fct_recode(
        .variable,
        "Hay meadow vs. Dry grassland" = "b_target_typedry_grassland",
        "Sand ratio: 0 vs. 25 %" = "b_sand_ratio25",
        "Sand ratio: 0 vs. 50 %" = "b_sand_ratio50",
        "Substrate depth: 30 vs. 15 cm" = "b_substrate_depth15",
        "Seed density: 4 vs. 8 g/m²" = "b_seed_density8"
      )
    ) %>%
    ggplot(aes(x = .value, y = .variable)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggdist::stat_halfeye(point_interval = "median_qi", .width = c(0.66, 0.95)) +
    scale_x_continuous(breaks = seq(-100, 400, .1)) +
    labs(x = expression(Delta ~ Persistence ~ "[# / 20]"), y = "") +
    theme_mb())

### Save ###

ggsave(
  here("outputs", "figures", "figure_3b_persistence_800dpi_24x8cm.tiff"),
  dpi = 800, width = 24, height = 8, units = "cm"
)
(graph_b <- p1 +
    labs(x = expression(Delta ~ Persistence)) +
    theme(
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    ))
