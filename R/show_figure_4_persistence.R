# Dike grassland field experiment
# Favourable conservation status ####
# Show figure

# Markus Bauer
# 2022-11-14



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)
library(ggrepel)

### Start ###
rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))

### Load data ###
sites <- read_csv(
  here("data", "processed", "data_processed_sites_temporal.csv"),
  col_names = TRUE, na = c("na", "NA", ""), col_types =
    cols(
      .default = "?",
      plot = "f",
      site = "f",
      sand_ratio = "f",
      substrate_depth = "f",
      target_type = col_factor(levels = c(
        "dry_grassland", "hay_meadow"
      )),
      seed_density = "f",
      exposition = col_factor(levels = c(
        "north", "south"
      )),
      survey_year = "c"
    )
  ) %>%
  ### Exclude data of seed mixtures
  filter(presabu == "presence") %>%
  mutate(
    survey_year_fct = factor(survey_year),
    id = factor(id),
    n = persistence
  ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, survey_year, n
  )

### * Model ####
load(file = here("data", "processed", "model_persistence_1.Rdata"))

model <- sites %>%
  add_epred_draws(m1, allow_new_levels = TRUE)

### * Functions ####
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



## 1 Expectations of predicted draws ###########################################

(p1 <- ggplot(data = sites) +
   geom_quasirandom(
     aes(y = n, x = sand_ratio, color = target_type),
     alpha = 0.2,
     dodge.width = 0.9,
     cex = .5
   ) +
   geom_hline(
     yintercept = c(25, 75),
     linetype = "dashed",
     linewidth = .3,
     color = "black"
   ) +
   geom_hline(
     yintercept = c(50),
     linetype = "solid",
     linewidth = .3,
     color = "black"
   ) +
   stat_pointinterval(
     aes(y = .epred, x = sand_ratio, color = target_type),
     data = model,
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
   scale_y_continuous(limits = c(0, 100), breaks = seq(-100, 400, 10)) +
   scale_color_manual(labels = c("Hay meadow", "Dry grassland"),
                      values = c("#00BFC4", "#F8766D")) +
   scale_fill_manual(labels = c("Hay meadow", "Dry grassland"),
                     values = c("#00BFC4", "#F8766D")) +
   labs(
     x = "Sand ratio [%]", fill = "", color = "",
     y = expression("Persistence [%]")
   ) +
   theme_mb())

### Save ###

ggsave(here("outputs", "figures",
            "figure_4_persistence_epred_800dpi_27x9cm.tiff"),
       dpi = 800, width = 27, height = 9, units = "cm")

p1 + theme(legend.position = "bottom")
ggsave(here("outputs", "figures",
            "figure_4_persistence_epred_800dpi_16.5x14cm.tiff"),
       dpi = 800, width = 16.5, height = 14, units = "cm")


## 2 Coefficients #############################################################

m1 %>%
  gather_draws(
    b_sand_ratio25, b_sand_ratio50, b_substrate_depth30,
    b_target_typehay_meadow, b_seed_density8,
    b_expositionsouth, b_survey_year_fct2019,
    b_survey_year_fct2020, b_survey_year_fct2021
  ) %>%
  mutate(
    .variable = as.factor(.variable),
    .variable = fct_relevel(
      .variable, "b_sand_ratio25", "b_sand_ratio50", "b_substrate_depth30",
      "b_target_typehay_meadow", "b_seed_density8",
      "b_expositionsouth", "b_survey_year_fct2019",
      "b_survey_year_fct2020", "b_survey_year_fct2021"
    ),
    .variable = fct_relevel(.variable, rev),
    .variable = fct_recode(
      .variable,
      "Dry grassland vs. hay meadow" = "b_target_typehay_meadow",
      "Sand ratio: 0 vs. 25 %" = "b_sand_ratio25",
      "Sand ratio: 0 vs. 50 %" = "b_sand_ratio50",
      "Substrate depth: 15 vs. 30 cm" = "b_substrate_depth30",
      "Seed density: 4 vs. 8 g/m²" = "b_seed_density8",
      "North vs. south" = "b_expositionsouth",
      "2018 vs. 2019" = "b_survey_year_fct2019",
      "2018 vs. 2020" = "b_survey_year_fct2020",
      "2018 vs. 2021" = "b_survey_year_fct2021"
    )
  ) %>%
  mutate(.value = ((1 - .value) * 100) - 100) %>%
  ggplot(aes(x = .value, y = .variable)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  stat_halfeye() +
  scale_x_continuous(breaks = seq(-100, 400, 10)) +
  labs(x = expression(Delta ~ Persistence ~ "[%]"), y = "") +
  theme_mb()

### Save ###

ggsave(here("outputs", "figures",
            "figure_4_persistence_coef_800dpi_27x9cm.tiff"),
       dpi = 800, width = 27, height = 9, units = "cm")
ggsave(here("outputs", "figures",
            "figure_4_persistence_coef_800dpi_16.5x14cm.tiff"),
       dpi = 800, width = 16.5, height = 14, units = "cm")
