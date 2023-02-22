# Dike grassland field experiment
# Community-weighted means of seed mixtures ####
# Show Appendix A3

# Markus Bauer
# 2023-02-08



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(naniar)
library(patchwork)

### Start ###
rm(list = ls())

data <- read_csv(
  here("data", "processed", "data_processed_sites_spatial.csv"),
  col_names = TRUE,
  na = c("", "NA", "na"),
  col_types =
    cols(
      target_type = "f"
      )
  ) %>%
  select(target_type, survey_year, starts_with("cwm")) %>%
  filter(survey_year == "seeded") %>%
  pivot_longer(starts_with("cwm"), names_to = "type", values_to = "value") %>%
  mutate(
    target_type = fct_recode(
      target_type,
      "Dry\ngrassland" = "dry_grassland",
      "Hay\nmeadow" = "hay_meadow"
    ))

### Functions ###
theme_mb <- function() {
  theme(
    panel.background = element_rect(fill = "transparent"),
    text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 10),
    axis.text = element_text(angle = 0, hjust = 0.5, size = 9,
                             color = "black"),
    axis.title = element_text(angle = 0, hjust = 0.5, size = 9,
                              color = "black"),
    axis.line = element_line(),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot #######################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(graph_a <- data %>%
    filter(type == "cwm_abu_sla") %>%
    ggplot(aes(x = target_type, y = value, color = target_type)) +
    geom_boxplot() +
    labs(x = "", y = expression("SLA [mm² mg"^-1*"]")) +
    scale_color_manual(breaks = c("Hay\nmeadow", "Dry\ngrassland"),
                       values = c("#00BFC4", "#F8766D")) +
    theme_mb())

(graph_b <- data %>%
    filter(type == "cwm_abu_seedmass") %>%
    ggplot(aes(x = target_type, y = value, color = target_type)) +
    geom_boxplot() +
    labs(x = "", y = expression("Seed mass [mg]")) +
    scale_color_manual(breaks = c("Hay\nmeadow", "Dry\ngrassland"),
                       values = c("#00BFC4", "#F8766D")) +
    theme_mb())

(graph_c <- data %>%
    filter(type == "cwm_abu_height") %>%
    ggplot(aes(x = target_type, y = value, color = target_type)) +
    geom_boxplot() +
    labs(x = "", y = expression("Canopy height [m]")) +
    scale_color_manual(breaks = c("Hay\nmeadow", "Dry\ngrassland"),
                       values = c("#00BFC4", "#F8766D")) +
    theme_mb())

graph_a + theme(legend.position = "none") +
  graph_b + theme(legend.position = "none") +
  graph_c + theme(legend.position = "none") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

### Save ###

ggsave(
  here("outputs", "figures", "figure_a3_800dpi_16.5x5cm_traits.tiff"),
  dpi = 800, width = 16.5, height = 5, units = "cm"
)
