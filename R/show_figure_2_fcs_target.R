# Dike grassland experiment
# Show_figure species composition ####
# Markus Bauer
# 2022-11-09


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
France <- read_csv("red_knot.csv")
france3_mbrms <- brms::brm(pop ~ I(year - 1976) + Location.of.population,
                           data = France, family = poisson(), chains = 3,
                           iter = 3000, warmup = 1000)
France_plot <- France %>%
  add_predicted_draws(france3_mbrms)
library(tidybayes)

(model_fit <- France %>%
    add_predicted_draws(france3_mbrms) %>%  # adding the posterior distribution
    ggplot(aes(x = year, y = pop)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = France, colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("Calidris canutus abundance\n") +  # latin name for red knot
    xlab("\nYear") +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = c(0.15, 0.85)))
sites <- read_csv("data_processed_sites.csv",
                  col_names = TRUE,
                  na = c("na", "NA", ""), col_types =
                    cols(
                      .default = "?",
                      id = "f",
                      plot = "f",
                      site = "f",
                      exposition = col_factor(levels = c("north", "south")),
                      sand_ratio = "f",
                      substrate_depth = col_factor(levels = c("30", "15")),
                      target_type = "f",
                      seed_density = "f"
                    )) %>%
  filter(survey_year != "seeded") %>%
  mutate(
    n = fcs_target,
    survey_year_fct = factor(survey_year),
    botanist_year = str_c(botanist, survey_year_fct, sep = "_"),
    botanist_year = str_replace_all(botanist_year, " ", "_")
  ) %>%
  select(
    id, plot, site, botanist_year, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, n
  )

### * Model ####
load(file = "model_fcs_3.Rdata")

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
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0, 0, 0, 0, "cm")
  )
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot ######################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Boxplots #################################################################

### a sandRatio x target Type -------------------------------------------------

test <- sites %>%
  add_predicted_draws(m3, allow_new_levels = TRUE)
  
ggplot(data = sites) +
  #geom_boxplot(aes(y = n, x = sand_ratio, color = target_type)) +
  geom_quasirandom(
    aes(y = n, x = sand_ratio, color = target_type),
    alpha = 0.2,
    dodge.width = 0.9,
    cex = .5
    ) +
  stat_pointinterval(
    aes(y = .prediction, x = sand_ratio, color = target_type),
    data = test,
    .width = c(0.66, 0.95),
    point_size = 2,
    position = "dodge"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = .3,
    color = "black"
  ) +
  facet_grid(
    exposition ~ survey_year_fct,
    labeller = as_labeller(
      c(south = "South", north = "North",
        "2018" = "2018", "2019" = "2019", "2020" = "2020", "2021" = "2021")
    )
  ) +
  scale_y_continuous(limits = c(-2.8, 1.9), breaks = seq(-100, 400, 1)) +
  scale_color_manual(labels = c("Hay meadow", "Dry grassland"),
                     values = c("#00BFC4", "#F8766D")) +
  scale_fill_manual(labels = c("Hay meadow", "Dry grassland"),
                    values = c("#00BFC4", "#F8766D")) +
  labs(
    x = "Sand ratio [%]", fill = "", color = "",
    y = expression(
      Favourable ~ Conservation ~ Status ~ "(FCS)"
    )
  ) +
  theme_mb()
  

### Save ###
ggsave(here("outputs", "figures",
            "figure_2_fcs_target_800dpi_27x9cm.tiff"),
       dpi = 800, width = 27, height = 9, units = "cm")
ggsave(here("outputs", "figures",
            "figure_2_fcs_target_800dpi_16.5x14cm.tiff"),
       dpi = 800, width = 16.5, height = 14, units = "cm")

### b sandRatio x substrateDepth -----------------------------------------------

(graph_b <- ggplot() +
   geom_quasirandom(
     aes(y = n, x = substrateDepth, color = sandRatio),
     data = sites,
     alpha = 0.5,
     dodge.width = 0.8,
     cex = .5
   ) +
   geom_hline(
     yintercept = 0,
     linetype = "dashed",
     size = .3,
     color = "grey70"
   ) +
   geom_boxplot(
     aes(y = n, x = substrateDepth, fill = sandRatio),
     data = sites,
     alpha = 0.5
   ) +
   facet_grid(
     exposition ~ surveyYear_fac,
     labeller = as_labeller(
       c(south = "South", north = "North",
         "2018" = "2018", "2019" = "2019", "2020" = "2020", "2021" = "2021")
     )
   ) +
   scale_y_continuous(limits = c(-2.8, 1.9), breaks = seq(-100, 400, 1)) +
   scale_color_manual(values = c("#990000", "#CC6600", "#FFFF00")) +
   scale_fill_manual(values = c("#990000", "#CC6600", "#FFFF00")) +
   labs(
     x = "Substrate depth [cm]", fill = "Sand ratio [%]",
     color = "Sand ratio [%]",
     y = expression(
       Favourable ~ Conservation ~ Status ~ "(FCS)"
     )
   ) +
   theme_mb())

### Save ###
ggsave(here("outputs", "figures",
            "figure_box_fcs_target_substrate_800dpi_27x9cm.tiff"),
       dpi = 800, width = 27, height = 9, units = "cm")
ggsave(here("outputs", "figures",
            "figure_box_fcs_target_substrate_800dpi_16.5x14cm.tiff"),
       dpi = 800, width = 16.5, height = 14, units = "cm")

## 2 Marginal effects ##########################################################

(graph_b <- ggeffects::ggemmeans(
  m,
  terms = c(
    "sandRatio", "targetType", "exposition", "surveyYear_fac"
  ),
  ci.lvl = 0.95,
  type = "fixed",
  typical = "mean",
  back.transform = TRUE,
  ppd = TRUE
) %>%
  mutate(facet = fct_relevel(facet, "north", "south")) %>%
  ggplot() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    size = .3,
    color = "grey70"
  ) +
  geom_point(
    aes(
      x = x, y = predicted, color = group
    ),
    position = position_dodge(.5)
  ) +
  geom_linerange(
    aes(
      x = x, color = group, ymin = conf.low, ymax = conf.high
    ),
    position = position_dodge(.5)
  ) +
  facet_grid(
    facet ~ panel,
    labeller = as_labeller(
      c(south = "South", north = "North",
        "2018" = "2018", "2019" = "2019", "2020" = "2020", "2021" = "2021")
    )
  ) +
  scale_y_continuous(limits = c(-2.1, 1.2), breaks = seq(-100, 400, 1)) +
  scale_color_manual(labels = c("Hay meadow", "Dry grassland"),
                     values = c("#00BFC4", "#F8766D")) +
  scale_fill_manual(labels = c("Hay meadow", "Dry grassland"),
                    values = c("#00BFC4", "#F8766D")) +
  labs(
    x = "Sand ratio [%]", color = "",
    y = expression(
      Favourable ~ Conservation ~ Status ~ "(FCS)"
    )
  ) +
  theme_mb())

### Save ###
ggsave(here("outputs", "figures",
            "figure_stat_fcs_target_800dpi_16.5x14cm.tiff"),
       dpi = 800, width = 16.5, height = 14, units = "cm")
