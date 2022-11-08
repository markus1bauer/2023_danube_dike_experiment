# Dike grassland experiment
# Show_figure species composition ####
# Markus Bauer
# 2022-04-27


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(ggbeeswarm)

### Start ###
#rm(list = setdiff(ls(), c("graph_a", "graph_b", "graph_c", "graph_d")))
setwd(here("data", "processed"))

### Load data ###
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
    survey_year_fct = factor(survey_year)
  ) %>%
  select(
    id, plot, site, exposition, sand_ratio, substrate_depth, target_type,
    seed_density, survey_year_fct, n
  )

### * Model ####
m <- m3

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
get_variables(m)[1:20]
m %>%
  spread_draws(b_target_typedry_grassland, b_expositionnorth, b_sand_ratio25,
               b_sand_ratio50, b_survey_year_fct2019, b_survey_year_fct2020,
               b_survey_year_fct2021) %>%
  summarise_draws() %>%
  ggplot(aes(x = mean, xmin = q5, xmax = q95, y = variable)) +
  geom_pointinterval()
m %>%
  gather_draws(b_Intercept, b_target_typedry_grassland, b_expositionnorth,
               b_sand_ratio25, b_sand_ratio50, b_survey_year_fct2019,
               b_survey_year_fct2020, b_survey_year_fct2021) %>%
  filter(.variable == "b_Intercept" | .variable == "b_target_typedry_grassland") %>%
  ggplot(aes(x = .value, y = .variable)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_mb()
m %>%
  modelr::data_grid(target_type) 
pred.rich <- posterior_epred(treat.eff.rich, newdata = newdata, 
                             scale = "response", 
                             re_formula = ~ (1|region.year), 
                             summary = FALSE
)
(graph_a <- ggplot() +
    geom_quasirandom(
      aes(y = n, x = sand_ratio, color = target_type),
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
    #geom_boxplot(
    #  aes(y = n, x = sand_ratio, fill = target_type),
    #  data = sites,
    #  alpha = 0.5
    #) +
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
    theme_mb())

### Save ###
ggsave(here("outputs", "figures",
            "figure_box_fcs_target_800dpi_27x9cm.tiff"),
       dpi = 800, width = 27, height = 9, units = "cm")
ggsave(here("outputs", "figures",
            "figure_box_fcs_target_800dpi_16.5x14cm.tiff"),
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
