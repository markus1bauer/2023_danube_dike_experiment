# Grassland experiment on dikes
# Establishment rates ####
# Markus Bauer
# 2022-10-20



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(gt)

### Start ###
rm(list = ls())
setwd(here("data", "processed"))

### Load data ###
data <- read_csv("data_processed_traits.csv",
                          col_names = TRUE, na = c("na", "NA"), col_types =
                            cols(
                              .default = "?",
                              name = "f"
                            )) %>%
  select(name, family, total_seeded_2018, starts_with("rate_")) %>%
  drop_na(total_seeded_2018) %>%
  arrange(match(family, c("Poaceae", "Fabaceae"))) %>%
  mutate(across(starts_with("rate"), ~ .x * 100))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot with gt ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



(table <- data %>%
   gt() %>%
   
   
   opt_table_lines("none") %>%
   tab_options(
     table.font.style = "Arial",
     table.font.size = px(11),
     table.font.color = "black",
     column_labels.font.weight = "bold",
     data_row.padding = px(4),
     table.align = "left",
     column_labels.border.top.style = "solid",
     table_body.border.top.style = "solid",
     table_body.border.bottom.style = "solid",
     column_labels.border.top.color = "black",
     table_body.border.top.color = "black",
     table_body.border.bottom.color = "black",
     column_labels.border.top.width = px(2),
     table_body.border.top.width = px(1),
     table_body.border.bottom.width = px(2)
   ) %>%
   tab_style(
     locations = cells_column_labels(),
     style = cell_borders(sides = "top", color = "black",
                          style = "solid", weight = px(1)) 
   ) %>%
   
   
   tab_style(
     locations = cells_column_labels(),
     style = cell_text(align = "center")
   ) %>%
   tab_style(
     locations = cells_body(),
     style = cell_text(align = "center")
   ) %>%
   tab_style(
     locations = cells_body(columns = "name"),
     style = cell_text(align = "left")
   ) %>%
   
   
   tab_spanner(
     label = "Establishment rate [%]",
     columns = starts_with("rate_")
   ) %>%
   
   
   cols_label(
     name = md("Species name"),
     family = md("Family"),
     total_seeded_2018 = md("Seeded plots"),
     rate_2018 = md("2018"),
     rate_2019 = md("2019"),
     rate_2020 = md("2020"),
     rate_2021 = md("2021")
   ) )#%>%
  
  ### 4 Highlight certain cells ####
tab_style(
  style = list(
    cell_fill(color = "grey"),
    cell_text(weight = "bold")
  ),
  locations = cells_body(
    columns = PC1,
    rows = PC1 >= 1.2 & PC1 < 2 | PC1 <= -1.16
  )
) %>%
  tab_style(
    style = list(
      cell_fill(color = "grey"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = PC2,
      rows = PC2 >= .91 & PC2 < 2 | PC2 <= -1
    )
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "grey"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      columns = PC3,
      rows = PC3 >= .9 & PC3 < 1.15 | PC3 <= -.88
    )
  ) %>%
  
  ### 5 Make subscripts ####
text_transform(
  locations = cells_body(
    columns = c(variables),
    rows = c(6, 9, 10)
  ),
  fn = function(x) {
    x <- c("3", "total", "concentration")
    text <- c("CaCO", "N", "N")
    glue::glue("{text}<sub>{x}</sub>")
  }
) %>%
  text_transform(
    locations = cells_body(
      columns = c(variables),
      rows = c(13)
    ),
    fn = function(x) {
      x <- c("2+")
      text <- c("Mg")
      glue::glue("{text}<sup>{x}</sup>")
    }
  ))


### Save ###
gtsave(table, here("outputs", "tables", "table_a1_establishment.png"))
