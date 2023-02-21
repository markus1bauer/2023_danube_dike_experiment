# Grassland experiment on dikes
# Appendix A4: Total species list of sown and derived species ####
# Markus Bauer
# 2023-02-21



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
library(gt)

### Start ###
rm(list = ls())

### Load data ###
species <- read_csv(
  here("data", "processed", "data_processed_species.csv"),
  col_names = TRUE,
  na = c("", "NA", "na"),
  col_types = cols(.default = "?")
  )



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot total species list ###################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



data <- species %>%
  mutate(across(where(is.numeric), ~1 * (. != 0))) %>% 
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = TRUE),
         .keep = "unused") %>%
  group_by(name) %>%
  summarise(sum = sum(sum))

(table <- data %>%
  mutate(name = str_replace_all(name, "_", " ")) %>%
  rename(Name = name, "Presence [# plots]" = sum) %>%
  gt() %>%
  
  opt_table_lines("none") %>% ### Set general options ###
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
    style = cell_borders(
      sides = "top", color = "black", style = "solid", weight = px(1)
      ) 
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
    locations = cells_body(columns = "Name"),
    style = cell_text(align = "left", style = "italic")
  ) %>%
  
  tab_source_note(source_note = md(
    "In total 274 species occured 2018–2021"
    )))

### Save ###
write_csv(data, here("outputs", "tables", "table_a4_specieslist.csv"))
gtsave(table, here("outputs", "tables", "table_a4_specieslist.html"))
