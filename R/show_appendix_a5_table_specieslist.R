# Grassland experiment on dikes
# Total species list of sown and derived species ####
# Show Appendix A5

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
traits <- read_csv(
  here("data", "processed", "data_processed_traits.csv"),
  col_names = TRUE,
  na = c("", "NA", "na"),
  col_types = cols(.default = "?")
  ) %>%
  select(name, target, seeded)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot total species list ###################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Number of established plots of Bromus hordeaceus ###
species %>%
  filter(name == "Bromus_hordeaceus") %>%
  mutate(across(where(is.numeric), ~1 * (. != 0))) %>%
  pivot_longer(cols = -name, names_to = "survey_year", values_to = "n") %>%
  mutate(
    survey_year = str_extract(
      survey_year, "[:digit:][:digit:][:digit:][:digit:]"
      )
    ) %>%
  group_by(name, survey_year) %>%
  summarise(total_established = sum(n, na.rm = TRUE))

### Plot the table ###
data <- species %>%
  mutate(across(where(is.numeric), ~1 * (. != 0))) %>%
  mutate(
    sum = rowSums(across(where(is.numeric)), na.rm = TRUE), .keep = "unused"
    ) %>%
  group_by(name) %>%
  summarise(sum = sum(sum)) %>%
  left_join(traits, by = "name") %>%
  select(name, seeded, target, sum)

(table <- data %>%
    mutate(name = str_replace_all(name, "_", " ")) %>%
    rename(
      Name = name,
      Seeded = seeded,
      "Target species" = target,
      "Presence [# plots]" = sum
    ) %>%
    gt() %>%
    sub_missing(columns = 2:3, missing_text = "") %>% # replace NA
    opt_table_lines("none") %>% ### Set general style options ###
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
    tab_source_note(source_note = md( ### Footnote ###
      "In total 274 species occured 2018–2021"
    )))

### Save ###

write_csv(data, here("outputs", "tables", "table_a5_specieslist.csv"))
gtsave(table, here("outputs", "tables", "table_a5_specieslist.html"))
