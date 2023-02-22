# Grassland experiment on dikes
# Establishment rates ####
# Show Appendix A3

# Markus Bauer
# 2023-02-21



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)

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
# B Plot with gt ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



data <- species %>%
  select(name, ends_with("seeded")) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  mutate(name = str_replace_all(name, "_", " ")) %>%
  pivot_longer(-name, names_to = "id", values_to = "seeded") %>%
  separate(
    id, c("plot"), sep = "_(?!.*_)",
    remove = TRUE, extra = "drop", fill = "warn", convert = FALSE
    ) %>%
  filter(seeded > 0) %>%
  select(plot, name, seeded) %>%
  arrange(plot, name) %>%
  mutate(seeded = round(seeded, digits = 1))


(table <- data %>%
    gt(
      groupname_col = "plot"
    ) %>%
    
    
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
      locations = cells_body(columns = "name"),
      style = cell_text(align = "left", style = "italic")
    ) %>%
    
    
    cols_label( ### Rename column names ###
      name = md("Species name"),
      seeded = md("Seed ratio of seed mixture [wt%]")
    )
  )

### Save ###

write_csv(data, here("outputs", "tables", "table_a2_seed_mixtures.csv"))
gtsave(table, here("outputs", "tables", "table_a2_seed_mixtures.html"))
