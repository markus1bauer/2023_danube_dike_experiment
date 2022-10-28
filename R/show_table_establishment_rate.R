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
  mutate(across(starts_with("rate"), ~ .x * 100),
         name = str_replace_all(name, "_", " "))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Plot with gt ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



range(data$rate_2018, data$rate_2019, data$rate_2020, data$rate_2021)
mean(data$rate_2018); sd(data$rate_2018)
mean(data$rate_2019); sd(data$rate_2019)
mean(data$rate_2020); sd(data$rate_2020)
mean(data$rate_2021); sd(data$rate_2021)


(table <- data %>%
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
      style = cell_text(align = "left", style = "italic")
    ) %>%
    
    
    tab_spanner( ### Create spanner column ###
      label = "Establishment rate [%]",
      columns = starts_with("rate_")
    ) %>%
    
    
    cols_label( ### Rename column names ###
      name = md("Species name"),
      family = md("Family"),
      total_seeded_2018 = md("Seeded plots"),
      rate_2018 = md("2018"),
      rate_2019 = md("2019"),
      rate_2020 = md("2020"),
      rate_2021 = md("2021")
    ) %>%
    
    
    tab_style( ### Highlight columns ###
      locations = cells_body(
        columns = rate_2018,
        rows = rate_2018 > 0
      ),
      style = list(
        cell_fill(color = "grey95"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2019,
        rows = rate_2019 > 0
      ),
      style = list(
        cell_fill(color = "grey95"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2020,
        rows = rate_2020 > 0
      ),
      style = list(
        cell_fill(color = "grey95"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2021,
        rows = rate_2021 > 0
      ),
      style = list(
        cell_fill(color = "grey95"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2018,
        rows = rate_2018 > 33
      ),
      style = list(
        cell_fill(color = "grey80"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2019,
        rows = rate_2019 > 33
      ),
      style = list(
        cell_fill(color = "grey80"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2020,
        rows = rate_2020 > 33
      ),
      style = list(
        cell_fill(color = "grey80"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2021,
        rows = rate_2021 > 33
      ),
      style = list(
        cell_fill(color = "grey80"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2018,
        rows = rate_2018 > 50
      ),
      style = list(
        cell_fill(color = "grey63"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2019,
        rows = rate_2019 > 50
      ),
      style = list(
        cell_fill(color = "grey63"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2020,
        rows = rate_2020 > 50
      ),
      style = list(
        cell_fill(color = "grey63"),
        cell_text(weight = "bold")
      )
    ) %>%
    tab_style(
      locations = cells_body(
        columns = rate_2021,
        rows = rate_2021 > 50
      ),
      style = list(
        cell_fill(color = "grey63"),
        cell_text(weight = "bold")
      )
    ))


### Save ###
gtsave(table, here("outputs", "tables", "table_a2_establishment.docx"))
