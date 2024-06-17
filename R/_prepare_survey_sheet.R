# Grassland experiment on dikes
# Species list of each plot for survey ####
# Print survey sheets

# Markus Bauer
# 2024-06-15



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



## After this script prepare in Excel the following steps:
# size 9 of text
# insert date in header
# insert persons in header
# insert page number in footnote
# print each plot on one page: Daten -> Teilergebnis -> sort by plot -> tick "plot" + tick "break page for each group"
# replicate row names on each page: Seitenlayout -> Drucktitel -> Weiederholungszeilen oben -> "$1:$1"

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
# B Plot species of each survey ###############################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


data <- species %>%
  pivot_longer(-name, names_to = "plot", values_to = "n") %>%
  separate(plot, c("block", "plot", "survey_year"), sep = "_", remove = FALSE) %>%
  pivot_wider(names_from = "survey_year", values_from = "n") %>%
  select(block, plot, seeded, "2019", "2020", "2021", name) %>%
  arrange(block, plot, name) %>%
  mutate(
    seeded = as.character(seeded),
    sum = rowSums(across(where(is.numeric))),
    current = ".5|2|4|10|20|30|40|50|60|70|80|90|97.5"
  ) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  group_by(plot)

write_excel_csv(data, here("outputs", "tables", "survey_sheet.csv"))
