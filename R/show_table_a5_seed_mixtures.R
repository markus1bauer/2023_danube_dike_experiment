### c Table of seedmixes -------------------------------------------------------

data <- species_experiment %>%
  select(name, ends_with("seeded")) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "seeded") %>%
  separate(id, c("plot"), sep = "_(?!.*_)",
           remove = TRUE, extra = "drop", fill = "warn", convert = FALSE) %>%
  filter(seeded > 0) %>%
  select(plot, name, seeded) %>%
  arrange(plot, name)
write_csv(data, here("outputs", "tables", "table_a2_seedmixes.csv"))