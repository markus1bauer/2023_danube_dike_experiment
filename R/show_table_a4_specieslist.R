
### Create list with species names and their frequency ###
specieslist <- species_experiment %>%
  mutate(across(where(is.numeric), ~1 * (. != 0))) %>% 
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = TRUE),
         .keep = "unused") %>%
  group_by(name) %>%
  summarise(sum = sum(sum))
write_csv(specieslist, here("outputs", "tables", "specieslist_20221102.csv"))