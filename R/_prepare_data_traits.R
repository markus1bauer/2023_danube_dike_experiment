#______________________________________________________________________________
## 7 CWM of Ellenberg ##########################################################


### a N value ------------------------------------------------------------------
data_traits <- traits %>%
  select(name, n) %>%
  filter(n > 0)
data_species <- semi_join(species, data_traits, by = "name") %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = TRUE)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id")
data_traits <- column_to_rownames(data_traits, "name")
### Calculate CWM ###
Nweighted <- dbFD(data_traits, data_species, w.abun = TRUE,
                  calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt")

### b F value ------------------------------------------------------------------
data_traits <- traits %>%
  select(name, f) %>%
  filter(f > 0)
data_species <- semi_join(species, data_traits, by = "name") %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = TRUE)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id")
data_traits <- column_to_rownames(data_traits, "name")
### Calculate CWM ###
Fweighted <- dbFD(data_traits, data_species, w.abun = TRUE,
                  calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt")

### c T value ------------------------------------------------------------------
data_traits <- traits %>%
  select(name, t) %>%
  filter(t > 0)
data_species <- semi_join(species, data_traits, by = "name") %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = TRUE)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id")
data_traits <- column_to_rownames(data_traits, "name")
### Calculate CWM ###
Tweighted <- dbFD(data_traits, data_species, w.abun = TRUE,
                  calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt")

### d implement in sites data set ----------------------------------------------
sites$cwmAbuN <- round(as.numeric(as.character(Nweighted$CWM$n)), 3)
sites$cwmAbuF <- round(as.numeric(as.character(Fweighted$CWM$f)), 3)
sites$cwmAbuT <- round(as.numeric(as.character(Tweighted$CWM$t)), 3)

rm(list = setdiff(ls(), c("sites", "species", "traits", "seedmixes")))



#______________________________________________________________________________
## 8 Functional plant traits ###################################################


### a LEDA data ---------------------------------------------------------------

data_sla <- data.table::fread(
  "data_raw_traitbase_leda_20210223_sla.txt",
  sep = ";",
  dec = ".",
  skip = 4,
  header = TRUE,
  select = c("SBS name", "single value [mm^2/mg]")
) %>%
  as_tibble() %>%
  rename(name = "SBS name",
         sla = "single value [mm^2/mg]") %>%
  mutate(name = str_replace_all(name, " ", "_"))

data_seedmass <- data.table::fread(
  "data_raw_traitbase_leda_20210223_seedmass.txt",
  sep = ";",
  dec = ".",
  skip = 3,
  header = TRUE,
  select = c("SBS name", "single value [mg]")
) %>%
  as_tibble() %>%
  rename(name = "SBS name",
         seedmass = "single value [mg]") %>%
  mutate(name = str_replace_all(name, " ", "_"),
         seedmass = round(seedmass, 3))

data_height <- data.table::fread(
  "data_raw_traitbase_leda_20210223_canopy_height.txt",
  sep = ";",
  dec = ".",
  skip = 3,
  header = TRUE,
  select = c("SBS name", "single value [m]")
) %>%
  as_tibble() %>%
  rename(name = "SBS name",
         height = "single value [m]") %>%
  mutate(name = str_replace_all(name, " ", "_"))
#### Join SLA, canopy height and seedmass of LEDA ###
data <- data_sla %>%
  full_join(data_seedmass, by = "name") %>%
  full_join(data_height, by = "name")

#### * Find synonyms ####
data <- data %>%
  mutate(name = as_factor(name),
         name = fct_recode(name,
                           Centaurea_stoebe = "Centaurea_stoebe_s.lat.",
                           Carex_filiformis = "Carex_tomentosa",
                           Carex_leporina = "Carex_ovalis",
                           Carex_praecox_ssp_praecox = "Carex_praecox",
                           Cerastium_fontanum_ssp_vulgare =
                             "Cerastium_fontanum",
                           Cerastium_fontanum_ssp_vulgare =
                             "Cerastium_fontanum_s._vulgare",
                           Cornus_controversa = "Cornus_sanguinea",
                           Erigeron_canadensis = "Conyza_canadensis",
                           Cota_tinctoria = "Anthemis_tinctoria",
                           Cyanus_segetum = "Centaurea_cyanus",
                           Euphorbia_verrucosa = "Euphorbia_brittingeri",
                           Helictotrichon_pubescens = "Avenula_pubescens",
                           Hypochaeris_radicata = "Hypochoeris_radicata",
                           Jacobaea_vulgaris = "Senecio_vulgaris",
                           Medicago_falcata = "Medicago_sativa_s._falcata",
                           Melilotus_albus = "Melilotus_alba",
                           Ononis_spinosa_ssp_procurrens = "Ononis_repens",
                           Persicaria_amphibia = "Polygonum_amphibium",
                           Persicaria_bistorta = "Polygonum_bistorta",
                           Persicaria_lapathifolia = "Polygonum_lapathifolium",
                           Persicaria_minor = "Polygonum_minus",
                           Persicaria_mitis = "Polygonum_mite",
                           Pilosella_caespitosa = "Hieracium_caespitosum",
                           Pilosella_officinarum = "Hieracium_pilosella",
                           Pilosella_piloselloides = "Hieracium_piloselloides",
                           Plantago_major = "Plantago_major",
                           Plantago_major_ssp_intermedia =
                             "Plantago_major_subsp._intermedia",
                           Rubus_fruticosus_agg = "Rubus_fruticosus",
                           Rubus_fruticosus_agg = "Rubus_fruticosus_ag._L.",
                           Securigera_varia = "Coronilla_varia",
                           Sedum_maximum = "Sedum_telephium_s._maximum",
                           "Silene_flos-cuculi" = "Lychnis_flos-cuculi",
                           Silene_latifolia_ssp_alba = "Silene_latifolia",
                           Taraxacum_campylodes = "Taraxacum_officinale",
                           Taraxacum_campylodes =
                             "Taraxacum_Sec._Ruderalia",
                           Tripleurospermum_maritimum =
                             "Matricaria_maritima",
                           Vicia_sativa_ssp_nigra =
                             "Vicia_sativa_s._nigra")) %>% 
  group_by(name) %>%
  ### summarise different rows of one species after renaming ###
  summarise(across(where(is.double), ~median(.x, na.rm = TRUE)))
traits$name[which(!(traits$name %in% data$name))]
#data %>%
#group_by(name) %>%
#summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
#filter(str_detect(name, "ovalis"))
traits <- traits %>%
  left_join(data, by = "name")

#### * check completeness of LEDA ####
test <- traits %>%
  select(name, sla, seedmass, height) %>%
  filter(!complete.cases(.))

rm(list = setdiff(ls(), c(
  "sites", "species", "traits", "seedmixes", "test"
)))

### b TRY data 1 --------------------------------------------------------------

data <- data.table::fread("data_raw_traitbase_try_20210306_13996.txt",
                          header = TRUE,
                          sep = "\t",
                          dec = ".",
                          quote = "") %>%
  as_tibble() %>%
  rename(name = "AccSpeciesName",
         value = "StdValue",
         trait = "TraitID") %>%
  select(name, value, trait) %>%
  mutate(name = str_replace_all(name, " ", "_")) %>%
  drop_na %>%
  mutate(trait = str_replace(trait, "26", "seedmass"),
         trait = str_replace(trait, "3106", "height"),
         trait = str_replace(trait, "3107", "height"),
         trait = str_replace(trait, "3115", "sla"),
         trait = str_replace(trait, "3116", "sla"),
         trait = str_replace(trait, "3117", "sla"))

##### * Find synonyms ####
data <- data %>%
  mutate(name = as_factor(name),
         name = fct_recode(name,
                           Carex_praecox_ssp_praecox =
                             "Carex_praecox",
                           Plantago_major_ssp_intermedia =
                             "Plantago_major_subsp._intermedia",
                           Ranunculus_serpens_ssp_nemorosus =
                             "Ranunculus_serpens_subsp._nemorosus"),
         trait = as_factor(trait)) %>%
  group_by(name, trait) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = TRUE))) %>%
  pivot_wider(names_from = "trait", values_from = "value")
test$name[which(!(test$name %in% data$name))]
#data %>%
#group_by(name) %>%
#summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
#filter(str_detect(name, "Equisetum"))
traits <- traits %>%
  left_join(data, by = "name") %>%
  mutate(sla = coalesce(sla.x, sla.y),
         seedmass = coalesce(seedmass.x, seedmass.y),
         height = coalesce(height.x, height.y),
         .keep = "unused")

### * check completeness of LEDA + TRY1 ####
(test <- traits %>%
   select(name, sla, seedmass, height) %>%
   filter(!complete.cases(sla, seedmass, height)))

### c TRY data 2 --------------------------------------------------------------

data <- data.table::fread("data_raw_traitbase_try_20210318_14157.txt",
                          header = TRUE,
                          sep = "\t",
                          dec = ".",
                          quote = "") %>%
  as_tibble() %>%
  rename(name = "AccSpeciesName",
         value = "StdValue",
         trait = "TraitID") %>%
  select(name, value, trait) %>%
  mutate(name = str_replace_all(name, " ", "_")) %>%
  drop_na %>%
  mutate(trait = str_replace(trait, "26", "seedmass"),
         trait = str_replace(trait, "3106", "height"),
         trait = str_replace(trait, "3107", "height"),
         trait = str_replace(trait, "3115", "sla"),
         trait = str_replace(trait, "3116", "sla"),
         trait = str_replace(trait, "3117", "sla"))

##### * Find Synonyms ####
data <- data %>%
  mutate(name = as_factor(name),
         name = fct_recode(name,
                           Plantago_major_ssp_intermedia =
                             "Plantago_major_subsp._intermedia",
                           Plantago_major =
                             "Plantago_major_subsp._major",
                           Cornus_controversa =
                             "Cornus_sanguinea",
                           Silene_latifolia_ssp_alba =
                             "Silene_latifolia_subsp._alba",
                           Silene_latifolia_ssp_alba =
                             "Silene_latifolia",
                           Silene_latifolia_ssp_alba =
                             "Silene_latifolia_ssp._alba"),
         trait = as_factor(trait)) %>%
  group_by(name, trait) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = TRUE))) %>%
  pivot_wider(names_from = "trait", values_from = "value")
test$name[which(!(test$name %in% data$name))]
#data %>%
#group_by(name) %>%
#summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
#filter(str_detect(name, "Silene"))
traits <- traits %>%
  left_join(data, by = "name") %>%
  mutate(sla = coalesce(sla.x, sla.y),
         seedmass = coalesce(seedmass.x, seedmass.y),
         height = coalesce(height.x, height.y),
         .keep = "unused")

### * check completeness of LEDA + TRY1 + TRY2 ####
(test <- traits %>%
   select(name, sla, seedmass, height) %>%
   filter(!complete.cases(sla, seedmass, height)))

### d GrooT data --------------------------------------------------------------

data <- read.csv("data_raw_traitbase_groot.csv", header = TRUE,
                 na.strings = c("", "NA", "na")) %>%
  filter(traitName == "Specific_root_length" |
           traitName == "Root_length_density_volume" |
           traitName == "Root_mass_fraction" |
           traitName == "Lateral_spread") %>%
  mutate(name = str_c(genusTNRS, speciesTNRS, sep = "_")) %>%
  rename(trait = traitName,
         value = medianSpecies) %>%
  select(name, trait, value) %>%
  pivot_wider(names_from = "trait", values_from = "value") %>%
  rename(rmf = Root_mass_fraction,
         rld = Root_length_density_volume,
         srl = Specific_root_length,
         lateral = Lateral_spread)

##### * Find Synonyms ####
data <- data %>%
  mutate(name = as_factor(name),
         name = fct_recode(name,
                           Cota_tinctoria = "Anthemis_tinctoria",
                           Carex_praecox_ssp_praecox = "Carex_praecox",
                           Carex_filiformis = "Carex_tomentosa",
                           Cerastium_fontanum_ssp_vulgare =
                             "Cerastium_fontanum",
                           Persicaria_bistorta = "Bistorta_officinalis",
                           Jacobaea_vulgaris = "Senecio_jacobaea",
                           "Silene_flos-cuculi" = "Lychnis_flos-cuculi",
                           Silene_latifolia_ssp_alba = "Silene_latifolia",
                           Vicia_sativa_ssp_nigra = "Vicia_sativa")) %>%
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = TRUE)))
traits$name[which(!(traits$name %in% data$name))]
#data %>%
#group_by(name) %>%
#summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
#filter(str_detect(name, "Centaurea_"))
traits <- traits %>%
  left_join(data, by = "name")

### * check completeness of Roots ####
(test <- traits %>%
   select(name, rmf, rld, srl, lateral) %>%
   filter(!complete.cases(rmf, rld, srl, lateral)))

### e check completeness of traits --------------------------------------------

test <- traits %>%
  select(name, group, t, n, f, sla, seedmass, height,
         rmf, rld, srl, lateral) %>%
  filter(!(str_detect(name, "_spec") | str_detect(name, "ceae"))) %>%
  filter(group != "tree" & group != "shrub") %>%
  select(-name, -group)
vis_miss(test, cluster = TRUE, sort_miss = TRUE)
miss_var_summary(test)
gg_miss_case(test, order_cases = TRUE)
(test <- traits %>%
    select(name, sla, seedmass, height, rmf, rld, srl, lateral) %>%
    filter(!complete.cases(sla, seedmass, height, rmf, rld, srl, lateral)))

rm(list = setdiff(ls(), c(
  "sites", "species", "traits", "seedmixes"
)))

### f prepare data frames for calculation of FD -------------------------------

traits <- traits %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(across(c(sla, seedmass, height, srl, rmf, rld),
                ~ round(.x, digits = 3)))
### gamma diversity herbs: 269 ###
(herbCount <- traits %>%
    filter(group != "tree" & group != "shrub") %>%
    count() %>%
    pull())
### gamma diversity of undefined species: 15 ###
(undefinedCount <- traits %>%
    filter((str_detect(name, "_spec") | str_detect(name, "ceae"))) %>%
    count() %>%
    pull())
data <- traits %>%
  filter(!(str_detect(name, "_spec") | str_detect(name, "ceae"))) %>%
  filter(group != "tree" & group != "shrub")
traits_lhs <- data %>%
  select(name, sla, seedmass, height) %>%
  drop_na()
traits_sla <- data %>%
  select(name, sla) %>%
  drop_na() 
traits_seedmass <- data %>%
  select(name, seedmass) %>%
  drop_na()
traits_height <- data %>%
  select(name, height) %>%
  drop_na()
traits_srl <- data %>%
  select(name, srl) %>%
  drop_na()
traits_rld <- data %>%
  select(name, rld) %>%
  drop_na()
traits_rmf <- data %>%
  select(name, rmf) %>%
  drop_na()
traits_all <- data %>%
  select(name, sla, seedmass, height, rmf) %>%
  drop_na()



#______________________________________________________________________________
## 9 CWM and FDis of functional plant traits #################################


### a LHS ----------------------------------------------------------------------

### Preparation ###
data_species <- semi_join(species, traits_lhs, by = "name")
data_traits <- semi_join(traits_lhs, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, names_to = "site", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = TRUE,
                          calc.FRic = FALSE, calc.FDiv = FALSE, corr = "cailliez")
### Integration to dataset ###
data <- data_diversityAbu$FDis %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  rename("fdisAbuLHS" = ".")
sites <- sites %>%
  left_join(data, by = "id")
length(traits_lhs$name) / (herbCount - undefinedCount)

### b SLA ----------------------------------------------------------------------

### Preparation ###
data_species <- semi_join(species, traits_sla, by = "name")
data_traits <- semi_join(traits_sla, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = TRUE,
                          calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt")
### Integration to dataset ###
data <- data_diversityAbu$FDis %>%
  as.data.frame() %>%
  add_column(data_diversityAbu$CWM$sla) %>%
  rownames_to_column("id") %>%
  rename("fdisAbuSla" = ".",
         "cwmAbuSla" = "data_diversityAbu$CWM$sla") %>%
  mutate(cwmAbuSla = exp(cwmAbuSla))
sites <- sites %>%
  left_join(data, by = "id")
### Describe
length(traits_sla$name) / (herbCount - undefinedCount)

### c Seed mass ----------------------------------------------------------------

### Preparation ###
data_species <- semi_join(species, traits_seedmass, by = "name")
data_traits <- semi_join(traits_seedmass, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, names_to = "site", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = TRUE,
                          calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt")
### Integration to dataset ###
data <- data_diversityAbu$FDis %>%
  as.data.frame() %>%
  add_column(data_diversityAbu$CWM$seedmass) %>%
  rownames_to_column("id") %>%
  rename("fdisAbuSeedmass" = ".",
         "cwmAbuSeedmass" = "data_diversityAbu$CWM$seedmass") %>%
  mutate(cwmAbuSeedmass = exp(cwmAbuSeedmass))
sites <- sites %>%
  left_join(data, by = "id")
### Describe ###
length(traits_seedmass$name) / (herbCount - undefinedCount)

### d Canopy height ------------------------------------------------------------

### Preparation ###
data_species <- semi_join(species, traits_height, by = "name")
data_traits <- semi_join(traits_height, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, names_to = "site", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = TRUE,
                          calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt")
### Integration to dataset ###
data <- data_diversityAbu$FDis %>%
  as.data.frame() %>%
  add_column(data_diversityAbu$CWM$height) %>%
  rownames_to_column("id") %>%
  rename("fdisAbuHeight" = ".",
         "cwmAbuHeight" = "data_diversityAbu$CWM$height") %>%
  mutate(cwmAbuHeight = exp(cwmAbuHeight))
sites <- sites %>%
  left_join(data, by = "id")
#Describe
length(traits_height$name) / (herbCount - undefinedCount)

### e Specific root length -----------------------------------------------------

### Preparation ###
data_species <- semi_join(species, traits_srl, by = "name")
data_traits <- semi_join(traits_srl, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, names_to = "site", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = TRUE,
                          calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt")
### Integration to dataset ###
data <- data_diversityAbu$FDis %>%
  as.data.frame() %>%
  add_column(data_diversityAbu$CWM$srl) %>%
  rownames_to_column("id") %>%
  rename("fdisAbuSrl" = ".",
         "cwmAbuSrl" = "data_diversityAbu$CWM$srl") %>%
  mutate(cwmAbuSrl = exp(cwmAbuSrl))
sites <- sites %>%
  left_join(data, by = "id")
### Describe ###
length(traits_srl$name) / (herbCount - undefinedCount)

### f Root mass fraction -------------------------------------------------------

### Prepartion ###
data_species <- semi_join(species, traits_rmf, by = "name")
data_traits <- semi_join(traits_rmf, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, names_to = "site", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = TRUE,
                          calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt")
### Integration to dataset ###
data <- data_diversityAbu$FDis %>%
  as.data.frame() %>%
  add_column(data_diversityAbu$CWM$rmf) %>%
  rownames_to_column("id") %>%
  rename("fdisAbuRmf" = ".",
         "cwmAbuRmf" = "data_diversityAbu$CWM$rmf") %>%
  mutate(cwmAbuRmf = exp(cwmAbuRmf))
sites <- sites %>%
  left_join(data, by = "id")
### Describe ###
length(traits_rmf$name) / (herbCount - undefinedCount)

### g All ----------------------------------------------------------------------

### Preparation ###
data_species <- semi_join(species, traits_all, by = "name")
data_traits <- semi_join(traits_all, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, names_to = "site", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = TRUE,
                          calc.FRic = FALSE, calc.FDiv = FALSE, corr = "cailliez")
### Integration to dataset ###
data <- data_diversityAbu$FDis %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  rename("fdisAbuAll" = ".")
sites <- sites %>%
  left_join(data, by = "id")
### Describe ###
length(traits_all$name) / (herbCount - undefinedCount)

### h Finalisation -------------------------------------------------------------

sites <- sites %>%
  mutate(across(c(fdisAbuLHS, fdisAbuSla, cwmAbuSla, fdisAbuSeedmass,
                  cwmAbuSeedmass, fdisAbuHeight, cwmAbuHeight,
                  fdisAbuSrl, cwmAbuSrl, fdisAbuRmf,
                  cwmAbuRmf, fdisAbuAll),
                ~ round(.x, digits = 3))) %>%
  select(id, tidyselect::peek_vars())

rm(list = setdiff(ls(), c(
  "sites", "species", "traits", "seedmixes"
)))
