# Prepare vegetation data for Danube experiment ####
# Markus Bauer


### Packages ###
library(here)
library(tidyverse)
library(vegan)
library(FD) #dbFD()
library(naniar) #are_na()

### Start ###
#installr::updateR(browse_news = F, install_R = T, copy_packages = T, copy_Rprofile.site = T, keep_old_packages = T, update_packages = T, start_new_R = F, quit_R = T, print_R_versions = T, GUI = F)
rm(list = ls())
setwd(here("data/raw"))


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Sites #####################################################################################

sites <- read_csv("data_raw_sites.csv", col_names = T, na = c("", "NA", "na"), col_types =
                    cols(
                      .default = "f",
                      surveyDate_seeded = col_date(format = "%Y-%m-%d"),
                      surveyDate_2018 = col_date(format = "%Y-%m-%d"),
                      surveyDate_2019 = col_date(format = "%Y-%m-%d"),
                      surveyDate_2020 = col_date(format = "%Y-%m-%d"),
                      surveyDate_2021 = col_date(format = "%Y-%m-%d"),
                      botanist_2018 = "c",
                      botanist_2019 = "c",
                      botanist_2020 = "c",
                      botanist_2021 = "c",
                      vegetationCov_2018 = "c",
                      vegetationCov_2019 = "c",
                      vegetationCov_2020 = "c",
                      vegetationCov_2021 = "c",
                      bioMass_2019 = "c"
                      )) %>%
  select(-starts_with("surveyDate")) %>%
  pivot_longer(starts_with("vegetationCov_") | starts_with("botanist_") | starts_with("bioMass_"), 
               names_to = c("x", "surveyYear"),
               names_sep = "_",
               values_to = "n") %>%
  pivot_wider(names_from = x, values_from = n) %>%
  mutate(plot = str_replace(plot, "-", "_"),
         plot = str_replace(plot, "L_", "L"),
         plot = str_replace(plot, "W_", "W"),
         id = str_c(plot, surveyYear, sep = "_"),
         plot = factor(plot),
         id = factor(id),
         vegetationCov = as.numeric(vegetationCov),
         surveyYear = as.numeric(surveyYear),
         bioMass = as.numeric(bioMass)) %>%
  filter(!(block == "C" & (surveyYear == "seeded" | 
                             surveyYear == "2018" | 
                             surveyYear == "2019" | 
                             surveyYear == "2020")))


## 2 Species #####################################################################################

species <- data.table::fread("data_raw_species_20211103.csv", 
                             sep = ",",
                             dec = ".",
                             skip = 0,
                             header = T,
                             na.strings = c("", "NA", "na"),
                             colClasses = list(
                               character = "name"
                             )) %>%
  ### Check that each species occurs at least one time ###
  group_by(name) %>%
  arrange(name) %>%
  select(name, all_of(sites$id)) %>%
  mutate(total = sum(c_across(starts_with("L") | starts_with("W") | starts_with("C")), na.rm = T),
         presence = if_else(total > 0, 1, 0)) %>%
  filter(presence == 1) %>% # filter only species which occur at least one time
  ungroup() %>%
  select(name, sort(tidyselect::peek_vars()), -total, -presence) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))

### Create list with species names and their frequency ###
specieslist <- species %>%
  mutate(across(where(is.numeric), ~1 * (. != 0))) %>%
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = T), .keep = "unused") %>%
  group_by(name) %>%
  summarise(sum = sum(sum))
#write_csv(specieslist, here("outputs/tables/specieslist_20211011.csv"))


## 3 Traits #####################################################################################

traits <- read_csv("data_raw_traits.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                     cols(
                       .default = "f",
                       name = "c",
                       l = "d",
                       t = "d",
                       k = "d",
                       f = "d",
                       r = "d",
                       n = "d",
                     )) %>%
  separate(name, c("genus", "species", "ssp", "subspecies"), "_", 
           remove = F, extra = "drop", fill = "right") %>%
  mutate(genus = str_sub(genus, 1, 4),
         species = str_sub(species, 1, 4),
         subspecies = str_sub(subspecies, 1, 4),
         name = as.character(name)) %>%
  unite(abb, genus, species, subspecies, sep = "", na.rm = T) %>%
  mutate(abb = str_replace(abb, "NA", ""),
         abb = as_factor(abb)) %>%
  select(-ssp, -synonym, -nomenclature, -legal, -l, -k, -fchange) %>%
  arrange(name)

### Check congruency of traits and species table ###
  traits[duplicated(traits$abb),]
  #traits$name[which(!(traits$name %in% species$name))]
  species$name[which(!(species$name %in% traits$name))]
traits <- semi_join(traits, species, by = "name")


## 4 Check data frames #####################################################################################

### Check typos ###
sites %>%
  filter(!str_detect(id, "_seeded$")) %>%
  janitor::tabyl(vegetationCov)
#sites %>% filter(vegetationCov == 17)
species %>%
  select(-name, -ends_with("_seeded")) %>%
  unlist() %>%
  janitor::tabyl()
species %>% # Check special typos
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  filter(value == 90)

### Compare vegetationCov and accumulatedCov ###
species %>%
    summarise(across(where(is.double), ~sum(.x, na.rm = T))) %>%
    pivot_longer(cols = everything(), names_to = "id", values_to = "value") %>%
    mutate(id = factor(id)) %>%
    full_join(sites, by = "id") %>% 
    mutate(diff = (value - vegetationCov)) %>%
    select(id, surveyYear, value, vegetationCov, diff) %>%
    filter(!str_detect(id, "_seeded$")) %>%
    filter(diff > 20 | diff < -5) %>%
    arrange(surveyYear, id, diff) %>%
    print(n = 100)

### Check plots over time ###
species %>%
  select(name, starts_with("L1_19"), -ends_with("_seeded")) %>%
  filter(if_any(starts_with("L"), ~ . > 0)) %>%
  print(n = 100)

### Check missing data ###
miss_var_summary(sites, order = T)
vis_miss(sites, cluster = F)
miss_var_summary(traits, order = T)
vis_miss(traits, cluster = F, sort_miss = T)

#sites[is.na(sites$id),]
rm(list = setdiff(ls(), c("sites", "species", "traits")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


## 1 Create simple variables #####################################################################################

sites <- sites %>%
  mutate(conf.low = c(1:length(id)),
         conf.high = c(1:length(id)))

traits <- traits %>%
  mutate(target = if_else(
    targetHerb == "yes" | targetGrass == "yes", "yes", "no"
  ))


## 2 Coverages #####################################################################################

cover <- left_join(species, traits, by = "name") %>%
  select(name, family, target, ruderalIndicator, seeded, starts_with("L"), starts_with("W")) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id)

### * graminoid, herb, and total coverage) ####
cover_total_and_graminoid <- cover %>%
  group_by(id, family) %>%
  summarise(total = sum(n, na.rm = T), .groups = "keep") %>%
  mutate(type = if_else(family == "Poaceae" | family == "Cyperaceae" | family == "Juncaceae", "graminoidCov", "herbCov")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = T), .groups = "keep") %>%
  spread(type, total) %>%
  mutate(accumulatedCov = graminoidCov + herbCov,
         accumulatedCov = round(accumulatedCov, 1)) %>%
  ungroup()

### * Target species' coverage ####
cover_target <- cover %>%
  filter(target == "yes") %>%
  summarise(targetCov = sum(n, na.rm = T)) %>%
  mutate(targetCov = round(targetCov, 1)) %>%
  ungroup()

### * Ruderal species' coverage ####
cover_ruderalIndicator <- cover %>%
  filter(ruderalIndicator == "yes") %>%
  summarise(ruderalCov = sum(n, na.rm = T)) %>%
  mutate(ruderalCov = round(ruderalCov, 1)) %>%
  ungroup()

### * Seeded species' coverage ####
cover_seeded <- species %>%
  select(-starts_with("C"), -contains("x0809")) %>%
  pivot_longer(-name, names_to = "id", values_to = "value", values_drop_na = T) %>%
  separate(id, c("plot", "surveyYear"), sep = "_(?!.*_)", 
           remove = F, extra = "merge", fill = "warn", convert = F) %>%
  pivot_wider(names_from = "surveyYear", values_from = "value") %>%
  group_by(plot, name) %>%
  summarise(across(where(is.double), ~sum(.x, na.rm = T)), .groups = "keep") %>%
  ungroup() %>%
  pivot_longer(starts_with("20"), names_to = "surveyYear", values_to = "value", 
               names_transform = list(surveyYear = as.factor)) %>%
  mutate(success = if_else(seeded > 0 & value > 0, value, 0)) %>%
  group_by(plot, surveyYear) %>%
  summarise(seededCov = sum(success), .groups = "keep") %>%
  ungroup() %>%
  unite(id, plot, surveyYear, sep = "_")

### * implement in sites data set ####
sites <- sites %>%
  left_join(cover_total_and_graminoid, by = "id") %>%
  left_join(cover_target, by = "id") %>%
  #right_join(cover_ruderalIndicator, by = "id") %>%
  left_join(cover_seeded, by = "id") %>%
  ### Calcute the ratio of target species richness of total species richness
  mutate(targetCovratio = targetCov / accumulatedCov,
         graminoidCovratio = graminoidCov / accumulatedCov,
         seededCovratio = seededCov / accumulatedCov,
         targetCovratio = round(targetCovratio, 3),
         graminoidCovratio = round(graminoidCovratio, 3),
         seededCovratio = round(seededCovratio, 3))

rm(list = setdiff(ls(), c("sites", "species", "traits")))


## 3 Alpha diversity #####################################################################################

### a Species richness ----------------------------------------------------------------------------------------------------

speciesRichness <- species %>%
  left_join(traits, by = "name") %>%
  select(name, rlg, rlb, target, targetHerb, ffh6510, ffh6210, starts_with("L"), starts_with("W")) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  group_by(id)

### * total species richness ####
speciesRichness_all <- speciesRichness %>%
  summarise(speciesRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * red list Germany (species richness) ####
speciesRichness_rlg <- speciesRichness %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  summarise(rlgRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * red list Bavaria (species richness) ####
speciesRichness_rlb <- speciesRichness %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  summarise(rlbRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * target species (species richness) ####
speciesRichness_target <- speciesRichness %>%
  filter(target == "yes") %>%
  summarise(targetRichness = sum(n, na.rm = T)) %>%
  ungroup()

### * seeded species (species richness) ####
speciesRichness_seeded <- species %>%
  select(-starts_with("C"), -contains("x0809")) %>%
  pivot_longer(-name, names_to = "id", values_to = "value", values_drop_na = T) %>%
  separate(id, c("plot", "surveyYear"), sep = "_(?!.*_)", 
           remove = F, extra = "merge", fill = "warn", convert = F) %>%
  pivot_wider(names_from = "surveyYear", values_from = "value") %>%
  group_by(plot, name) %>%
  summarise(across(where(is.double), ~sum(.x, na.rm = T)), .groups = "keep") %>%
  ungroup() %>%
  pivot_longer(starts_with("20"), names_to = "surveyYear", values_to = "value", 
               names_transform = list(surveyYear = as.factor)) %>%
  mutate(successCov = if_else(seeded > 0 & value > 0, 1, 0)) %>%
  group_by(plot, surveyYear) %>%
  summarise(seededRichness = sum(successCov), .groups = "keep") %>%
  ungroup() %>%
  unite(id, plot, surveyYear, sep = "_")

### * ffh6510 species (species richness) ####
speciesRichness_ffh6510 <- speciesRichness %>%
  filter(ffh6510 == "yes") %>%
  summarise(ffh6510Richness = sum(n, na.rm = T)) %>%
  ungroup()

### * ffh6210 species (species richness) ####
speciesRichness_ffh6210 <- speciesRichness %>%
  filter(ffh6210 == "yes") %>%
  summarise(ffh6210Richness = sum(n, na.rm = T)) %>%
  ungroup()

### * implement in sites data set ####
sites <- sites %>%
  left_join(speciesRichness_all, by = "id") %>%
  left_join(speciesRichness_rlg, by = "id") %>%
  left_join(speciesRichness_rlb, by = "id") %>%
  left_join(speciesRichness_target, by = "id") %>%
  left_join(speciesRichness_seeded, by = "id") %>%
  left_join(speciesRichness_ffh6510, by = "id") %>%
  left_join(speciesRichness_ffh6210, by = "id") %>%
  ### Create targetRichratio ###
  mutate(targetRichratio = targetRichness / speciesRichness,
         seededRichratio = seededRichness / speciesRichness,
         targetRichratio = round(targetRichratio, 3),
         seededRichratio = round(seededRichratio, 3))

### b Species eveness and shannon ----------------------------------------------------------------------
data <- species  %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id") %>%
  diversity(index = "shannon") %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "id") %>%
  mutate(id = factor(id)) %>%
  rename(shannon = value)
sites <- sites %>%
  left_join(data, by = "id") %>%
  mutate(eveness = shannon / log(speciesRichness))

rm(list = ls(pattern = "[^species|traits|sites]"))


## 5 CWM of Ellenberg #####################################################################################

### a N value -------------------------------------------------------------------------------------------
data_traits <- traits %>%
  select(name, n) %>%
  filter(n > 0)
data_species <- semi_join(species, data_traits, by = "name") %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = T)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id")
data_traits <- column_to_rownames(data_traits, "name")
### Calculate CWM ###
Nweighted <- dbFD(data_traits, data_species, w.abun = T,
                        calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### b F value -------------------------------------------------------------------------------------------
data_traits <- traits %>%
  select(name, f) %>%
  filter(f > 0) 
data_species <- semi_join(species, data_traits, by = "name") %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = T)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id")
data_traits <- column_to_rownames(data_traits, "name")
### Calculate CWM ###
Fweighted <- dbFD(data_traits, data_species, w.abun = T,
                  calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### c T value -------------------------------------------------------------------------------------------
data_traits <- traits %>%
  select(name, t) %>%
  filter(t > 0) 
data_species <- semi_join(species, data_traits, by = "name") %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = T)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id")
data_traits <- column_to_rownames(data_traits, "name")
### Calculate CWM ###
Tweighted <- dbFD(data_traits, data_species, w.abun = T,
                  calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### d implement in sites data set -------------------------------------------------------------------------------------------
sites$cwmAbuN <- round(as.numeric(as.character(Nweighted$CWM$n)), 3)
sites$cwmAbuF <- round(as.numeric(as.character(Fweighted$CWM$f)), 3)
sites$cwmAbuT <- round(as.numeric(as.character(Tweighted$CWM$t)), 3)

rm(list = setdiff(ls(), c("sites", "species", "traits")))


## 6 Functional plant traits #####################################################################################

### a LEDA data #####
data_sla <- data.table::fread("data_raw_traitbase_leda_20210223_sla.txt", 
                             sep = ";",
                             dec = ".",
                             skip = 4,
                             header = T,
                             select = c("SBS name", "single value [mm^2/mg]")
                             ) %>%
  as_tibble() %>%
  rename(name = "SBS name", 
         sla = "single value [mm^2/mg]") %>%
  mutate(name = str_replace_all(name, " ", "_"))

data_seedmass <- data.table::fread("data_raw_traitbase_leda_20210223_seedmass.txt", 
                            sep = ";",
                            dec = ".",
                            skip = 3,
                            header = T,
                            select = c("SBS name", "single value [mg]"),
                            ) %>%
  as_tibble() %>%
  rename(name = "SBS name",
         seedmass = "single value [mg]") %>%
  mutate(name = str_replace_all(name, " ", "_"),
         seedmass = round(seedmass, 3))

data_height <- data.table::fread("data_raw_traitbase_leda_20210223_canopy_height.txt", 
                           sep = ";",
                           dec = ".",
                           skip = 3,
                           header = T,
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
                           Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum",
                           Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum_s._vulgare",
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
                           Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia",
                           Rubus_fruticosus_agg = "Rubus_fruticosus",
                           Rubus_fruticosus_agg = "Rubus_fruticosus_ag._L.",
                           Securigera_varia = "Coronilla_varia",
                           Sedum_maximum = "Sedum_telephium_s._maximum",
                           "Silene_flos-cuculi" = "Lychnis_flos-cuculi",
                           Silene_latifolia_ssp_alba = "Silene_latifolia",
                           Taraxacum_campylodes = "Taraxacum_officinale",
                           Taraxacum_campylodes = "Taraxacum_Sec._Ruderalia",
                           Tripleurospermum_maritimum = "Matricaria_maritima",
                           Vicia_sativa_ssp_nigra= "Vicia_sativa_s._nigra")) %>% 
  group_by(name) %>%
### summarise different rows of one species after renaming ###  
  summarise(across(where(is.double), ~median(.x, na.rm = T)))
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

rm(list = setdiff(ls(), c("sites", "species", "traits", "test")))


### b TRY data 1 ####
data <- data.table::fread("data_raw_traitbase_try_20210306_13996.txt", 
              header = T, 
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
                           Carex_praecox_ssp_praecox = "Carex_praecox",
                           Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia",
                           Ranunculus_serpens_ssp_nemorosus = "Ranunculus_serpens_subsp._nemorosus"),
         trait = as_factor(trait)) %>%
  group_by(name, trait) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
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

### c TRY data 2 ####
data <- data.table::fread("data_raw_traitbase_try_20210318_14157.txt", 
                          header = T, 
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
                           Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia",
                           Plantago_major = "Plantago_major_subsp._major",
                           Cornus_controversa = "Cornus_sanguinea",
                           Silene_latifolia_ssp_alba = "Silene_latifolia_subsp._alba",
                           Silene_latifolia_ssp_alba = "Silene_latifolia",
                           Silene_latifolia_ssp_alba = "Silene_latifolia_ssp._alba"),
         trait = as_factor(trait)) %>%
  group_by(name, trait) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
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

### d GrooT data ####
data <- read.csv("data_raw_traitbase_groot.csv", header = T, na.strings = c("", "NA", "na")) %>%
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
                           Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum",
                           Persicaria_bistorta = "Bistorta_officinalis",
                           Jacobaea_vulgaris = "Senecio_jacobaea",
                           "Silene_flos-cuculi" = "Lychnis_flos-cuculi",
                           Silene_latifolia_ssp_alba = "Silene_latifolia",
                           Vicia_sativa_ssp_nigra = "Vicia_sativa")) %>%
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T)))
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

### e check completeness of traits ####
test <- traits %>%
  select(name, group, t, n, f, sla, seedmass, height, rmf, rld, srl, lateral) %>%
  filter(!(str_detect(name, "_spec") | str_detect(name, "ceae"))) %>%
  filter(group != "tree" & group != "shrub") %>%
  select(-name, -group)
vis_miss(test, cluster = T, sort_miss = T)
miss_var_summary(test)
gg_miss_case(test, order_cases = T)
(test <- traits %>%
    select(name, sla, seedmass, height, rmf, rld, srl, lateral) %>%
    filter(!complete.cases(sla, seedmass, height, rmf, rld, srl, lateral)))

rm(list = ls(pattern = "[^species|traits|sites]"))

### f prepare data frames ####
species <- species %>%
  mutate(name = as.character(name)) %>%
  arrange(name)
traits <- traits %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(across(c(sla, seedmass, height, srl, rmf, rld), 
                ~ round(.x, digits = 3)))
(herbCount <- traits %>% #gamma diversity herbs: 269
    filter(group != "tree" & group != "shrub") %>%
    count() %>%
    pull())
(undefinedCount <- traits %>% #gamma diversity of undefined species: 15
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


## 7 CWM and FDis of functional plant traits #####################################################################################

### a LHS -------------------------------------------------------------------------------------------
data_species <- semi_join(species, traits_lhs, by = "name")
data_traits <- semi_join(traits_lhs, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = T,
               calc.FRic = F, calc.FDiv = F, corr = "cailliez")
sites$fdisAbuLHS <- data_diversityAbu$FDis
length(traits_lhs$name) / (herbCount - undefinedCount)

### b SLA -------------------------------------------------------------------------------------------
data_species <- semi_join(species, traits_sla, by = "name")
data_traits <- semi_join(traits_sla, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = T, 
                   calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSla <- data_diversityAbu$FDis
sites$cwmAbuSla <- data_diversityAbu$CWM$sla %>%
  as.character() %>% 
  as.numeric() %>% 
  exp()
length(traits_sla$name) / (herbCount - undefinedCount)

### c Seed mass -------------------------------------------------------------------------------------------
data_species <- semi_join(species, traits_seedmass, by = "name")
data_traits <- semi_join(traits_seedmass, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSeedmass <- data_diversityAbu$FDis
sites$cwmAbuSeedmass <- data_diversityAbu$CWM$seedmass %>%
  as.character() %>% 
  as.numeric() %>% 
  exp()
length(traits_seedmass$name) / (herbCount - undefinedCount)

### d Canopy height -------------------------------------------------------------------------------------------
data_species <- semi_join(species, traits_height, by = "name")
data_traits <- semi_join(traits_height, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuHeight <- data_diversityAbu$FDis
sites$cwmAbuHeight <- data_diversityAbu$CWM$height %>%
  as.character() %>% 
  as.numeric() %>% 
  exp()
length(traits_height$name) / (herbCount - undefinedCount)

### e Specific root length -------------------------------------------------------------------------------------------
data_species <- semi_join(species, traits_srl, by = "name")
data_traits <- semi_join(traits_srl, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSrl <- data_diversityAbu$FDis
sites$cwmAbuSrl <- data_diversityAbu$CWM$srl %>%
  as.character() %>% 
  as.numeric() %>% 
  exp()
length(traits_srl$name) / (herbCount - undefinedCount)

### f Root mass fraction -------------------------------------------------------------------------------------------
data_species <- semi_join(species, traits_rmf, by = "name")
data_traits <- semi_join(traits_rmf, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuRmf <- data_diversityAbu$FDis
sites$cwmAbuRmf <- data_diversityAbu$CWM$rmf %>%
  as.character() %>% 
  as.numeric() %>% 
  exp()
length(traits_rmf$name) / (herbCount - undefinedCount)

### g All -------------------------------------------------------------------------------------------
data_species <- semi_join(species, traits_all, by = "name")
data_traits <- semi_join(traits_all, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
data_diversityAbu <- dbFD(log_data_traits, data_species, w.abun = T,
                      calc.FRic = F, calc.FDiv = F, corr = "cailliez")
sites$fdisAbuAll <- data_diversityAbu$FDis
length(traits_all$name) / (herbCount - undefinedCount)

### h Finalisation -------------------------------------------------------------------------------------------
sites <- sites %>%
  mutate(across(c(fdisAbuLHS, fdisAbuSla, cwmAbuSla, fdisAbuSeedmass, cwmAbuSeedmass, fdisAbuHeight, cwmAbuHeight, fdisAbuSrl, cwmAbuSrl, fdisAbuRmf, cwmAbuRmf, fdisAbuAll), 
                ~ round(.x, digits = 3)))

rm(list = setdiff(ls(), c("sites", "species", "traits")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


write_csv(sites, here("data/processed/data_processed_sites.csv"))
write_csv(traits, here("data/processed/data_processed_traits.csv"))
write_csv(species, here("data/processed/data_processed_species.csv"))

