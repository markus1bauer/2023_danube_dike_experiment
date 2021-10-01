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


### 1 Sites #####################################################################################

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


### 2 Species #####################################################################################

species <- data.table::fread("data_raw_species.csv", 
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
  mutate(total = sum(c_across(starts_with("L") | starts_with("W") | starts_with("C")), na.rm = T),
         presence = if_else(total > 0, 1, 0),
         name = factor(name)) %>%
  filter(presence == 1) %>% # filter species
  ungroup() %>%
  select(name, sort(tidyselect::peek_vars()), -total, -presence) %>%
  na_if(0)

### Create list with species names and their frequency ###
specieslist <- species %>%
  mutate(across(where(is.numeric), ~1 * (. != 0))) %>%
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = T), .keep = "unused") %>%
  group_by(name) %>%
  summarise(sum = sum(sum))
#write_csv(specieslist, here("data/raw/specieslist_2022xxxx.csv"))


### 3 Traits #####################################################################################

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
         name = factor(name)) %>%
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


### 4 Check data frames #####################################################################################

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
    summarise(across(where(is.double),~sum(.x, na.rm = T))) %>%
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
  select(name, starts_with("W5_24")) %>%
  filter(if_any(starts_with("W"), ~ . > 0)) %>%
  print(n = 100)

### Check missing data ###
miss_var_summary(sites, order = T)
vis_miss(sites, cluster = F)
miss_var_summary(traits, order = T)
vis_miss(traits, cluster = F, sort_miss = T)

#sites[is.na(sites$id),]
rm(list=setdiff(ls(), c("sites", "species", "traits")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Create simple variables #####################################################################################

sites <- sites %>%
  mutate(conf.low = c(1:length(id)),
         conf.high = c(1:length(id)))

traits <- traits %>%
  mutate(target = if_else(
    targetHerb == "yes" | targetGrass == "yes", "yes", "no"
  ))


### 2 Coverages #####################################################################################

cover <- left_join(species, traits, by = "name") %>%
  select(name, family, target, ruderalIndicator, seeded, starts_with("L"), starts_with("W")) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id)

### * graminoid, herb, and total coverage) ####
cover_total_and_graminoid <- cover %>%
  group_by(id, family) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  mutate(type = if_else(family == "Poaceae" | family == "Cyperaceae" | family == "Juncaceae", "graminoidCov", "herbCov")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = T)) %>%
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
  summarise(across(where(is.double), ~sum(.x, na.rm = T))) %>%
  ungroup() %>%
  pivot_longer(starts_with("20"), names_to = "surveyYear", values_to = "value", 
               names_transform = list(surveyYear = as.factor)) %>%
  mutate(success = if_else(seeded > 0 & value > 0, 1, 0)) %>%
  group_by(plot, surveyYear) %>%
  summarise(seededRichness = sum(success)) %>%
  ungroup() %>%
  unite(id, plot, surveyYear, sep = "_")

### * implement in sites data set ####
sites <- sites %>%
  full_join(cover_total_and_graminoid, by = "id") %>%
  full_join(cover_target, by = "id") %>%
  #right_join(cover_ruderalIndicator, by = "id") %>%
  full_join(cover_seeded, by = "id") %>%
  ### Calcute the ratio of target species richness of total species richness
  mutate(targetCovratio = targetCov / accumulatedCov,
         graminoidCovratio = graminoidCov / accumulatedCov,
         seededCovratio = seededCov / accumulatedCov,
         targetCovratio = round(targetCovratio, 3),
         graminoidCovratio = round(graminoidCovratio, 3),
         seededCovratio = round(seededCovratio, 3))

rm(list = setdiff(ls(), c("sites", "species", "traits")))

### 3 Species richness #####################################################################################

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
  summarise(across(where(is.double), ~sum(.x, na.rm = T))) %>%
  ungroup() %>%
  pivot_longer(starts_with("20"), names_to = "surveyYear", values_to = "value", 
               names_transform = list(surveyYear = as.factor)) %>%
  mutate(successCov = if_else(seeded > 0 & value > 0, value, 0)) %>%
  group_by(plot, surveyYear) %>%
  summarise(seededCov = sum(successCov)) %>%
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
  right_join(speciesRichness_all, by = "id") %>%
  right_join(speciesRichness_rlg, by = "id") %>%
  right_join(speciesRichness_rlb, by = "id") %>%
  right_join(speciesRichness_target, by = "id") %>%
  right_join(speciesRichness_seeded, by = "id") %>%
  right_join(speciesRichness_ffh6510, by = "id") %>%
  right_join(speciesRichness_ffh6210, by = "id") %>%
  ### Create targetRichratio ###
  mutate(targetRichratio = targetRichness / speciesRichness,
         seededRichratio = seededRichness / speciesRichness,
         targetRichratio = round(targetRichratio, 3),
         seededRichratio = round(seededRichratio, 3))

rm(list = setdiff(ls(), c("sites", "species", "traits")))


### 5 TBI #####################################################################################

tbi <- species %>%
  pivot_longer(-name, "id", "n") %>%
  pivot_wider(id, name) %>%
  left_join(sites, by = "id") %>%
  select(id, exposition, surveyYear, Acer_campestre:Vulpia_myuros) %>%
  filter(surveyYear == 2018 | surveyYear == 2020) %>%
#tbiStart <- sites %>% select(id, plot, surveyYear) %>% filter(surveyYear == 2018)
#tbiEnd <- sites %>% select(id, plot, surveyYear) %>% filter(surveyYear == 2020)
#anti_join(tbiStart, tbiEnd, by = "plot")
#anti_join(tbiEnd, tbiStart, by = "plot")
#rm(tbiStart, tbiEnd)
  column_to_rownames(var = "id") %>%
  select(-surveyYear)

### a Abundance data --------------------------------------------------------------------------------------------
tbi2 <- select(tbi, -exposition)
tbiAbu <- tbi2[,colSums(tbi2) > 0]

### b Presence-absence data--------------------------------------------------------------------------------------------
tbiPa <- tbiAbu
tbiPa[tbiPa > 0] = 1

### c Abundance Exposition--------------------------------------------------------------------------------------------
#tbi2 <- tbi %>%
  #filter(exposition == "north") %>%
 # select(-exposition)
#tbiAbuN <- tbi2[,colSums(tbi2) > 0]

tbiAbu <- rownames_to_column(tbiAbu, var = "id")
tbiPa <- rownames_to_column(tbiPa, var = "id")
#tbiAbuN <- rownames_to_column(tbiAbuN, var = "id")
rm(tbi, tbi2)


### 6 CWM of Ellenberg #####################################################################################

### a N value -------------------------------------------------------------------------------------------
Ntraits <- traits %>%
  select(name, n) %>%
  filter(n > 0)
Nspecies <- semi_join(species, Ntraits, by = "name") %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = T)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id")
Ntraits <- column_to_rownames(Ntraits, "name")
### Calculate CWM ###
Nweighted <- dbFD(Ntraits, Nspecies, w.abun = T,
                        calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### b F value -------------------------------------------------------------------------------------------
Ftraits <- traits %>%
  select(name, f) %>%
  filter(f > 0) 
Fspecies <- semi_join(species, Ftraits, by = "name") %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = T)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id")
Ftraits <- column_to_rownames(Ftraits, "name")
### Calculate CWM ###
Fweighted <- dbFD(Ftraits, Fspecies, w.abun = T,
                  calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### c T value -------------------------------------------------------------------------------------------
Ttraits <- traits %>%
  select(name, t) %>%
  filter(t > 0) 
Tspecies <- semi_join(species, Ttraits, by = "name") %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  group_by(id) %>%
  mutate(sum = sum(value, na.rm = T)) %>%
  filter(sum > 0) %>%
  select(-sum) %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id")
Ttraits <- column_to_rownames(Ttraits, "name")
### Calculate CWM ###
Tweighted <- dbFD(Ttraits, Tspecies, w.abun = T,
                  calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### d implement in sites data set -------------------------------------------------------------------------------------------
sites$cwmAbuN <- round(as.numeric(as.character(Nweighted$CWM$n)), 3)
sites$cwmAbuF <- round(as.numeric(as.character(Fweighted$CWM$f)), 3)
sites$cwmAbuT <- round(as.numeric(as.character(Tweighted$CWM$t)), 3)
rm(list=setdiff(ls(), c("sites", "species", "traits", "tbiPa", "tbiAbu")))


### 7 Functional plant traits #####################################################################################

#### * read and join LEDA data #### 
dataSLA <- data.table::fread("data_raw_LEDA_20210223_sla.txt", 
                             sep = ";",
                             dec = ".",
                             skip = 3,
                             header = T,
                             select = c("SBS name", "single value [mm^2/mg]"),
                             ) %>%
  as_tibble() %>%
  rename(name = "SBS name") %>%
  rename(sla = "single value [mm^2/mg]") %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_")))

dataSM <- data.table::fread("data_raw_LEDA_20210223_seedmass.txt", 
                            sep = ";",
                            dec = ".",
                            skip = 3,
                            header = T,
                            select = c("SBS name", "single value [mg]"),
                            ) %>%
  as_tibble() %>%
  rename(name = "SBS name") %>%
  rename(seedmass = "single value [mg]") %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_")))

dataH <- data.table::fread("data_raw_LEDA_20210223_canopy_height.txt", 
                           sep = ";",
                           dec = ".",
                           skip = 3,
                           header = T,
                           select = c("SBS name", "single value [m]")
                             ) %>%
  as_tibble() %>%
  rename(name = "SBS name") %>%
  rename(height = "single value [m]") %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_")))
### Join SLA, canopy height and seedmass of LEDA ###
data <- dataSLA %>%
  full_join(dataSM, by = "name") %>%
  full_join(dataH, by = "name")
rm(dataSLA, dataSM, dataH)
##### Find synonyms ###
#traits$name[which(!(traits$name %in% data$name))]
#data %>%
  #group_by(name) %>%
  #summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  #filter(str_detect(name, "incana"))
traits <- data %>%
  mutate(name = fct_recode(name, Centaurea_stoebe = "Centaurea_stoebe_s.lat.")) %>%
  mutate(name = fct_recode(name, Carex_praecox_ssp_praecox = "Carex_praecox")) %>%
  mutate(name = fct_recode(name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum")) %>%
  mutate(name = fct_recode(name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum_s._vulgare")) %>%
  mutate(name = fct_recode(name, Cota_tinctoria = "Anthemis_tinctoria")) %>%
  mutate(name = fct_recode(name, Cyanus_segetum = "Centaurea_cyanus")) %>%
  mutate(name = fct_recode(name, Euphorbia_verrucosa = "Euphorbia_brittingeri")) %>%
  mutate(name = fct_recode(name, Helictotrichon_pubescens = "Avenula_pubescens")) %>%
  mutate(name = fct_recode(name, Hypochaeris_radicata = "Hypochoeris_radicata")) %>%
  mutate(name = fct_recode(name, Jacobaea_vulgaris = "Senecio_vulgaris")) %>%
  mutate(name = fct_recode(name, Medicago_falcata = "Medicago_sativa_s._falcata")) %>%
  mutate(name = fct_recode(name, Ononis_spinosa_ssp_procurrens = "Ononis_repens")) %>%
  mutate(name = fct_recode(name, Persicaria_amphibia = "Polygonum_amphibium")) %>%
  mutate(name = fct_recode(name, Pilosella_caespitosa = "Hieracium_caespitosum")) %>%
  mutate(name = fct_recode(name, Pilosella_officinarum = "Hieracium_pilosella")) %>%
  mutate(name = fct_recode(name, Pilosella_piloselloides = "Hieracium_piloselloides")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia")) %>%
  mutate(name = fct_recode(name, Rubus_fruticosus_agg = "Rubus_fruticosus")) %>%
  mutate(name = fct_recode(name, Rubus_fruticosus_agg = "Rubus_fruticosus_ag._L.")) %>%
  mutate(name = fct_recode(name, Securigera_varia = "Coronilla_varia")) %>%
  mutate(name = fct_recode(name, "Silene_flos-cuculi" = "Lychnis_flos-cuculi")) %>%
  mutate(name = fct_recode(name, Taraxacum_campylodes = "Taraxacum_officinale")) %>%
  mutate(name = fct_recode(name, Taraxacum_campylodes = "Taraxacum_Sec._Ruderalia")) %>%
  mutate(name = fct_recode(name, Tripleurospermum_maritimum = "Matricaria_maritima")) %>%
  mutate(name = fct_recode(name, Vicia_sativa_ssp_nigra= "Vicia_sativa_s._nigra")) %>%  group_by(name) %>%
### summarise different rows of one species after renaming ###  
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  right_join(traits, by = "name")
### check completeness of LEDA ###
#(test2 <- traits %>%
  #select(name, sla, seedmass, height) %>%
  #filter(!complete.cases(.)))

### * read TRY data 1 ####
data <- data.table::fread("data_raw_TRY_20210306_13996.txt", 
              header = T, 
              sep = "\t", 
              dec = ".", 
              quote = "") %>%
  as_tibble() %>%
  rename(name = "AccSpeciesName") %>%
  rename(value = "StdValue") %>%
  rename(trait = "TraitID") %>%
  select(name, value, trait) %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_"))) %>%
  drop_na %>%
  mutate(trait = str_replace(trait, "26", "seedmass")) %>%
  mutate(trait = str_replace(trait, "3106", "height")) %>%
  mutate(trait = str_replace(trait, "3107", "height")) %>%
  mutate(trait = str_replace(trait, "3115", "sla")) %>%
  mutate(trait = str_replace(trait, "3116", "sla")) %>%
  mutate(trait = str_replace(trait, "3117", "sla"))
##### Find synonyms ###
#test2$name[which(!(test2$name %in% data$name))]
#data %>%
  #group_by(name) %>%
  #summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  #filter(str_detect(name, "Equisetum"))
traits <- data %>%
  mutate(name = fct_recode(name, Carex_praecox_ssp_praecox = "Carex_praecox")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia")) %>%
  mutate(name = fct_recode(name, Ranunculus_serpens_ssp_nemorosus = "Ranunculus_serpens_subsp._nemorosus")) %>%
  mutate(trait = as_factor(trait)) %>%
  group_by(name, trait) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  pivot_wider(names_from = "trait", values_from = "value") %>%
  right_join(traits, by = "name") %>%
  mutate(sla = coalesce(sla.x, sla.y), .keep = "unused") %>%
  mutate(seedmass = coalesce(seedmass.x, seedmass.y), .keep = "unused") %>%
  mutate(height = coalesce(height.x, height.y), .keep = "unused")
### check completeness of LEDA + TRY1 ###
#(test2 <- traits %>%
#select(name, sla, seedmass, height) %>%
#filter(!complete.cases(sla, seedmass, height)))

### * read TRY data 2 ####
data <- data.table::fread("data_raw_TRY_20210318_14157.txt", 
                          header = T, 
                          sep = "\t", 
                          dec = ".", 
                          quote = "") %>%
  as_tibble() %>%
  rename(name = "AccSpeciesName") %>%
  rename(value = "StdValue") %>%
  rename(trait = "TraitID") %>%
  select(name, value, trait) %>%
  mutate(name = as_factor(str_replace_all(name, " ", "_"))) %>%
  drop_na %>%
  mutate(trait = str_replace(trait, "26", "seedmass")) %>%
  mutate(trait = str_replace(trait, "3106", "height")) %>%
  mutate(trait = str_replace(trait, "3107", "height")) %>%
  mutate(trait = str_replace(trait, "3115", "sla")) %>%
  mutate(trait = str_replace(trait, "3116", "sla")) %>%
  mutate(trait = str_replace(trait, "3117", "sla"))
##### Find Synonyms ###
#test2$name[which(!(test2$name %in% data$name))]
#data %>%
#group_by(name) %>%
#summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
#filter(str_detect(name, "Silene"))
traits <- data %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_major = "Plantago_major_subsp._major")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_major = "Plantago_major")) %>%
  mutate(name = fct_recode(name, Cornus_controversa = "Cornus_sanguinea")) %>%
  mutate(name = fct_recode(name, Silene_latifolia_ssp_alba = "Silene_latifolia_subsp._alba")) %>%
  mutate(name = fct_recode(name, Silene_latifolia_ssp_alba = "Silene_latifolia")) %>%
  mutate(name = fct_recode(name, Silene_latifolia_ssp_alba = "Silene_latifolia_ssp._alba")) %>%
  mutate(trait = as_factor(trait)) %>%
  group_by(name, trait) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  pivot_wider(names_from = "trait", values_from = "value") %>%
  right_join(traits, by = "name") %>%
  mutate(sla = coalesce(sla.x, sla.y), .keep = "unused") %>%
  mutate(seedmass = coalesce(seedmass.x, seedmass.y), .keep = "unused") %>%
  mutate(height = coalesce(height.x, height.y), .keep = "unused") 
### check completeness of LEDA + TRY1 + TRY2 ###
#(test2 <- traits %>%
#select(name, sla, seedmass, height) %>%
#filter(!complete.cases(sla, seedmass, height)))

### * read GrooT data ####
data <- read.csv("data_raw_GrooT.csv", header = T, na.strings = c("", "NA")) %>%
  filter(traitName == "Specific_root_length" |
           traitName == "Root_length_density_volume" |
           traitName == "Root_mass_fraction" |
           traitName == "Lateral_spread") %>%
  mutate(name = as_factor(str_c(genusTNRS, speciesTNRS, sep = "_"))) %>%
  rename(trait = traitName, value = medianSpecies) %>%
  select(name, trait, value) %>%
  pivot_wider(names_from = "trait", values_from = "value") %>%
  rename(rmf = Root_mass_fraction, rld = Root_length_density_volume, srl = Specific_root_length, lateral = Lateral_spread)
##### Find Synonyms ###
#traits$name[which(!(traits$name %in% data$name))]
#data %>%
  #group_by(name) %>%
  #summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  #filter(str_detect(name, "Centaurea_"))
traits <- data %>%
  mutate(name = fct_recode(name, Cota_tinctoria = "Anthemis_tinctoria")) %>%
  mutate(name = fct_recode(name, Carex_praecox_ssp_praecox = "Carex_praecox")) %>%
  mutate(name = fct_recode(name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum")) %>%
  mutate(name = fct_recode(name, Persicaria_bistorta = "Bistorta_officinalis")) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_major = "Plantago_major")) %>%
  mutate(name = fct_recode(name, Jacobaea_vulgaris = "Senecio_jacobaea")) %>%
  mutate(name = fct_recode(name, "Silene_flos-cuculi" = "Lychnis_flos-cuculi")) %>%
  mutate(name = fct_recode(name, Silene_latifolia_ssp_alba = "Silene_latifolia")) %>%
  mutate(name = fct_recode(name, Vicia_sativa_ssp_nigra = "Vicia_sativa")) %>%
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = T))) %>%
  right_join(traits, by = "name")
### check completeness of Roots ###
(test2 <- traits %>%
   select(name, rmf, rld, srl, lateral) %>%
   filter(!complete.cases(rmf, rld, srl, lateral)))

### * check completeness of traits ####
test <- traits %>%
  select(t, n, f, sla, seedmass, height, rmf, rld, srl, lateral)
vis_miss(test, cluster = F, sort_miss = T)
#gg_miss_var(test)
gg_miss_case(test, order_cases = F)
(test2 <- traits %>%
  select(name, sla, seedmass, height, rmf, rld, srl, lateral) %>%
  filter(!complete.cases(sla, seedmass, height, rmf, rld, srl, lateral)))
rm(test, test2, data)

### * prepare data frames ####
species <- species %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(name = as_factor(name))
traits <- traits %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(name = as_factor(name))
herbCount <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  left_join(species, by = "name") %>%
  count() %>%
  pull()
undefinedSpeciesCount <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  filter(str_detect(name, "_spec|Öhrchen")) %>%
  left_join(species, by = "name") %>%
  count() %>%
  pull()
traitsLHS <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, sla, seedmass, height) %>%
  drop_na()
traitsSLA <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, sla) %>%
  drop_na()  
traitsSM <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, seedmass) %>%
  drop_na()
traitsH <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, height) %>%
  drop_na()
traitsSRL <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, srl) %>%
  drop_na()
traitsRLD <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, rld) %>%
  drop_na()
traitsRMF <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, rmf) %>%
  drop_na()
traitsAll <- traits %>%
  filter(group != "tree" & group != "shrub") %>%
  select(name, sla, seedmass, height, srl) %>%
  drop_na()


### 8 CWM and FDis of functional plant traits #####################################################################################

### a LHS -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsLHS, by = "name")
Ttraits <- semi_join(traitsLHS, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T,
               calc.FRic = F, calc.FDiv = F, corr = "cailliez")
sites$fdisAbuLHS <- TdiversityAbu$FDis
length(traitsLHS$name) / (herbCount - undefinedSpeciesCount)

### b SLA -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsSLA, by = "name")
Ttraits <- semi_join(traitsSLA, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                   calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSla <- TdiversityAbu$FDis
sites$cwmAbuSla <- TdiversityAbu$CWM$sla %>%
  as.character() %>% as.numeric() %>% exp()
length(traitsSLA$name) / (herbCount - undefinedSpeciesCount)

### c Seed mass -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsSM, by = "name")
Ttraits <- semi_join(traitsSM, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSeedmass <- TdiversityAbu$FDis
sites$cwmAbuSeedmass <- TdiversityAbu$CWM$seedmass %>%
  as.character() %>% as.numeric() %>% exp()
length(traitsSM$name) / (herbCount - undefinedSpeciesCount)

### d Canopy height -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsH, by = "name")
Ttraits <- semi_join(traitsH, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuHeight <- TdiversityAbu$FDis
sites$cwmAbuHeight <- TdiversityAbu$CWM$height %>%
  as.character() %>% as.numeric() %>% exp()
length(traitsH$name) / (herbCount - undefinedSpeciesCount)

### e Specific root length -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsSRL, by = "name")
Ttraits <- semi_join(traitsSRL, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuSrl <- TdiversityAbu$FDis
sites$cwmAbuSrl <- TdiversityAbu$CWM$srl %>%
  as.character() %>% as.numeric() %>% exp()
length(traitsSRL$name) / (herbCount - undefinedSpeciesCount)

### f Root mass fraction -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsRMF, by = "name")
Ttraits <- semi_join(traitsRMF, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T, 
                      calc.FRic = F, calc.FDiv = F, corr = "sqrt");
sites$fdisAbuRmf <- TdiversityAbu$FDis
sites$cwmAbuRmf <- TdiversityAbu$CWM$rmf %>%
  as.character() %>% as.numeric() %>% exp()
length(traitsRMF$name) / (herbCount - undefinedSpeciesCount)

### g All -------------------------------------------------------------------------------------------
Tspecies <- semi_join(species, traitsAll, by = "name")
Ttraits <- semi_join(traitsAll, Tspecies, by = "name")
Tspecies <- Tspecies %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ttraits <- column_to_rownames(Ttraits, "name")
log_Ttraits <- log(Ttraits)
TdiversityAbu <- dbFD(log_Ttraits, Tspecies, w.abun = T,
                      calc.FRic = F, calc.FDiv = F, corr = "cailliez")
sites$fdisAbuAll <- TdiversityAbu$FDis
length(traitsAll$name) / (herbCount - undefinedSpeciesCount)
rm(list=setdiff(ls(), c("sites", "species", "traits", "tbiPa", "tbiAbu")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


write_csv(sites, here("data/processed/data_processed_sites.csv"))
write_csv(species, here("data/processed/data_processed_species.csv"))
write_csv(traits, here("data/processed/data_processed_traits.csv"))
write_csv(tbiAbu, here("data/processed/data_processed_tbiAbu.csv"))
write_csv(tbiPa, here("data/processed/data_processed_tbiPa.csv"))
