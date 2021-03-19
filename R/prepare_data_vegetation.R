# Prepare vegetation data for Danube experiment ####
# Markus Bauer


### Packages ###
library(tidyverse)
library(vegan)
library(FD) #dbFD()
library(naniar) #are_na()

### Start ###
#installr::updateR(browse_news = F, install_R = T, copy_packages = T, copy_Rprofile.site = T, keep_old_packages = T, update_packages = T, start_new_R = F, quit_R = T, print_R_versions = T, GUI = F)
rm(list = ls())
setwd("Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/raw")



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Species #####################################################################################

speciesS <- read_csv2("data_raw_species_seeded.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                          cols(
                            .default = col_double(),
                            side = col_factor(),
                            id = col_factor()
                          )) #%>%
  mutate(name = str_replace(name, " ", "_")) %>%
  rename_with( ~ paste0("L_", .x), -name) %>%
  rename_with(~ paste0(.x, "_2018"), -name) %>%
  separate(name, c("name", "descriptor"), " ", extra = "drop", fill = "right") %>%
  mutate(name = as_factor(str_replace(name, "\\.", ""))) %>%
  select(-descriptor)

### a Load species tables of all years -------------------------------------------------------------------------------------------
speciesL18 <- read_csv2("data_raw_species_2018_land.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                       cols(
                         .default = col_double(),
                         name = col_factor()
                       )) %>%
  mutate(name = str_replace(name, " ", "_")) %>%
  rename_with( ~ paste0("L_", .x), -name) %>%
  rename_with(~ paste0(.x, "_2018"), -name) %>%
  separate(name, c("name", "descriptor"), " ", extra = "drop", fill = "right") %>%
  mutate(name = as_factor(str_replace(name, "\\.", ""))) %>%
  select(-descriptor)
speciesW18 <- read_csv2("data_raw_species_2018_water.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                          cols(
                            .default = col_double(),
                            name = col_factor()
                          )) %>%
  mutate(name = str_replace(name, " ", "_")) %>%
  rename_with( ~ paste0("W_", .x), -name) %>%
  rename_with(~ paste0(.x, "_2018"), -name) %>%
  separate(name, c("name", "descriptor"), " ", extra = "drop", fill = "right") %>%
  mutate(name = as_factor(str_replace(name, "\\.", ""))) %>%
  select(-descriptor) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_intermedia = "Plantago_major"))
speciesL19 <- read_csv2("data_raw_species_2019_land.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                          cols(
                            .default = col_double(),
                            name = col_factor()
                          )) %>%
  mutate(name = str_replace(name, " ", "_")) %>%
  rename_with( ~ paste0("L_", .x), -name) %>%
  rename_with(~ paste0(.x, "_2019"), -name) %>%
  separate(name, c("name", "descriptor"), " ", extra = "drop", fill = "right") %>%
  mutate(name = as_factor(str_replace(name, "\\.", ""))) %>%
  select(-descriptor) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_major = "Plantago_major"))
speciesW19 <- read_csv2("data_raw_species_2019_water.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                          cols(
                            .default = col_double(),
                            name = col_factor()
                          )) %>%
  mutate(name = str_replace(name, " ", "_")) %>%
  rename_with( ~ paste0("W_", .x), -name) %>%
  rename_with(~ paste0(.x, "_2019"), -name) %>%
  separate(name, c("name", "descriptor"), " ", extra = "drop", fill = "right") %>%
  mutate(name = as_factor(str_replace(name, "\\.", ""))) %>%
  select(-descriptor) %>%
  mutate(name = fct_recode(name, Plantago_major_ssp_major = "Plantago_major"))
speciesL20 <- read_csv2("data_raw_species_2020_land.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                          cols(
                            .default = col_double(),
                            name = col_factor()
                          )) %>%
  mutate(name = str_replace(name, " ", "_")) %>%
  rename_with( ~ paste0("L_", .x), -name) %>%
  rename_with(~ paste0(.x, "_2020"), -name) %>%
  separate(name, c("name", "descriptor"), " ", extra = "drop", fill = "right") %>%
  mutate(name = as_factor(str_replace(name, "\\.", ""))) %>%
  select(-descriptor)
speciesW20 <- read_csv2("data_raw_species_2020_water.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                          cols(
                            .default = col_double(),
                            name = col_factor()
                          )) %>%
  mutate(name = str_replace(name, " ", "_")) %>%
  rename_with( ~ paste0("W_", .x), -name) %>%
  rename_with(~ paste0(.x, "_2020"), -name) %>%
  mutate(name = str_replace(name, "Plantago_major ssp. major", "Plantago_major_ssp_major")) %>%
  mutate(name = str_replace(name, "Plantago_major ssp. intermedia", "Plantago_major_ssp_intermedia")) %>%
  separate(name, c("name", "descriptor"), " ", extra = "drop", fill = "right") %>%
  mutate(name = as_factor(str_replace(name, "\\.", ""))) %>%
  select(-descriptor)

### b Join species tables -------------------------------------------------------------------------------------------
speciesL <- speciesL18 %>%
  full_join(speciesL19, by = "name") %>%
  full_join(speciesL20, by = "name") 
speciesW <- speciesW18 %>%
  full_join(speciesW19, by = "name") %>%
  full_join(speciesW20, by = "name") 
species <- full_join(speciesL, speciesW) %>%
  drop_na(name) %>%
### Rename species ###
  mutate(name = fct_recode(name, Carex_praecox_ssp_praecox = "Carex_praecox")) %>%
  mutate(name = fct_recode(name, Cerastium_fontanum_ssp_vulgare = "Cerastium_fontanum")) %>%
  mutate(name = fct_recode(name, Cerastium_fontanum_ssp_vulgare = "Cerastium_holosteoides")) %>%
  mutate(name = fct_recode(name, Cornus_controversa = "Cornus_sanguinea")) %>%
  mutate(name = fct_recode(name, Cota_tinctoria = "Anthemis_tinctoris")) %>%
  mutate(name = fct_recode(name, Erigeron_canadensis = "Conyza_canadensis")) %>%
  mutate(name = fct_recode(name, Fallopia_convolvulus = "Polygonum_convolvulus")) %>%
  mutate(name = fct_recode(name, Galium_mollugo = "Galium_album")) %>%
  mutate(name = fct_recode(name, Jacobaea_vulgaris = "Senecio_jacobaea")) %>%
  mutate(name = fct_recode(name, Persicaria_amphibia = "Polygonum_amphibium")) %>%
  mutate(name = fct_recode(name, Persicaria_bistorta = "Bistorta_officinalis")) %>%
  mutate(name = fct_recode(name, Persicaria_minor = "Polygonum_minus")) %>%
  mutate(name = fct_recode(name, Populus_alba = "Popolus_alba")) %>%
  mutate(name = fct_recode(name, Silene_latifolia_ssp_alba = "Silene_latifolia")) %>%
  mutate(name = fct_recode(name, Silene_latifolia_ssp_alba = "Silene_alba")) %>%
  mutate(name = fct_recode(name, "Silene_flos-cuculi" = "Lychnis_flos-cuculi")) %>%
  mutate(name = fct_recode(name, Taraxacum_campylodes = "Taraxacum_officinale")) %>%
  mutate(name = fct_recode(name, Vicia_sativa_ssp_nigra = "Vicia_angustifolia")) %>%
### combine created duplicates ###
  group_by(name) %>%
  summarise(across(where(is.double), ~sum(.x, na.rm = T))) %>%
### Implement counted inidviduals as 0.5 % in normal columns ###
  mutate(across(where(is.numeric) & contains("i"), ~ 0.5 * (. > 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "n") %>%
  mutate(id  =  str_replace(id, "i", "")) %>%
  group_by(name, id) %>%
  summarise(max = max(n, na.rm = T)) %>%
  mutate(max = ifelse(is.infinite(max), NA, max)) %>%
  pivot_wider(names_from = "id", values_from = "max") %>%
  mutate(name = as_factor(name)) %>%
### Check that each species occurs at least one time ###
  group_by(name) %>%
  mutate(total = sum(c_across(starts_with("L") | starts_with("W")), na.rm = T)) %>%
  mutate(presence = if_else(total > 0, 1, 0)) %>%
  filter(presence == 1) %>%
  ungroup() %>%
  select(-total, -presence) %>%
  rename_with(~ str_replace(.x, "_1_201", "_01_201")) %>%
  rename_with(~ str_replace(.x, "_1_2020", "_01_2020")) %>%
  rename_with(~ str_replace(.x, "_2_201", "_02_201")) %>%
  rename_with(~ str_replace(.x, "_2_2020", "_02_2020")) %>%
  rename_with(~ str_replace(.x, "_3_201", "_03_201")) %>%
  rename_with(~ str_replace(.x, "_3_2020", "_03_2020")) %>%
  rename_with(~ str_replace(.x, "_4_201", "_04_201")) %>%
  rename_with(~ str_replace(.x, "_4_2020", "_04_2020")) %>%
  rename_with(~ str_replace(.x, "_5_201", "_05_201")) %>%
  rename_with(~ str_replace(.x, "_5_2020", "_05_2020")) %>%
  rename_with(~ str_replace(.x, "_6_201", "_06_201")) %>%
  rename_with(~ str_replace(.x, "_6_2020", "_06_2020")) %>%
  rename_with(~ str_replace(.x, "_7_201", "_07_201")) %>%
  rename_with(~ str_replace(.x, "_7_2020", "_07_2020")) %>%
  rename_with(~ str_replace(.x, "_8_201", "_08_201")) %>%
  rename_with(~ str_replace(.x, "_8_2020", "_08_2020")) %>%
  rename_with(~ str_replace(.x, "_9_201", "_09_201")) %>%
  rename_with(~ str_replace(.x, "_9_2020", "_09_2020")) %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(name = as_factor(name))
rm(list = setdiff(ls(), c("species")))


### 2 Traits #####################################################################################

traits <- read_csv2("data_raw_traits.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                      cols(
                        .default = col_factor(),
                        name = col_factor(),
                        l = col_double(),
                        t = col_double(),
                        k = col_double(),
                        f = col_double(),
                        r = col_double(),
                        n = col_double()
                      )) %>%
  separate(name, c("genus", "species", "ssp", "subspecies"), "_", remove = F, extra = "drop", fill = "right") %>%
  mutate(genus = str_sub(genus, 1, 4)) %>%
  mutate(species = str_sub(species, 1, 4)) %>%
  mutate(subspecies = str_sub(subspecies, 1, 4)) %>%
  unite(abb, genus, species, subspecies, sep = "") %>%
  mutate(abb = str_replace(abb, "NA", "")) %>%
  mutate(abb = as_factor(abb)) %>%
  select(-ssp, - descriptor, -synonym, -nomenclature, -legal, -l, -k, -fchange) %>%
#traits[duplicated(traits$abb),]
#Check congruency of traits and species table
#traits$name[which(!(traits$name %in% species$name))]
#species$name[which(!(species$name %in% traits$name))]
  semi_join(species, by = "name")


### 3 Sites #####################################################################################

sites <- read_csv2("data_raw_sites.csv", col_names = T, na = c("", "NA", "na"), col_types = 
                     cols(
                       .default = col_factor(),
                       vegetationCov_2018 = col_double(),
                       vegetationCov_2019 = col_double(),
                       vegetationCov_2020 = col_double(),
                       surveyDate_2018 = col_date(),
                       surveyDate_2019 = col_date(),
                       surveyDate_2020 = col_date()
                       )) %>%
  select(-starts_with("surveyDate"), -starts_with("botanist")) %>%
  pivot_longer(starts_with("vegetationCov"), names_to = "surveyYear", values_to ="vegetationCov") %>%
  mutate(plot = str_replace(plot, "-", "_")) %>%
  mutate(surveyYear = str_replace(surveyYear, "vegetationCov_", "")) %>%
  mutate(id = as_factor(str_c(plot, surveyYear, sep = "_")), .keep = "all") %>%
  mutate(plot = as_factor(plot))
### Remove plots with no species ###
ids <- as_factor(colnames(species[-1]))
#sites$id[which(!(sites$id %in% ids))]
#ids[which(!(ids %in% sites$id))]
sites <- sites %>%
  filter(id %in% ids)
### Remove plots which are not part of the sites tibble ##
species <- species %>%
  select(name, all_of(sites$id))
rm(ids)
### Check missing values ###
#miss_var_summary(species, order = T)
#miss_var_summary(sites, order = T)
#vis_miss(species, cluster = F)
#vis_miss(sites, cluster = F)
#sites[is.na(sites$vegetationCov_2018),]



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Create simple variables #####################################################################################

sites <- sites %>%
  mutate(conf.low = c(1:length(id))) %>%
  mutate(conf.high = c(1:length(id)))

traits <- traits %>%
  mutate(target = if_else(
    targetHerb == "yes" | targetGrass == "yes", "yes", "no"
  ))


### 2 Coverages #####################################################################################

### a Graminoid covratio; graminoid, herb, and total coverage)-------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$family) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  mutate(type = if_else(type == "Poaceae" | type == "Cyperaceae" | type == "Juncaceae", "graminoidCov", "herbCov")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = T)) %>%
  spread(type, total) %>%
  mutate(graminoidCovratio = graminoidCov / (graminoidCov + herbCov)) %>%
  mutate(accumulatedCov = graminoidCov + herbCov) %>%
  ungroup() %>%
  mutate(graminoidCovratio = round(graminoidCovratio, 3), .keep = "unused") %>%
  mutate(accumulatedCov = round(accumulatedCov, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### b Target species' coverage -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$target) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(targetCov = round(total, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### c Ruderal species' coverage -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$ruderalExperiment) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(ruderalCov = round(total, 3), .keep = "unused") %>%
  right_join(sites, by = "id")

### d Seeded species' coverage -------------------------------------------------------------------------------------------
sites <- species %>%
  mutate(type = traits$seeded) %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id, type) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(type == "yes") %>%
  select(-type) %>%
  ungroup() %>%
  mutate(seededCov = round(total, 3), .keep = "unused") %>%
  right_join(sites, by = "id") %>%
### Create ruderalCovratio and seededCovratio ###
  mutate(ruderalCovratio = ruderalCov / accumulatedCov) %>%
  mutate(seededCovratio = seededCov / accumulatedCov)


### 3 Species richness #####################################################################################

specRich <- species %>%
  left_join(traits, by = "name") %>%
  select(rlg, rlb, target, targetHerb, ffh6510, ffh6210, starts_with("L"), starts_with("W"))

### a total species richness -------------------------------------------------------------------------------------------
specRich_all <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  group_by(id) %>%
  summarise(speciesRichness = sum(total, na.rm = T)) %>%
  ungroup()

### b red list Germany (species richness) -------------------------------------------------------------------------------------------
specRich_rlg <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id, rlg) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  group_by(id) %>%
  summarise(rlgRichness = sum(total, na.rm = T)) %>%
  ungroup()

### c red list Bavaria (species richness) -------------------------------------------------------------------------------------------
specRich_rlb <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id, rlb) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  group_by(id) %>%
  summarise(rlbRichness = sum(total, na.rm = T)) %>%
  ungroup()

### d target species (species richness) -------------------------------------------------------------------------------------------
specRich_target <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id, target) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(target == "yes") %>%
  group_by(id) %>%
  summarise(targetRichness = sum(total, na.rm = T)) %>%
  ungroup()

### g ffh6510 species (species richness) -------------------------------------------------------------------------------------------
specRich_ffh6510 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id, ffh6510) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(ffh6510 == "yes") %>%
  group_by(id) %>%
  summarise(ffh6510Richness = sum(total, na.rm = T)) %>%
  ungroup()

### h ffh6210 species (species richness) -------------------------------------------------------------------------------------------
specRich_ffh6210 <- specRich %>%
  pivot_longer(names_to = "id", values_to = "n", cols = starts_with("L") | starts_with("W")) %>%
  group_by(id, ffh6210) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  summarise(total = sum(n, na.rm = T)) %>%
  filter(ffh6210 == "yes") %>%
  group_by(id) %>%
  summarise(ffh6210Richness = sum(total, na.rm = T)) %>%
  ungroup()

### m implement in sites data set -------------------------------------------------------------------------------------------
sites <- sites %>%
  right_join(specRich_all, by = "id") %>%
  right_join(specRich_rlg, by = "id") %>%
  right_join(specRich_rlb, by = "id") %>%
  right_join(specRich_target, by = "id") %>%
  right_join(specRich_ffh6510, by = "id") %>%
  right_join(specRich_ffh6210, by = "id") %>%
### Create targetRichratio ###
  mutate(targetRichratio = targetRichness / speciesRichness) %>%
  mutate(targetRichratio = round(targetRichratio, 3))
rm(list=setdiff(ls(), c("sites", "species", "traits")))


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
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ntraits <- column_to_rownames(Ntraits, "name")
### Calculate CWM ###
Nweighted <- dbFD(Ntraits, Nspecies, w.abun = T,
                        calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### b F value -------------------------------------------------------------------------------------------
Ftraits <- traits %>%
  select(name, f) %>%
  filter(f > 0) 
Fspecies <- semi_join(species, Ftraits, by = "name") %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
Ftraits <- column_to_rownames(Ftraits, "name")
### Calculate CWM ###
Fweighted <- dbFD(Ftraits, Fspecies, w.abun = T,
                  calc.FRic = F, calc.FDiv = F, corr = "sqrt")

### c T value -------------------------------------------------------------------------------------------
Ttraits <- traits %>%
  select(name, t) %>%
  filter(t > 0) 
Tspecies <- semi_join(species, Ttraits, by = "name") %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) %>%
  column_to_rownames("site")
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
### check completeness of Roots ####
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
#traits %>%
  #filter(group != "tree" & group != "shrub") %>%
  #left_join(species, by = "name") %>%
  #count() #241 species remain
traits <- traits %>%
  mutate(name = as.character(name)) %>%
  arrange(name) %>%
  mutate(name = as_factor(name))
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
224 / (241 - 12) #97.8%

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
226 / (241 - 12) #98.7%

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
226 / (241 - 12) #98.7%

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
228 / (241 - 12) #99.6%

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
137 / (241 - 12) #59.8% available

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
157 / (241 - 12) #68.6% available

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
137 / (241 - 12) #59.8% available
rm(list=setdiff(ls(), c("sites", "species", "traits", "tbiPa", "tbiAbu")))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ##############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


setwd("Z:/Documents/0_Uni/2022_Donaudeiche/3_Aufnahmen_und_Ergebnisse/2022_Danube_dike_experiment/data/processed")
write_csv2(sites, "data_processed_sites.csv")
write_csv2(species, "data_processed_species.csv")
write_csv2(traits, "data_processed_traits.csv")
write_csv2(tbiAbu, "data_processed_tbiAbu.csv")
write_csv2(tbiPa, "data_processed_tbiPa.csv")
