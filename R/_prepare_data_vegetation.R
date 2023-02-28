# Grassland experiment on dikes
# Prepare species, sites, and traits data ####

# Markus Bauer
# 2023-02-18



# Content #####################################################################

# A Load data +++++++++++++++++++++++++++++++++++++
## 1 Sites and sPlotOpen
## 2 Species and sPlotOpen
## 3 Traits
## 5 Check data frames

# B Create variables ++++++++++++++++++++++++++++++
## 1 Create simple variables
## 2 Coverages
## 3 Alpha diversity
## 4 Reference sites
## 5 NMDS
## 6 TBI: Temporal beta diversity index
## 7 Functional traits

# C Save processed data +++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
suppressPackageStartupMessages(library(lubridate))
library(naniar)
library(vegan)
suppressPackageStartupMessages(library(adespatial))
suppressPackageStartupMessages(library(FD))
#remotes::install_github(file.path("inbo", "checklist"))
#checklist::setup_project(path = here())

### Start ###
rm(list = ls())
#installr::updateR(browse_news = FALSE, install_R = TRUE, copy_packages = TRUE, copy_Rprofile.site = TRUE, keep_old_packages = TRUE, update_packages = TRUE, start_new_R = TRUE, quit_R = TRUE, print_R_versions = TRUE, GUI = FALSE)
checklist::check_source()
#renv::status()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#_______________________________________________________________________________
## 1 Sites ####################################################################


sites_experiment <- read_csv(
  here("data", "raw", "data_raw_sites.csv"),
  col_names = TRUE,
  na = c("", "NA", "na"),
  col_types =
    cols(
      .default = "?",
      survey_date.seeded = col_date(format = "%Y-%m-%d"),
      survey_date.2018 = col_date(format = "%Y-%m-%d"),
      survey_date.2019 = col_date(format = "%Y-%m-%d"),
      survey_date.2020 = col_date(format = "%Y-%m-%d"),
      survey_date.2021 = col_date(format = "%Y-%m-%d"),
      botanist.2018 = "c",
      botanist.2019 = "c",
      botanist.2020 = "c",
      botanist.2021 = "c",
      vegetation_cover.2018 = "d",
      vegetation_cover.2019 = "d",
      vegetation_cover.2020 = "d",
      vegetation_cover.2021 = "d",
      biomass.2019 = "d"
      )
  ) %>%
  pivot_longer(
    starts_with("vegetation_cover") |
      starts_with("botanist") |
      starts_with("biomass") |
      starts_with("survey_date"),
    names_to = c("x", "survey_year"),
    names_sep = "\\.",
    values_to = "n",
    values_transform = list(n = as.character)
  ) %>%
  pivot_wider(names_from = "x", values_from = "n") %>%
  mutate(
    plot = str_replace(plot, "-", "_"),
    plot = str_replace(plot, "L_", "L"),
    plot = str_replace(plot, "W_", "W"),
    id = str_c(plot, survey_year, sep = "_"),
    plot = factor(plot),
    id = factor(id),
    vegetation_cover = as.numeric(vegetation_cover),
    biomass = as.numeric(biomass)
  ) %>%
  filter(
    !(site == "C" & (
      survey_year == "seeded" |
        survey_year == "2018" |
        survey_year == "2019" |
        survey_year == "2020" |
        survey_year == "2021"
    ))
  )

### Sabatini et al. (2021) Global Ecol Biogeogr:
### https://doi.org/10.1111/geb.13346
sites_splot <- read_delim(
  here("data", "raw", "sabatini_etal_2021", "sPlotOpen_header.txt"),
  col_names = TRUE, na = c("", "NA", "na"),
  col_types = cols(
    .default = "?",
    Cover_algae_layer = "d"
  )
)

### Bauer et al. (2022) Zenodo:
### https://doi.org/10.5281/zenodo.6334100
sites_bauer <- read_csv(
  here("data", "raw", "bauer_etal_2022", "data_sites_bauer_etal_2022.csv"),
  col_names = TRUE,
  na = c("", "NA", "na"),
  col_types = cols(
    .default = "?"
  )
)




#_______________________________________________________________________________
## 2 Species ###################################################################


species_experiment <- data.table::fread(
  here("data", "raw", "data_raw_species_20211112.csv"),
  sep = ",",
  dec = ".",
  skip = 0,
  header = TRUE,
  na.strings = c("", "NA", "na"),
  colClasses = list(
    character = "name"
    )
  ) %>%
  ### Check that each species occurs at least one time ###
  group_by(name) %>%
  arrange(name) %>%
  select(name, all_of(sites_experiment$id)) %>%
  mutate(
    total = sum(c_across(
      starts_with("L") | starts_with("W")),
      na.rm = TRUE),
    presence = if_else(total > 0, 1, 0)
    ) %>%
  # filter only species which occur at least one time:
  filter(presence == 1) %>%
  ungroup() %>%
  select(name, sort(tidyselect::peek_vars()), -total, -presence) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))

### Sabatini et al. (2021) Global Ecol Biogeogr:
### https://doi.org/10.1111/geb.13346
species_splot <- read_delim(
  here("data", "raw", "sabatini_etal_2021", "sPlotOpen_DT.txt"),
  col_names = TRUE, na = c("", "NA", "na"), col_types =
    cols(
      .default = "?"
      )
  ) %>%
  filter(Abundance_scale == "CoverPerc")

### Bauer et al. (2022) Zenodo:
### https://doi.org/10.5281/zenodo.6334100
species_bauer <- read_csv(
  here("data", "raw", "bauer_etal_2022",
       "data_species_bauer_etal_2022.csv"),
  col_names = TRUE, na = c("", "NA", "na"),
  col_types = cols(
    .default = "?"
  )
)



#_______________________________________________________________________________
## 3 Traits ####################################################################


traits <- read_csv(
  here("data", "raw", "data_raw_traits.csv"),
  col_names = TRUE,
  na = c("", "NA", "na"),
  col_types =
    cols(
      .default = "c",
      name = "c",
      l = "d",
      t = "d",
      k = "d",
      f = "d",
      r = "d",
      n = "d",
    )
) %>%
  separate(name, c("genus", "species", "ssp", "subspecies"), "_",
           remove = FALSE, extra = "drop", fill = "right") %>%
  mutate(genus = str_sub(genus, 1, 4),
         species = str_sub(species, 1, 4),
         subspecies = str_sub(subspecies, 1, 4),
         name = as.character(name)) %>%
  unite(abb, genus, species, subspecies, sep = "", na.rm = TRUE) %>%
  mutate(abb = str_replace(abb, "NA", ""),
         abb = as_factor(abb)) %>%
  select(-ssp, -synonym, -nomenclature, -legal, -l, -k, -fchange) %>%
  arrange(name)

### Check congruency of traits and species table ###
anti_join(traits, species_experiment, by = "name") %>% select(name)
anti_join(species_experiment, traits, by = "name") %>% select(name)
data <- traits %>%
  mutate(name = str_replace_all(name, "_", " ")) %>%
  select(abb, name)
TNRS::TNRS(
  taxonomic_names = data,
  sources = c("wfo", "tropicos", "wcvp"),
  classification = "wfo",
  mode = "resolve"
) %>%
  select(ID, Name_submitted, Overall_score, Source, Accepted_name,
         Accepted_name_rank, Accepted_name_author, Accepted_family,
         Accepted_name_url, Taxonomic_status) %>%
  filter(
    Overall_score < 1 | Source != "wfo" | Taxonomic_status != "Accepted"
  ) %>%
  arrange(Source)


### Combine with species_experiment table ###
traits <- traits %>%
  semi_join(species_experiment, by = "name")

rm(list = setdiff(ls(), c(
  "sites_experiment", "sites_splot", "sites_bauer",
  "species_experiment", "species_splot", "species_bauer",
  "traits"
  )))



#_______________________________________________________________________________
## 4 Check data frames #########################################################


### Check typos ###
sites_experiment %>%
  filter(!str_detect(id, "_seeded$")) %>%
  janitor::tabyl(vegetation_cover)
#sites %>% filter(vegetation_cover == 17)
species_experiment %>%
  select(-name, -ends_with("_seeded")) %>%
  unlist() %>%
  janitor::tabyl()
species_experiment %>% # Check special typos
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  filter(value == 90)

### Compare vegetation_cover and accumulated_cover ###
species_experiment %>%
    summarise(across(where(is.double), ~sum(.x, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "id", values_to = "value") %>%
    mutate(id = factor(id)) %>%
    full_join(sites_experiment, by = "id") %>%
    mutate(diff = (value - vegetation_cover)) %>%
    select(id, survey_year, vegetation_cover, value, diff) %>%
    filter(!str_detect(id, "_seeded$")) %>%
    filter(diff > 20 | diff < -5) %>%
    arrange(survey_year, id, diff) %>%
    print(n = 100)

### Check plots over time ###
species_experiment %>%
  select(name, starts_with("L1_19"), -ends_with("_seeded")) %>%
  filter(if_any(starts_with("L"), ~ . > 0)) %>%
  print(n = 100)

### Check missing data ###
miss_var_summary(sites_experiment, order = TRUE)
vis_miss(sites_experiment, cluster = FALSE)
miss_var_summary(traits, order = TRUE)
vis_miss(traits, cluster = FALSE, sort_miss = TRUE)

rm(list = setdiff(ls(), c(
  "sites_experiment", "sites_splot", "sites_bauer",
  "species_experiment", "species_splot", "species_bauer",
  "traits"
)))



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ###########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#_______________________________________________________________________________
## 1 Create variables #########################################################


### a Simple variables --------------------------------------------------------

traits <- traits %>%
  mutate(
    target = if_else(
      ffh6510 == "1" | ffh6210 == "1" | biotope_target == "1" |
        table_31 == "1" | table_34 == "1" | table_35 == "1" | table_36 == "1" |
        mapped_2011 == "1" | special_target == "1",
    "1", "0"
    )
  )

sites <- sites_experiment %>%
  mutate(
    surveyDate = as_date(survey_date),
    seeding_date = if_else(
      exposition == "north", ymd("20180413"), ymd("20180427")
    ),
    age = interval(seeding_date, survey_date) %/% days(1),
    block = str_c(site, exposition, sep = "_"),
    block = factor(block)
  )


### b Establishment of species -------------------------------------------------

data <- species_experiment %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "n") %>%
  left_join(sites_experiment, by = "id") %>%
  mutate(n = if_else(n > 0 & survey_year == "seeded", 1, n)) %>%
  select(plot, name, survey_year, n) %>%
  pivot_wider(names_from = "survey_year", values_from = "n") %>%
  pivot_longer(-c(plot, name, seeded),
               names_to = "survey_year", values_to = "n") %>%
  mutate(n = if_else(seeded == 1 & n > 0, 1, 0)) %>%
  group_by(name, survey_year) %>%
  summarise(
    total_established = sum(n, na.rm = TRUE),
    total_seeded = sum(seeded, na.rm = TRUE),
    .groups = "keep"
    ) %>%
  filter(total_seeded > 0) %>%
  mutate(
    rate = total_established / total_seeded,
    rate = round(rate, digits = 2),
    seeded = "1"
  ) %>%
  pivot_wider(names_from = "survey_year",
              values_from = c("rate", "total_seeded", "total_established")) %>%
  select(-total_seeded_2019, -total_seeded_2020, -total_seeded_2021)

traits <- traits %>%
  left_join(data, by = "name")




#_______________________________________________________________________________
## 2 Coverages #################################################################


cover <- species_experiment %>%
  left_join(traits, by = "name") %>%
  select(name, family, target, seeded,
         starts_with("L"), starts_with("W"), starts_with("C")) %>%
  pivot_longer(
    names_to = "id",
    values_to = "n",
    cols = starts_with("L") | starts_with("W") | starts_with("C")
    ) %>%
  group_by(id)

#### Graminoid, herb, and total coverage) ###
cover_total_and_graminoid <- cover %>%
  group_by(id, family) %>%
  summarise(total = sum(n, na.rm = TRUE), .groups = "keep") %>%
  mutate(
    type = if_else(
      family == "Poaceae" |
        family == "Cyperaceae" |
        family == "Juncaceae",
      "graminoid_cover", "herb_cover")
  ) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = TRUE), .groups = "keep") %>%
  spread(type, total) %>%
  mutate(
    accumulated_cover = graminoid_cover + herb_cover,
    accumulated_cover = round(accumulated_cover, 1)
  ) %>%
  ungroup()

#### Target specis' coverage ###
cover_target <- cover %>%
  filter(target == "1") %>%
  summarise(target_cover = sum(n, na.rm = TRUE)) %>%
  mutate(target_cover = round(target_cover, 1)) %>%
  ungroup()

#### Seeded species' coverage ###
cover_seeded <- species_experiment %>%
  select(-contains("x0809")) %>%
  # Make two columns out of column id:
  pivot_longer(-name, names_to = "id", values_to = "value",
               values_drop_na = TRUE) %>%
  separate(id, c("plot", "survey_year"), sep = "_(?!.*_)",
           remove = FALSE, extra = "merge", fill = "warn", convert = FALSE) %>%
  # Summarise for plot:
  pivot_wider(names_from = "survey_year", values_from = "value") %>%
  group_by(plot, name) %>%
  summarise(across(where(is.double), ~sum(.x, na.rm = TRUE)),
            .groups = "keep") %>%
  ungroup() %>%
  # Combine with seedmixes:
  pivot_longer(starts_with("20"),
               names_to = "survey_year", values_to = "value",
               names_transform = list(survey_year = as.factor)) %>%
  mutate(success = if_else(seeded > 0 & value > 0, value, 0)) %>%
  group_by(plot, survey_year) %>%
  summarise(seeded_cover = sum(success, na.rm = TRUE), .groups = "keep") %>%
  ungroup() %>%
  unite(id, plot, survey_year, sep = "_")

#### Implement in sites data set ###
sites_experiment <- sites_experiment %>%
  left_join(cover_total_and_graminoid, by = "id") %>%
  left_join(cover_target, by = "id") %>%
  left_join(cover_seeded, by = "id") %>%
  ### Calculate the ratio of target cover of total cover
  mutate(
    target_cover_ratio = target_cover / accumulated_cover,
    graminoid_cover_ratio = graminoid_cover / accumulated_cover,
    seeded_cover_ratio = seeded_cover / accumulated_cover,
    target_cover_ratio = round(target_cover_ratio, 3),
    graminoid_cover_ratio = round(graminoid_cover_ratio, 3),
    seeded_cover_ratio = round(seeded_cover_ratio, 3)
    )

rm(list = setdiff(ls(), c(
  "sites_experiment", "sites_splot", "sites_bauer",
  "species_experiment", "species_splot", "species_bauer",
  "traits"
)))



#_______________________________________________________________________________
## 3 Alpha diversity ##########################################################


### a Species richness ---------------------------------------------------------

species_richness <- species_experiment %>%
  left_join(traits, by = "name") %>%
  select(name, rlg, rlb, target, ffh6510, ffh6210,
         starts_with("L"), starts_with("W"), starts_with("C")) %>%
  pivot_longer(
    names_to = "id", values_to = "n",
    cols = starts_with("L") | starts_with("W") | starts_with("C")
    ) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  group_by(id)

#### Total species richness ###
species_richness_all <- species_richness %>%
  summarise(species_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### Red list Germany (species richness) ###
species_richness_rlg <- species_richness %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  summarise(rlg_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### Red list Bavaria (species richness) ###
species_richness_rlb <- species_richness %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  summarise(rlb_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### Target species (species richness) ###
species_richness_target <- species_richness %>%
  filter(target == "1") %>%
  summarise(target_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### Seeded species (species richness) ###
species_richness_seeded <- species_experiment %>%
  select(-contains("x0809")) %>%
  # Make two columns out of column id:
  pivot_longer(-name, names_to = "id", values_to = "value",
               values_drop_na = TRUE) %>%
  separate(id, c("plot", "survey_year"), sep = "_(?!.*_)",
           remove = FALSE, extra = "merge", fill = "warn", convert = FALSE) %>%
  # Summarise for plot:
  pivot_wider(names_from = "survey_year", values_from = "value") %>%
  group_by(plot, name) %>%
  summarise(across(where(is.double), ~sum(.x, na.rm = TRUE)),
            .groups = "keep") %>%
  ungroup() %>%
  # Combine with seedmixes:
  pivot_longer(starts_with("20"), names_to = "survey_year", values_to = "value",
               names_transform = list(survey_year = as.factor)) %>%
  mutate(success_cover = if_else(seeded > 0 & value > 0, 1, 0)) %>%
  group_by(plot, survey_year) %>%
  summarise(seeded_richness = sum(success_cover, na.rm = TRUE),
            .groups = "keep") %>%
  ungroup() %>%
  unite(id, plot, survey_year, sep = "_")

### FFH6510 species (species richness) ###
species_richness_ffh6510 <- species_richness %>%
  filter(ffh6510 == "1") %>%
  summarise(ffh6510_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### FFH6210 species (species richness) ###
species_richness_ffh6210 <- species_richness %>%
  filter(ffh6210 == "1") %>%
  summarise(ffh6210_richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

#### Implement in sites data set ###
sites_experiment <- sites_experiment %>%
  left_join(species_richness_all, by = "id") %>%
  left_join(species_richness_rlg, by = "id") %>%
  left_join(species_richness_rlb, by = "id") %>%
  left_join(species_richness_target, by = "id") %>%
  left_join(species_richness_seeded, by = "id") %>%
  left_join(species_richness_ffh6510, by = "id") %>%
  left_join(species_richness_ffh6210, by = "id") %>%

  ### * Calculate Favourable Conservation Status (FCS) ####
  ### Helm et al. 2014 Divers Distrib
  ### https://doi.org/10.1111/ddi.12285

  mutate(
    fcs_target = log(
      (target_richness + 1) / (species_richness - target_richness + 1)
      ),
    fcs_seeded = log(
      (seeded_richness + 1) / (species_richness - seeded_richness + 1)
      ),
    fcs_target = round(fcs_target, 3),
    fcs_seeded = round(fcs_seeded, 3)
    )

### b Species evenness and shannon ---------------------------------------------

data <- species_experiment  %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("id") %>%
  diversity(index = "shannon") %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "id") %>%
  mutate(id = factor(id)) %>%
  rename(shannon = value)
sites_experiment <- sites_experiment %>%
  left_join(data, by = "id") %>%
  mutate(
    evenness = shannon / log(species_richness),
    shannon = round(shannon, digits = 3),
    evenness = round(evenness, digits = 3)
    )

rm(list = setdiff(ls(), c(
  "sites_experiment", "sites_splot", "sites_bauer",
  "species_experiment", "species_splot", "species_bauer",
  "traits"
)))



#_______________________________________________________________________________
## 4 Reference sites ##########################################################


### a ESy: EUNIS expert vegetation classification system ----------------------

#### Start ###
### Bruelheide et al. 2021 Appl Veg Sci
### https://doi.org/10.1111/avsc.12562

expertfile <- "EUNIS-ESy-2020-06-08.txt" ### file of 2021 is not working

obs <- species_experiment %>%
  pivot_longer(cols = -name,
               names_to = "RELEVE_NR",
               values_to = "Cover_Perc") %>%
  rename(TaxonName = "name") %>%
  mutate(
    TaxonName = str_replace_all(TaxonName, "_", " "),
    TaxonName = str_replace_all(TaxonName, "ssp", "subsp."),
    TaxonName = as.factor(TaxonName),
    TaxonName = fct_recode(
      TaxonName,
      "Cerastium fontanum" = "Cerastium fontanum subsp. vulgare",
      "Silene latifolia" = "Silene latifolia subsp. alba"
    )
  ) %>%
  data.table::as.data.table()

header <- sites_experiment %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  rename(
    RELEVE_NR = id
  ) %>%
  mutate(
    "Altitude (m)" = 313,
    Latitude = sf::st_coordinates(.)[, 2],
    Longitude = sf::st_coordinates(.)[, 1],
    Country = "Germany",
    Coast_EEA = "N_COAST",
    Dunes_Bohn = "N_DUNES",
    Ecoreg = 686,
    dataset = "Danube_experiment"
  ) %>%
  select(RELEVE_NR, "Altitude (m)", Latitude, Longitude, Country,
         Coast_EEA, Dunes_Bohn, Ecoreg, dataset) %>%
  sf::st_drop_geometry()

setwd(here("R", "esy"))
source(here("R", "esy", "code", "prep.R"))

#### Step 1 and 2: Load and parse the expert file ###
source(here("R", "esy", "code", "step1and2_load-and-parse-the-expert-file.R"))

#### Step 3: Create a numerical plot x membership condition matrix  ###
plot.cond <- array(
  0,
  c(length(unique(obs$RELEVE_NR)), length(conditions)),
  dimnames = list(
    as.character(unique(obs$RELEVE_NR)),
    conditions
  )
)

### Step 4: Aggregate taxon levels ###
source(here("R", "esy", "code", "step4_aggregate-taxon-levels.R"))

(data <- obs %>%
    group_by(TaxonName) %>%
    slice(1) %>%
    anti_join(AGG, by = c("TaxonName" = "ind")))

#### Step 5: Solve the membership conditions ###
mc <- 1
source(here("R", "esy", "code",
            "step3and5_extract-and-solve-membership-conditions.R"))

table(result.classification)
eval.EUNIS(which(result.classification == "V12")[4], "V12")

sites_experiment <- sites_experiment %>%
  mutate(
    esy = result.classification,
    esy = if_else(id == "W1_18_2021", "R", esy),
    esy = if_else(id == "L6_01_2018", "V", esy),
    esy = if_else(id == "L6_04_2018", "V", esy),
    esy = if_else(id == "W1_18_2019", "R", esy),
    esy = if_else(id == "L6_13_2021", "V", esy),
    esy = if_else(id == "W3_06_2019", "R", esy),
    esy = if_else(id == "W5_06_2019", "R", esy),
    esy = if_else(id == "W5_24_2019", "R", esy),
    esy = if_else(id == "W6_12_2019", "R", esy),
    esy = if_else(id == "L2_14_seeded", "R", esy),
    esy = if_else(id == "W2_11_2019", "R", esy),
    esy = if_else(id == "W3_07_seeded", "R", esy),
    esy = if_else(id == "W3_24_2019", "R", esy),
    esy = if_else(id == "W4_17_2019", "R", esy),
    esy = if_else(id == "L2_07_2018", "V", esy),
    esy = if_else(id == "L2_12_2019", "V", esy),
    esy = if_else(id == "L2_12_2020", "V", esy),
    esy = if_else(id == "L3_12_2018", "V", esy),
    esy = if_else(id == "L4_03_2018", "V", esy),
    esy = if_else(id == "L5_06_2018", "V", esy),
    esy = if_else(id == "L5_09_2019", "V", esy),
    esy = if_else(id == "L6_18_2018", "R", esy),
    esy = if_else(id == "W1_21_2018", "V", esy),
    esy = if_else(id == "W3_03_2018", "V", esy),
    esy = if_else(id == "W4_11_2018", "V", esy),
    esy = if_else(id == "W4_24_2018", "V", esy)
    )
table(sites_experiment$esy)

rm(list = setdiff(ls(), c(
  "sites_experiment", "sites_splot", "sites_bauer",
  "species_experiment", "species_splot", "species_bauer",
  "traits", "result.classification"
)))


### b sPlotOpen ---------------------------------------------------------------

### Sabatini et al. (2021) Global Ecol Biogeogr
### https://doi.org/10.1111/geb.13346

data_sites <- sites_splot %>%
  filter(
    # Hay meadow: EUNIS2007 code E2.2
    # Chytry et al. 2020 Appl Veg Sci
    # https://doi.org/10.1111/avsc.12519
    (ESY == "E22" |
       # Dry grassland: EUNIS2007 code E1.2a
       ESY == "E12a") &
      Releve_area >= 1 &
      Releve_area <= 25 &
      Longitude > 10.89845 & # West: Augsburg
      Longitude < 13.46434 & # East: Passau
      Latitude > 	47.85298 & # South: Rosenheim
      Latitude < 49.45095 & # North: Nuernberg
      Elevation < 700
  ) %>%
  rename_with(tolower) %>%
  rename(id = plotobservationid, survey_year = date_of_recording,
         plot_size = releve_area, reference = country) %>%
  mutate(
    id = paste0("S", id),
    reference = str_replace(reference, "Germany", "reference"),
    survey_year = year(survey_year),
    source = if_else(
      givd_id == "EU-DE-014",
      "Jandt & Bruelheide (2012) Biodivers Ecol https://doi.org/10.7809/b-e.00146",
      "other"
      ),
    longitude = longitude * 10^5,
    latitude = latitude * 10^5
  ) %>%
  select(id, givd_id, source, longitude, latitude, elevation, plot_size,
         survey_year, reference, esy) %>%
  mutate(
    survey_year = as.character(survey_year),
    source = "Sabatini et al. (2021) Global Ecol Biogeogr https://doi.org/10.1111/geb.13346",
    reference = if_else(
      esy == "E12a", "positive_reference", if_else(
        esy == "E22", "positive_reference", "other"
      )
    ),
    target_type = if_else(
      esy == "E12a", "dry_grassland", if_else(
        esy == "E22", "hay_meadow", "other"
      )
    ),
    esy = if_else(
      esy == "E12a", "R1A", if_else(
        esy == "E22", "R22", "other"
      )
    ),
    exposition = "other"
  )
sites_splot <- data_sites

data_species <- species_splot %>%
  rename(id = PlotObservationID, name = Species,
         abundance = Original_abundance) %>%
  mutate(id = paste0("S", id)) %>%
  semi_join(data_sites, by = "id") %>%
  select(id, name, abundance) %>%
  pivot_wider(names_from = "id",
              values_from = "abundance",
              values_fn = sum) %>%
  mutate(
    name = str_replace(name, " ", "_"),
    name = str_replace(name, "Helianthemum_ovatum", "Helianthemum_nummularium"),
    name = str_replace(name, "Galium_album", "Galium_mollugo"),
    name = str_replace(name, "Taraxacum", "Taraxacum_campylodes"),
    name = str_replace(
      name, "Cerastium_fontanum", "Cerastium_fontanum_ssp_vulgare"
    ),
    name = str_replace(name, "Leucanthemum_ircutianum", "Leucanthemum_vulgare"),
    name = str_replace(name, "Tragopogon_orientalis", "Tragopogon_pratensis"),
    name = factor(name)
  ) %>%
  group_by(name) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)))

### Check species name congruency ###
data <- data_species %>%
  anti_join(traits, by = "name") %>%
  select(name) %>%
  print(n = 50)

species_splot <- data_species


### c Our own regional dike surveys -------------------------------------------

### Bauer et al. (2022) Zenodo:
### https://doi.org/10.5281/zenodo.6334100

data_sites <- sites_bauer %>%
  mutate(
    reference = "reference",
    plot_size = 25
    ) %>%
  select(id, longitude, latitude, plot_size, survey_year, reference,
         esy, exposition, orientation, plot_age, species_richness) %>%
  filter(exposition == "south" | exposition == "north") %>%
  mutate(
    survey_year = as.character(survey_year),
    source = "Bauer et al. (2023) EcoEvoRxiv https://doi.org/10.32942/X2959J",
    reference = if_else(
      esy == "R1A", "positive_reference", if_else(
        esy == "R22", "positive_reference", if_else(
          esy == "R", "Grassland", if_else(
            esy == "?", "no", if_else(
              esy == "+", "no", if_else(
                esy == "R21", "Grassland", if_else(
                  esy == "V38", "negative_reference", "other"
                )
              )
            )
          )
        )
      )
    ),
    target_type = if_else(
      esy == "R1A", "dry_grassland", if_else(
        esy == "R22", "hay_meadow", "other"
      )
    )
  ) %>%
  filter(reference != "no" & reference != "Grassland")
sites_bauer <- data_sites

data_species <- species_bauer %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  semi_join(data_sites, by = "id") %>%
  pivot_longer(cols = -id, names_to = "name", values_to = "value") %>%
  pivot_wider(names_from = "id", values_from = "value")

### Check species name congruency ###
data <- data_species %>%
  anti_join(traits, by = "name") %>%
  select(name) %>%
  print(n = 50)
  
species_bauer <- data_species

rm(list = setdiff(ls(), c(
  "sites_experiment", "sites_splot", "sites_bauer",
  "species_experiment", "species_splot", "species_bauer",
  "traits", "results.classification"
)))



#_______________________________________________________________________________
## 5 NMDS ######################################################################


### Create reference variable ###
data_experiment <- sites_experiment %>%
  mutate(
    reference = survey_year,
    source = "experiment"
  )

data_sites <- data_experiment %>%
  bind_rows(sites_splot, sites_bauer) %>%
  select(
    id, plot, site, esy, reference,
    exposition, sand_ratio, substrate_depth, target_type, seed_density,
    survey_year, longitude, latitude, elevation, plot_size, botanist,
    givd_id, source
  ) %>%
  arrange(id)
#### Prepare species data ###
### Calculate rare species (< 0.5% accumulated cover in all plots)
rare <- species_experiment %>%
  full_join(species_splot, by = "name") %>%
  full_join(species_bauer, by = "name") %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  group_by(name) %>%
  summarise(total_cover_species = sum(value)) %>%
  filter(total_cover_species < 0.5)

data_species <- species_experiment %>%
  full_join(species_splot, by = "name") %>%
  full_join(species_bauer, by = "name") %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = -name, names_to = "id", values_to = "value") %>%
  filter(!(name %in% rare$name)) %>% # Use 'rare' to exclude rare species
  pivot_wider(names_from = "name", values_from = "value") %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id")

anti_join(data_sites, data_species, by = "id")

data_sites <- data_sites %>%
  semi_join(data_species, by = "id")

data_species <- data_species %>%
  column_to_rownames("id")

### * Model ####
set.seed(12)
(ordi <- metaMDS(
  data_species,
  distance = "bray",
  binary = TRUE,
  na.rm = TRUE,
  try = 50,
  previous.best = TRUE
  ))
#save(ordi, file = here("outputs", "models", "model_nmds.Rdata"))
base::load(here("outputs", "models", "model_nmds.Rdata"))
ordi
stressplot(ordi)


### Integrate NMDS in dataset ###
sites_nmds <- data_sites %>%
  mutate(NMDS1 = ordi$points[, 1], NMDS2 = ordi$points[, 2])

### Calculate centroids of positive_reference ###
sites_nmds %>%
  filter(reference == "positive_reference") %>%
  group_by(exposition, target_type) %>%
  summarise(mean_reference = mean(NMDS1), sd_reference = sd(NMDS1))

sites_nmds <- sites_nmds %>%
  mutate(
    mean_reference = if_else(
      exposition == "north" & target_type == "dry_grassland", 0.644, if_else(
        exposition == "north" & target_type == "hay_meadow", 0.491, if_else(
          exposition == "south" & target_type == "dry_grassland", 0.512,
          if_else(
            exposition == "south" & target_type == "hay_meadow", 0.478, NA_real_
          )
        )
      )
    ),
    sd_reference = if_else(
      exposition == "north" & target_type == "dry_grassland", 0.279, if_else(
        exposition == "north" & target_type == "hay_meadow", 0.217, if_else(
          exposition == "south" & target_type == "dry_grassland", 0.337,
          if_else(
            exposition == "south" & target_type == "hay_meadow", 0.186, NA_real_
          )
        )
      )
    ),
    recovery_completeness = NMDS1 - mean_reference,
    recovery_completeness = round(recovery_completeness, digits = 3)
  )

rm(list = setdiff(ls(), c(
  "sites_experiment", "sites_splot", "sites_bauer", "sites_nmds",
  "species_experiment", "species_splot", "species_bauer",
  "traits", "results.classification"
)))



#_______________________________________________________________________________
## 6 Temporal beta diversity ###################################################


### Prepare data ###
data_sites <- sites_experiment %>%
  # Choose only plots which were surveyed in each year:
  filter(accumulated_cover > 0) %>%
  add_count(plot) %>%
  filter(n == max(n)) %>%
  select(id, plot)
data_species <- species_experiment %>%
  select(where(~ !all(is.na(.x)))) %>%
  mutate(across(where(is.numeric), ~ replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  separate(id, c("plot", "year"), sep = "_(?!.*_)",
           remove = FALSE, extra = "merge", fill = "warn", convert = FALSE) %>%
  arrange(id) %>%
  semi_join(data_sites, by = "id") %>%
  select(plot, year, tidyselect::peek_vars(), -id)

### Separate each year in several tibbles ###
for (i in unique(data_species$year)) {
  nam <- paste("species", i, sep = "_")
  
  assign(
    nam,
    data_species %>%
      filter(year == i) %>%
      select(-year) %>%
      column_to_rownames(var = "plot")
  )
}

### a Calculate TBI Presence --------------------------------------------------

### Temporal Beta diversity Index
### Legendre (2019) Ecol Evol
### https://doi.org/10.1002/ece3.4984

#### * seedmixes vs. 2018 ####
res18 <- TBI(species_seeded, species_2018,
             method = "sorensen",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res18$BCD.summary # B = .571, C = .304, D = .875 (.652 vs. .347)
res18$t.test_B.C
tbi18 <- res18$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2018")
#### Test plot
plot(res18, type = "BC")

#### * seedmixes vs. 2019 ####
res19 <- TBI(species_seeded, species_2019,
             method = "sorensen",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res19$BCD.summary # B = .377, C = .357, D = .735 (.513 vs. .486)
res19$t.test_B.C
tbi19 <- res19$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2019")
#### Test plot
plot(res19, type = "BC")

#### * seedmixes vs. 2020 ####
res20 <- TBI(species_seeded, species_2020,
             method = "sorensen",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res20$BCD.summary # B = .324, C = .392, D = .717 (.452 vs. .547)
res20$t.test_B.C
tbi20 <- res20$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2020")
#### Test plot
plot(res20, type = "BC")

#### * seedmixes vs. 2021 ####
res21 <- TBI(species_seeded, species_2021,
             method = "sorensen",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res21$BCD.summary # B = .349, C = .391, D = .740 (.471 vs. .528)
res21$t.test_B.C
tbi21 <- res21$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2021")
#### Test plot
plot(res21, type = "BC")

#### * Combine datasets ####
data_presence <- bind_rows(tbi18, tbi19, tbi20, tbi21) %>%
  mutate(presabu = "presence")

rm(list = setdiff(ls(), c(
  "sites_experiment", "species_experiment", "traits",
  "sites_bauer", "sites_splot", "species_bauer", "species_splot",
  "species_seeded", "species_2018", "species_2019", "species_2020",
  "species_2021", "data_presence", "data_sites", "data_species"
  )))

### b Calculate TBI Abundance -------------------------------------------------

#### * seedmixes vs. 2018 ####
res18 <- TBI(species_seeded, species_2018,
             method = "%difference",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res18$BCD.summary # B = .657, C = .315, D = .973 (.675 vs. .324)
res18$t.test_B.C
tbi18 <- res18$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2018")
#### Test plot
plot(res18, type = "BC")

#### * seedmixes vs. 2019 ####
res19 <- TBI(species_seeded, species_2019,
             method = "%difference",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res19$BCD.summary # B = .522, C = .383, D = .906 (.576 vs. .423)
res19$t.test_B.C
tbi19 <- res19$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2019")
#### Test plot
plot(res19, type = "BC")

#### * seedmixes vs. 2020 ####
res20 <- TBI(species_seeded, species_2020,
             method = "%difference",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res20$BCD.summary # B = .488, C = .382, D = .870 (.560 vs. .439)
res20$t.test_B.C
tbi20 <- res20$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2020")
#### Test plot
plot(res20, type = "BC")

#### * seedmixes vs. 2021 ####
res21 <- TBI(species_seeded, species_2021,
             method = "%difference",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res21$BCD.summary # B = .473, C = .390, D = .863 (.548 vs. .451)
res21$t.test_B.C
tbi21 <- res21$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2021")
#### Test plot
plot(res21, type = "BC")

#### * Combine datasets ####
data_abundance <- bind_rows(tbi18, tbi19, tbi20, tbi21) %>%
  mutate(presabu = "abundance")
plot <- data_sites %>%
  mutate(plot = factor(plot)) %>%
  filter(str_detect(id, "2018")) %>%
  pull(plot)

#### * Combine all TBIs ####
data <- data_presence %>%
  add_row(data_abundance) %>%
  mutate(
    plot = rep(plot, length(data_abundance$comparison) * 2 / 288),
    id = str_c(plot, comparison, sep = "_")
  ) %>%
  rename(
    B = "B/(2A+B+C)", C = "C/(2A+B+C)", D = "D=(B+C)/(2A+B+C)",
    change = Change
  ) %>%
  mutate(
    change = C - B,
    across(c(B, C, D, change), ~ round(.x, digits = 4))
  ) %>%
  select(id, B, C, D, presabu)

sites_temporal <- sites_experiment %>%
  filter(survey_year != "seeded") %>%
  left_join(data, by = "id") %>%
  select(
    id, plot, site, longitude, latitude, elevation, plot_size, exposition,
    orientation, sand_ratio, substrate_depth, target_type, seed_density,
    survey_year, botanist, B, C, D, presabu
    ) %>%
  mutate(persistence = (1 - B) * 100)

rm(list = setdiff(ls(), c(
  "sites_experiment", "species_experiment", "traits", "sites_temporal",
  "sites_bauer", "sites_splot", "species_bauer", "species_splot"
  )))



#_______________________________________________________________________________
## 7 Functional plant traits ###################################################


### a LEDA database ------------------------------------------------------------

### Kleyer et al. (2008) J Ecol
### https://doi.org/10.1111/j.1365-2745.2008.01430.x

data_sla <- data.table::fread(
  here("data", "raw", "leda_database",
       "data_raw_traitbase_leda_20210223_sla.txt"),
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
  here("data", "raw", "leda_database",
       "data_raw_traitbase_leda_20210223_seedmass.txt"),
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
  here("data", "raw", "leda_database",
       "data_raw_traitbase_leda_20210223_canopy_height.txt"),
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

### Check names of 'traits' which are not present in trait database ###
traits %>%
  filter(!str_detect(name, "_spec") & !str_detect(name, "aceae")) %>%
  filter(group != "tree" & group != "shrub") %>%
  anti_join(data, by = "name") %>%
  select(name)

### Search in 'name' for certain names ###
data %>%
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = TRUE))) %>%
  filter(str_detect(name, "ovalis"))

data <- data %>%
  mutate(name = as_factor(name),
         name = fct_recode(
           name,
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
           Vicia_sativa_ssp_nigra = "Vicia_sativa_s._nigra"
           )
         ) %>%
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = TRUE)))

### Check names of 'traits' which are not present in LEDA database ###
traits %>%
  filter(!str_detect(name, "_spec") & !str_detect(name, "aceae")) %>%
  filter(group != "tree" & group != "shrub") %>%
  anti_join(data, by = "name") %>%
  select(name)

traits <- traits %>%
  left_join(data, by = "name")

#### * Check completeness of LEDA ####

traits %>%
  filter(!str_detect(name, "_spec") & !str_detect(name, "aceae")) %>%
  filter(group != "tree" & group != "shrub") %>%
  select(sla, seedmass, height) %>%
  miss_var_summary(order = TRUE)
traits %>%
  filter(!str_detect(name, "_spec") & !str_detect(name, "aceae")) %>%
  filter(group != "tree" & group != "shrub") %>%
  select(sla, seedmass, height) %>%
  vis_miss(cluster = FALSE, sort_miss = TRUE)
(uncomplete_cases <- traits %>%
  select(name, sla, seedmass, height) %>%
  filter(!complete.cases(.)))


### b TRY database -------------------------------------------------------------

### Kattge et al. (2020) Glob Change Biol
### https://doi.org/10.1111/gcb.14904

data <- data.table::fread(
  here("data", "raw", "try_database",
       "data_raw_traitbase_try_20210306_13996.txt"),
  header = TRUE,
  sep = "\t",
  dec = ".",
  quote = ""
  ) %>%
  as_tibble() %>%
  bind_rows(
    data.table::fread(
      here("data", "raw", "try_database",
           "data_raw_traitbase_try_20210318_14157.txt"),
      header = TRUE,
      sep = "\t",
      dec = ".",
      quote = ""
    ) %>%
      as_tibble()
  ) %>%
  rename(
    name = "AccSpeciesName",
    value = "StdValue",
    trait = "TraitID"
  ) %>%
  select(name, value, trait, Reference) %>%
  mutate(name = str_replace_all(name, " ", "_")) %>%
  drop_na() %>%
  mutate(
    trait = str_replace(trait, "26", "seedmass"),
    trait = str_replace(trait, "3106", "height"),
    trait = str_replace(trait, "3107", "height"),
    trait = str_replace(trait, "3115", "sla"),
    trait = str_replace(trait, "3116", "sla"),
    trait = str_replace(trait, "3117", "sla")
  )

##### * Find synonyms ####

### Check names of 'uncomplete_cases' which are not present in TRY database ###
uncomplete_cases %>%
  filter(!str_detect(name, "_spec") & !str_detect(name, "aceae")) %>%
  anti_join(data, by = "name")

### Search in 'name' for certain names ###
data %>%
  group_by(name) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = TRUE))) %>%
  filter(str_detect(name, "Equisetum"))

data <- data %>%
  mutate(
    name = as_factor(name),
    name = fct_recode(
      name,
      Carex_praecox_ssp_praecox = "Carex_praecox",
      Cornus_controversa = "Cornus_sanguinea",
      Plantago_major_ssp_intermedia = "Plantago_major_subsp._intermedia",
      Plantago_major = "Plantago_major_subsp._major",
      Ranunculus_serpens_ssp_nemorosus = "Ranunculus_serpens_subsp._nemorosus",
      Silene_latifolia_ssp_alba = "Silene_latifolia_subsp._alba",
      Silene_latifolia_ssp_alba = "Silene_latifolia",
      Silene_latifolia_ssp_alba = "Silene_latifolia_ssp._alba"
    ),
    trait = as_factor(trait)
  )

### Check names of 'uncomplete_cases' which are not present in TRY database ###
uncomplete_cases %>%
  anti_join(data, by = "name") %>%
  filter(!str_detect(name, "_spec") & !str_detect(name, "aceae"))

### Show references ###
data %>%
  group_by(Reference) %>%
  count(sort = TRUE) %>%
  print(n = 10)

data <- data %>%
  group_by(name, trait) %>%
  summarise(across(where(is.double), ~median(.x, na.rm = TRUE))) %>%
  pivot_wider(names_from = "trait", values_from = "value")

traits <- traits %>%
  left_join(data, by = "name") %>%
  mutate(
    sla = coalesce(sla.x, sla.y),
    seedmass = coalesce(seedmass.x, seedmass.y),
    height = coalesce(height.x, height.y),
    .keep = "unused",
    name = as.character(name)
  ) %>%
  mutate(
    across(c(sla, seedmass, height), ~ round(.x, digits = 3))
    ) %>%
  arrange(name)

### * Check completeness of LEDA + TRY ####

traits %>%
  filter(!str_detect(name, "_spec") & !str_detect(name, "aceae")) %>%
  filter(group != "tree" & group != "shrub") %>%
  select(sla, seedmass, height) %>%
  miss_var_summary(order = FALSE)
traits %>%
  filter(!str_detect(name, "_spec") & !str_detect(name, "aceae")) %>%
  filter(group != "tree" & group != "shrub") %>%
  select(sla, seedmass, height) %>%
  vis_miss(cluster = TRUE, sort_miss = FALSE)
(uncomplete_cases <- traits %>%
    select(name, sla, seedmass, height) %>%
    filter(!complete.cases(sla, seedmass, height)) %>%
    filter(!str_detect(name, "_spec") & !str_detect(name, "aceae")))


### c Prepare data frames for calculation of functional diversity --------------

data <- traits %>%
  filter(!(str_detect(name, "_spec") | str_detect(name, "aceae"))) %>%
  filter(group != "tree" & group != "shrub")
data_lhs <- data %>%
  select(name, sla, seedmass, height) %>%
  drop_na()
data_sla <- data %>%
  select(name, sla) %>%
  drop_na()
data_seedmass <- data %>%
  select(name, seedmass) %>%
  drop_na()
data_height <- data %>%
  select(name, height) %>%
  drop_na()


### d Leaf-height-seed (LHS) functional diversity ------------------------------

### Preparation ###
data_species <- semi_join(species_experiment, data_lhs, by = "name")
data_traits <- semi_join(data_lhs, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, names_to = "site", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversity <- dbFD(
  log_data_traits, data_species,
  w.abun = TRUE, calc.FRic = FALSE, calc.FDiv = FALSE, corr = "cailliez"
)
### Integration to dataset ###
data <- data_diversity$FDis %>%
  as.data.frame() %>%
  rownames_to_column("id") %>%
  rename("fdis_abu_lhs" = ".") %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 3)))
sites_experiment <- sites_experiment %>%
  left_join(data, by = "id")


### e Specific leaf area (SLA) -------------------------------------------------

### Preparation ###
data_species <- semi_join(species_experiment, data_sla, by = "name")
data_traits <- semi_join(data_sla, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversity <- dbFD(
  log_data_traits, data_species,
  w.abun = TRUE, calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt"
)
### Integration to dataset ###
data <- data_diversity$FDis %>%
  as.data.frame() %>%
  add_column(data_diversity$CWM$sla) %>%
  rownames_to_column("id") %>%
  rename("fdis_abu_sla" = ".",
         "cwm_abu_sla" = "data_diversity$CWM$sla") %>%
  mutate(cwm_abu_sla = exp(cwm_abu_sla)) %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 3))) %>%
  select(id, cwm_abu_sla)
sites_experiment <- sites_experiment %>%
  left_join(data, by = "id")


### f Seed mass ----------------------------------------------------------------

### Preparation ###
data_species <- semi_join(species_experiment, data_seedmass, by = "name")
data_traits <- semi_join(data_seedmass, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, names_to = "site", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversity <- dbFD(
  log_data_traits, data_species,
  w.abun = TRUE, calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt"
)
### Integration to dataset ###
data <- data_diversity$FDis %>%
  as.data.frame() %>%
  add_column(data_diversity$CWM$seedmass) %>%
  rownames_to_column("id") %>%
  rename("fdis_abu_seedmass" = ".",
         "cwm_abu_seedmass" = "data_diversity$CWM$seedmass") %>%
  mutate(cwm_abu_seedmass = exp(cwm_abu_seedmass)) %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 3))) %>%
  select(id, cwm_abu_seedmass)
sites_experiment <- sites_experiment %>%
  left_join(data, by = "id")


### g Canopy height ------------------------------------------------------------

### Preparation ###
data_species <- semi_join(species_experiment, data_height, by = "name")
data_traits <- semi_join(data_height, data_species, by = "name")
data_species <- data_species %>%
  pivot_longer(-name, names_to = "site", values_to = "value") %>%
  pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames("site")
data_traits <- column_to_rownames(data_traits, "name")
log_data_traits <- log(data_traits)
### Calculation ###
data_diversity <- dbFD(
  log_data_traits, data_species,
  w.abun = TRUE, calc.FRic = FALSE, calc.FDiv = FALSE, corr = "sqrt"
)
### Integration to dataset ###
data <- data_diversity$FDis %>%
  as.data.frame() %>%
  add_column(data_diversity$CWM$height) %>%
  rownames_to_column("id") %>%
  rename("fdis_abu_height" = ".",
         "cwm_abu_height" = "data_diversity$CWM$height") %>%
  mutate(cwm_abu_height = exp(cwm_abu_height)) %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 3))) %>%
  select(id, cwm_abu_height)
sites_experiment <- sites_experiment %>%
  left_join(data, by = "id") %>%
  select(id, tidyselect::peek_vars())



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#### Data of experiment ###
write_csv(
  sites_experiment,
  here("data", "processed", "data_processed_sites_spatial.csv")
  )
write_csv(
  sites_temporal,
  here("data", "processed", "data_processed_sites_temporal.csv")
)
write_csv(
  sites_nmds,
  here("data", "processed", "data_processed_sites_nmds.csv")
)
write_csv(
  species_experiment,
  here("data", "processed", "data_processed_species.csv")
  )
write_csv(
  traits,
  here("data", "processed", "data_processed_traits.csv")
  )
