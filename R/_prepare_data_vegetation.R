# Grassland experiment on dikes
# Prepare species, sites, and traits data ####
# Markus Bauer
# 2022-05-24



# Content #####################################################################

# A Load data +++++++++++++++++++++++++++++++++++++
## 1 Sites and sPlotOpen
## 2 Species and sPlotOpen
## 3 Traits
## 4 Temperature and precipitation
## 5 Check data frames

# B Create variables ++++++++++++++++++++++++++++++
## 1 Create simple variables
## 2 Coverages
## 3 Alpha diversity
## 4 sPlotOpen
## 5 TBI: Temporal beta diversity index

# C Save processed data +++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Packages ###
library(here)
library(tidyverse)
suppressPackageStartupMessages(library(lubridate))
library(naniar) #are_na
library(vegan)
suppressPackageStartupMessages(library(adespatial))
suppressPackageStartupMessages(library(FD))
#remotes::install_github(file.path("inbo", "checklist"))

### Start ###
rm(list = ls())
setwd(here("data", "raw"))
#installr::updateR(browse_news = FALSE, install_R = TRUE, copy_packages = TRUE, copy_Rprofile.site = TRUE, keep_old_packages = TRUE, update_packages = TRUE, start_new_R = TRUE, quit_R = TRUE, print_R_versions = TRUE, GUI = FALSE)
#checklist::setup_source()
#checklist::check_source()
#renv::status()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Load data ##################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#______________________________________________________________________________
## 1 Sites ####################################################################


sites_experiment <- read_csv("data_raw_sites.csv", col_names = TRUE,
                  na = c("", "NA", "na"),
                  col_types =
                    cols(
                      .default = "f",
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
                      )) %>%
  pivot_longer(
    starts_with("vegetation_cover") |
      starts_with("botanist") |
      starts_with("biomass") |
      starts_with("survey_date"),
    names_to = c("x", "survey_year"),
    names_sep = "\\.",
    values_to = "n",
    values_transform = list (n = as.character)
  ) %>%
  pivot_wider(names_from = "x", values_from = "n") %>%
  mutate(plot = str_replace(plot, "-", "_"),
         plot = str_replace(plot, "L_", "L"),
         plot = str_replace(plot, "W_", "W"),
         id = str_c(plot, survey_year, sep = "_"),
         plot = factor(plot),
         id = factor(id),
         vegetation_cover = as.numeric(vegetation_cover),
         biomass = as.numeric(biomass)) %>%
  filter(!(site == "C" & (survey_year == "seeded" |
                             survey_year == "2018" |
                             survey_year == "2019" |
                             survey_year == "2020" |
                            survey_year == "2021")))

### Sabatini et al. (2021) Global Ecol Biogeogr:
### https://doi.org/10.1111/geb.13346
sites_splot <- read_delim(here("data", "raw", "sabatini_etal_2021",
                               "sPlotOpen_header.txt"),
                          col_names = TRUE, na = c("", "NA", "na"),
                          col_types = cols(
                            .default = "?",
                            Cover_algae_layer = "d"
                          ))

### Bauer et al. (2022) Zenodo:
### https://doi.org/10.5281/zenodo.6334100
sites_bauer <- read_csv(here("data", "raw", "bauer_etal_2022",
                               "data_sites_bauer_etal_2022.csv"),
                          col_names = TRUE, na = c("", "NA", "na"),
                          col_types = cols(
                            .default = "?"
                          ))




#______________________________________________________________________________
## 2 Species ###################################################################


species_experiment <- data.table::fread("data_raw_species_20211112.csv",
                             sep = ",",
                             dec = ".",
                             skip = 0,
                             header = TRUE,
                             na.strings = c("", "NA", "na"),
                             colClasses = list(
                               character = "name"
                             )) %>%
  ### Check that each species occurs at least one time ###
  group_by(name) %>%
  arrange(name) %>%
  select(name, all_of(sites_experiment$id)) %>%
  mutate(total = sum(c_across(
    starts_with("L") | starts_with("W")),
    na.rm = TRUE),
    presence = if_else(total > 0, 1, 0)) %>%
  # filter only species which occur at least one time:
  filter(presence == 1) %>%
  ungroup() %>%
  select(name, sort(tidyselect::peek_vars()), -total, -presence) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0)))

### Sabatini et al. (2021) Global Ecol Biogeogr:
### https://doi.org/10.1111/geb.13346
species_splot <- read_delim(here("data", "raw", "sabatini_etal_2021",
                                 "sPlotOpen_DT.txt"),
                            col_names = TRUE, na = c("", "NA", "na"), col_types =
                              cols(
                                .default = "?"
                              )) %>%
  filter(Abundance_scale == "CoverPerc")

### Bauer et al. (2022) Zenodo:
### https://doi.org/10.5281/zenodo.6334100
species_bauer <- read_csv(here("data", "raw", "bauer_etal_2022",
                             "data_species_bauer_etal_2022.csv"),
                        col_names = TRUE, na = c("", "NA", "na"),
                        col_types = cols(
                          .default = "?"
                        ))

### Create list with species names and their frequency ###
specieslist <- species_experiment %>%
  mutate(across(where(is.numeric), ~1 * (. != 0))) %>%
  mutate(sum = rowSums(across(where(is.numeric)), na.rm = TRUE),
         .keep = "unused") %>%
  group_by(name) %>%
  summarise(sum = sum(sum))
#write_csv(specieslist, here("outputs", "tables", "specieslist_20221102.csv"))



#______________________________________________________________________________
## 3 Seedmixes #################################################################


seedmixes <- data.table::fread("data_raw_species_20211112.csv",
                             sep = ",",
                             dec = ".",
                             skip = 0,
                             header = TRUE,
                             na.strings = c("", "NA", "na"),
                             colClasses = list(
                               character = "name"
                             )) %>%
  arrange(name) %>%
  select(name, ends_with("seeded")) %>%
  filter(name %in% species_experiment$name) %>%
  select(name, sort(tidyselect::peek_vars())) %>%
  mutate(across(where(is.numeric), ~replace(., is.na(.), 0))) %>%
  pivot_longer(-name, names_to = "id", values_to = "seeded") %>%
  separate(id, c("plot"), sep = "_(?!.*_)",
           remove = TRUE, extra = "drop", fill = "warn", convert = FALSE)

 

#______________________________________________________________________________
## 4 Traits ####################################################################


traits <- read_csv("data_raw_traits.csv", col_names = TRUE,
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
                     )) %>%
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



#______________________________________________________________________________
## 5 Check data frames #########################################################


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

rm(list = setdiff(ls(), c("sites_experiment", "sites_splot", "sites_bauer",
                          "species_experiment", "species_splot", "sites_bauer",
                          "traits", "seedmixes")))




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Create variables ###########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#______________________________________________________________________________
## 1 Create variables #########################################################


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
  mutate(surveyDate = as_date(survey_date),
         seeding_date = if_else(
           exposition == "north", ymd("20180413"), ymd("20180427")
           ),
         age = interval(seeding_date, survey_date) %/% days(1),
         block = str_c(site, exposition, sep = "_"),
         block = factor(block)
         )

### Establishment ###
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
  mutate(rate = total_established / total_seeded,
         rate = round(rate, digits = 2),
         seeded = "1") %>%
  pivot_wider(names_from = "survey_year",
              values_from = c("rate", "total_seeded", "total_established")) %>%
  select(-total_seeded_2019, -total_seeded_2020, -total_seeded_2021)

### Check size of species pool for seeding ###
traits %>%
  select(seeded_hay_meadow, seeded_dry_grassland) %>%
  mutate(across(where(is.character), ~as.numeric(.x))) %>%
  summarise(across(where(is.numeric), ~sum(.x)))

traits <- traits %>%
  left_join(data, by = "name")



#______________________________________________________________________________
## 2 Coverages #################################################################


cover <- species_experiment %>%
  left_join(traits, by = "name") %>%
  select(name, family, target, seeded,
         starts_with("L"), starts_with("W"), starts_with("C")) %>%
  pivot_longer(
    names_to = "id",
    values_to = "n",
    cols = starts_with("L") |starts_with("W") | starts_with("C")
    ) %>%
  group_by(id)

#### Graminoid, herb, and total coverage) ###
cover_total_and_graminoid <- cover %>%
  group_by(id, family) %>%
  summarise(total = sum(n, na.rm = TRUE), .groups = "keep") %>%
  mutate(type = if_else(
    family == "Poaceae" |
      family == "Cyperaceae" |
      family == "Juncaceae",
    "graminoid_cover", "herb_cover")) %>%
  group_by(id, type) %>%
  summarise(total = sum(total, na.rm = TRUE), .groups = "keep") %>%
  spread(type, total) %>%
  mutate(accumulated_cover = graminoid_cover + herb_cover,
         accumulated_cover = round(accumulated_cover, 1)) %>%
  ungroup()

#### Target specis' coverage ###
cover_target <- cover %>%
  filter(target == "yes") %>%
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
  summarise(seededCov = sum(success, na.rm = TRUE), .groups = "keep") %>%
  ungroup() %>%
  unite(id, plot, survey_year, sep = "_")

#### Implement in sites data set ###
sites_experiment <- sites_experiment %>%
  left_join(cover_total_and_graminoid, by = "id") %>%
  left_join(cover_target, by = "id") %>%
  left_join(cover_seeded, by = "id") %>%
  ### Calcute the ratio of target richness of total species richness
  mutate(
    target_cover_ratio = target_cover / accumulated_cover,
    graminoid_cover_ratio = graminoid_cover / accumulated_cover,
    seeded_cover_ratio = seeded_cover / accumulated_cover,
    target_cover_ratio = round(target_cover_ratio, 3),
    graminoid_cover_ratio = round(graminoid_cover_ratio, 3),
    seeded_cover_ratio = round(seeded_cover_ratio, 3)
    )

rm(list = setdiff(ls(), c("sites_experiment", "sites_splot", "sites_bauer",
                          "species_experiment", "species_splot", "sites_bauer",
                          "traits", "seedmixes")))



#______________________________________________________________________________
## 3 Alpha diversity ##########################################################


### a Species richness ---------------------------------------------------------

speciesRichness <- species %>%
  left_join(traits, by = "name") %>%
  select(name, rlg, rlb, target, ffh6510, ffh6210,
         starts_with("L"), starts_with("W"), starts_with("C")) %>%
  pivot_longer(
    names_to = "id", values_to = "n",
    cols = starts_with("L") | starts_with("W") | starts_with("C")
    ) %>%
  mutate(n = if_else(n > 0, 1, 0)) %>%
  group_by(id)

### * total species richness ####
speciesRichness_all <- speciesRichness %>%
  summarise(speciesRichness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * red list Germany (species richness) ####
speciesRichness_rlg <- speciesRichness %>%
  filter(rlg == "1" | rlg == "2" | rlg == "3" | rlg == "V") %>%
  summarise(rlgRichness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * red list Bavaria (species richness) ####
speciesRichness_rlb <- speciesRichness %>%
  filter(rlb == "1" | rlb == "2" | rlb == "3" | rlb == "V") %>%
  summarise(rlbRichness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * target species (species richness) ####
speciesRichness_target <- speciesRichness %>%
  filter(target == "yes") %>%
  summarise(targetRichness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * seeded species (species richness) ####
speciesRichness_seeded <- species %>%
  select(-starts_with("C"), -contains("x0809")) %>%
  # Make two columns out of column id:
  pivot_longer(-name, names_to = "id", values_to = "value",
               values_drop_na = TRUE) %>%
  separate(id, c("plot", "surveyYear"), sep = "_(?!.*_)",
           remove = FALSE, extra = "merge", fill = "warn", convert = FALSE) %>%
  # Summarise for plot:
  pivot_wider(names_from = "surveyYear", values_from = "value") %>%
  group_by(plot, name) %>%
  summarise(across(where(is.double), ~sum(.x, na.rm = TRUE)),
            .groups = "keep") %>%
  ungroup() %>%
  # Combine with seedmixes:
  pivot_longer(starts_with("20"), names_to = "surveyYear", values_to = "value",
               names_transform = list(surveyYear = as.factor)) %>%
  mutate(successCov = if_else(seeded > 0 & value > 0, 1, 0)) %>%
  group_by(plot, surveyYear) %>%
  summarise(seededRichness = sum(successCov, na.rm = TRUE),
            .groups = "keep") %>%
  ungroup() %>%
  unite(id, plot, surveyYear, sep = "_")

### * ffh6510 species (species richness) ####
speciesRichness_ffh6510 <- speciesRichness %>%
  filter(ffh6510 == "yes") %>%
  summarise(ffh6510Richness = sum(n, na.rm = TRUE)) %>%
  ungroup()

### * ffh6210 species (species richness) ####
speciesRichness_ffh6210 <- speciesRichness %>%
  filter(ffh6210 == "yes") %>%
  summarise(ffh6210Richness = sum(n, na.rm = TRUE)) %>%
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

  ### * Calculate Favourable Conservation Status (FCS) ####
  ### Helm et al. 2014 Divers Distrib
  ### https://doi.org/10.1111/ddi.12285

  mutate(
    fcs_target = log(
      (targetRichness + 1) / (speciesRichness - targetRichness + 1)
      ),
    fcs_seeded = log(
      (seededRichness + 1 ) / (speciesRichness - seededRichness + 1)
      ),
    fcs_target = round(fcs_target, 3),
    fcs_seeded = round(fcs_seeded, 3)
    )

### b Species eveness and shannon ----------------------------------------------

data <- species  %>%
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

rm(list = setdiff(ls(), c("sites", "species", "traits", "seedmixes")))



#______________________________________________________________________________
## 4 Reference sites ##########################################################


### a ESy: EUNIS expert vegetation classification system ----------------------

#### Start ###
### Bruelheide et al. 2021 Appl Veg Sci
### https://doi.org/10.1111/avsc.12562

expertfile <- "EUNIS-ESy-2020-06-08.txt" ### file of 2021 is not working

obs <- species %>%
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

header <- sites %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 31468) %>%
  sf::st_transform(4326) %>%
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
    dataset = "Danube_dikes"
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
eval.EUNIS(which(result.classification == "V39")[1], "V39")

sites_dikes <- sites_dikes %>%
  mutate(
    esy = result.classification,
    esy = if_else(id == "X05_m_2021", "R1A", esy),
    esy = if_else(id == "X62_m_2019", "R", esy),
    esy = if_else(id == "X67_o_2021", "R", esy)
  )
table(sites_dikes$esy)
rm(list = setdiff(ls(), c("sites_dikes", "sites_splot",
                          "species_dikes", "species_splot",
                          "traits", "pca_soil", "pca_construction_year",
                          "pca_survey_year", "result.classification")))

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
      Releve_area >= 10 &
      Releve_area <= 40 &
      Longitude > 10.89845 & # West: Augsburg
      Longitude < 13.46434 & # East: Passau
      Latitude > 	47.85298 & # South: Rosenheim
      Latitude < 49.45095 & # North: Nuernberg
      Elevation < 700
  ) %>%
  rename_with(tolower) %>%
  rename(id = plotobservationid, survey_year = date_of_recording,
         plotSize = releve_area, reference = country) %>%
  mutate(
    id = paste0("X", id),
    reference = str_replace(reference, "Germany", "reference"),
    survey_year = year(survey_year),
    givd_database = if_else(
      givd_id == "EU-DE-014",
      "Jandt & Bruelheide (2012) https://doi.org/10.7809/b-e.00146",
      "other"
    )
  ) %>%
  select(id, givd_id, longitude, latitude, elevation, plotSize, survey_year,
         reference, esy)
sites_splot <- data_sites

data_species <- species_splot %>%
  rename(id = PlotObservationID, name = Species,
         abundance = Original_abundance) %>%
  mutate(id = paste0("X", id)) %>%
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
species_splot <- data_species

### Check species name congruency ###
data <- anti_join(species_splot, traits, by = "name") %>%
  select(name) %>%
  print(n = 50)

rm(list = setdiff(ls(), c("sites_dikes", "sites_splot",
                          "species_dikes", "species_splot",
                          "traits", "pca_soil", "pca_construction_year",
                          "pca_survey_year")))



#______________________________________________________________________________
## 5 Temporal beta diversity ###################################################


### Prepare data ###
data_sites <- sites %>%
  # Choose only plots which were surveyed in each year:
  filter(accumulatedCov > 0) %>%
  add_count(plot) %>%
  filter(n == max(n)) %>%
  select(id, plot)
data_species <- species %>%
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
  
  assign(nam, data_species %>%
           filter(year == i) %>%
           select(-year) %>%
           column_to_rownames(var = "plot"))
}

### a Calculate TBI Presence --------------------------------------------------

### Legendre (2019) Ecol Evol
### https://doi.org/10.1002/ece3.4984

#### * seedmixes vs. 2018 ####
res18 <- TBI(species_seeded, species_2018,
             method = "sorensen",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res18$BCD.summary # B = .566, C = .308, D = .874 (.647 vs. .352)
res18$t.test_B.C # p.perm = 1e-04
tbi18 <- res18$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2018")
#### Test plot
plot(res18, type = "BC")

#### * seedmixes vs. 2019 ####
res19 <- TBI(species_seeded, species_2019,
             method = "sorensen",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res19$BCD.summary # B = .372, C = .361, D = .733 (.507 vs. .492)
res19$t.test_B.C # p.perm = .422
tbi19 <- res19$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2019")
#### Test plot
plot(res19, type = "BC")

#### * seedmixes vs. 2020 ####
res20 <- TBI(species_seeded, species_2020,
             method = "sorensen",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res20$BCD.summary # B = .318, C = .396, D = .715 (.445 vs. .554)
res20$t.test_B.C # p.perm = 1e-04
tbi20 <- res20$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2020")
#### Test plot
plot(res20, type = "BC")

#### * seedmixes vs. 2021 ####
res21 <- TBI(species_seeded, species_2021,
             method = "sorensen",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res21$BCD.summary # B = .343, C = .394, D = .738 (.465 vs. .534)
res21$t.test_B.C # p.perm = 1e-04
tbi21 <- res21$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2021")
#### Test plot
plot(res21, type = "BC")

#### * Combine datasets ####
data_presence <- bind_rows(tbi18, tbi19, tbi20, tbi21) %>%
  mutate(presabu = "presence")

rm(list = setdiff(ls(), c(
  "sites", "species", "traits", "seedmixes", "species_seeded",
  "species_2018", "species_2019", "species_2020", "species_2021",
  "data_presence", "data_sites", "data_species"
)))

### b Calculate TBI Abundance -------------------------------------------------

#### * seedmixes vs. 2018 ####
res18 <- TBI(species_seeded, species_2018,
             method = "%difference",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res18$BCD.summary # B = .654, C = .318, D = .973 (.672 vs. .327)
res18$t.test_B.C # p.perm = 1e-04
tbi18 <- res18$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2018")
#### Test plot
plot(res18, type = "BC")

#### * seedmixes vs. 2019 ####
res19 <- TBI(species_seeded, species_2019,
             method = "%difference",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res19$BCD.summary # B = .519, C = .386, D = .905 (.573 vs. .426)
res19$t.test_B.C # p.perm = .1e-04
tbi19 <- res19$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2019")
#### Test plot
plot(res19, type = "BC")

#### * seedmixes vs. 2020 ####
res20 <- TBI(species_seeded, species_2020,
             method = "%difference",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res20$BCD.summary # B = .484, C = .385, D = .870 (.557 vs. .442)
res20$t.test_B.C # p.perm = 1e-04
tbi20 <- res20$BCD.mat %>%
  as_tibble() %>%
  mutate(comparison = "2020")
#### Test plot
plot(res20, type = "BC")

#### * seedmixes vs. 2021 ####
res21 <- TBI(species_seeded, species_2021,
             method = "%difference",
             nperm = 9999, test.t.perm = TRUE, clock = TRUE)
res21$BCD.summary # B = .470, C = .392, D = .862 (.544 vs. .455)
res21$t.test_B.C # p.perm = 1e-04
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
sites <- sites %>%
  left_join(data, by = "id")

rm(list = setdiff(ls(), c(
  "sites", "species", "traits", "seedmixes"
)))



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



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# C Save processed data ########################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


write_csv(sites, here("data", "processed", "data_processed_sites.csv"))
write_csv(traits, here("data", "processed", "data_processed_traits.csv"))
write_csv(species, here("data", "processed", "data_processed_species.csv"))
