# Grassland experiment on dikes
# Prepare project metadata ####

# Markus Bauer
# 2023-02-28


### Packages ###
library(here)
library(tidyverse)
library(EML)
library(emld)
#remotes::install_github("ropenscilabs/emldown", build = FALSE)
library(emldown)
#remotes::install_github("EDIorg/EMLassemblyline")
library(EMLassemblyline)

### Start ###
rm(list = ls())



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Collect metadata ##########################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### 1 Methods and units #######################################################


methods_file <- here("data", "text", "methods.odt")
methods <- set_methods(methods_file)

EMLassemblyline::view_unit_dictionary()
# List of standard units, which should be used in metadata file



### 2 Raw data ################################################################

### a data_raw_species  -------------------------------------------------------

### b data_raw_traits  --------------------------------------------------------

### c data_raw_sites  ---------------------------------------------------------

### 3 Processed data ##########################################################

### a data_processed_species  -------------------------------------------------

### b data_processed_traits  --------------------------------------------------

### c data_processed_sites  ---------------------------------------------------

### 4 Put data table together #################################################



### 5 Contact #################################################################


address <- list(
  deliveryPoint = "Emil-Ramann-Strasse 6",
  city = "Freising",
  administrativeArea = "Bayern",
  postalCode = "85354",
  country = "Germany"
  )

creator <- eml$creator(
  individualName = eml$individualName(
    givenName = "Markus",
    surName = "Bauer"
  ),
  positionName = "PhD student",
  organizationName = "Technical University of Munich",
  address = address,
  electronicMailAddress = "markusbauer@mailbox.org",
  phone = "0049-152-56391781",
  id = "https://orcid.org/0000-0001-5372-4174"
)

associatedParty <- list(
  eml$associatedParty(
    individualName = eml$individualName(
      givenName = "Jakob",
      surName = "Huber"
    ),
    role = "Researcher",
    organizationName = "Technical University of Munich",
    electronicMailAddress = "jakob.huber@posteo.de"
  ),
  eml$associatedParty(
    individualName = eml$individualName(
      givenName = "Johannes",
      surName = "Kollmann"
    ),
    role = "Professor",
    organizationName = "Technical University of Munich",
    address = address,
    electronicMailAddress = "johannes.kollmann@.tum.de",
    phone = "0049-8161-714144",
    id = "https://orcid.org/0000-0002-4990-3636"
  )
)

contact <- list(
  individualName = creator$individualName,
  electronicMailAddress = creator$electronicMailAddress,
  address = address,
  organizationName = "Technical University of Munich",
  onlineUrl = "DOI address to the database"
)



### 6 Temporal and spatial coverage ###########################################


geographicDescription <- "Danube dike near Deggendorf"

coverage <- set_coverage(
  begin = "2018-06-01", end = "2021-07-31",
  sci_names = list(list(
    Subdivision = "Spermatophytina"
  )),
  geographicDescription = geographicDescription,
  west = 12.87704, east = 12.89251,
  north = 48.84048, south = 48.83844,
  altitudeMin = 313, altitudeMaximum = 315,
  altitudeUnits = "meter"
)



### 7 Description #############################################################


pubDate <- "2023"

title <- "Danube old dikes"

abstract <- "Not written yet"

keywordSet <- list(
  list(
    keywordThesaurus = "LTER controlled vocabulary",
    keyword = list(
      "rivers",
      "vegetation dynamics",
      "restoration"
      )
  ),
  list(
    keywordThesaurus = "own vocabulary",
    keyword = list(
      "beta diversity",
      "temperate grassland",
      "dike"
      )
  )
)

intellectualRights <-
  "CC-BY-4.0: https://creativecommons.org/licenses/by/4.0/deed.en"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B finalize EML ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



dataset <- list(
  title = title,
  pubDate = pubDate,
  creator = creator,
  associatedParty = associatedParty,
  intellectualRights = intellectualRights,
  abstract = abstract,
  keywordSet = keywordSet,
  coverage = coverage,
  contact = contact,
  methods = methods#,
  #dataTable = dataTable,
  #additonalMetadata = list(metadata = list(unitList = unitList))
  )

eml <- list(
  packageId = uuid::UUIDgenerate(),
  system = "uuid", # type of identifier
  dataset = dataset
  )

write_eml(eml, here("METADATA.xml"))
eml_validate(here("METADATA.xml"))

#render_eml(file = here("METADATA.xml"), outfile = "METADATA.html",
#           open = TRUE, publish_mode = TRUE)
