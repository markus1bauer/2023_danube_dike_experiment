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



## 1 Methods and units #######################################################


#methods_file <- here("data", "text", "methods.odt")
#methods <- set_methods(methods_file)

#EMLassemblyline::view_unit_dictionary()
# List of standard units, which should be used in metadata file



## 2 Raw data ################################################################

## 3 Processed data ##########################################################

## 4 Put data table together #################################################



## 5 Contact #################################################################


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
      givenName = "Jakob Kaspar",
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
  onlineUrl = "https://www3.ls.tum.de/roek/mitarbeiter-in/prof-dr-johannes-kollmann/"
)



## 6 Temporal and spatial coverage ###########################################


geographic_description <- "Danube dike near Deggendorf"

coverage <- set_coverage(
  begin = "2018-06-01", end = "2021-07-31",
  sci_names = list(list(
    Subdivision = "Spermatophytina"
  )),
  geographicDescription = geographic_description,
  west = 12.87704, east = 12.89251,
  north = 48.84048, south = 48.83844,
  altitudeMin = 313, altitudeMaximum = 315,
  altitudeUnits = "meter"
)



## 7 Description #############################################################


title <- "Fit by design: Developing substrate-specific seed mixtures for functional dike grasslands"

pubDate <- "2023"

alternate_identifier <- "https://doi.org/10.5281/zenodo.7713396"

abstract <- "A multifactorial experiment was set up on a dike of the Danube 
River in SE Germany. The aim of the experiment was to find suitable seed-
substrate combinations. Four treatments were tested on both dike slopes: nutrient 
reduction of agricultural soil by sand admixture (0, 25, 50%), substrate depths 
(15, 30 cm), seed density (4, 8 g / m²), and seed mixture type (mesic hay 
meadow, semi-dry calcareous grasslands. The 48 treatment combinations were 
replicated six times resulting in 288 plots of a size of 2 x 3 m². The 
experiment was established in April 2018 and observed until 2021. The 
vegetation surveys according to Braun-Blanquet were conducted using the Londo 
scale."

keyword_set <- list(
  list(
    keywordThesaurus = "LTER controlled vocabulary",
    keyword = list(
      "droughts",
      "grasslands",
      "meadows",
      "monitoring",
      "permanent plots",
      "plant communities",
      "restoration",
      "rivers",
      "soil samples",
      "species composition",
      "vegetation"
      )
  ),
  list(
    keywordThesaurus = "own vocabulary",
    keyword = list(
      "experiment",
      "dike",
      "temperate grassland"
      )
  )
)

license <- list(
  licenseName = "CC-BY-4.0",
  url = "https://creativecommons.org/licenses/by/4.0/deed.en"
)

short_name <- "Experiment on Danube dike"

language <- "English"



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B finalize EML ##############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



dataset <- list(
  title = title,
  shortName = short_name,
  pubDate = pubDate,
  creator = creator,
  associatedParty = associatedParty,
  licensed = license,
  alternateIdentifier = alternate_identifier,
  abstract = abstract,
  keywordSet = keyword_set,
  coverage = coverage,
  language = language,
  contact = contact#,
  #methods = methods,
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

#render_eml(
#  file = here("METADATA.xml"), outfile = "METADATA.html",
#  open = TRUE, publish_mode = TRUE
#  )
