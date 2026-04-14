# ============================================================
# Title: Data Import and Preprocessing Pipeline
# Author: Dr. Carlo Andrea Cossu
# Description:
# Loads raw microplastics dataset, performs cleaning,
# feature engineering, and builds analysis-ready datasets.
# ============================================================

# ---- 1. Load packages ----
library(data.table)
library(readxl)
library(dplyr)
library(RSQLite)

# ---- 2. Set working directory ----
wd <- "project-root-directory"
setwd(wd)

# ---- 3. Load external scripts ----
source("App/data_prep/packages.R")
source("App/data_prep/shapefiles.R")

# ---- 4. Load raw dataset ----
mps <- as.data.table(
  read_xlsx("data/raw/MPs_results.xlsx")
)

# ---- 5. Basic transformations ----
mps[, m_feret := as.numeric(m_feret)]
mps[, concentration := mps_nr / weight]

# ---- 6. Summary metadata ----
animals_n <- mps[, .(n_animals = uniqueN(animal_id)), by = animal_species]

grams_n <- mps[
  , .(n_grams = round(sum(unique(weight,
       by = c("animal_species", "sample_type"))), 1)),
  by = c("animal_species", "sample_type")
]

mps <- merge(mps, animals_n, by = "animal_species", all.x = TRUE)
mps <- merge(mps, grams_n, by = c("animal_species", "sample_type"), all.x = TRUE)

mps[, animal_species_label := paste0(animal_species, " (", n_grams, ")")]

# ---- 7. Connection with database of original project ----
database_name <- "PhD_db_240113.sqlite"

data_retrieve <- function(string, string_tail, conditions = NULL) {
  wd <- "project-root-directory"
  setwd(wd)

  mydb <- dbConnect(RSQLite::SQLite(), file.path(getwd(), database_name))

  if (is.null(conditions)) {
    data <- as.data.table(dbGetQuery(mydb, string))
  } else {
    data <- as.data.table(dbGetQuery(mydb, paste0(string, " ", string_tail),
                                     params = conditions))
  }

  dbDisconnect(mydb)
  data
}

# ---- 8. Load database components ----
source("App/data_prep/database_parts/database_animals.R")

class(mps$animal_id) <- "integer"

# ---- 9. Merge spatial data ----
mps <- merge(
  mps,
  unique(animal_samples[, .(animal_id, latitude, longitude)], by = "animal_id"),
  by = "animal_id",
  all.x = TRUE
)

# ---- 10. MP type aggregation ----
mps <- mps %>%
  group_by(sample_id, mps_type) %>%
  mutate(
    count_mps_type = n(),
    concentration_by_mps_type = count_mps_type / weight
  ) %>%
  ungroup()

setDT(mps)

# ---- 11. Handle missing MP types ----
mps[mps_type == "-", `:=`(
  count_mps_type = 0,
  concentration_by_mps_type = 0
)]

# ---- 12. Create analysis datasets ----
mps_sample <- unique(mps, by = "sample_id")
mps_by_type <- unique(mps, by = c("sample_id", "mps_type"))
