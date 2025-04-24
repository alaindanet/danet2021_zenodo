#

library(tidyverse)
library(magrittr)
library(here)
library(furrr)

# Get raw data
dir.create(here("data-raw"))
zenodo_zip_file <- here("data-raw", "danet2021zenodo.zip")
download.file(
  "https://zenodo.org/api/records/5095656/files-archive",
  destfile = zenodo_zip_file
)
unzip(zenodo_zip_file, exdir = here("data-raw"))
list.files(here("data-raw"))

# Load data
community_data <- read_csv(here("data-raw", "community_data.csv")) %>%
  mutate(
    opcod = as.integer(opcod),
    nind = as.integer(nind)
  )
station <- read_csv(here("data-raw", "station_basin.csv"),
  col_types = list(col_integer(), col_character(), col_double(), col_double()))
fishing_protocol <- read_csv(here("data-raw", "fishing_protocol.csv")) %>%
  mutate(across(c(opcod, station, nb_sp, nb_ind, year), as.integer))
load("class_network.rda")

# Associate fish sampling id (opcod) to station id (station), protocol, etc...
min_data_sampling <- fishing_protocol %>%
  select(station, opcod, protocol, date, year)
community_data %<>%
  left_join(min_data_sampling) %>%
  filter(!is.na(station))
class_network %<>%
  left_join(min_data_sampling) %>%
  filter(!is.na(station))


# Compute species cv by station
## Using custom function from my post-doc
source("https://raw.githubusercontent.com/alaindanet/fishcom/refs/heads/master/R/synchrony.R")
sync_cv <- get_sync_cv_mat(
  com_analysis = community_data %>% select(opcod, species, biomass),
  op_analysis = min_data_sampling,
  presence_threshold = 0.1
)

## Only station and species_cv
species_cv <- sync_cv %>%
  select(station, species_cv = avg_sp) %>%
  mutate(species_cv = map(species_cv, ~enframe(.x, name = "species", value = "cv"))) %>%
  unnest(cols = species_cv)

# Compute network metrics
var_chr <- "sp_class"
## Compute obs_troph_lvl
network_analysis <- class_network %>%
  mutate(
    network = future_map(network, igraph::graph_from_data_frame, directed = TRUE),
    network = future_map(network, igraph::as_adjacency_matrix, sparse = FALSE),
    troph = future_map(network, NetIndices::TrophInd)
  )

network_analysis %>%
  select(opcod, troph) %>%
  unnest(troph) %>%
  filter(is.na(TL))

## get obs_troph_level
network_analysis %<>%
  mutate(
    obs_troph_level = map(troph, function (x) {
      out <- tibble(
        !!sym(var_chr) := row.names(x),
        obs_troph_level = x$TL
      )
      return(out)
})
  )
## Put obs_troph_level in composition
network_analysis %<>%
  mutate(
  composition = furrr::future_map2(composition, obs_troph_level, function (compo, troph_group, var2join){

    if ("obs_troph_level" %in% names(compo)) {
      compo %<>% dplyr::select(-obs_troph_level)
    }
    left_join(compo, troph_group, by = var2join)
}, var2join = var_chr))

## Compute trophic level by species by taking their avg trophic level weighted
##by their biomass:
composition <-
  network_analysis %>%
  select(-obs_troph_level) %>%
  unnest(composition) %>%
  mutate(species = str_extract(sp_class, "[A-Z]{3}")) %>%
  group_by(opcod, species) %>%
  summarise(
    # Average species trophic level by biomass:
    troph_level_std = sum(bm_std * obs_troph_level) / sum(bm_std),
    troph_level = sum(biomass * obs_troph_level) / sum(biomass),
    bm_std = sum(bm_std),
    biomass = sum(biomass)
    ) %>%
  ungroup()

composition %<>%
  group_by(opcod) %>%
  nest(composition = !opcod)

network_analysis %<>%
  select(-composition) %>%
  left_join(., composition, by = "opcod")

# Check that everything is ok for Ismael
## You should have to play with:
#  1. the temporal cv of each species by station
species_cv %>%
  filter(station == 122)
#  2. The average trophic level of each species by sampling event (you might
#  need to summarise by station to match the temporal CV)
network_analysis %>%
  select(opcod, station, date, composition) %>%
  unnest(composition) %>%
  filter(station == 122)
# 3. You might want to play with community data directly to compute
# self-regulation and so
community_data %>%
  filter(station == 122)
