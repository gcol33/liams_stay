# Prepare Austrian subset from ASAAS data for Liam's exercise
# Run this ONCE before the exercise to create the filtered data files

library(tidyverse)
library(data.table)

# Paths
asaas_path <- "J:/Phd Local/Gilles_paper2/Data/ASAAS/"
output_path <- "../data/"

cat("Loading header data...\n")
header <- fread(paste0(asaas_path, "01_HEADER_ASAAS_ALL.csv"))

cat("Filtering to Austria (ISO_EVA == 'AT')...\n
")
# Filter to Austria only
austria_header <- header[ISO_EVA == "AT"]
cat(paste("Austrian plots:", nrow(austria_header), "\n"))

# Keep only essential columns for the exercise
austria_header <- austria_header[, .(
  PlotObservationID,
  Year,
  Longitude,
  Latitude,
  Eunis_lvl1,
  Eunis_lvl2,
  EVA_Country
)]

cat("Saving austria_header.csv...\n")
fwrite(austria_header, paste0(output_path, "austria_header.csv"))

# Get the plot IDs we need
austria_plot_ids <- austria_header$PlotObservationID

cat("Loading species data (this may take a while)...\n")
species <- fread(paste0(asaas_path, "02_PLOT_DATA_ASAAS_ALL.csv"))

cat("Filtering species to Austrian plots...\n")
austria_species <- species[PlotObservationID %in% austria_plot_ids]
cat(paste("Austrian species records:", nrow(austria_species), "\n"))

# Keep only essential columns
austria_species <- austria_species[, .(
  PlotObservationID,
  WFO_TAXON,
  WFO_FAMILY,
  STATUS
)]

# Filter to native and neo only (remove NA or other statuses)
austria_species <- austria_species[STATUS %in% c("native", "neo")]

cat("Saving austria_species.csv...\n")
fwrite(austria_species, paste0(output_path, "austria_species.csv"))

# Summary
cat("\n=== Summary ===\n")
cat(paste("Plots:", nrow(austria_header), "\n"))
cat(paste("Species records:", nrow(austria_species), "\n"))
cat(paste("Unique species:", n_distinct(austria_species$WFO_TAXON), "\n"))

cat("\nStatus breakdown:\n")
print(austria_species[, .N, by = STATUS])

cat("\nDone! Files saved to:", output_path, "\n")
