# Make Austria Halo Data using areaOfEffect package
# This script reads the full ASAAS dataset and extracts only the HALO points
# (points outside Austria but within the Area of Effect)

library(areaOfEffect)
library(data.table)
library(sf)

# Paths
data_dir <- "J:/Phd Local/Gilles_paper2/Data/ASAAS"
output_dir <- "C:/Users/Gilles Colling/Documents/liams stay/data"

# Read headers using data.table for speed (39M rows in species file)
cat("Reading header data...\n")
header <- fread(file.path(data_dir, "01_HEADER_ASAAS_ALL.csv"))
cat(sprintf("  %d plots loaded\n", nrow(header)))

# Filter to only valid coordinates (not dropped)
header <- header[DROP_COORD == FALSE]
cat(sprintf("  %d plots with valid coordinates\n", nrow(header)))

# Convert to sf object first (aoe expects sf or data.frame, not data.table)
# Use PlotObservationID as row names so aoe() uses them as point_id
cat("Converting to sf object...\n")
header_df <- as.data.frame(header)
row.names(header_df) <- header_df$PlotObservationID
header_sf <- st_as_sf(
  header_df,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE  # Keep coordinate columns
)

# Use areaOfEffect to classify points as core (inside Austria) or halo
cat("Classifying points using areaOfEffect...\n")
result <- aoe(
  points = header_sf,
  support = "Austria",
  mask = "land"
)

cat(sprintf("  %d points in AoE\n", nrow(result)))
cat("  Summary:\n")
print(table(result$aoe_class))

# Create lookup table for aoe_class by PlotObservationID
aoe_lookup <- data.table(
  PlotObservationID = as.numeric(result$point_id),
  aoe_class = result$aoe_class
)

# Get all IDs in AoE (both core and halo)
all_ids_num <- as.numeric(result$point_id)

# Filter header to AoE points and add aoe_class
header_aoe <- header[PlotObservationID %in% all_ids_num,
                     .(PlotObservationID, Year, Longitude, Latitude,
                       Eunis_lvl1, Eunis_lvl2)]
header_aoe <- merge(header_aoe, aoe_lookup, by = "PlotObservationID")

cat(sprintf("\nWriting %d header rows to austria_header_halo.csv...\n", nrow(header_aoe)))
cat(sprintf("  Core: %d, Halo: %d\n", sum(header_aoe$aoe_class == "core"), sum(header_aoe$aoe_class == "halo")))
fwrite(header_aoe, file.path(output_dir, "austria_header_halo.csv"))

# Now read and filter species data
cat("Reading species data (this may take a while)...\n")
species <- fread(file.path(data_dir, "02_PLOT_DATA_ASAAS_ALL.csv"))
cat(sprintf("  %d species records loaded\n", nrow(species)))

# Filter to AoE plots and add aoe_class
species_aoe <- species[PlotObservationID %in% all_ids_num,
                       .(PlotObservationID, WFO_TAXON, WFO_FAMILY, STATUS)]
species_aoe <- merge(species_aoe, aoe_lookup, by = "PlotObservationID")

cat(sprintf("Writing %d species rows to austria_species_halo.csv...\n", nrow(species_aoe)))
cat(sprintf("  Core: %d, Halo: %d\n", sum(species_aoe$aoe_class == "core"), sum(species_aoe$aoe_class == "halo")))
fwrite(species_aoe, file.path(output_dir, "austria_species_halo.csv"))

cat("\nDone!\n")
cat(sprintf("Header: %s\n", file.path(output_dir, "austria_header_halo.csv")))
cat(sprintf("Species: %s\n", file.path(output_dir, "austria_species_halo.csv")))
