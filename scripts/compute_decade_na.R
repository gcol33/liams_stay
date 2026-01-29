## Compute species accumulation per decade, split by native/alien
## Output: decade_na_results.rds

library(spacc)
library(data.table)

comp_dir  <- file.path(dirname(getwd()), "comp")
data_path <- file.path(dirname(getwd()), "data")

# --- Load data ---------------------------------------------------------------
header  <- fread(file.path(data_path, "austria_header.csv"))
species <- fread(file.path(data_path, "austria_species.csv"))

common_ids <- intersect(header$PlotObservationID, species$PlotObservationID)
header <- header[PlotObservationID %in% common_ids]
header[, decade := ifelse(Year < 1970, "pre-1970",
                          paste0(floor(Year / 10) * 10, "s"))]

mats    <- readRDS(file.path(comp_dir, "matrices.rds"))
native_mat <- mats$native; alien_mat <- mats$alien
rm(mats); invisible(gc())

coords <- data.frame(x = header$Longitude, y = header$Latitude)

# --- Decades -----------------------------------------------------------------
decades <- c("pre-1970", "1970s", "1980s", "1990s", "2000s", "2010s")

n_seeds <- 50
cat("=== Computing SAC per decade (native vs alien) ===\n")

decade_na_results <- list()
for (d in decades) {
  idx <- which(header$decade == d)
  n   <- length(idx)
  cat(sprintf("  %s: %d sites\n", d, n))

  sub_native <- native_mat[idx, , drop = FALSE]
  sub_alien  <- alien_mat[idx, , drop = FALSE]
  sub_coords <- coords[idx, ]

  # Remove species with 0 occurrences in this decade
  nat_present <- colSums(sub_native) > 0
  ali_present <- colSums(sub_alien) > 0
  sub_native <- sub_native[, nat_present, drop = FALSE]
  sub_alien  <- sub_alien[, ali_present, drop = FALSE]

  cat("    native ... ")
  nat_res <- spacc(sub_native, sub_coords,
                   n_seeds = min(n_seeds, n), method = "knn",
                   distance = "haversine", seed = 42)

  cat("done\n    alien ... ")
  ali_res <- spacc(sub_alien, sub_coords,
                   n_seeds = min(n_seeds, n), method = "knn",
                   distance = "haversine", seed = 42)

  decade_na_results[[d]] <- list(
    native = nat_res,
    alien  = ali_res,
    n_sites = n
  )
  cat("done\n")
}

saveRDS(decade_na_results, file.path(comp_dir, "decade_na_results.rds"))
cat("\nSaved decade_na_results.rds\n")
