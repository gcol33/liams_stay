# Compute beta diversity decomposition for native and alien species separately
library(spacc)
library(data.table)

comp_dir  <- normalizePath("../comp", winslash = "/")
data_path <- normalizePath("../data", winslash = "/")

# Load matrices
mats <- readRDS(file.path(comp_dir, "matrices.rds"))
native_mat <- mats[["native"]]
alien_mat  <- mats[["alien"]]

# Load header
header <- fread(file.path(data_path, "austria_header.csv"))
cat("Header rows:", nrow(header), "\n")

# Keep only plots present in both
common_ids <- intersect(header[["PlotObservationID"]], as.integer(rownames(native_mat)))
header <- header[PlotObservationID %in% common_ids]

# Decade assignment
header[, decade := ifelse(Year < 1970, "pre-1970",
                          paste0(floor(Year / 10) * 10, "s"))]

# 1990s subset
idx_90s <- which(header[["decade"]] == "1990s")
coords  <- header[, .(x = Longitude, y = Latitude)]

cat("1990s plots:", length(idx_90s), "\n")
cat("Native species:", ncol(native_mat), "\n")
cat("Alien species:", ncol(alien_mat), "\n")

# Compute beta for native
cat("Computing native beta...\n")
beta_native <- spaccBeta(native_mat[idx_90s, , drop = FALSE], coords[idx_90s, ],
                         n_seeds = 30, index = "sorensen",
                         distance = "haversine", seed = 42)
cat("Native beta done\n")

# Compute beta for alien
cat("Computing alien beta...\n")
beta_alien <- spaccBeta(alien_mat[idx_90s, , drop = FALSE], coords[idx_90s, ],
                        n_seeds = 30, index = "sorensen",
                        distance = "haversine", seed = 42)
cat("Alien beta done\n")

# Save
out_path <- file.path(comp_dir, "beta_native_alien.rds")
saveRDS(list(native = beta_native, alien = beta_alien), out_path)
cat("Saved:", out_path, "\n")
