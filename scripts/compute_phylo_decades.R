## Recompute phylogenetic diversity (MPD + MNTD) per decade
## - phylo_decades.rds       : pooled (all species), per decade
## - phylo_na_decades.rds    : native vs alien, per decade

library(spacc)
library(ape)
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
all_mat <- mats$all; native_mat <- mats$native; alien_mat <- mats$alien
rm(mats); invisible(gc())

tree <- readRDS(file.path(data_path, "phylo_tree.rds"))
coords <- data.frame(x = header$Longitude, y = header$Latitude)

# Match species to tree
in_tree        <- colnames(all_mat) %in% tree$tip.label
native_in_tree <- colnames(native_mat) %in% tree$tip.label
alien_in_tree  <- colnames(alien_mat) %in% tree$tip.label

cat(sprintf("Pooled: %d / %d species in tree (%.0f%%)\n",
            sum(in_tree), ncol(all_mat), mean(in_tree) * 100))
cat(sprintf("Native: %d / %d in tree\n", sum(native_in_tree), ncol(native_mat)))
cat(sprintf("Alien:  %d / %d in tree\n", sum(alien_in_tree), ncol(alien_mat)))

# Prune tree to matched species (union of all)
all_in_tree <- unique(c(colnames(all_mat)[in_tree],
                         colnames(native_mat)[native_in_tree],
                         colnames(alien_mat)[alien_in_tree]))
tree <- keep.tip(tree, all_in_tree)

# --- Decades -----------------------------------------------------------------
decades <- c("pre-1970", "1970s", "1980s", "1990s", "2000s", "2010s")

n_seeds <- 10
cat("\n=== Computing phylo per decade (pooled) ===\n")

phylo_decades <- list()
for (d in decades) {
  idx <- which(header$decade == d)
  n   <- length(idx)
  cat(sprintf("  %s: %d sites ... ", d, n))

  sub_mat <- all_mat[idx, in_tree, drop = FALSE]
  sub_coords <- coords[idx, ]

  # Remove species with 0 occurrences in this decade
  present <- colSums(sub_mat) > 0
  sub_mat <- sub_mat[, present, drop = FALSE]
  sub_tree <- keep.tip(tree, colnames(sub_mat))

  res <- spaccPhylo(sub_mat, sub_coords, tree = sub_tree,
                     metric = c("mpd", "mntd"),
                     n_seeds = n_seeds, method = "knn",
                     distance = "haversine", seed = 42)

  phylo_decades[[d]] <- list(
    mpd_mean  = colMeans(res$curves$mpd),
    mntd_mean = colMeans(res$curves$mntd),
    n_sites   = n
  )
  cat("done\n")
}

saveRDS(phylo_decades, file.path(comp_dir, "phylo_decades.rds"))
cat("Saved phylo_decades.rds\n")

# --- Native vs Alien per decade ----------------------------------------------
cat("\n=== Computing phylo per decade (native vs alien) ===\n")

phylo_na_decades <- list()
for (d in decades) {
  idx <- which(header$decade == d)
  n   <- length(idx)
  cat(sprintf("  %s: %d sites\n", d, n))

  # Native
  sub_nat <- native_mat[idx, native_in_tree, drop = FALSE]
  nat_present <- colSums(sub_nat) > 0
  sub_nat <- sub_nat[, nat_present, drop = FALSE]
  nat_tree <- keep.tip(tree, colnames(sub_nat))

  cat("    native ... ")
  nat_res <- spaccPhylo(sub_nat, coords[idx, ], tree = nat_tree,
                         metric = c("mpd", "mntd"),
                         n_seeds = n_seeds, method = "knn",
                         distance = "haversine", seed = 42)

  # Alien
  sub_ali <- alien_mat[idx, alien_in_tree, drop = FALSE]
  ali_present <- colSums(sub_ali) > 0
  sub_ali <- sub_ali[, ali_present, drop = FALSE]
  ali_tree <- keep.tip(tree, colnames(sub_ali))

  cat("done\n    alien ... ")
  ali_res <- spaccPhylo(sub_ali, coords[idx, ], tree = ali_tree,
                         metric = c("mpd", "mntd"),
                         n_seeds = n_seeds, method = "knn",
                         distance = "haversine", seed = 42)

  phylo_na_decades[[d]] <- list(
    native = list(
      mpd_mean  = colMeans(nat_res$curves$mpd),
      mntd_mean = colMeans(nat_res$curves$mntd),
      n_sites   = n,
      n_species = ncol(sub_nat)
    ),
    alien = list(
      mpd_mean  = colMeans(ali_res$curves$mpd),
      mntd_mean = colMeans(ali_res$curves$mntd),
      n_sites   = n,
      n_species = ncol(sub_ali)
    )
  )
  cat("done\n")
}

saveRDS(phylo_na_decades, file.path(comp_dir, "phylo_na_decades.rds"))
cat("Saved phylo_na_decades.rds\n\nAll done!\n")
