library(data.table)
h <- fread("data/austria_header.csv")
s <- fread("data/austria_species.csv")
ci <- intersect(h$PlotObservationID, s$PlotObservationID)
h <- h[PlotObservationID %in% ci]
cat("Header rows:", nrow(h), "\n")

mats <- readRDS("comp/matrices.rds")
cat("Matrix rows:", nrow(mats$all), "\n")

rn <- as.integer(rownames(mats$all))
cat("First 5 rownames:", paste(head(rn, 5), collapse = ", "), "\n")
cat("First 5 header IDs:", paste(head(h$PlotObservationID, 5), collapse = ", "), "\n")
cat("All match:", all(rn == h$PlotObservationID), "\n")
