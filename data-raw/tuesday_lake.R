# Tuesday Lake data for the vignette "tuesday-lake" (inst/extdata/tuesday_lake)
#
# Source: the communities TL84 and TL86 of the R package cheddar (version
# 0.1-639; Copyright (c) 2014-2016 Lawrence Hudson, Dan Reuman and Rob
# Emerson; BSD 2-Clause licence, reproduced in inst/COPYRIGHTS), from
# Jonsson T, Cohen JE, Carpenter SR (2005). Food webs, body size, and species
# abundance in ecological community description. Advances in Ecological
# Research 36, 1-84, and Cohen JE, Jonsson T, Carpenter SR (2003). Ecological
# community description using the food web, species abundance, and body size.
# PNAS 100, 1781-1786. Tuesday Lake (Michigan, USA), pelagic community sampled
# in 1984 and 1986. No mass of Tuesday Lake is in the reference data of
# bmassa.
#
# Files written (CSV, UTF-8):
#   nodes.csv  one row per year and node: year, node (name in cheddar),
#              category (producer, invertebrate, vert.ecto), family (as in
#              cheddar), M_kg (mean wet mass of an individual, kg; per algal
#              unit in the phytoplankton, a cell or a colony, from the volume
#              with a density of 1 kg/L) and N_m3 (individuals per m3)
#   taxa.csv   the classification of the taxon units of the vignette in the
#              Catalogue of Life, from standardize_taxa() with the names
#              cleaned as in the vignette (needs internet)
#
# Use, from the package root, in a UTF-8 locale and with internet:
#   LC_ALL=C.UTF-8 Rscript data-raw/tuesday_lake.R [folder with TL84.RData and TL86.RData]
# Without a folder, the data are taken from the installed package cheddar.
# Packages: pkgload, curl and jsonlite (for standardize_taxa()).

suppressMessages(pkgload::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
out <- file.path("inst", "extdata", "tuesday_lake")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

load_web <- function(w) {
  e <- new.env()
  if (length(args)) load(file.path(args[1], paste0(w, ".RData")), envir = e)
  else utils::data(list = w, package = "cheddar", envir = e)
  get(w, envir = e)$nodes
}
blank <- function(x) ifelse(is.na(x) | trimws(x) == "", NA_character_, trimws(x))
nodes <- do.call(rbind, lapply(c(TL84 = 1984L, TL86 = 1986L), function(y) {
  n <- load_web(if (y == 1984L) "TL84" else "TL86")
  data.frame(year = y, node = n$node, category = n$category, family = blank(n$family),
             M_kg = n$M, N_m3 = n$N, stringsAsFactors = FALSE)
}))
rownames(nodes) <- NULL
utils::write.csv(nodes, file.path(out, "nodes.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# The cleaning of the names, as in the vignette.
clean_name <- function(x) {
  x <- sub("-spiny", "", x, fixed = TRUE)
  x <- sub(" \\((solitary|colonial)\\)$", "", x)
  x <- sub(" sp\\. ?[0-9]*$", "", x)
  x[grepl("^Unclassified", x)] <- NA_character_
  # A misspelling of Cyclops varicans rubellus Lilljeborg, 1901, a synonym of
  # Microcyclops rubellus in the GBIF Backbone Taxonomy.
  x[x == "Cyclops varians rubellus"] <- "Microcyclops rubellus"
  x
}
units <- nodes[!duplicated(nodes$node), c("node", "category", "family")]
units$taxon_id <- units$node
units$broad_group <- c(producer = "phytoplankton", invertebrate = "zooplankton",
                       vert.ecto = "fish")[units$category]
units$scientific_name <- clean_name(units$node)
named <- units[!is.na(units$scientific_name), c("taxon_id", "broad_group", "scientific_name",
                                                "family")]
taxa <- standardize_taxa(named, quiet = TRUE)
utils::write.csv(taxa, file.path(out, "taxa.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("%d nodes (%d taxon units, %d named); COL release %s\n", nrow(nodes), nrow(units),
            nrow(taxa), unique(taxa$col_release)[1]))
print(taxa[!is.na(taxa$note), c("taxon_id", "col_name", "note")])
