# Example data set of bmassa: example_survey.
#
# Run from the package root:  Rscript data-raw/examples.R
# Packages used only here: usethis.
#
# These tables were constructed for illustration (taxa typical of a shallow
# eutrophic Mediterranean lake; abundances and lengths chosen to be
# plausible). They are NOT field data and must not be used as such.

# The input tables of estimate_biomass(), for two sites of the
# same lake: the taxon units (a taxon, or a taxon in one life stage and sex, with
# its counting unit and, if measured, its mean length), the abundances of each
# unit by site, and the lengths measured in a particular site (not field data). The
# classification follows the Catalogue of Life (COL26.6 XR), with the orders
# that the COL XR leaves empty (Cyclopoida, Hygrophila) from other
# classifications, and Squalius placed in Leuciscidae, as in other
# classifications. `col_id` is the identifier of each taxon in the COL26.6 XR,
# as returned by standardize_taxa() on 2026-11-02 (all exact matches).
lineage <- function(...) {
  x <- c(...)
  stopifnot(length(x) == 7L)
  x
}
tx <- rbind(
  lineage("Chromista", "Ochrophyta", "Bacillariophyceae", "Thalassiosirales", "Stephanodiscaceae", "Cyclotella", "Cyclotella meneghiniana"),
  lineage("Bacillati", "Cyanobacteriota", "Cyanophyceae", "Chroococcales", "Microcystaceae", "Microcystis", "Microcystis aeruginosa"),
  lineage("Chromista", "Cryptista", "Cryptophyceae", "Cryptomonadales", "Cryptomonadaceae", "Cryptomonas", NA),
  lineage("Bacillati", "Cyanobacteriota", "Cyanophyceae", "Oscillatoriales", "Microcoleaceae", "Planktothrix", "Planktothrix agardhii"),
  lineage("Plantae", "Chlorophyta", "Chlorophyceae", "Sphaeropleales", "Scenedesmaceae", "Scenedesmus", "Scenedesmus quadricauda"),
  lineage("Animalia", "Rotifera", "Eurotatoria", "Ploima", "Brachionidae", "Keratella", "Keratella cochlearis"),
  lineage("Animalia", "Rotifera", "Eurotatoria", "Ploima", "Brachionidae", "Brachionus", "Brachionus calyciflorus"),
  lineage("Animalia", "Arthropoda", "Branchiopoda", "Diplostraca", "Daphnidae", "Daphnia", "Daphnia magna"),
  lineage("Animalia", "Arthropoda", "Copepoda", "Cyclopoida", "Cyclopidae", "Acanthocyclops", "Acanthocyclops robustus"),
  lineage("Animalia", "Arthropoda", "Insecta", "Diptera", "Chironomidae", "Chironomus", "Chironomus plumosus"),
  lineage("Animalia", "Arthropoda", "Insecta", "Ephemeroptera", "Baetidae", "Cloeon", "Cloeon dipterum"),
  lineage("Animalia", "Arthropoda", "Insecta", "Odonata", "Coenagrionidae", "Ischnura", "Ischnura elegans"),
  lineage("Animalia", "Arthropoda", "Malacostraca", "Amphipoda", "Gammaridae", "Gammarus", NA),
  lineage("Animalia", "Mollusca", "Gastropoda", "Hygrophila", "Physidae", "Physella", "Physella acuta"),
  lineage("Animalia", "Annelida", "Clitellata", "Tubificida", "Naididae", "Limnodrilus", "Limnodrilus hoffmeisteri"),
  lineage("Animalia", "Chordata", "Teleostei", "Cypriniformes", "Cyprinidae", "Luciobarbus", "Luciobarbus graellsii"),
  lineage("Animalia", "Chordata", "Teleostei", "Cypriniformes", "Leuciscidae", "Squalius", "Squalius pyrenaicus"),
  lineage("Animalia", "Chordata", "Teleostei", "Cyprinodontiformes", "Poeciliidae", "Gambusia", "Gambusia holbrooki"),
  lineage("Animalia", "Chordata", "Amphibia", "Anura", "Ranidae", "Pelophylax", "Pelophylax perezi")
)
colnames(tx) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
# The taxon units: the frog Pelophylax perezi twice, as adults (t19) and as
# larvae (t20).
tx <- rbind(tx, tx[nrow(tx), ])
units <- data.frame(
  taxon_id = sprintf("t%02d", seq_len(nrow(tx))),
  broad_group = c(rep("phytoplankton", 5), rep("zooplankton", 4), rep("macroinvertebrate", 6),
                  rep("fish", 3), "amphibian", "amphibian"),
  scientific_name = ifelse(is.na(tx[, "species"]), tx[, "genus"], tx[, "species"]),
  rank = ifelse(is.na(tx[, "species"]), "genus", "species"),
  tx,
  col_id = c("MYH9J", "CP3G2", "3WNK", "4JK5W", "SWLMR", "3R5TG", "5X2DW", "6CCSV", "VKDSB",
             "8BGDB", "5ZC6V", "6N9HG", "7P34G", "4HGZ8", "6QDFZ", "3WCKH", "6ZDS2", "3F7VR",
             "76CHM", "76CHM"),
  life_stage = c(rep(NA, 8), "adult", "larva", NA, "larva", NA, NA, NA, NA, NA, NA, "adult",
                 "larva"),
  sex = NA_character_,
  count_unit = c(rep("cell", 5), rep("individual", 15)),
  # Mean lengths of the units in both sites (NA: not measured).
  length_mean = c(rep(NA, 7), 2.2, NA, 12, 6, 16, NA, NA, NA, 15, NA, NA, 58, NA),
  length_unit = c(rep(NA, 7), "mm", NA, "mm", "mm", "mm", NA, NA, NA, "cm", NA, NA, "mm", NA),
  length_type = c(rep(NA, 7), "body_length", NA, "body_length", "body_length", "body_length",
                  NA, NA, NA, "fork_length", NA, NA, "snout_vent_length", NA),
  length_sd = c(rep(NA, 9), 3, NA, 3, NA, NA, NA, 4, NA, NA, NA, NA),
  length_cv = c(rep(NA, 7), 0.25, rep(NA, 10), 0.15, NA),
  stringsAsFactors = FALSE
)
example_taxa <- units
ab_one <- data.frame(
  taxon_id = sprintf("t%02d", 1:20),
  per = c(rep("mL", 5), rep("L", 4), rep("m2", 6), rep("100 m2", 5)),
  stringsAsFactors = FALSE
)
example_abundance <- data.frame(
  site = rep(c("L1", "L2"), each = nrow(ab_one)),
  rbind(ab_one, ab_one),
  abundance = c(
    # L1: cells per mL; individuals per L, per m2 and per 100 m2
    1800, 25000, 300, 9000, 600, 350, 120, 4, 25, 420, 55, 12, 80, 25, 300, 3, 12, 40, 4, 30,
    # L2: no Scenedesmus seen (0); Daphnia not counted (NA)
    750, 41000, 120, 16000, 0, 410, 260, NA, 12, 260, 30, 6, 150, 40, 520, 5, 8, 65, 2, 55
  ),
  stringsAsFactors = FALSE
)
example_abundance <- example_abundance[, c("site", "taxon_id", "abundance", "per")]
# Lengths measured in each site: the chub Squalius pyrenaicus (t17).
example_length_table <- data.frame(
  taxon_id = c("t17", "t17"),
  site = c("L1", "L2"),
  length_mean = c(110, 85),
  length_unit = "mm",
  length_type = "fork_length",
  length_sd = c(20, 15),
  stringsAsFactors = FALSE
)
example_survey <- list(abundance = example_abundance, taxa = example_taxa,
                       lengths = example_length_table)

usethis::use_data(example_survey, overwrite = TRUE)
