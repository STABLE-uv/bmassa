# Raw data of the azoxystrobin mesocosm experiment (inst/extdata/azoxystrobin)
#
# Source: the data workbook of Amador P, Vega C, Navarro Pacheco NI,
# Moratalla-López J, Palacios J, Crettaz Minaglia MC, López I, Díaz M, Rico A
# (2024). Effects of the fungicide azoxystrobin in two habitats
# representative of Mediterranean coastal wetlands: a mesocosm experiment.
# Aquatic Toxicology 267, 106828. https://doi.org/10.1016/j.aquatox.2023.106828
# File `Amador_et_al_2024.xlsx`, supplied by the first author (2026-11-07),
# MD5 86e35bdab33b30c446653bb75561b138. Licence of the data in the package:
# CC BY 4.0 (the authors' choice), with the citation of the article.
#
# The experiment: 24 outdoor mesocosms (140 cm in diameter, 80 cm of water,
# 1230 L) of two habitats of the Albufera Natural Park (Valencia, Spain),
# eutrophic lake (code W in the workbook) and water spring (S), each with
# controls and three concentrations of azoxystrobin (2, 20 and 200 µg/L,
# applied twice, 14 days apart), three replicates each; plankton and
# macroinvertebrates sampled on days -7, 7, 21, 35 and 56 relative to the
# first application. The habitat of each code follows the conductivity of the
# sheet `fq` (W, the lake with added salt, about 1380 µS/cm; S about 860).
#
# Files written (CSV, UTF-8), raw as in the workbook, one row per sample:
#   design.csv              one row per mesocosm: tank, treatment, habitat,
#                           concentration_ug_L, replicate, volume_L, depth_m,
#                           area_m2
#   phytoplankton.csv       Day, Treatment, Concentration and cells/L of each
#                           taxon (units of the metadata of the workbook)
#   zooplankton.csv         the same, individuals/L
#   macroinvertebrates.csv  the same, individuals in the two traps (30 x 10 cm
#                           each) of each mesocosm and date, together
# The workbook has no tank identifiers: within each day, treatment and
# concentration the three rows are the three replicates, always in the same
# order (between dates and sheets), which `design.csv` numbers 1 to 3.
#
# Changes to the workbook, all in the macroinvertebrates: the column labelled
# "Palaemonetes zariquieyi" is the shrimp Dugastella valentina (as in the
# article), and is named so; the Spanish label "Dytiscidae (larva pequeña)"
# becomes "Dytiscidae (small larva)"; an unlabelled last column (filled only
# in the water spring mesocosms, with a single individual on day -7 in a
# control) is left out. Names are otherwise as in the workbook (with their
# misspellings), without leading or trailing spaces.
#
# taxa.csv: the classification of the taxon units of the vignette
# "azoxystrobin" in the Catalogue of Life, from standardize_taxa() (needs
# internet; run with BMASSA_AZOXY_TAXA=true), with the names cleaned as in the
# vignette (the code below is that of the vignette, which checks that its
# names are those of this file).
#
# Packages: readxl (to read the workbook). Run from the package root in a UTF-8
# locale (LC_ALL=C.UTF-8 Rscript data-raw/azoxystrobin.R).

xlsx <- Sys.getenv("BMASSA_AZOXY_XLSX",
                   "/mnt/user-data/uploads/Azoxy_fw/data/raw/Amador_et_al_2024.xlsx")
stopifnot(unname(tools::md5sum(xlsx)) == "86e35bdab33b30c446653bb75561b138")
out <- file.path("inst", "extdata", "azoxystrobin")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_sheet <- function(sheet) {
  d <- as.data.frame(readxl::read_excel(xlsx, sheet = sheet, .name_repair = "minimal"))
  names(d) <- trimws(names(d))
  d <- d[, nzchar(names(d)), drop = FALSE]
  stopifnot(identical(names(d)[1:3], c("Day", "Treatment", "Concentration")),
            !anyNA(d), !anyDuplicated(names(d)))
  d
}
phyto <- read_sheet("phyto")
zoo <- read_sheet("zoo")
macro <- read_sheet("minvertebrates")
names(macro)[names(macro) == "Palaemonetes zariquieyi"] <- "Dugastella valentina"
names(macro)[names(macro) == "Dytiscidae (larva pequeña)"] <- "Dytiscidae (small larva)"
stopifnot(all(c("Dugastella valentina", "Dytiscidae (small larva)") %in% names(macro)))

# The same samples, in the same order, in the three sheets.
key <- function(d) paste(d$Day, d$Treatment, d$Concentration)
stopifnot(identical(key(phyto), key(zoo)), identical(key(phyto), key(macro)),
          all(table(key(phyto)) == 3L))

# Habitat of each code, from the conductivity.
fq <- as.data.frame(readxl::read_excel(xlsx, sheet = "fq"))
ec <- tapply(fq$EC, fq$Treatment, mean)
stopifnot(ec[["W"]] > ec[["S"]])

design <- expand.grid(replicate = 1:3, concentration_ug_L = c(0, 2, 20, 200),
                      treatment = c("W", "S"), stringsAsFactors = FALSE)
design <- data.frame(
  tank = sprintf("%s-%g-%d", design$treatment, design$concentration_ug_L, design$replicate),
  treatment = design$treatment,
  habitat = ifelse(design$treatment == "W", "eutrophic_lake", "water_spring"),
  concentration_ug_L = design$concentration_ug_L, replicate = design$replicate,
  volume_L = 1230, depth_m = 0.8, area_m2 = round(pi * 0.7^2, 3), stringsAsFactors = FALSE)

write_raw <- function(d, file) {
  utils::write.csv(d, file.path(out, file), row.names = FALSE, fileEncoding = "UTF-8")
}
write_raw(design, "design.csv")
write_raw(phyto, "phytoplankton.csv")
write_raw(zoo, "zooplankton.csv")
write_raw(macro, "macroinvertebrates.csv")
cat(sprintf("%s: %d rows, %d taxa\n", c("phytoplankton", "zooplankton", "macroinvertebrates"),
            c(nrow(phyto), nrow(zoo), nrow(macro)), c(ncol(phyto), ncol(zoo), ncol(macro)) - 3L),
    sep = "")

# Classification of the taxon units of the vignette (needs internet).
if (identical(Sys.getenv("BMASSA_AZOXY_TAXA"), "true")) {
  suppressMessages(pkgload::load_all(".", quiet = TRUE))
  # From here to standardize_taxa(), the code of the vignette.
  to_long <- function(d, broad_group) {
    taxa_cols <- setdiff(names(d), c("Day", "Treatment", "Concentration"))
    replicate <- stats::ave(seq_len(nrow(d)), d$Day, d$Treatment, d$Concentration, FUN = seq_along)
    k <- length(taxa_cols)
    data.frame(day = rep(d$Day, k), treatment = rep(d$Treatment, k),
               concentration = rep(d$Concentration, k), replicate = rep(replicate, k),
               broad_group = broad_group, raw_name = rep(taxa_cols, each = nrow(d)),
               abundance = unlist(d[taxa_cols], use.names = FALSE))
  }
  long <- rbind(to_long(phyto, "phytoplankton"), to_long(zoo, "zooplankton"),
                to_long(macro, "macroinvertebrate"))
  units <- unique(long[c("broad_group", "raw_name")])
  rownames(units) <- NULL
  units$taxon_id <- gsub("^_+|_+$", "", gsub("[^a-z0-9]+", "_", tolower(units$raw_name)))
  nm <- units$raw_name
  units$life_stage <- ifelse(grepl("\\(adult\\)$", nm), "adult",
                      ifelse(grepl("larva", nm), "larva",
                      ifelse(grepl(" pupa$", nm), "pupa",
                      ifelse(nm == "Nauplii", "nauplius", NA_character_))))
  nm <- sub(" \\(adult\\)$| larvae$| \\(small larva\\)$| pupa$| \\(lineal\\)$", "", nm)
  nm <- sub("(_| )spp?\\.?( [0-9]+)?( \\(small\\))?$", "", nm)
  fixes <- c(Dyctiosphaerium = "Dictyosphaerium", Ankystrodesmus = "Ankistrodesmus",
             Chroccoccus = "Chroococcus", Pseudonitzschia = "Pseudo-nitzschia",
             "Coccal < 5 um" = "Chlorophyta", Chidorus = "Chydorus", Leidigia = "Leydigia",
             Ciclopoida = "Cyclopoida", Nauplii = "Copepoda",
             "Brachionus havanensis" = "Brachionus havanaensis",
             Orthocladineae = "Orthocladiinae", "Physella Acuta" = "Physella acuta",
             Hydracarina = "Trombidiformes")
  fix <- nm %in% names(fixes)
  nm[fix] <- fixes[nm[fix]]
  units$scientific_name <- nm
  units <- units[!units$raw_name %in% c("Vorticella", "Coccal colony", "Pleurogyra"), ]
  hints <- data.frame(
    raw_name = c("Zygoptera", "Anisoptera", "Ciclopoida", "Calanoida", "Hirudinea"),
    phylum = c("Arthropoda", "Arthropoda", "Arthropoda", "Arthropoda", "Annelida"),
    class = c("Insecta", "Insecta", "Copepoda", "Copepoda", "Clitellata"),
    order = c("Odonata", "Odonata", "Cyclopoida", "Calanoida", NA))
  k <- match(units$raw_name, hints$raw_name)
  for (col in c("phylum", "class", "order")) units[[col]] <- hints[[col]][k]
  stopifnot(!anyDuplicated(units$taxon_id))
  taxa <- standardize_taxa(units[c("taxon_id", "broad_group", "scientific_name", "phylum",
                                   "class", "order")], quiet = TRUE)
  utils::write.csv(taxa, file.path(out, "taxa.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("taxa.csv: %d taxon units; COL release %s\n", nrow(taxa), unique(taxa$col_release)[1]))
  print(taxa[!is.na(taxa$note), c("taxon_id", "col_name", "note")])
}
