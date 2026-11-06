# Taxonomic normalization of the reference names with the Catalogue of Life
# eXtended Release (COL XR), through the GBIF name-matching API v2 (version 0.2.0)
#
# Builds, in data-raw/derived/:
#   colxr_matches.csv    one row per query sent to the API (name and the
#                        classification hints sent with it) and the parsed
#                        answer; it is also the cache: queries already in it
#                        are not sent again
#   taxa_colxr.csv       one row per source and name of the reference tables,
#                        with the COL taxon chosen, its classification and flags
#   colxr_version.csv    version of the COL XR used by the matching service
#                        when the queries were sent
#   taxa_colxr_log.md    counts, checks and the names to review
# from the names of the package data of 0.1.0 (typical_mass, length_mass, with
# the lineages of `taxonomy`) and of the tables in data-raw/derived/, with the
# classification given by each source as hints, and the corrections in
# data-raw/taxa_manual.csv.
#
# Internet (api.gbif.org) is needed only for queries not yet in
# colxr_matches.csv, or for all of them with BMASSA_COLXR_REFRESH = "true".
# Run from the package root, after the other scripts of data-raw/, in a UTF-8
# locale:  Rscript data-raw/taxonomy_colxr.R
# Packages used only here: curl, jsonlite, readxl.

source("data-raw/sources.R")

out_dir <- "data-raw/derived"
reg <- read_sources()
checklist_key <- "7ddf754f-d193-4cc9-b351-99906754a03b" # COL XR in GBIF
api <- "https://api.gbif.org/v2/species/match"
cache_file <- file.path(out_dir, "colxr_matches.csv")
refresh <- tolower(Sys.getenv("BMASSA_COLXR_REFRESH", "false")) == "true"
log_lines <- c("# Taxonomic normalization with the Catalogue of Life (COL XR)", "",
               paste0("Built by `data-raw/taxonomy_colxr.R` on ", format(Sys.Date()), " (",
                      R.version.string, ")."), "")
say <- function(...) log_lines <<- c(log_lines, sprintf(...))
na_if_empty <- function(v) {
  v <- trimws(as.character(v))
  v[!is.na(v) & (!nzchar(v) | v %in% c("-", "--", "?"))] <- NA_character_
  v
}
most_frequent <- function(v) {
  v <- v[!is.na(v)]
  if (length(v)) names(sort(table(v), decreasing = TRUE))[1] else NA_character_
}

# Phyla of animals as named in the COL, and names used by the sources for them.
animal_phyla <- c("Acanthocephala", "Annelida", "Arthropoda", "Brachiopoda", "Bryozoa",
                  "Chaetognatha", "Chordata", "Cnidaria", "Ctenophora", "Cycliophora",
                  "Echinodermata", "Entoprocta", "Gastrotricha", "Gnathostomulida",
                  "Hemichordata", "Kinorhyncha", "Loricifera", "Micrognathozoa", "Mollusca",
                  "Nematoda", "Nematomorpha", "Nemertea", "Onychophora", "Phoronida",
                  "Placozoa", "Platyhelminthes", "Porifera", "Priapulida", "Rotifera",
                  "Sipuncula", "Tardigrada", "Xenacoelomorpha")
phylum_alias <- c(Crustacea = "Arthropoda", Insecta = "Arthropoda", Chelicerata = "Arthropoda",
                  Osteichthyes = "Chordata", Chondrichthyes = "Chordata", Agnatha = "Chordata",
                  Tunicata = "Chordata", Cephalochordata = "Chordata", Nemata = "Nematoda",
                  Plathelminthes = "Platyhelminthes", Priapula = "Priapulida")
col_kingdoms <- c("Animalia", "Plantae", "Chromista", "Protozoa", "Fungi")

# 1. Names of the reference tables, with the classification of their source ----

hint_cols <- c("kingdom", "phylum", "class", "order", "family")
names_of <- function(source_id, taxon, kingdom = NA, phylum = NA, class = NA, order = NA,
                     family = NA, domain) {
  x <- data.frame(source_id = source_id, taxon = na_if_empty(taxon),
                  kingdom = na_if_empty(kingdom), phylum = na_if_empty(phylum),
                  class = na_if_empty(class), order = na_if_empty(order),
                  family = na_if_empty(family), domain = domain, stringsAsFactors = FALSE)
  x <- x[!is.na(x$taxon), ]
  # One hint set per source and name: the most frequent value of each rank.
  do.call(rbind, lapply(split(x, paste(x$source_id, x$taxon), drop = TRUE), function(d) {
    d[1, hint_cols] <- lapply(d[, hint_cols, drop = FALSE], most_frequent)
    d[1, ]
  }))
}
nm <- list()

# Package data of 0.1.0: lineages compiled in `taxonomy` (Kremer et al. 2014,
# Moreau et al. 2021, van Leeuwen et al. 2025, Dumont et al. 1975, WoRMS).
pkg <- new.env()
for (f in c("typical_mass", "length_mass", "taxonomy")) {
  # Tables of 0.1.0 (data-raw/references.R).
  assign(f, read.csv(file.path(out_dir, paste0(f, "_0.1.0.csv")), stringsAsFactors = FALSE,
                     encoding = "UTF-8"), envir = pkg)
}
old <- unique(rbind(pkg$typical_mass[, c("source", "taxon")], pkg$length_mass[, c("source", "taxon")]))
lin <- pkg$taxonomy[match(old$taxon, pkg$taxonomy$name), ]
nm$old <- names_of(old$source, old$taxon, lin$kingdom, lin$phylum, lin$class, lin$order,
                   lin$family, domain = ifelse(lin$kingdom %in% "Animalia", "animal", "phytoplankton"))

d_tr <- read.csv(file.path(out_dir, "typical_records.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
d_lf <- read.csv(file.path(out_dir, "lwr_fitted.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
d_lp <- read.csv(file.path(out_dir, "lwr_published.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
d_cv <- read.csv(file.path(out_dir, "conversions.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
taxa_of <- function(id) {
  unique(c(d_tr$taxon[d_tr$source_id == id], d_lf$taxon[d_lf$source_id == id],
           d_lp$taxon[d_lp$source_id == id], d_cv$taxon[d_cv$source_id == id]))
}

# Hrycik et al. (2024): phylum from the group of SpeciesList.csv (typical
# masses) or from Table 3 (conversion factors).
if (use_source("hrycik2024", reg)) {
  group_phylum <- c(Amphipoda = "Arthropoda", Chironomidae = "Arthropoda", Dreissena = "Mollusca",
                    Gastropoda = "Mollusca", Hirudinea = "Annelida", Hydrozoa = "Cnidaria",
                    Isopoda = "Arthropoda", Mysis = "Arthropoda", Nemertea = "Nemertea",
                    Oligochaeta = "Annelida", OtherInsecta = "Arthropoda", Insecta = "Arthropoda",
                    Lumbriculidae = "Annelida", Platyhelminthes = "Platyhelminthes",
                    Polychaeta = "Annelida", Sphaeriidae = "Mollusca")
  sl <- read.csv(source_files("hrycik2024", reg)[3], stringsAsFactors = FALSE)
  tr_h <- d_tr[d_tr$source_id == "hrycik2024", c("taxon", "name_as_given")]
  ph <- unname(group_phylum[sl$Group[match(tr_h$name_as_given, sl$Name)]])
  ph[is.na(ph)] <- unname(group_phylum[tr_h$taxon[is.na(ph)]])
  t3 <- read.csv("data-raw/hrycik2024_table3.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
  cv_h <- d_cv[d_cv$source_id == "hrycik2024", c("taxon", "note")]
  row <- as.integer(sub("^Table 3, rows? ([0-9]+).*$", "\\1", cv_h$note))
  nm$hrycik <- names_of("hrycik2024", c(tr_h$taxon, cv_h$taxon), "Animalia",
                        c(ph, t3$phylum[match(row, t3$row)]), domain = "animal")
}

# CADSES: classification of each record (phylum, class, order, family), used
# only when the phylum is an animal phylum: in the data, "Oligochaeta" is
# classified as a plant genus (Tracheophyta, Asterales, Compositae).
for (id in c("cadses_boxcorer", "cadses_leafpack")) {
  if (!use_source(id, reg)) next
  x <- read.table(source_files(id, reg)[1], sep = ";", header = TRUE, quote = "", comment.char = "",
                  na.strings = "NA", fileEncoding = "UTF-8-BOM", colClasses = "character")
  x <- x[x$scientificName %in% taxa_of(id), ]
  ok <- x$phylum %in% animal_phyla
  for (k in c("phylum", "class", "order", "family")) {
    if (!k %in% names(x)) x[[k]] <- NA_character_
    x[[k]][!ok] <- NA_character_
  }
  nm[[id]] <- names_of(id, x$scientificName, "Animalia", x$phylum, x$class, x$order, x$family,
                       domain = "animal")
}

# Fish: Arranz et al. (2022) give the family; Koning et al. (2024), Bagley
# et al. (2022) and NEON only the species.
if (use_source("arranz2022", reg)) {
  x <- read.delim(source_files("arranz2022", reg)[1], stringsAsFactors = FALSE)
  nm$arranz <- names_of("arranz2022", x$Scientific_name, "Animalia", "Chordata",
                        family = x$Family, domain = "animal")
}
for (id in c("koning2024", "bagley2022", "neon_fish")) {
  if (use_source(id, reg)) nm[[id]] <- names_of(id, taxa_of(id), "Animalia", "Chordata", domain = "animal")
}
if (use_source("coughlan2021", reg)) {
  nm$coughlan <- names_of("coughlan2021", taxa_of("coughlan2021"), "Animalia", "Mollusca",
                          "Bivalvia", domain = "animal")
}

# Amphibians: order and family from Santini et al. (2018), Table S1
# ("Urodela" there; Caudata in the COL, where Urodela is a genus).
if (use_source("santini2018", reg)) {
  s1 <- suppressMessages(as.data.frame(readxl::read_excel(source_files("santini2018", reg)[1],
                                                          sheet = "Data")))
  s1$Order[s1$Order %in% "Urodela"] <- "Caudata"
  s1$Species <- sub("^(.+?)\\1+$", "\\1", s1$Species, perl = TRUE)
  amph <- function(id) {
    tx <- taxa_of(id)
    o <- vapply(tx, function(t) most_frequent(s1$Order[s1$Species == t]), "")
    f <- vapply(tx, function(t) most_frequent(s1$Family[s1$Species == t]), "")
    names_of(id, tx, "Animalia", "Chordata", "Amphibia", o, f, domain = "animal")
  }
  nm$santini <- amph("santini2018")
  if (use_source("trochet2014", reg)) nm$trochet <- amph("trochet2014")
}

# Rimet and Druart (2018): their classification; kingdom "Bacteria" and the
# phyla of algae are not passed (other names in the COL).
if (use_source("rimet2018", reg)) {
  x <- suppressMessages(as.data.frame(readxl::read_excel(source_files("rimet2018", reg)[1], sheet = 1)))
  nm$rimet <- names_of("rimet2018", x[["Genus + species name"]], x$Kingdom, NA, x$Class,
                       x$Order, x$Family, domain = "phytoplankton")
}

# Brey et al. (2010): lineage of each record; the major taxon gives the
# phylum of animals and the kingdom of algae and plants.
if (use_source("brey2010", reg)) {
  tmp <- tempfile("brey")
  dir.create(tmp)
  utils::unzip(source_files("brey2010", reg)[1], files = "Conversion03.xls", exdir = tmp)
  raw <- as.data.frame(readxl::read_excel(file.path(tmp, "Conversion03.xls"), sheet = "Data",
                                          col_names = FALSE, col_types = "text",
                                          .name_repair = "minimal"))[-(1:5), ]
  sp <- na_if_empty(raw[[3]])
  taxon <- ifelse(is.na(sp), na_if_empty(raw[[4]]), sp)
  major <- na_if_empty(raw[[9]])
  phy <- ifelse(major %in% names(phylum_alias), unname(phylum_alias[major]), major)
  animal <- na_if_empty(raw[[10]]) %in% "Animalia"
  plant_major <- c("Rhodophycota", "Chlorophycota", "Prasinophyta", "Charophyta", "Magnoliophyta",
                   "Pteridophyta", "Equisetophyta")
  chromist_major <- c("Phaeophycophyta", "Bacillariophyta", "Pyrrophycophyta", "Haptophyta",
                      "Cryptophycophyta", "Chrysophyta", "Ciliophora")
  kingdom <- ifelse(animal, "Animalia", ifelse(major %in% plant_major, "Plantae",
                                               ifelse(major %in% chromist_major, "Chromista", NA)))
  keep <- taxon %in% taxa_of("brey2010")
  nm$brey <- names_of("brey2010", taxon[keep], kingdom[keep],
                      ifelse(animal & phy %in% animal_phyla, phy, NA)[keep],
                      na_if_empty(raw[[8]])[keep], na_if_empty(raw[[7]])[keep],
                      na_if_empty(raw[[5]])[keep],
                      domain = ifelse(animal, "animal", "other")[keep])
}

# Energy densities and mass ratios from articles (data-raw/energy_literature.csv):
# phylum and class from the group of each animal.
lit_groups <- data.frame(
  taxon = c("Brachionus", "Geocrinia vitellina", "Diplomystes viedmensis", "Chilina gibbosa",
            "Anguilla anguilla", "Rana temporaria", "Bufo bufo", "Brachionus plicatilis",
            "Brachionus calyciflorus",
            # Méthot et al. (2012), Table 2
            "Oligochaeta", "Ancylidae", "Bithyniidae", "Hydrobiidae", "Planorbidae", "Asellidae",
            "Gammaridae", "Caenidae", "Coenagrionidae", "Corixidae", "Leptoceridae", "Pyralidae",
            "Chironomidae", "Tipulidae",
            # Cianciara (1980)
            "Cloeon dipterum"),
  phylum = c("Rotifera", "Chordata", "Chordata", "Mollusca", "Chordata", "Chordata", "Chordata",
             "Rotifera", "Rotifera",
             "Annelida", rep("Mollusca", 4), rep("Arthropoda", 9),
             "Arthropoda"),
  class = c(NA, "Amphibia", NA, "Gastropoda", NA, "Amphibia", "Amphibia", NA, NA,
            "Clitellata", rep("Gastropoda", 4), rep("Malacostraca", 2), rep("Insecta", 7),
            "Insecta"),
  stringsAsFactors = FALSE)
for (id in c("yufera1997", "mitchell2001", "ciancio2006", "nelson1997", "theilacker1971",
             "theilacker1984", "jensen2006", "methot2012", "cianciara1980")) {
  if (!use_source(id, reg) || !length(taxa_of(id))) next
  tx <- taxa_of(id)
  g <- lit_groups[match(tx, lit_groups$taxon), ]
  if (anyNA(g$phylum)) stop("Group unknown for: ", paste(tx[is.na(g$phylum)], collapse = ", "))
  nm[[id]] <- names_of(id, tx, "Animalia", g$phylum, g$class, domain = "animal")
}

# Cummins and Wuycheck (1971): phylum from the first level of their
# systematic position and family from the level ending in -idae; Protozoa are
# not animals.
if (use_source("cummins1971", reg) && file.exists("data-raw/cummins1971_table3.csv")) {
  cw <- read.csv("data-raw/cummins1971_table3.csv", stringsAsFactors = FALSE, encoding = "UTF-8",
                 colClasses = "character", na.strings = "")
  levels_cw <- strsplit(gsub(" *\\([^)]*\\)", "", cw$higher_taxa), " > ", fixed = TRUE)
  first <- vapply(levels_cw, function(v) v[1], character(1))
  cw_alias <- c(Coelenterata = "Cnidaria", Vertebrata = "Chordata", Pisces = "Chordata",
                Pelecypoda = "Mollusca", Bryozoa = "Bryozoa", Ectoprocta = "Bryozoa",
                Brachyopoda = "Brachiopoda")
  ph <- ifelse(first %in% animal_phyla, first,
               ifelse(first %in% names(phylum_alias), unname(phylum_alias[first]),
                      ifelse(first %in% names(cw_alias), unname(cw_alias[first]), NA)))
  fam <- vapply(levels_cw, function(v) {
    f <- v[grepl("idae$", v)]
    if (length(f)) f[length(f)] else NA_character_
  }, character(1))
  protozoa <- first %in% "Protozoa"
  nm$cummins <- names_of("cummins1971", cw$taxon, ifelse(protozoa, NA, "Animalia"), ph,
                         family = fam, domain = ifelse(protozoa, "other", "animal"))
}

# Evans-White et al. (2005): phylum, class and family of each taxon, written
# in the transcription (data-raw/evanswhite2005_appendix.csv).
if (use_source("evanswhite2005", reg) && file.exists("data-raw/evanswhite2005_appendix.csv")) {
  ew <- read.csv("data-raw/evanswhite2005_appendix.csv", stringsAsFactors = FALSE, encoding = "UTF-8",
                 colClasses = "character", na.strings = "")
  nm$evanswhite <- names_of("evanswhite2005", ew$taxon, "Animalia", ew$phylum, ew$class,
                            family = ew$family, domain = "animal")
}

# Zooplankton transcriptions (Culver et al. 1985, Bottrell et al. 1976,
# McCauley 1984, Ejsmont-Karabin 1998): phylum and class written in each
# transcription (columns `taxon`, `phylum`, `class` of data-raw/<source>_*.csv).
for (id in c("culver1985", "bottrell1976", "mccauley1984", "ejsmontkarabin1998")) {
  if (!use_source(id, reg) || !length(taxa_of(id))) next
  files <- list.files("data-raw", pattern = paste0("^", id, "_.*\\.csv$"), full.names = TRUE)
  z <- do.call(rbind, lapply(files, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE, encoding = "UTF-8", colClasses = "character",
                  na.strings = c("", "NA"))
    if (all(c("taxon", "phylum", "class") %in% names(d))) d[, c("taxon", "phylum", "class")]
  }))
  z <- z[!is.na(z$taxon) & !duplicated(z$taxon), ]
  tx <- taxa_of(id)
  g <- z[match(tx, z$taxon), ]
  if (anyNA(g$phylum)) stop("Group unknown for: ", paste(tx[is.na(g$phylum)], collapse = ", "))
  nm[[id]] <- names_of(id, tx, "Animalia", g$phylum, g$class, domain = "animal")
}

# Hébert et al. (2016): cladocerans and copepods (group at the start of the
# note of each typical record).
if (use_source("hebert2016", reg) && length(taxa_of("hebert2016"))) {
  hb <- d_tr[d_tr$source_id == "hebert2016", c("taxon", "note")]
  grp <- sub(";.*$", "", hb$note)
  stopifnot(all(grp %in% c("Cladocera", "Copepoda")))
  nm$hebert <- names_of("hebert2016", hb$taxon, "Animalia", "Arthropoda",
                        ifelse(grp == "Cladocera", "Branchiopoda", "Copepoda"), domain = "animal")
}

nm <- do.call(rbind, nm)
rownames(nm) <- NULL
# Hints passed to the API: kingdoms of the COL only, phyla only for animals.
nm$kingdom[!nm$kingdom %in% col_kingdoms] <- NA
nm$phylum[!(nm$phylum %in% animal_phyla)] <- NA
missing_pairs <- setdiff(
  unique(c(paste(old$source, old$taxon), paste(d_tr$source_id, d_tr$taxon),
           paste(d_lf$source_id, d_lf$taxon), paste(d_lp$source_id, d_lp$taxon),
           paste(d_cv$source_id, d_cv$taxon))),
  paste(nm$source_id, nm$taxon))
if (length(missing_pairs)) stop("Names without hints: ", paste(head(missing_pairs, 10), collapse = "; "))
say("## Names")
say("")
say("- %d names of %d sources (%d distinct names).", nrow(nm), length(unique(nm$source_id)),
    length(unique(nm$taxon)))
say("")

# 2. Query names and the rank they should have ------------------------------------
# The names are sent as given by the source, after transliterating letters
# with diacritics, dropping "cf." and "aff.", reducing "sp. 1", "sp1" and
# similar to "sp.", and dropping the descriptions of morphotypes that follow
# the name in Rimet and Druart (single cell, colony of 10 um, 4 cells, plein,
# vide...). The form of the name says what rank it should get: genus for
# "Genus sp.", species for binomials, the species or below for names with
# var., f. or three words, the rank of the suffix (-idae, -aceae, -ini...),
# and genus or higher for other single words.

clean_query <- function(x) {
  y <- iconv(x, "UTF-8", "ASCII//TRANSLIT")
  y[is.na(y)] <- x[is.na(y)]
  y <- gsub("\\s+", " ", trimws(y))
  y <- sub("^(.+?)\\1+$", "\\1", y, perl = TRUE)
  y <- gsub("(^| )(cf|aff)\\.? ", "\\1", y)
  y <- sub("^([A-Z][a-z-]+) \\([A-Z][a-z-]+\\) ", "\\1 ", y)
  y <- sub("^([A-Z][a-z-]+) (sp|spp)(\\.|[0-9]| |$).*$", "\\1 sp.", y)
  y <- sub(" \\((single cell|[0-9]+ cells?|petit|grande?|large)\\)$", "", y)
  y <- sub(" (single cell|colony .*|[0-9]+ cells?|filament .*|plein|vide|petit morphotype|tab\\. .*)$", "", y)
  y
}
suffix_rank <- c(idae = "FAMILY", aceae = "FAMILY", inae = "SUBFAMILY", oideae = "SUBFAMILY",
                 ini = "TRIBE", eae = "TRIBE", oidea = "SUPERFAMILY", ales = "ORDER",
                 ineae = "SUBORDER", phyceae = "CLASS", opsida = "CLASS", mycetes = "CLASS",
                 phyta = "PHYLUM", mycota = "PHYLUM")
name_form <- function(q) {
  single <- grepl("^[A-Z][a-z-]+$", q)
  suf <- vapply(q, function(n) {
    hit <- names(suffix_rank)[endsWith(n, names(suffix_rank))]
    if (length(hit)) suffix_rank[[hit[which.max(nchar(hit))]]] else NA_character_
  }, "")
  ifelse(grepl("^[A-Z][a-z-]+ sp\\.$", q), "genus_sp",
         ifelse(grepl("^[A-Z][a-z-]+ [a-z-]+$", q), "species",
                ifelse(grepl("^[A-Z][a-z-]+ [a-z-]+ (var\\. |f\\. |fo\\. |subsp\\. |ssp\\. )?[a-z-]+$", q),
                       "infraspecific",
                       ifelse(single & !is.na(suf), tolower(suf),
                              ifelse(single, "genus_or_higher", "other")))))
}
nm$query_name <- clean_query(nm$taxon)
nm$form <- name_form(nm$query_name)
nm$genus_part <- sub(" .*$", "", nm$query_name)
nm$species_part <- ifelse(nm$form %in% c("species", "infraspecific"),
                          sub("^([A-Z][a-z-]+ [a-z-]+).*$", "\\1", nm$query_name), NA_character_)

# 3. Queries, with the cache -----------------------------------------------------

query_cols <- c("query_name", hint_cols, "taxon_rank")
query_key <- function(q) {
  k <- paste(q$query_name, q$kingdom, q$phylum, q$class, q$order, q$family, sep = "|")
  ifelse(is.na(q$taxon_rank), k, paste0(k, "|rank=", q$taxon_rank))
}
query_url <- function(q) {
  p <- c(scientificName = q$query_name, checklistKey = checklist_key, kingdom = q$kingdom,
         phylum = q$phylum, class = q$class, order = q$order, family = q$family,
         taxonRank = q$taxon_rank)
  p <- p[!is.na(p)]
  paste0(api, "?", paste(names(p), vapply(p, utils::URLencode, "", reserved = TRUE),
                         sep = "=", collapse = "&"))
}
answer_cols <- c("match_type", "confidence", "match_note", "issues", "usage_key", "usage_name",
                 "usage_authorship", "usage_rank", "usage_status", "accepted_key", "accepted_name",
                 "accepted_rank", "classification", "classification_keys", "queried_on")
parse_answer <- function(txt) {
  r <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
  g <- function(x, f) if (is.null(x) || is.null(x[[f]])) NA_character_ else as.character(x[[f]])
  d <- r$diagnostics
  u <- r$usage
  a <- r$acceptedUsage
  cl <- r$classification
  data.frame(match_type = g(d, "matchType"), confidence = g(d, "confidence"),
             match_note = g(d, "note"),
             issues = if (length(d$issues)) paste(unlist(d$issues), collapse = ";") else NA_character_,
             usage_key = g(u, "key"), usage_name = g(u, "canonicalName"),
             usage_authorship = g(u, "authorship"), usage_rank = g(u, "rank"),
             usage_status = g(u, "status"), accepted_key = g(a, "key"),
             accepted_name = g(a, "canonicalName"), accepted_rank = g(a, "rank"),
             classification = if (length(cl)) paste(vapply(cl, function(z) paste0(z$rank, ":", z$name), ""), collapse = "|") else NA_character_,
             classification_keys = if (length(cl)) paste(vapply(cl, function(z) paste0(z$rank, ":", z$key), ""), collapse = "|") else NA_character_,
             stringsAsFactors = FALSE)
}
fetch <- function(urls, pool_size = 6, tries = 3) {
  body <- rep(NA_character_, length(urls))
  status <- rep(NA_integer_, length(urls))
  todo <- seq_along(urls)
  for (attempt in seq_len(tries)) {
    if (!length(todo)) break
    if (attempt > 1) Sys.sleep(5 * attempt)
    pool <- curl::new_pool(total_con = pool_size, host_con = pool_size)
    for (j in todo) {
      local({
        k <- j
        h <- curl::new_handle(useragent = "bmassa data-raw script (R)", timeout = 60)
        curl::curl_fetch_multi(urls[k], handle = h, pool = pool,
          done = function(res) {
            status[k] <<- res$status_code
            body[k] <<- rawToChar(res$content)
          },
          fail = function(msg) status[k] <<- -1L)
      })
    }
    curl::multi_run(pool = pool)
    todo <- which(is.na(status) | status != 200L)
  }
  if (length(todo)) stop(length(todo), " queries failed (HTTP status ",
                          paste(unique(status[todo]), collapse = ", "), ")")
  Encoding(body) <- "UTF-8"
  body
}
cache_cols <- c("key", query_cols, answer_cols)
cache <- if (file.exists(cache_file) && !refresh) {
  read.csv(cache_file, stringsAsFactors = FALSE, encoding = "UTF-8", colClasses = "character")
} else NULL
if (!is.null(cache)) {
  if (!"taxon_rank" %in% names(cache)) cache$taxon_rank <- NA_character_
  cache <- cache[, cache_cols]
}
n_sent <- 0
# Version of the COL XR of the matching service now, and of the cache. The
# cache holds the answers of one version: when the service has moved to
# another one, new queries are refused unless all the queries are sent again
# (BMASSA_COLXR_REFRESH = "true"), so that the names of the reference are never
# matched against two versions.
version_file <- file.path(out_dir, "colxr_version.csv")
meta <- tryCatch(jsonlite::fromJSON(paste0(api, "/metadata?checklistKey=", checklist_key)),
                 error = function(e) NULL)
now <- if (is.null(meta)) NULL else data.frame(
  dataset_alias = meta$mainIndex$datasetAlias, dataset_title = meta$mainIndex$datasetTitle,
  clb_dataset_key = meta$mainIndex$clbDatasetKey, name_usages = meta$mainIndex$nameUsageCount,
  index_built = substr(meta$created, 1, 10), checked_on = format(Sys.Date()))
cached_version <- if (!is.null(cache) && file.exists(version_file)) {
  read.csv(version_file, stringsAsFactors = FALSE)
} else NULL
run_queries <- function(q) {
  # Sends the queries that are not in the cache, in chunks of 500, and saves
  # the cache after each chunk.
  if (!"taxon_rank" %in% names(q)) q$taxon_rank <- NA_character_
  q <- q[, query_cols]
  q$key <- query_key(q)
  q <- q[!duplicated(q$key), ]
  new <- q[!(q$key %in% cache$key), ]
  if (nrow(new)) {
    if (is.null(now)) stop("The metadata of the matching service could not be read.")
    if (!is.null(cached_version) && cached_version$clb_dataset_key != now$clb_dataset_key) {
      stop(nrow(new), " names are not in the cache, which holds answers of the ",
           cached_version$dataset_alias, ", but the service now uses the ", now$dataset_alias,
           ": run with BMASSA_COLXR_REFRESH = \"true\" to send all the queries again.")
    }
    message(nrow(new), " queries to send")
    chunks <- split(seq_len(nrow(new)), ceiling(seq_len(nrow(new)) / 500))
    for (ch in chunks) {
      b <- fetch(vapply(ch, function(i) query_url(new[i, ]), ""))
      parsed <- do.call(rbind, lapply(b, parse_answer))
      add <- cbind(new[ch, c("key", query_cols)], parsed, queried_on = format(Sys.Date()),
                   stringsAsFactors = FALSE)
      add[] <- lapply(add, as.character)
      cache <<- rbind(cache, add[, cache_cols])
      write.csv(cache[order(cache$key, method = "radix"), ], cache_file, row.names = FALSE,
                fileEncoding = "UTF-8")
      message("  ", nrow(cache), " queries in the cache")
    }
  }
  n_sent <<- n_sent + nrow(new)
  invisible(nrow(new))
}
answer <- function(q) {
  if (!"taxon_rank" %in% names(q)) q$taxon_rank <- NA_character_
  k <- match(query_key(q), cache$key)
  if (anyNA(k)) stop("Queries missing from the cache: ", paste(q$query_name[is.na(k)], collapse = ", "))
  a <- cache[k, answer_cols]
  rownames(a) <- NULL
  a
}


nm$taxon_rank <- NA_character_
run_queries(nm)
nm <- cbind(nm, answer(nm))
nm$attempt <- "hints"

# 4. Checks of an answer and second queries ---------------------------------------
# An answer is accepted when it is an exact match (or a spelling variant of
# the same name, below), at the rank the name should have, and its
# classification is consistent with the source: Animalia for animals, not
# Animalia for phytoplankton and the algae of Brey et al., and the phylum
# given by the source for animals.
# The hints only break ties between homonyms, so a genus that shares its
# name with its nominotypical subgenus (Chironomus, Daphnia...), or one whose
# hints are outdated, comes back placed at a higher rank. Those names are
# sent again with the rank they should have (taxonRank), first with the hints,
# then with the kingdom (and phylum) only, and then alone, and the first
# accepted answer is taken. "Genus sp." is sent as the genus.

good_type <- c("EXACT", "FUZZY", "VARIANT")
rank_name <- function(cls, rank) {
  vapply(strsplit(ifelse(is.na(cls), "", cls), "|", fixed = TRUE), function(p) {
    hit <- p[startsWith(p, paste0(rank, ":"))]
    if (length(hit)) sub("^[A-Z_]+:", "", hit[1]) else NA_character_
  }, "")
}
in_classification <- function(cls, nms) {
  mapply(function(c1, n1) is.na(n1) || (!is.na(c1) && n1 %in% sub("^[A-Z_]+:", "", strsplit(c1, "|", fixed = TRUE)[[1]])),
         cls, nms, USE.NAMES = FALSE)
}
consistent <- function(x) {
  k <- rank_name(x$classification, "KINGDOM")
  ok_dom <- ifelse(x$domain == "animal", k %in% "Animalia",
                   !(k %in% "Animalia") & !is.na(x$classification))
  ok_dom & (x$domain != "animal" | in_classification(x$classification, x$phylum))
}
# Same name, or a spelling variant: a genus that differs in at most two
# letters (one for a name of one word) and epithets that differ in one letter
# or only in their ending (-us, -a, -um, -is, -e, -i, -y). Epithets that
# differ more (minuta and minima, contorta and contexta) are other species.
same_name <- function(a, b) {
  mapply(function(a1, b1) {
    if (is.na(a1) || is.na(b1)) return(FALSE)
    if (a1 == b1) return(TRUE)
    wa <- strsplit(a1, " ")[[1]]
    wb <- strsplit(b1, " ")[[1]]
    if (length(wa) != length(wb)) return(FALSE)
    if (length(wa) == 1) return(nchar(a1) >= 5 && adist(a1, b1) <= 1)
    if (adist(wa[1], wb[1]) > 2) return(FALSE)
    stem <- function(w) sub("(us|um|is|a|e|i|y)$", "", w)
    all(mapply(function(x, y) adist(x, y) <= 1 || stem(x) == stem(y), wa[-1], wb[-1]))
  }, a, b, USE.NAMES = FALSE)
}
# Ambiguous synonyms and misapplied names are followed only to the same
# species (Acartia clausi to Acartia clausii, Macrocyclops albidus to its
# nominotypical subspecies); otherwise the name goes to its genus.
# A one-word name may differ in one letter only when the answer is a genus
# (Dasyoma and Dasyomma), not for higher taxa (Hirudinida, an order, and
# Hirudinidae, a family).
same_taxon <- function(target, a) {
  same_name(target, a$usage_name) &
    (target %in% a$usage_name | grepl(" ", target) | a$usage_rank %in% "GENUS")
}
first_words <- function(x) sub("^([A-Z][a-z-]+( [a-z-]+)?).*$", "\\1", x)
status_ok <- function(a, target) {
  odd <- a$usage_status %in% c("AMBIGUOUS_SYNONYM", "MISAPPLIED")
  !odd | same_name(first_words(target), first_words(a$accepted_name))
}
expected_rank <- function(form) {
  ifelse(form %in% c("genus_sp", "genus_or_higher"), "GENUS",
         ifelse(form %in% c("species", "infraspecific"), "SPECIES",
                ifelse(form %in% tolower(suffix_rank), toupper(form), NA_character_)))
}
# A spelling variant of the genus is accepted only when the family (or, if
# the source gives none, the class) of the source is in its classification
# (not Geminella for Gesinella). In the second queries, a genus without family
# in the COL XR is not accepted for a one-word name: the COL XR has
# placeholder genera named after higher taxa (Oligochaeta, Brachyura).
accept <- function(x, a, target, second = FALSE) {
  # x: rows of nm; a: answers; target: name the answer should carry.
  g_t <- sub(" .*$", "", target)
  g_a <- sub(" .*$", "", a$usage_name)
  other_genus <- !is.na(g_a) & g_t != g_a
  hint_ok <- ifelse(!is.na(x$family), in_classification(a$classification, x$family),
                    in_classification(a$classification, x$class))
  placeholder <- second & !grepl(" ", target) & a$usage_rank %in% "GENUS" &
    is.na(rank_name(a$classification, "FAMILY"))
  a$match_type %in% good_type & consistent(cbind(x[, c("domain", "phylum")], a)) &
    same_taxon(target, a) & status_ok(a, target) & (!other_genus | hint_ok) & !placeholder
}
nm$expected_rank <- expected_rank(nm$form)
nm$target <- ifelse(nm$form == "genus_sp", nm$genus_part,
                    ifelse(nm$form == "infraspecific", nm$query_name, nm$query_name))
ok_first <- accept(nm, nm[, answer_cols], nm$target) |
  (nm$form == "infraspecific" & nm$match_type %in% c(good_type, "HIGHERRANK") &
     same_name(nm$species_part, nm$usage_name) & consistent(nm) & status_ok(nm, nm$species_part)) |
  (nm$form == "genus_sp" & same_taxon(nm$genus_part, nm) & consistent(nm) &
     status_ok(nm, nm$genus_part))
retry <- which(!ok_first & !is.na(nm$expected_rank))
second <- list(
  rank_hints = function(q) q,
  rank_kingdom = function(q) { q[, c("class", "order", "family")] <- NA; q },
  rank_only = function(q) { q[, hint_cols] <- NA; q })
if (length(retry)) {
  base <- nm[retry, c("query_name", hint_cols)]
  base$taxon_rank <- nm$expected_rank[retry]
  base$query_name <- ifelse(nm$form[retry] == "genus_sp", nm$genus_part[retry],
                            ifelse(nm$form[retry] == "infraspecific", nm$species_part[retry],
                                   nm$query_name[retry]))
  tgt <- base$query_name
  for (s in names(second)) run_queries(second[[s]](base))
  done <- rep(FALSE, length(retry))
  for (s in names(second)) {
    a <- answer(second[[s]](base))
    ok <- !done & accept(nm[retry, ], a, tgt, second = TRUE)
    nm[retry[ok], answer_cols] <- a[ok, answer_cols]
    nm$attempt[retry[ok]] <- s
    done <- done | ok
  }
  # Species (and varieties) still without an accepted answer go to their
  # genus, sent with the rank GENUS.
  gen <- retry[!done & nm$form[retry] %in% c("species", "infraspecific", "genus_sp")]
  if (length(gen)) {
    gbase <- nm[gen, c("query_name", hint_cols)]
    gbase$query_name <- nm$genus_part[gen]
    gbase$taxon_rank <- "GENUS"
    for (s in names(second)) run_queries(second[[s]](gbase))
    gdone <- rep(FALSE, length(gen))
    for (s in names(second)) {
      a <- answer(second[[s]](gbase))
      ok <- !gdone & accept(nm[gen, ], a, gbase$query_name, second = TRUE)
      nm[gen[ok], answer_cols] <- a[ok, answer_cols]
      nm$attempt[gen[ok]] <- paste0("genus_", s)
      gdone <- gdone | ok
    }
  }
}

# Manual corrections (data-raw/taxa_manual.csv), by source and name ("*" for
# any source): "query" sends another name or other hints and takes the answer
# as it is; "accept" keeps the answer although it was flagged; "drop" leaves
# the name without a COL taxon.
man_file <- "data-raw/taxa_manual.csv"
man <- if (file.exists(man_file)) read.csv(man_file, stringsAsFactors = FALSE, encoding = "UTF-8",
                                           na.strings = "") else NULL
nm$manual <- NA_character_
nm$manual_action <- NA_character_
if (!is.null(man) && nrow(man)) {
  stopifnot(all(man$action %in% c("query", "accept", "drop")))
  hits <- lapply(seq_len(nrow(man)), function(i) {
    which(nm$taxon == man$taxon[i] & (man$source_id[i] == "*" | nm$source_id == man$source_id[i]))
  })
  unused <- man$taxon[lengths(hits) == 0]
  if (length(unused)) stop("Rows of taxa_manual.csv that match no name: ", paste(unused, collapse = ", "))
  q <- man[man$action == "query", ]
  if (nrow(q)) run_queries(q)
  for (i in seq_len(nrow(man))) {
    k <- hits[[i]]
    nm$manual[k] <- man$reason[i]
    nm$manual_action[k] <- man$action[i]
    if (man$action[i] == "drop") {
      nm[k, answer_cols] <- NA
    } else if (man$action[i] == "query") {
      a <- answer(man[i, ])
      nm[k, answer_cols] <- a[rep(1, length(k)), answer_cols]
      nm$query_name[k] <- man$query_name[i]
      nm$target[k] <- sub(" sp\\.$", "", man$query_name[i])
      nm$attempt[k] <- "manual"
    }
  }
}

# 5. Accepted taxon, classification and flags ---------------------------------------

syn <- !is.na(nm$accepted_key)
nm$col_id <- ifelse(syn, nm$accepted_key, nm$usage_key)
nm$col_name <- ifelse(syn, nm$accepted_name, nm$usage_name)
nm$col_rank <- tolower(ifelse(syn, nm$accepted_rank, nm$usage_rank))
for (r in c("DOMAIN", "KINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")) {
  nm[[paste0("col_", tolower(r))]] <- rank_name(nm$classification, r)
}
rank_order <- c("domain", "realm", "kingdom", "subkingdom", "infrakingdom", "phylum", "subphylum",
                "infraphylum", "parvphylum", "gigaclass", "megaclass", "superclass", "class",
                "subclass", "infraclass", "subterclass", "superorder", "order", "suborder",
                "infraorder", "parvorder", "superfamily", "epifamily", "family", "subfamily",
                "supertribe", "tribe", "subtribe", "genus", "subgenus", "section_botany",
                "species_aggregate", "species", "subspecies", "variety", "subvariety", "form",
                "infraspecific_name")
rank_pos <- function(r) match(tolower(r), rank_order)
manual_ok <- nm$manual_action %in% c("query", "accept")
nm$consistent <- !is.na(nm$classification) & consistent(nm)
own_genus <- nm$form %in% c("genus_sp", "species", "infraspecific") & nm$usage_rank %in% "GENUS" &
  same_taxon(nm$genus_part, nm)
own <- nm$match_type %in% good_type & same_taxon(nm$target, nm) & !own_genus
own_species <- nm$form %in% "infraspecific" & nm$usage_rank %in% "SPECIES" &
  same_name(nm$species_part, nm$usage_name)
nm$placement <- ifelse(nm$manual_action %in% "drop", "dropped",
                ifelse(is.na(nm$usage_key), "unmatched",
                ifelse(nm$form == "genus_sp" & own_genus, "genus (identified to genus)",
                ifelse(nm$form == "genus_sp" & same_taxon(nm$genus_part, nm), "own rank",
                ifelse(own & nm$match_type == "EXACT", "own rank",
                ifelse(own, "own rank (spelling variant)",
                ifelse(own_species, "species (variety or form not in the COL)",
                ifelse(own_genus, "genus (species not in the COL)",
                       paste0("higher: ", tolower(nm$usage_rank), " ", nm$usage_name)))))))))
flag <- function(cond, text) ifelse(cond, text, NA_character_)
fl <- cbind(
  flag(nm$placement == "unmatched", paste0("no match", ifelse(is.na(nm$match_note), "", paste0(": ", nm$match_note)))),
  flag(startsWith(nm$placement, "higher") & nm$form != "other", "placed at a higher rank than the name"),
  flag(!is.na(nm$classification) & !nm$consistent, "classification inconsistent with the source (homonym?)"),
  flag(nm$placement == "own rank (spelling variant)", paste0("spelling variant: ", nm$usage_name)),
  flag(syn & rank_pos(nm$usage_rank) <= rank_pos("genus") &
         rank_pos(nm$accepted_rank) < rank_pos(nm$usage_rank),
       paste0("synonym of a taxon of higher rank (", tolower(nm$accepted_rank), " ", nm$accepted_name, ")")),
  flag(!status_ok(nm, nm$target), tolower(nm$usage_status)),
  flag(nm$col_kingdom %in% "Fungi", "fungus"))
nm$flag <- apply(fl, 1, function(v) if (all(is.na(v))) NA_character_ else paste(v[!is.na(v)], collapse = "; "))
review_flags <- "no match|higher rank than the name|inconsistent|synonym of a taxon of higher rank|ambiguous_synonym|misapplied|fungus"
nm$review <- !is.na(nm$flag) & grepl(review_flags, nm$flag) & !manual_ok & !(nm$manual_action %in% "drop")
nm$use <- !is.na(nm$col_id) & (manual_ok | (nm$consistent & !nm$review))
# Version of the COL XR of the matching service (read before the queries),
# saved with the cache: it is rewritten when queries are sent and read from
# the file otherwise.
if (!is.null(now) && (n_sent > 0 || !file.exists(version_file))) {
  write.csv(now, version_file, row.names = FALSE, fileEncoding = "UTF-8")
}
ver <- if (file.exists(version_file)) read.csv(version_file, stringsAsFactors = FALSE) else NULL
nm$col_release <- if (is.null(ver)) NA_character_ else ver$dataset_alias

out <- nm[, c("source_id", "taxon", "domain", hint_cols, "query_name", "form", "attempt",
              "match_type", "confidence", "match_note", "usage_key", "usage_name",
              "usage_authorship", "usage_rank", "usage_status", "col_id", "col_name", "col_rank",
              "col_domain", "col_kingdom", "col_phylum", "col_class", "col_order", "col_family",
              "col_genus", "col_species", "classification", "placement", "use", "review", "flag",
              "manual", "col_release")]
names(out)[names(out) %in% hint_cols] <- paste0("hint_", hint_cols)
out <- out[order(out$source_id, out$taxon, method = "radix"), ]
write.csv(out, file.path(out_dir, "taxa_colxr.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 6. Log ---------------------------------------------------------------------------

say("## COL XR")
say("")
if (!is.null(ver)) {
  say("- Matching index of the GBIF API v2 when the queries were sent: %s (%s, ChecklistBank dataset %s), %s name usages, index built %s (checked %s).",
      ver$dataset_alias, ver$dataset_title, ver$clb_dataset_key,
      format(ver$name_usages, big.mark = " "), ver$index_built, ver$checked_on)
}
if (!is.null(now) && !is.null(ver) && now$clb_dataset_key != ver$clb_dataset_key) {
  say("- The service now uses %s: run with BMASSA_COLXR_REFRESH = \"true\" to send all the queries again.",
      now$dataset_alias)
}
say("- Queries sent in this run: %d; in the cache: %d.", n_sent, nrow(cache))
say("")
say("## Results")
say("")
pl <- sub("^higher: .*$", "higher rank", out$placement)
lev <- c("own rank", "own rank (spelling variant)", "genus (identified to genus)",
         "genus (species not in the COL)", "species (variety or form not in the COL)",
         "higher rank", "unmatched", "dropped")
tab <- table(out$source_id, factor(pl, lev))
say("| Source | %s | Used | To review |", paste(lev, collapse = " | "))
say("|---|%s---|---|", paste(rep("---|", length(lev)), collapse = ""))
for (s in rownames(tab)) {
  say("| %s | %s | %d | %d |", s, paste(tab[s, ], collapse = " | "), sum(out$use[out$source_id == s]),
      sum(out$review[out$source_id == s]))
}
say("")
say("- %d of %d names used; synonyms replaced by their accepted name: %d; matched by the second queries: %d (%s); manual corrections: %d.",
    sum(out$use), nrow(out), sum(out$usage_status %in% c("SYNONYM", "AMBIGUOUS_SYNONYM", "MISAPPLIED") & out$use),
    sum(out$attempt %in% c(names(second), paste0("genus_", names(second)))),
    paste(sprintf("%s %d", c(names(second), paste0("genus_", names(second))),
                  vapply(c(names(second), paste0("genus_", names(second))), function(s) sum(out$attempt == s), 0L)),
          collapse = ", "),
    sum(!is.na(out$manual)))
agree <- vapply(c("phylum", "class", "order", "family"), function(r) {
  h <- out[[paste0("hint_", r)]]
  k <- !is.na(h) & out$placement == "own rank"
  c(n = sum(k), same = sum(in_classification(out$classification[k], h[k])))
}, c(n = 0, same = 0))
say("- Rank given by the source found in the COL classification of the names matched at their own rank: %s.",
    paste(sprintf("%s %d of %d", colnames(agree), agree["same", ], agree["n", ]), collapse = ", "))
say("")
# Gaps of the COL XR classification, for the model (phase 3).
gs <- out[out$use & out$col_rank %in% c("genus", "species"), ]
gap <- vapply(c("col_phylum", "col_class", "col_order", "col_family"), function(r) sum(is.na(gs[[r]])), 0L)
say("- Genera and species used without some main rank in their COL classification: %s (of %d); the model will have to fill these gaps.",
    paste(sprintf("%s %d", sub("col_", "", names(gap)), gap), collapse = ", "), nrow(gs))
no_order <- sort(table(gs$col_class[is.na(gs$col_order)]), decreasing = TRUE)
say("- Classes with genera or species without order: %s.",
    paste(sprintf("%s %d", names(head(no_order, 8)), as.integer(head(no_order, 8))), collapse = ", "))
if (!is.null(man) && nrow(man)) {
  say("- Manual corrections (`taxa_manual.csv`): %s.",
      paste(sprintf("%s %d", names(table(man$action)), as.integer(table(man$action))), collapse = ", "))
}
say("")
variants <- out[out$placement == "own rank (spelling variant)", ]
say("## Spelling variants accepted (%d)", nrow(variants))
say("")
say("%s.", paste(sprintf("%s -> %s", variants$query_name, variants$usage_name), collapse = "; "))
say("")
rev <- out[out$review, ]
say("## Names to review (%d)", nrow(rev))
say("")
say("| Source | Name | COL | Flag |")
say("|---|---|---|---|")
for (i in seq_len(nrow(rev))) {
  say("| %s | %s | %s | %s |", rev$source_id[i], rev$taxon[i],
      ifelse(is.na(rev$col_name[i]), "-", paste0(rev$col_name[i], " (", rev$col_rank[i], ")")),
      rev$flag[i])
}
writeLines(log_lines, file.path(out_dir, "taxa_colxr_log.md"), useBytes = TRUE)
cat(head(log_lines, 45), sep = "\n")
