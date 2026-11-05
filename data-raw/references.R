# Builds the reference tables of bmassa 0.1.0, which are no longer
# distributed with the package (removed in 0.2.0) but are inputs of the
# scripts of version 0.2.0 (taxonomy_colxr.R and model_data.R):
#   data-raw/derived/taxonomy_0.1.0.csv, typical_mass_0.1.0.csv and
#   length_mass_0.1.0.csv (numbers written with 17 significant digits, so
#   that they are read back exactly), and kremer_colonial_genera.csv (the
#   colonial genera of Kremer et al. 2014, read by model_tables.R)
#
# Run from the package root:  Rscript data-raw/references.R
# Packages used only here: readxl.
#
# Source files are read from the folder given by the environment variable
# BMASSA_SOURCES (default "data-raw/sources"); they are not distributed with
# the package. Where to obtain them, their licences and MD5 checksums are in
# data-raw/README.md. Small tables transcribed by hand from the papers are
# kept in data-raw/*.csv (see the header of each section).

source("data-raw/names.R") # normalize_taxon_names(), rank_levels, is_binomial()

`%||%` <- function(x, y) if (is.null(x)) y else x

src <- Sys.getenv("BMASSA_SOURCES", "data-raw/sources")
expected_md5 <- c(
  "taxa_table_030614.csv" = "147cbba2ba0fd96d33dd0eb9ad32a0e4",
  "bvd_species_ag_030614.csv" = "b3638fd46f5d1518b8206bf96ea6bbaf",
  "bvd_genus_ag_030614.csv" = "3a6f1dc5c9be57c1d4db5731205cba50",
  "bvd_species_raw_030614.csv" = "1e3e815263e9d5a65dc5f91506ca30a3",
  "bvd_genus_raw_030614.csv" = "9d8889b1d03b470bc60c812f7177f758",
  "Database_invertebrates.txt" = "b14b453a45a39b91552f03485c743f03",
  "Database_references.txt" = "3c52052294fbeeebfb5a2c1b88bb2680",
  "vanLeeuwen2025_TableS2.xlsx" = "eeb61729f6fb47065b292871c9c8b442"
)
md5 <- tools::md5sum(file.path(src, names(expected_md5)))
if (!identical(unname(md5), unname(expected_md5))) {
  stop("Source files differ from the expected versions:\n",
       paste(names(expected_md5)[md5 != expected_md5], collapse = "\n"))
}

lineage_cols <- rev(rank_levels) # kingdom ... species

# Empty lineage table with the standard columns.
empty_lineage <- function(n) {
  as.data.frame(setNames(rep(list(rep(NA_character_, n)), length(lineage_cols)),
                         lineage_cols), stringsAsFactors = FALSE)
}

blank_to_na <- function(x) {
  x <- as.character(x)
  x <- gsub("^\\s+|\\s+$", "", x)
  x[!is.na(x) & !nzchar(x)] <- NA_character_
  x
}

# =============================================================================
# 1. Kremer et al. (2014), Ecology 95: 2984, doi:10.1890/14-0603.1
#    Ecological Archives E095-257. Biovolumes in log10(um3).
#    Metadata: "Copyright restrictions: None"; the paper must be cited.
# =============================================================================
sp_raw <- read.csv(file.path(src, "bvd_species_raw_030614.csv"), stringsAsFactors = FALSE)
ge_raw <- read.csv(file.path(src, "bvd_genus_raw_030614.csv"), stringsAsFactors = FALSE)
sp_ag <- read.csv(file.path(src, "bvd_species_ag_030614.csv"), stringsAsFactors = FALSE)
ge_ag <- read.csv(file.path(src, "bvd_genus_ag_030614.csv"), stringsAsFactors = FALSE)
kr_taxa <- read.csv(file.path(src, "taxa_table_030614.csv"), stringsAsFactors = FALSE)

# Verification: the published taxon means are means of log10 biovolumes of all
# the raw records of each accepted name (species, varieties and forms).
chk <- tapply(sp_raw$nu.biovol, sp_raw$accepted.taxa.name, mean)
stopifnot(max(abs(sp_ag$mn_nu.biovol - chk[sp_ag$accepted.taxa.name])) < 1e-12)
chk_n <- tapply(sp_raw$nu.biovol, sp_raw$accepted.taxa.name, length)
stopifnot(all(sp_ag$n == chk_n[sp_ag$accepted.taxa.name]))
# The genus product pools all the records of each genus (genus raw file:
# records identified to species, variety, form or genus only).
chk_g <- tapply(ge_raw$nu.biovol, ge_raw$genus, mean)
stopifnot(max(abs(ge_ag$mn_nu.biovol - chk_g[ge_ag$accepted.taxa.name])) < 1e-6)
stopifnot(all(ge_ag$n == tapply(ge_raw$nu.biovol, ge_raw$genus, length)[ge_ag$accepted.taxa.name]))

# Per-cell biovolumes come from records that report the number of cells per
# natural unit. In colonial or filamentous genera some data sources report a
# whole colony or filament as one "cell" (in species with both kinds of
# records, the per-cell volume of their one-cell records is up to 30 times
# that of their multi-cell records), so in those genera only records with more
# than one cell per natural unit are used. A genus is treated as colonial when
# at least 5% of its records with an explicit cell count report more than one
# cell per natural unit; below that, multi-cell records are occasional
# (e.g. 3 of 3212 in Cyclotella) and all records with an explicit cell count
# are used.
colonial_min_fraction <- 0.05
explicit_ge <- !is.na(ge_raw$cells.per.nu) & !is.na(ge_raw$cell.biovol)
multi_fraction <- tapply(ge_raw$cells.per.nu[explicit_ge] > 1, ge_raw$genus[explicit_ge], mean)
colonial_genera <- sort(names(multi_fraction)[multi_fraction >= colonial_min_fraction])
occasional_genera <- sort(names(multi_fraction)[multi_fraction > 0 &
                                                  multi_fraction < colonial_min_fraction])

# Summaries of linear biovolumes (um3) from log10 values.
summarise_biovolume <- function(log10_values, sources) {
  v <- 10^log10_values
  c(
    n = length(v),
    value = mean(v),
    geometric_mean = 10^mean(log10_values),
    median = stats::median(v),
    value_min = min(v),
    value_max = max(v),
    sd_log10 = if (length(v) > 1L) stats::sd(log10_values) else NA_real_,
    n_sources = length(unique(sources))
  )
}

# Records of one product (species: all records of each accepted name; genus:
# all records of each genus, as in the genus product of Kremer et al.).
kremer_records <- function(raw, product, count_unit) {
  if (count_unit == "natural_unit") {
    use <- raw
    use$lv <- use$nu.biovol
  } else {
    explicit <- !is.na(raw$cells.per.nu) & !is.na(raw$cell.biovol)
    colonial <- raw$genus %in% colonial_genera
    keep <- explicit & (!colonial | raw$cells.per.nu > 1)
    use <- raw[keep, , drop = FALSE]
    use$lv <- use$cell.biovol
  }
  key <- if (product == "genus") use$genus else use$accepted.taxa.name
  groups <- split(seq_len(nrow(use)), key)
  groups <- groups[order(names(groups), method = "radix")]
  stats <- t(vapply(groups, function(i) summarise_biovolume(use$lv[i], use$data.source[i]),
                    numeric(8)))
  out <- data.frame(taxon = names(groups), stats, stringsAsFactors = FALSE, row.names = NULL)
  out$count_unit <- count_unit
  out$product <- product
  out$colonial_evidence <- vapply(strsplit(out$taxon, " "), `[`, character(1), 1) %in% colonial_genera
  out
}

kr <- rbind(
  kremer_records(sp_raw, "species", "natural_unit"),
  kremer_records(sp_raw, "species", "cell"),
  kremer_records(ge_raw, "genus", "natural_unit"),
  kremer_records(ge_raw, "genus", "cell")
)
# The rows built reproduce the published products (natural units).
kr_nu <- kr[kr$count_unit == "natural_unit", ]
pub <- rbind(data.frame(taxon = sp_ag$accepted.taxa.name, product = "species",
                        mn = sp_ag$mn_nu.biovol, n = sp_ag$n),
             data.frame(taxon = ge_ag$accepted.taxa.name, product = "genus",
                        mn = ge_ag$mn_nu.biovol, n = ge_ag$n))
m <- match(paste(pub$product, pub$taxon), paste(kr_nu$product, kr_nu$taxon))
stopifnot(!anyNA(m), nrow(kr_nu) == nrow(pub))
stopifnot(all(kr_nu$n[m] == pub$n))
stopifnot(max(abs(log10(kr_nu$geometric_mean[m]) - pub$mn)) < 1e-6)

kr$note <- ifelse(
  kr$count_unit == "natural_unit" & kr$colonial_evidence,
  "colonial or filamentous genus: natural units (cells, colonies, filaments) differ among the data sources",
  ifelse(kr$count_unit == "cell" & kr$colonial_evidence,
         "per-cell biovolume from records with more than one cell per natural unit (colonial genus)",
         ifelse(kr$count_unit == "cell",
                "per-cell biovolume from records that report the number of cells per natural unit",
                NA_character_))
)
kr$note <- ifelse(kr$product == "genus",
                  paste0("genus value pooling all records of the genus (Kremer et al. genus product)",
                         ifelse(is.na(kr$note), "", paste0("; ", kr$note))),
                  kr$note)
kr$record_id <- sprintf("kremer2014_%s_%04d", ifelse(kr$count_unit == "cell", "cell", "nu"),
                        ave(seq_len(nrow(kr)), kr$count_unit, FUN = seq_along))
kremer_typical <- data.frame(
  record_id = kr$record_id,
  taxon = kr$taxon,
  group = "phytoplankton",
  count_unit = kr$count_unit,
  life_stage = NA_character_,
  currency = "biovolume",
  unit = "um3",
  value = kr$value,
  geometric_mean = kr$geometric_mean,
  median = kr$median,
  value_min = kr$value_min,
  value_max = kr$value_max,
  sd_log10 = kr$sd_log10,
  n = as.integer(kr$n),
  source = "kremer2014",
  source_detail = sprintf("Kremer et al. 2014, %s product (%d data source%s)", kr$product,
                          as.integer(kr$n_sources), ifelse(kr$n_sources > 1, "s", "")),
  pooled = kr$product == "genus",
  note = kr$note,
  stringsAsFactors = FALSE
)

# Lineage of the accepted names (ITIS-based classification used by Kremer et al.).
kr_lin <- rbind(
  # Four species are labelled "Variety" in the species product (their records
  # were identified to variety); their accepted names are binomials.
  data.frame(name = sp_ag$accepted.taxa.name,
             rank = ifelse(is_binomial(sp_ag$accepted.taxa.name), "species", tolower(sp_ag$level)),
             kingdom = sp_ag$kingdom, phylum = ifelse(is.na(sp_ag$phylum), sp_ag$division, sp_ag$phylum),
             class = sp_ag$class, subclass = sp_ag$subclass, order = sp_ag$order,
             family = sp_ag$family, genus = sp_ag$genus,
             species = sp_ag$accepted.taxa.name, stringsAsFactors = FALSE),
  data.frame(name = ge_ag$accepted.taxa.name, rank = "genus",
             kingdom = ge_ag$kingdom, phylum = ifelse(is.na(ge_ag$phylum), ge_ag$division, ge_ag$phylum),
             class = ge_ag$class, subclass = ge_ag$subclass, order = ge_ag$order,
             family = ge_ag$family, genus = ge_ag$genus, species = NA_character_,
             stringsAsFactors = FALSE)
)
kr_lin <- kr_lin[kr_lin$rank %in% c("species", "genus"), , drop = FALSE]
kr_lin$source <- "kremer2014"

# Synonyms: names as reported by the data sources -> accepted names (ITIS 2014).
kr_syn <- unique(kr_taxa[, c("original.taxa.name", "accepted.taxa.name")])
kr_syn <- kr_syn[!is.na(kr_syn$original.taxa.name) & !is.na(kr_syn$accepted.taxa.name), ]

# =============================================================================
# 2. Moreau et al. (2021), Biodiversity Data Journal 9: e70214,
#    doi:10.3897/BDJ.9.e70214. Data: Dryad doi:10.5061/dryad.j3tx95xfg (CC0).
#    Units: lengths mm, WW and DW mg per individual.
# =============================================================================
mo <- read.delim(file.path(src, "Database_invertebrates.txt"), check.names = FALSE,
                 colClasses = "character", na.strings = "", fileEncoding = "latin1")
mo <- mo[, nzchar(names(mo)), drop = FALSE]
for (j in names(mo)) mo[[j]] <- blank_to_na(mo[[j]])
mo_refs <- read.delim(file.path(src, "Database_references.txt"), colClasses = "character",
                      fileEncoding = "CP1252", check.names = FALSE)
names(mo_refs) <- c("code", "reference")

# Documented corrections of the classification (see lineage_manual.csv).
mo$Genus[mo$Genus %in% "Pleuroxusv"] <- "Pleuroxus"
mo$Order[mo$Genus %in% c("Acanthocyclops", "Mesocyclops")] <- "Cyclopoida"
mo$Family[mo$Genus %in% "Moina"] <- "Moinidae"
mo$Family[mo$Family %in% "Coenargionidae"] <- "Coenagrionidae"
mo$Family[mo$Family %in% "Odonthoceridae"] <- "Odontoceridae"

mo_names <- normalize_taxon_names(mo$Taxa)
mo$taxon <- mo_names$taxon_clean
mo$taxon[mo$Taxa %in% "Pleuroxusv aduncus"] <- "Pleuroxus aduncus"

mo_lin <- data.frame(
  name = mo$taxon, rank = NA_character_, kingdom = "Animalia", phylum = mo$Phylum,
  subphylum = mo$SubPhylum, class = mo$Class, subclass = mo$SubClass, order = mo$Order,
  suborder = mo$SubOrder, family = mo$Family, subfamily = NA_character_, genus = mo$Genus,
  species = ifelse(is_binomial(mo$taxon), mo$taxon, NA_character_),
  stringsAsFactors = FALSE
)
mo_lin$source <- "moreau2021"

# Short citation of a Moreau reference code: "First author year".
short_citation <- function(code) {
  ref <- mo_refs$reference[match(code, mo_refs$code)]
  first <- sub("^([^ ,]+).*$", "\\1", ref)
  year <- sub("^[^(]*\\(([0-9]{4}[a-z]?)\\).*$", "\\1", ref)
  n_auth <- lengths(regmatches(sub("\\(.*$", "", ref), gregexpr(",", sub("\\(.*$", "", ref))))
  ifelse(is.na(ref), NA_character_,
         paste0(first, ifelse(n_auth >= 1, " et al. ", " "), year))
}

# Dry mass values. Point values ("12.3"), means with SD ("3.47+/-0.98") or
# with range ("12.1 [1.3-26.6]") are used; ranges ("1.1-101.0") are not a
# central value and are left out. Mollusc values are used only when the shell
# was removed (note "SR"); values with shell ("WS") or without information on
# the shell are left out.
parse_moreau_value <- function(x) {
  out <- data.frame(value = rep(NA_real_, length(x)), value_min = NA_real_, value_max = NA_real_,
                    kind = NA_character_, shell = NA_character_, stringsAsFactors = FALSE)
  for (i in seq_along(x)) {
    s <- x[i]
    if (is.na(s)) next
    if (grepl("\\(WS\\)", s)) out$shell[i] <- "with_shell"
    if (grepl("\\(SR\\)", s)) out$shell[i] <- "shell_removed"
    s <- gsub("\\((WS|SR)\\)", "", s)
    s <- gsub("\\s+", "", s)
    if (grepl("^[0-9.]+$", s)) {
      out$value[i] <- as.numeric(s); out$kind[i] <- "point"
    } else if (grepl("^[0-9.]+\\+/-[0-9.]+$", s)) {
      out$value[i] <- as.numeric(sub("\\+/-.*$", "", s)); out$kind[i] <- "mean_sd"
    } else if (grepl("^[0-9.]+\\[[0-9.]+-[0-9.]+\\]$", s)) {
      out$value[i] <- as.numeric(sub("\\[.*$", "", s))
      rng <- as.numeric(strsplit(sub("^.*\\[(.*)\\]$", "\\1", s), "-")[[1]])
      out$value_min[i] <- rng[1]; out$value_max[i] <- rng[2]; out$kind[i] <- "mean_range"
    } else if (grepl("^[0-9.]+-[0-9.]+$", s)) {
      out$kind[i] <- "range"
    } else {
      out$kind[i] <- "unparsed"
    }
  }
  out
}

dw <- cbind(mo, parse_moreau_value(mo[["DW (value)"]]))
dw <- dw[!is.na(mo[["DW (value)"]]), , drop = FALSE]
stopifnot(!any(dw$kind == "unparsed"))
mollusc <- dw$Phylum %in% "Mollusca"
excluded_moreau <- data.frame(
  taxon = dw$taxon, value = dw[["DW (value)"]], reference = dw[["DW (reference)"]],
  reason = ifelse(dw$kind == "range", "range, not a central value",
           ifelse(dw[["DW (reference)"]] == "41", "duplicate of Dumont et al. 1975 (used directly)",
           ifelse(mollusc & !(dw$shell %in% "shell_removed"), "mollusc value with shell or without information on the shell",
           ifelse(is.na(dw$taxon), "name not usable", NA_character_)))),
  stringsAsFactors = FALSE
)
keep <- is.na(excluded_moreau$reason)
excluded_moreau <- excluded_moreau[!keep, , drop = FALSE]
dw <- dw[keep, , drop = FALSE]

stage_map <- c(Adult = "adult", Larvae = "larva", `Larvae I` = "larva", `Larvae II` = "larva",
               `Larvae III` = "larva", `Larvae IV` = "larva", `Larvae V` = "larva",
               `Larvae VI` = "larva", Nymph = "nymph", Pupae = "pupa")
stopifnot(all(is.na(dw$LifeStage) | dw$LifeStage %in% names(stage_map)))
moreau_typical <- data.frame(
  record_id = sprintf("moreau2021_id%04d", as.integer(dw$ID)),
  taxon = dw$taxon,
  group = NA_character_,
  count_unit = "individual",
  life_stage = unname(stage_map[dw$LifeStage]),
  currency = "dry_mass",
  unit = "ug",
  value = dw$value * 1000,                          # mg -> ug
  geometric_mean = NA_real_,
  median = NA_real_,
  value_min = dw$value_min * 1000,
  value_max = dw$value_max * 1000,
  sd_log10 = NA_real_,
  n = NA_integer_,
  source = "moreau2021",
  source_detail = paste0("Moreau et al. 2021, ref. ", dw[["DW (reference)"]], ": ",
                         short_citation(dw[["DW (reference)"]])),
  note = ifelse(dw$shell %in% "shell_removed", "shell removed",
         ifelse(grepl("^Larvae [IV]+$", dw$LifeStage), paste("instar", sub("^Larvae ", "", dw$LifeStage)),
                NA_character_)),
  stringsAsFactors = FALSE
)

# =============================================================================
# 3. Dumont et al. (1975), Oecologia 19: 75-97, doi:10.1007/BF00377592
#    Transcribed from the page images (data-raw/dumont1975_*.csv); the
#    equations were checked against the representative values of Tables 4
#    and 6 (tests/testthat/fixtures/dumont1975_representative.csv).
# =============================================================================
du_t <- read.csv("data-raw/dumont1975_typical.csv", stringsAsFactors = FALSE,
                 colClasses = c(life_stage = "character", origin = "character"))
du_t$life_stage <- blank_to_na(du_t$life_stage)
dumont_typical <- data.frame(
  record_id = du_t$record_key,
  taxon = du_t$taxon,
  group = "zooplankton",
  count_unit = "individual",
  life_stage = du_t$life_stage,
  currency = "dry_mass",
  unit = "ug",
  value = du_t$dw_ug,
  geometric_mean = NA_real_,
  median = NA_real_,
  value_min = NA_real_,
  value_max = NA_real_,
  sd_log10 = NA_real_,
  n = as.integer(du_t$n_weighed),
  source = "dumont1975",
  source_detail = paste0("Dumont et al. 1975, ", du_t$table,
                         ifelse(is.na(du_t$mean_length_mm), "",
                                sprintf(" (mean length %s mm)", format(du_t$mean_length_mm)))),
  note = blank_to_na(du_t$note),
  stringsAsFactors = FALSE
)

du_l <- read.csv("data-raw/dumont1975_length_mass.csv", stringsAsFactors = FALSE,
                 colClasses = c(life_stage = "character", origin = "character",
                                exclude_reason = "character"))
du_l$life_stage <- blank_to_na(du_l$life_stage)
excl <- !is.na(blank_to_na(du_l$exclude_reason))
du_excluded <- du_l[excl, , drop = FALSE]
du_l <- du_l[!excl, , drop = FALSE]
# W[ug] = a L[um]^b  ->  W[mg] = a 1000^b L[mm]^b / 1000 ;  W[ug] = a L[mm]^b -> /1000.
a_mm_mg <- du_l$a_printed *
  ifelse(du_l$length_unit == "um", 1000^du_l$b_printed, 1) /
  ifelse(du_l$mass_unit == "ug", 1000, 1)
# Origin of the animals (codes of Tables 2 and 4-8): a-g in Table 8, h-l in
# the footnote of Table 2 (continued), p. 79. Turkish sites: Aci Gol is a
# saline lake in Anatolia (Asia); the location of Kanle Cesme is not given, so
# its continent is left undetermined. Pooled equations take the origin of the
# animals pooled when it is single.
dumont_origin <- data.frame(
  code = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)",
           "(l)"),
  continent = c("Europe", "Europe", "Europe", "Asia", "Europe", "Asia", "Africa", "Europe",
                NA, "Africa", "Europe", "Europe"),
  place = c("Lake Donk, Belgium", "University Pond, Gent, Belgium", "Sambre River, Belgium",
            "Pond at Patan, Nepal", "Wellemeersen, Belgium", "Aci Gol, Turkey",
            "Aguelman Sidi Ali, Morocco", "Putse Moer, Belgium",
            "Kanle Cesme, Turkey (continent not determined)", "Daiet Iffer, Morocco",
            "Pond at Teralfene, Belgium", "Lake Donk, Belgium (1967)"),
  stringsAsFactors = FALSE
)
du_code <- blank_to_na(du_l$origin)
du_code[du_l$record_key %in% "dumont1975_t2_07"] <- "(a)" # pooled Daphnia of origin (a)
du_code[du_l$record_key %in% "dumont1975_f7_02"] <- "(a)" # Euchlanis: Lake Donk (Table 8)
oi <- match(du_code, dumont_origin$code)
du_continent <- dumont_origin$continent[oi]
du_place <- ifelse(!is.na(oi), dumont_origin$place[oi], NA_character_)
du_place[du_l$record_key %in% "dumont1975_t2_06"] <- "Sambre River, Belgium" # Table 2 label
# Pooled copepod equations (Table 6): origins of the species pooled (Table 7).
du_place[du_l$record_key %in% sprintf("dumont1975_t6_%02d", 1:4)] <-
  "pooled from Belgium and Morocco (Table 7)"
du_place[du_l$record_key %in% "dumont1975_t6_05"] <- "pooled from Belgium and Turkey (Table 7)"
stopifnot(!anyNA(du_place))
dumont_length <- data.frame(
  equation_id = du_l$record_key,
  taxon = du_l$taxon,
  group = "zooplankton",
  life_stage = du_l$life_stage,
  length_measure = "body_length",
  a = a_mm_mg,
  b = du_l$b_printed,
  length_min = du_l$length_min_mm,
  length_max = du_l$length_max_mm,
  n = NA_real_,
  continent = du_continent,
  region_detail = du_place,
  source = "dumont1975",
  source_detail = paste0("Dumont et al. 1975, ", du_l$table),
  original_model = sprintf("W[%s] = %s * L[%s]^%s", du_l$mass_unit, format(du_l$a_printed),
                           du_l$length_unit, format(du_l$b_printed)),
  note = paste0(du_l$note, ifelse(is.na(blank_to_na(du_l$origin)), "",
                                  paste0("; origin ", du_l$origin))),
  stringsAsFactors = FALSE
)

# =============================================================================
# 4. van Leeuwen et al. (2025), Freshwater Biology 70: e70056,
#    doi:10.1111/fwb.70056 (CC BY 4.0). Table S2 (Figshare
#    doi:10.6084/m9.figshare.29085596): DM [mg] = a BL[mm]^b.
# =============================================================================
vl <- as.data.frame(readxl::read_excel(file.path(src, "vanLeeuwen2025_TableS2.xlsx"),
                                       sheet = "data"))
names(vl) <- gsub(" ", "_", names(vl))
for (j in c("order", "family", "genus", "species", "tax_resolution", "continent",
            "region_detail", "source_publication", "relation", "notes")) {
  vl[[j]] <- blank_to_na(vl[[j]])
}
vl$row <- seq_len(nrow(vl))
vl$note <- NA_character_
add_note <- function(i, txt) {
  vl$note[i] <<- ifelse(is.na(vl$note[i]), txt, paste(vl$note[i], txt, sep = "; "))
}
# Consistency of the standardised coefficients with the original ones.
ln_rows <- vl$relation %in% "linear model natural log"
log10_rows <- vl$relation %in% "linear model log10 base"
rel_dev <- rep(NA_real_, nrow(vl))
rel_dev[ln_rows] <- vl$a[ln_rows] / exp(vl$orginal_ln_a[ln_rows]) - 1
rel_dev[log10_rows] <- vl$a[log10_rows] / 10^vl$orginal_ln_a[log10_rows] - 1
vl_coef_check <- summary(abs(rel_dev))
# The standardised a is rounded in a few rows (e.g. 1e-04 for exp(-9.33)):
# when the original log-scale intercept and slope are given, the equation is
# recomputed from them, intercept and slope together. In 12 rows (Giustini et
# al. 2008) the standardised b of Table S2 is 3.5 while the original slopes
# range from 2.7 to 3.7; the original pair is used and the row is noted.
orig_ok <- (ln_rows | log10_rows) & !is.na(vl$orginal_ln_a) & !is.na(vl$orginal_b)
vl$a_used <- vl$a
vl$b_used <- vl$b
vl$a_used[orig_ok & ln_rows] <- exp(vl$orginal_ln_a[orig_ok & ln_rows])
vl$a_used[orig_ok & log10_rows] <- 10^vl$orginal_ln_a[orig_ok & log10_rows]
vl$b_used[orig_ok] <- vl$orginal_b[orig_ok]
b_diff <- orig_ok & abs(vl$b - vl$orginal_b) > 0.0005
vl_b_check <- vl[b_diff, c("row", "source_publication", "b", "orginal_b")]

add_note(which(b_diff), sprintf("b = %s in Table S2 differs from the original slope %s; the original intercept and slope are used",
                                format(vl$b[b_diff]), format(vl$orginal_b[b_diff])))

# Subfamilies written in the family column.
sub_rows <- vl$family %in% c("Chironominae", "Orthocladiinae", "Tanypodinae")
vl$subfamily <- NA_character_
vl$subfamily[sub_rows] <- vl$family[sub_rows]
vl$family[sub_rows] <- "Chironomidae"
partial <- vl$family %in% "Chironomidae (non-Tanypodinae)"
vl$family[partial] <- "Chironomidae"
add_note(which(partial), "Chironomidae excluding Tanypodinae")
# Genus and species that disagree: the species name is kept.
bad_gen <- !is.na(vl$species) & !is.na(vl$genus) & !startsWith(remove_diacritics(vl$species), remove_diacritics(vl$genus))
add_note(which(bad_gen), sprintf("genus '%s' in the source disagrees with species name", vl$genus[bad_gen]))
vl$genus[bad_gen] <- sub(" .*$", "", vl$species[bad_gen])
# Life stage from the notes.
vl$life_stage <- ifelse(vl$notes %in% "larvae", "larva", ifelse(vl$notes %in% "pupae", "pupa", NA_character_))
other_notes <- !is.na(vl$notes) & !vl$notes %in% c("larvae", "pupae")
add_note(which(other_notes), paste("source note:", vl$notes[other_notes]))
# Taxon of each equation at its stated resolution.
vl_name_raw <- ifelse(vl$tax_resolution %in% "species", vl$species,
               ifelse(vl$tax_resolution %in% "genus", vl$genus,
               ifelse(vl$tax_resolution %in% "subfamily", vl$subfamily,
               ifelse(vl$tax_resolution %in% "family", vl$family,
               ifelse(vl$tax_resolution %in% "order", vl$order, NA_character_)))))
# Rows without resolution, or whose stated resolution has no name (e.g.
# "species" with only the genus given): the finest name given.
no_res <- is.na(vl$tax_resolution) | is.na(vl_name_raw)
add_note(which(!is.na(vl$tax_resolution) & is.na(vl_name_raw)),
         sprintf("resolution '%s' in the source but no name at that rank", vl$tax_resolution[!is.na(vl$tax_resolution) & is.na(vl_name_raw)]))
vl_name_raw[no_res] <- ifelse(!is.na(vl$species[no_res]), vl$species[no_res],
                       ifelse(!is.na(vl$genus[no_res]), vl$genus[no_res], vl$family[no_res]))
vl_nn <- normalize_taxon_names(vl_name_raw)
vl$taxon <- vl_nn$taxon_clean
has_nn <- !is.na(vl_nn$name_note)
add_note(which(has_nn), vl_nn$name_note[has_nn])
# The single row labelled "class" names an order (Trichoptera): ambiguous, left out.
vl_excluded <- vl[vl$tax_resolution %in% "class" | is.na(vl$taxon), , drop = FALSE]
vl <- vl[!(vl$tax_resolution %in% "class" | is.na(vl$taxon)), , drop = FALSE]

vanleeuwen_length <- data.frame(
  equation_id = sprintf("vanleeuwen2025_row%03d", vl$row),
  taxon = vl$taxon,
  group = "macroinvertebrate",
  life_stage = vl$life_stage,
  length_measure = "body_length",
  a = vl$a_used,
  b = vl$b_used,
  length_min = vl$minimum_size,
  length_max = vl$maximum_size,
  n = vl$n,
  continent = vl$continent,
  region_detail = vl$region_detail,
  source = "vanleeuwen2025",
  source_detail = paste0("van Leeuwen et al. 2025, Table S2: ", vl$source_publication),
  original_model = ifelse(ln_rows[vl$row], sprintf("ln(DM) = %s + %s ln(BL)", format(vl$orginal_ln_a), format(vl$orginal_b)),
                   ifelse(log10_rows[vl$row], sprintf("log10(DM) = %s + %s log10(BL)", format(vl$orginal_ln_a), format(vl$orginal_b)),
                          sprintf("DM = %s BL^%s", format(vl$original_a), format(vl$orginal_b)))),
  note = vl$note,
  stringsAsFactors = FALSE
)
vl_lin <- data.frame(
  name = vl$taxon, rank = NA_character_, kingdom = "Animalia", phylum = "Arthropoda",
  subphylum = "Hexapoda", class = "Insecta", subclass = "Pterygota", order = vl$order,
  suborder = NA_character_, family = vl$family, subfamily = vl$subfamily,
  genus = ifelse(is.na(vl$genus), NA_character_, remove_diacritics(vl$genus)),
  species = ifelse(is_binomial(vl$taxon), vl$taxon, NA_character_),
  stringsAsFactors = FALSE
)
vl_lin$source <- "vanleeuwen2025"

# =============================================================================
# 5. Taxonomy: one lineage per name, built from parent links.
#    Precedence when sources disagree: manual (documented corrections and
#    WoRMS) > Moreau et al. 2021 > van Leeuwen et al. 2025 > Kremer et al. 2014.
# =============================================================================
manual <- read.csv("data-raw/lineage_manual.csv", stringsAsFactors = FALSE, na.strings = "")
manual$source <- "manual"
manual$species <- NA_character_

du_names <- unique(c(du_t$taxon, du_l$taxon))
du_lin <- empty_lineage(length(du_names))
du_lin$genus <- genus_of_binomial(du_names)
du_lin$species <- ifelse(is_binomial(du_names), du_names, NA_character_)
du_lin <- data.frame(name = du_names, rank = NA_character_, du_lin, source = "dumont1975",
                     stringsAsFactors = FALSE)

all_lin <- rbind(
  manual[, c("name", "rank", lineage_cols, "source")],
  du_lin[, c("name", "rank", lineage_cols, "source")],
  mo_lin[, c("name", "rank", lineage_cols, "source")],
  vl_lin[, c("name", "rank", lineage_cols, "source")],
  transform(kr_lin, subphylum = NA_character_, suborder = NA_character_,
            subfamily = NA_character_)[, c("name", "rank", lineage_cols, "source")]
)
precedence <- c(manual = 1L, moreau2021 = 2L, vanleeuwen2025 = 3L, dumont1975 = 4L,
                kremer2014 = 5L)

# Parent links: every non-missing rank of a lineage row points to the next
# non-missing rank above it.
edges <- list()
for (i in seq_len(nrow(all_lin))) {
  ranks_present <- lineage_cols[!is.na(unlist(all_lin[i, lineage_cols]))]
  vals <- unlist(all_lin[i, ranks_present])
  if (length(vals) == 0L) next
  for (k in seq_along(vals)) {
    edges[[length(edges) + 1L]] <- c(
      name = vals[[k]], rank = ranks_present[k],
      parent = if (k > 1L) vals[[k - 1L]] else NA_character_,
      source = all_lin$source[i]
    )
  }
}
edges <- as.data.frame(do.call(rbind, edges), stringsAsFactors = FALSE)
edges$prec <- precedence[edges$source]
edges <- edges[order(edges$name, edges$prec, method = "radix"), ]

# Choose rank and parent per name: best source; within it, the most frequent.
choose <- function(d) {
  d <- d[d$prec == min(d$prec), , drop = FALSE]
  key <- paste(d$rank, d$parent, sep = "\r")
  tab <- table(key)
  best <- names(tab)[tab == max(tab)][1]
  d[match(best, key), c("name", "rank", "parent", "source")]
}
by_name <- split(edges, edges$name)
nodes <- do.call(rbind, lapply(by_name, choose))
# The genus of a binomial is its first word, whatever a source says.
bin_nodes <- nodes$rank == "species" & is_binomial(nodes$name)
gen_of <- genus_of_binomial(nodes$name)
wrong_genus <- bin_nodes & !is.na(nodes$parent) & nodes$parent != gen_of & gen_of %in% nodes$name
genus_overrides <- data.frame(name = nodes$name[wrong_genus], source_parent = nodes$parent[wrong_genus],
                              used_parent = gen_of[wrong_genus], stringsAsFactors = FALSE)
nodes$parent[wrong_genus] <- gen_of[wrong_genus]

# Walk up the parents.
idx <- match(nodes$parent, nodes$name)
taxonomy <- empty_lineage(nrow(nodes))
for (i in seq_len(nrow(nodes))) {
  j <- i
  steps <- 0L
  while (!is.na(j)) {
    r <- nodes$rank[j]
    if (is.na(taxonomy[i, r])) taxonomy[i, r] <- nodes$name[j]
    j <- idx[j]
    steps <- steps + 1L
    if (steps > 20L) stop("cycle in the taxonomy at ", nodes$name[i])
  }
}
taxonomy <- data.frame(name = nodes$name, rank = nodes$rank, accepted_name = nodes$name,
                       taxonomy, lineage_source = nodes$source, stringsAsFactors = FALSE)
stopifnot(all(taxonomy$rank %in% rank_levels))

# Real conflicts: a source gives, for one of its names, an ancestor that
# differs from the final lineage at the same rank.
fin <- taxonomy[match(all_lin$name, taxonomy$name), lineage_cols]
conflicts <- list()
for (r in lineage_cols) {
  bad <- !is.na(all_lin[[r]]) & !is.na(fin[[r]]) & all_lin[[r]] != fin[[r]]
  if (any(bad)) {
    conflicts[[r]] <- unique(data.frame(name = all_lin$name[bad], rank = r,
                                        source_value = all_lin[[r]][bad],
                                        used_value = fin[[r]][bad],
                                        source = all_lin$source[bad], stringsAsFactors = FALSE))
  }
}
conflicts <- if (length(conflicts) > 0L) do.call(rbind, conflicts) else NULL

# Synonyms from Kremer et al. (names reported by the sources -> accepted).
syn <- kr_syn
syn$name <- normalize_taxon_names(syn$original.taxa.name)$taxon_clean
syn <- syn[!is.na(syn$name) & syn$name != syn$accepted.taxa.name &
           !syn$name %in% taxonomy$name & syn$accepted.taxa.name %in% taxonomy$name, ]
syn <- syn[!duplicated(syn$name) & !syn$name %in% syn$name[duplicated(syn$name)], ]
acc <- taxonomy[match(syn$accepted.taxa.name, taxonomy$name), ]
syn_rows <- acc
syn_rows$name <- syn$name
syn_rows$accepted_name <- syn$accepted.taxa.name
syn_rows$lineage_source <- "kremer2014 (synonym)"
taxonomy <- rbind(taxonomy, syn_rows)
taxonomy <- taxonomy[order(taxonomy$name, method = "radix"), ]
rownames(taxonomy) <- NULL
stopifnot(!anyDuplicated(taxonomy$name))

# =============================================================================
# 6. Final tables
# =============================================================================
dumont_typical$pooled <- FALSE
moreau_typical$pooled <- FALSE
typical_mass <- rbind(kremer_typical, dumont_typical, moreau_typical)
lin <- taxonomy[match(typical_mass$taxon, taxonomy$name), ]
stopifnot(!anyNA(lin$name))
typical_mass$rank <- lin$rank
zoo <- lin$suborder %in% "Cladocera" | lin$subclass %in% "Copepoda" | lin$phylum %in% "Rotifera"
typical_mass$group[is.na(typical_mass$group)] <- ifelse(zoo, "zooplankton", "macroinvertebrate")[is.na(typical_mass$group)]
typical_mass <- typical_mass[, c("record_id", "taxon", "rank", "group", "count_unit", "life_stage",
                                 "currency", "unit", "value", "geometric_mean", "median",
                                 "value_min", "value_max", "sd_log10", "n", "source",
                                 "source_detail", "pooled", "note")]
rownames(typical_mass) <- NULL

length_mass <- rbind(dumont_length, vanleeuwen_length)
lin <- taxonomy[match(length_mass$taxon, taxonomy$name), ]
stopifnot(!anyNA(lin$name))
length_mass$rank <- lin$rank
# Screening for unit or transcription errors: the prediction of each equation
# at the middle of its calibration range (5 mm when the range is unknown) is
# compared with the median prediction of the other equations of the same
# family (at least 3) at that length. 99% of the equations are within a factor
# of 15; those beyond a factor of 20 are left out as probable errors.
screen_factor <- 20
lm_family <- lin$family
lm_L <- ifelse(!is.na(length_mass$length_min) & !is.na(length_mass$length_max),
               (length_mass$length_min + length_mass$length_max) / 2, 5)
lm_ratio <- rep(NA_real_, nrow(length_mass))
for (i in seq_len(nrow(length_mass))) {
  others <- which(lm_family == lm_family[i] & !is.na(lm_family) & seq_along(lm_family) != i)
  if (length(others) >= 3L) {
    lm_ratio[i] <- length_mass$a[i] * lm_L[i]^length_mass$b[i] /
      stats::median(length_mass$a[others] * lm_L[i]^length_mass$b[others])
  }
}
lm_outlier <- !is.na(lm_ratio) & (lm_ratio > screen_factor | lm_ratio < 1 / screen_factor)
lm_screened_out <- data.frame(equation_id = length_mass$equation_id, taxon = length_mass$taxon,
                              family = lm_family, length_mm = lm_L,
                              ratio_to_family_median = signif(lm_ratio, 3),
                              source_detail = length_mass$source_detail,
                              stringsAsFactors = FALSE)[lm_outlier, , drop = FALSE]
length_mass <- length_mass[!lm_outlier, , drop = FALSE]
length_mass <- length_mass[, c("equation_id", "taxon", "rank", "group", "life_stage",
                               "length_measure", "a", "b", "length_min", "length_max", "n",
                               "continent", "region_detail", "source", "source_detail",
                               "original_model", "note")]
rownames(length_mass) <- NULL


# Only ASCII in the tables (as when they were distributed).
for (obj in c("taxonomy", "typical_mass", "length_mass")) {
  d <- get(obj)
  for (j in names(d)) if (is.character(d[[j]])) d[[j]] <- remove_diacritics(d[[j]])
  bad <- vapply(d, function(v) is.character(v) && any(grepl("[^ -~]", v, perl = TRUE)), logical(1))
  if (any(bad)) {
    for (j in names(d)[bad]) d[[j]] <- iconv(d[[j]], "UTF-8", "ASCII", sub = "?")
    message(obj, ": non-ASCII characters replaced in ", paste(names(d)[bad], collapse = ", "))
  }
  assign(obj, d)
}

# =============================================================================
# 7. Log and save
# =============================================================================
log_file <- "data-raw/references_log.md"
sink(log_file)
cat("# Build log of the reference tables\n\n")
cat("Built with R", R.version.string, "on", format(Sys.Date()), "\n\n")
cat("## Kremer et al. 2014\n\n")
cat("- Published log10 means reproduced from the raw files (max. difference < 1e-12).\n")
cat("- Species and genus rows (natural units) reproduce the published species and genus\n")
cat("  products: same n, log10 geometric means within 1e-6.\n")
cat("- Colonial genera (>=", 100 * colonial_min_fraction, "% of records with an explicit cell count\n")
cat("  report more than one cell per natural unit; per-cell values only from those records):",
    length(colonial_genera), "\n")
cat(strwrap(paste(colonial_genera, collapse = ", "), width = 78, prefix = "  "), sep = "\n")
cat("- Genera with occasional multi-cell records (all records with an explicit cell count used):",
    length(occasional_genera), "\n")
cat(strwrap(paste(occasional_genera, collapse = ", "), width = 78, prefix = "  "), sep = "\n")
print(table(kremer_typical$count_unit, ifelse(kremer_typical$pooled, "genus (pooled)", "species")))
cat("\n## Moreau et al. 2021: dry mass records left out\n\n")
print(table(excluded_moreau$reason))
cat("\n## Dumont et al. 1975: equations left out\n\n")
print(du_excluded[, c("record_key", "taxon", "exclude_reason")], row.names = FALSE)
cat("\n## van Leeuwen et al. 2025\n\n")
cat("- Relative deviation of standardised a from the original log-scale intercept:\n")
print(vl_coef_check)
cat("- Rows left out:\n")
print(vl_excluded[, c("row", "order", "family", "tax_resolution", "source_publication")], row.names = FALSE)
cat("- Rows whose standardised b differs from the original slope (original pair used):\n")
print(vl_b_check, row.names = FALSE)
cat("\n## Length-mass equations left out by the screening (factor >", screen_factor,
    "from the median of their family)\n\n")
print(lm_screened_out, row.names = FALSE)
cat("\n## Species whose genus in the source differs from the first word of the binomial\n\n")
print(genus_overrides, row.names = FALSE)
cat("\n## Taxonomy conflicts resolved by precedence\n\n")
if (!is.null(conflicts)) print(conflicts, row.names = FALSE) else cat("none\n")
cat("\n## Sizes\n\n")
cat("taxonomy:", nrow(taxonomy), "names (", sum(taxonomy$name != taxonomy$accepted_name), "synonyms)\n")
cat("typical_mass:", nrow(typical_mass), "records\n")
print(table(typical_mass$source, typical_mass$count_unit))
cat("length_mass:", nrow(length_mass), "equations\n")
print(table(length_mass$source))
sink()

# Exact numbers: doubles with 17 significant digits; text quoted.
write_csv_exact <- function(x, file) {
  qcols <- which(vapply(x, is.character, logical(1)))
  for (j in names(x)) {
    if (is.double(x[[j]])) {
      x[[j]] <- ifelse(is.na(x[[j]]), NA_character_, sprintf("%.17g", x[[j]]))
    }
  }
  write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8",
            quote = if (length(qcols)) qcols else FALSE)
}
for (obj in c("taxonomy", "typical_mass", "length_mass")) {
  write_csv_exact(get(obj), file.path("data-raw/derived", paste0(obj, "_0.1.0.csv")))
}
# Colonial genera of Kremer et al. (2014), read by model_tables.R: the package
# warns that the natural units of these genera in the reference are those of
# its sources (often single cells or small groups), not whole colonies.
write.csv(data.frame(genus = colonial_genera), "data-raw/derived/kremer_colonial_genera.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
