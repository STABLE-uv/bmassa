# Data of the hierarchical models, by domain, in their native currency (version 0.2.0)
#
# Builds, in data-raw/derived/:
#   model_typical.csv   one row per typical mass (or biovolume) used to fit the
#                       models without lengths, with its lineage and log value
#   model_lwr.csv       one row per length-mass equation used to fit the
#                       models with lengths, with its lineage, ln a, b and range
#   model_data_log.md   counts and what was left out, and why
# from the package data of 0.1.0 (typical_mass, length_mass), the tables of
# data-raw/derived/ (typical_records.csv, lwr_fitted.csv), the hierarchy
# (ref_taxa.csv) and the conversion models (conversion_*.csv). Run
# hierarchy.R and fit_conversions.R first.
#
# Models and native currencies (decided 2026-11-02): phytoplankton, biovolume
# (um3) per cell and per natural unit; invertebrates, dry mass (mg; shell-free
# in molluscs); fish and amphibians, wet mass (g); tadpoles (from version
# 0.6.0, domain "amphibian_larva"), dry mass (mg, without the gut), with a
# model of their typical total length (mm) for tadpoles without lengths. Records in another currency
# are converted with the conversion models (point estimate of the factor of
# the taxon; its variance is kept in `conversion_var` but not added to the
# fit): the wet masses of Hrycik et al. (2024), and the dry masses with shell
# of the CADSES molluscs. CADSES barnacles (dry mass with shell) are left out:
# there are no shell-free factors for them.
#
# Records not identified to species get their own nodes below the lowest rank
# they have ("Chironomus sp." of a source is a species of its own of
# Chironomus), so that every record has a full path from kingdom to species.
#
# Run from the package root, in a UTF-8 locale:  Rscript data-raw/model_data.R

source("data-raw/sources.R")

out_dir <- "data-raw/derived"
reg <- read_sources()
log_lines <- c("# Data of the hierarchical models", "",
               paste0("Built by `data-raw/model_data.R` on ", format(Sys.Date()), " (",
                      R.version.string, ")."), "")
say <- function(...) log_lines <<- c(log_lines, sprintf(...))
main <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

rt <- read.csv(file.path(out_dir, "ref_taxa.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
rt_key <- paste(rt$source_id, rt$taxon)
pkg <- new.env()
for (f in c("typical_mass", "length_mass")) {
  # Tables of 0.1.0 (data-raw/references.R).
  assign(f, read.csv(file.path(out_dir, paste0(f, "_0.1.0.csv")), stringsAsFactors = FALSE,
                     encoding = "UTF-8"), envir = pkg)
}
tr <- read.csv(file.path(out_dir, "typical_records.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
lf <- read.csv(file.path(out_dir, "lwr_fitted.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
Sys.unsetenv("BMASSA_FIT_ONLY_FUNCTIONS")
with_lineage <- function(d) add_lineage(d, rt)
stopifnot(!anyDuplicated(rt_key))
cn <- lapply(c(nodes = "conversion_nodes.csv", components = "conversion_components.csv",
               fixed = "conversion_fixed.csv"),
             function(f) read.csv(file.path(out_dir, f), stringsAsFactors = FALSE, encoding = "UTF-8"))
# ln of the conversion factor for rows with a lineage: deepest node of the
# conversion model, with the variance for a new population of the taxon.
convert_ln <- function(d, model, habitat) {
  tabs <- lapply(cn, function(x) x[x$model == model, , drop = FALSE])
  stopifnot(nrow(tabs$nodes) > 0)
  d$habitat <- habitat
  predict_tables(tabs, d, "typical", model_levels(tabs))
}
stage_group <- function(s) {
  s <- tolower(s)
  ifelse(is.na(s) | s == "", "unspecified",
  ifelse(s %in% c("adult", "female", "male", "female_ovigerous"), "adult",
  ifelse(s %in% c("larva", "nymph"), "larva",
  ifelse(s %in% c("pupa"), "pupa",
  ifelse(s %in% c("copepodite", "juvenile", "nauplius"), "juvenile", "unspecified")))))
}
# Stages of the copepods: fixed contrasts of copepodites and nauplii against
# the adults (records without a stage are taken as adults). Records of
# several stages pooled are not expected among the typical masses.
copepod_contrast <- function(s) {
  s <- tolower(s)
  ok <- is.na(s) | s %in% c("adult", "female", "male", "female_ovigerous", "copepodite", "nauplius")
  if (!all(ok)) stop("Copepod records of unexpected stages: ", paste(unique(s[!ok]), collapse = ", "))
  ifelse(s %in% "copepodite", "copepod_copepodite",
         ifelse(s %in% "nauplius", "copepod_nauplius", NA_character_))
}

# 1. Typical masses -------------------------------------------------------------------

tm <- pkg$typical_mass
typ <- list()
# Kremer et al. (2014): taxon means of raw records; the genus rows built from
# all records of a genus ("pooled") repeat the species and are left out. The
# value used is the geometric mean (the authors' means are of log10 values).
k <- tm$source == "kremer2014" & !tm$pooled
typ$kremer <- data.frame(record_id = tm$record_id[k], source_id = "kremer2014", taxon = tm$taxon[k],
                         model = ifelse(tm$count_unit[k] == "cell", "phytoplankton_cell", "phytoplankton_natural_unit"),
                         life_stage = NA_character_, stage = "unspecified", value = ifelse(is.na(tm$geometric_mean[k]), tm$value[k],
                                                               tm$geometric_mean[k]),
                         unit = "um3", n = tm$n[k], stringsAsFactors = FALSE)
k <- tr$source_id == "rimet2018" & tr$use_in_model
typ$rimet <- data.frame(record_id = tr$record_id[k], source_id = "rimet2018", taxon = tr$taxon[k],
                        model = ifelse(tr$count_unit[k] == "cell", "phytoplankton_cell", "phytoplankton_natural_unit"),
                        life_stage = NA_character_, stage = "unspecified", value = tr$mean[k], unit = "um3", n = tr$n[k],
                        stringsAsFactors = FALSE)
# Moreau et al. (2021) and Dumont et al. (1975): dry mass (ug in 0.1.0, here mg).
k <- tm$source %in% c("moreau2021", "dumont1975")
typ$invert <- data.frame(record_id = tm$record_id[k], source_id = tm$source[k], taxon = tm$taxon[k],
                         model = "invertebrate_individual", life_stage = tm$life_stage[k],
                         stage = stage_group(tm$life_stage[k]), value = tm$value[k] / 1000, unit = "mg", n = tm$n[k], stringsAsFactors = FALSE)
# Amphibians: adults (Santini et al. 2018 means, Trochet et al. 2014), wet mass in g.
k <- tr$source_id %in% c("santini2018", "trochet2014") & tr$use_in_model & tr$life_stage == "adult"
typ$amph <- data.frame(record_id = tr$record_id[k], source_id = tr$source_id[k], taxon = tr$taxon[k],
                       model = "amphibian_adult", life_stage = "adult", stage = "adult",
                       value = tr$mean[k], unit = "g",
                       n = tr$n[k], stringsAsFactors = FALSE)
# Hrycik et al. (2024): mean individual wet mass (formalin) of Great Lakes
# benthos, in mg, converted to dry mass (freshwater): soft bodies with the
# model wet_mass_to_dry_mass, molluscs (weighed with the shell) with
# wet_mass+shell_to_dry_mass (shell-free dry mass).
k <- tr$source_id == "hrycik2024" & tr$use_in_model
stopifnot(all(tr$currency[k] == "wet_mass"), all(tr$unit[k] == "mg"))
typ$hrycik <- data.frame(record_id = tr$record_id[k], source_id = "hrycik2024", taxon = tr$taxon[k],
                         model = "invertebrate_individual", life_stage = tr$life_stage[k],
                         stage = stage_group(tr$life_stage[k]), value = tr$mean[k], unit = "mg",
                         n = tr$n[k],
                         converted_from = ifelse(tr$shell_included[k], "wet_mass+shell", "wet_mass"),
                         stringsAsFactors = FALSE)
# Zooplankton (phase 4k; data-raw/typical_records.R): dry masses of
# crustaceans (Culver et al. 1985; Hébert et al. 2016, from typical lengths)
# and rotifers, and wet masses of rotifers from their volumes (Ejsmont-Karabin
# 1998; Bottrell et al. 1976, Table IX), converted to dry mass with the
# convention of US EPA (2016) used by convert_currency(): 10 %, 3.9 % in
# Asplanchna, without uncertainty.
k <- tr$source_id %in% c("culver1985", "bottrell1976", "ejsmontkarabin1998", "hebert2016") & tr$use_in_model
stopifnot(all(tr$unit[k] == "mg"), all(tr$currency[k] %in% c("dry_mass", "wet_mass")))
typ$zoo <- data.frame(record_id = tr$record_id[k], source_id = tr$source_id[k], taxon = tr$taxon[k],
                      model = "invertebrate_individual", life_stage = tr$life_stage[k],
                      stage = stage_group(tr$life_stage[k]), value = tr$mean[k], unit = "mg",
                      n = tr$n[k],
                      converted_from = ifelse(tr$currency[k] == "wet_mass", "wet_mass (rotifer)", NA),
                      stringsAsFactors = FALSE)
# Fish without lengths (decision of 2026-11-05; data-raw/typical_records.R):
# mean wet mass (g) of the populations of Arranz et al. (2022), of the
# species in the catch of Koning et al. (2024) and of the species in each
# sampling of NEON (decision of 2026-11-06). The samplings of a species at a
# NEON site are one population (`population`, completed with its species
# below, so that two names of one species are one population): the model has
# a random effect of population, and the residual is the variation between
# samplings of a population (data-raw/fit_models.R).
k <- tr$source_id %in% c("arranz2022", "koning2024", "neon_fish") & tr$use_in_model &
  tr$count_unit == "individual"
stopifnot(all(tr$currency[k] == "wet_mass"), all(tr$unit[k] == "g"))
typ$fish <- data.frame(record_id = tr$record_id[k], source_id = tr$source_id[k], taxon = tr$taxon[k],
                       model = "fish_individual", life_stage = NA_character_, stage = "unspecified",
                       value = tr$mean[k], unit = "g", n = tr$n[k],
                       population = ifelse(is.na(tr$sampling[k]), NA_character_,
                                           paste(tr$source_id[k], tr$site[k], sep = "|")),
                       stringsAsFactors = FALSE)
# Typical total lengths of tadpoles (decision of 2026-11-09;
# data-raw/typical_records.R): the mean total length (mm) of the tadpoles of
# each sample (species, wetland and date) of Knapp et al. (2021), before the
# metamorphic climax. The samples of a species in a wetland are one population
# (random effect), and the residual is the variation between its dates. The
# package predicts from this model the mean length of tadpoles without
# lengths, and their mass from it with the length-mass model of tadpoles.
k <- tr$source_id == "knapp2021" & tr$use_in_model & tr$length_type %in% "total_length"
stopifnot(all(tr$length_unit[k] == "mm"), all(tr$life_stage[k] == "larva"))
typ$tadpole_length <- data.frame(record_id = tr$record_id[k], source_id = "knapp2021", taxon = tr$taxon[k],
                                 model = "amphibian_larva_length", life_stage = "larva", stage = "larva",
                                 value = tr$length_mean[k], unit = "mm", n = tr$n[k],
                                 population = paste(tr$source_id[k], tr$site[k], sep = "|"),
                                 stringsAsFactors = FALSE)
for (j in setdiff(names(typ), c("hrycik", "zoo"))) typ[[j]]$converted_from <- NA_character_
for (j in setdiff(names(typ), c("fish", "tadpole_length"))) typ[[j]]$population <- NA_character_
typ <- do.call(rbind, typ)
rownames(typ) <- NULL
typ <- typ[is.finite(typ$value) & typ$value > 0, ]
typ <- with_lineage(typ)
typ$population <- ifelse(is.na(typ$population), NA_character_, paste(typ$population, typ$species, sep = "|"))
# Hrycik et al. (2024) give no life stage: their benthic samples hold the
# aquatic stages, so their insects are taken as larvae (larva group: larvae
# and nymphs); the other animals stay "unspecified".
hi <- typ$source_id == "hrycik2024" & typ$class == "Insecta"
typ$stage[hi] <- "larva"
# Stage groups are random effects in insects (larva, pupa, adult,
# unspecified), the class with records of several stages of the same taxa;
# the other invertebrates form one group ("non_insect"), so that the stage
# effects are not confounded with the taxa that have no stage. Copepods are
# in that group too, and their copepodites and nauplii have fixed contrasts
# (`contrast`) against the adults: most records of nauplii are of an order,
# and as random effects the nodes of their own taxa took up most of their
# difference (fit of 2026-11-04: nauplii only 2.5 times lighter than adults,
# against about 20 times in the data).
typ$contrast <- NA_character_
cop <- typ$model == "invertebrate_individual" & typ$class %in% "Copepoda"
typ$contrast[cop] <- copepod_contrast(typ$life_stage[cop])
# Method of the masses of rotifers (decision of 2026-11-05, after the test with
# Tuesday Lake, dev/validation/tuesday_lake.md): direct weighings (Dumont et
# al. 1975, Table 8; Bottrell et al. 1976, Table VIII, weighings compiled from
# several studies) are several times the dry masses from calculated volumes
# (Bottrell et al. 1976, Table IX; Ejsmont-Karabin 1998) with the convention
# of 10 % (3.9 % in Asplanchna). The weighings get a fixed contrast, and the
# model predicts for volumes, the usual method of field studies.
rw <- typ$model == "invertebrate_individual" & typ$phylum %in% "Rotifera" &
  typ$source_id %in% c("dumont1975", "bottrell1976") & is.na(typ$converted_from)
typ$contrast[rw] <- "rotifer_weighing"
# Kind of value of the benthic invertebrates of Moreau et al. (2021) (decision
# of 2026-11-05): masses of individual prey compiled from studies of the diet
# of birds and fish and of calorimetry (Driver et al. 1974; Driver 1981;
# Kaminski and Prince 1981; Boreham 1994; Pizzul et al. 2009, among others),
# often of large or adult animals, against the means of all the individuals of
# field samples of Hrycik et al. (2024), the quantity a biomass needs. For the
# same taxa they are about six times heavier
# (dev/validation/source_effect_invertebrates.md). They get a fixed contrast,
# and the model predicts for population means.
mb <- typ$model == "invertebrate_individual" & typ$source_id %in% "moreau2021" &
  !typ$phylum %in% "Rotifera" & !typ$class %in% c("Branchiopoda", "Copepoda", "Ostracoda") &
  is.na(typ$contrast)
typ$contrast[mb] <- "moreau_benthos"
ni <- typ$model == "invertebrate_individual" & typ$class != "Insecta"
typ$stage[ni] <- "non_insect"
typ$conversion_ln <- NA_real_
typ$conversion_var <- NA_real_
typ$conversion_level <- NA_character_
for (from in c("wet_mass", "wet_mass+shell")) {
  i <- which(typ$converted_from %in% from)
  p <- convert_ln(typ[i, ], paste0("invertebrate ", from, "_to_dry_mass"), "freshwater")
  typ$conversion_ln[i] <- p$eta
  typ$conversion_var[i] <- p$v_total
  typ$conversion_level[i] <- p$matched_level
}
rot <- which(typ$converted_from %in% "wet_mass (rotifer)")
stopifnot(all(typ$phylum[rot] == "Rotifera"))
typ$conversion_ln[rot] <- log(ifelse(typ$genus[rot] %in% "Asplanchna", 0.039, 0.1))
typ$conversion_var[rot] <- 0
typ$conversion_level[rot] <- "convention"
conv <- !is.na(typ$converted_from)
stopifnot(all(is.finite(typ$conversion_ln[conv])))
typ$value[conv] <- typ$value[conv] * exp(typ$conversion_ln[conv])
mismatch <- typ$domain != sub("_.*$", "", typ$model)
say("## Typical masses")
say("")
say("- %d records (%d without a COL taxon left out; %d whose taxon is outside the domain of the model left out%s).",
    nrow(typ) - sum(mismatch), attr(typ, "lost"), sum(mismatch),
    if (any(mismatch)) paste0(": ", paste(unique(paste(typ$source_id[mismatch], typ$taxon[mismatch])), collapse = ", ")) else "")
cl <- table(typ$converted_from[conv], factor(typ$conversion_level[conv], levels = c("root", main)))
cl <- cl[, colSums(cl) > 0, drop = FALSE]
say("- Hrycik et al. (2024): stage of the %d insect records set to larva (benthic samples; the data give no stage).", sum(hi))
say("- Hrycik et al. (2024): %d wet masses converted to dry mass (%s); level of the factor used: %s; median factor %.3f (soft bodies) and %.3f (molluscs, shell-free dry mass / wet mass with shell).",
    sum(conv), paste(sprintf("%s %d", rownames(cl), rowSums(cl)), collapse = ", "),
    paste(sprintf("%s %d", colnames(cl), colSums(cl)), collapse = ", "),
    exp(stats::median(typ$conversion_ln[typ$converted_from %in% "wet_mass"])),
    exp(stats::median(typ$conversion_ln[typ$converted_from %in% "wet_mass+shell"])))
typ <- typ[!mismatch, ]
typ$y <- log(typ$value)
tab <- table(typ$model, typ$source_id)
say("")
say("| Model | %s |", paste(colnames(tab), collapse = " | "))
say("|---|%s", paste(rep("---|", ncol(tab)), collapse = ""))
for (i in seq_len(nrow(tab))) say("| %s | %s |", rownames(tab)[i], paste(tab[i, ], collapse = " | "))
say("")
for (m in unique(typ$model)) {
  s <- typ[typ$model == m, ]
  say("- %s: %d records of %d species nodes (%d with more than one record), %d genera, %d families; stages: %s.",
      m, nrow(s), length(unique(s$species)), sum(table(s$species) > 1), length(unique(s$genus)),
      length(unique(s$family)), paste(sprintf("%s %d", names(table(s$stage)), as.integer(table(s$stage))), collapse = ", "))
  if (any(!is.na(s$contrast))) {
    say("  Fixed contrasts: %s (copepods: %d records, those of adults or without a stage, %d, the baseline; rotifers: %d records, the masses from volumes, %d, the baseline).",
        paste(sprintf("%s %d", names(table(s$contrast)), as.integer(table(s$contrast))), collapse = ", "),
        sum(s$class %in% "Copepoda"), sum(s$class %in% "Copepoda" & is.na(s$contrast)),
        sum(s$phylum %in% "Rotifera"), sum(s$phylum %in% "Rotifera" & is.na(s$contrast)))
  }
  if (any(!is.na(s$population))) {
    np <- table(s$population)
    say("  Populations sampled more than once (%s): %d records of %d populations (species and site), %d to %d samplings each (median %d); the other records are one population each.",
        paste(unique(s$source_id[!is.na(s$population)]), collapse = ", "), sum(np), length(np),
        min(np), max(np), as.integer(stats::median(np)))
  }
}
say("")
typ <- typ[, c("model", "record_id", "source_id", "taxon", "col_id", "domain", "lowest_rank", main,
               "life_stage", "stage", "contrast", "population", "value", "unit", "y", "n", "converted_from",
               "conversion_ln", "conversion_var", "conversion_level")]
write.csv(typ, file.path(out_dir, "model_typical.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 2. Length-mass equations ------------------------------------------------------------

# One model per domain and length type: an equation is only used with lengths
# of its own type. Lengths in mm; masses in mg
# (invertebrates, dry) or g (fish and amphibians, wet).
lm0 <- pkg$length_mass
stopifnot(all(lm0$length_measure == "body_length"))
# In Table S2 of van Leeuwen et al. (2025) the row of Ecdyonurus dispar
# (Baumgartner & Rothaupt 2003) has n = 3.33, minimum size 13.84 and no
# maximum: one column is shifted, and the range is read as 3.33-13.84 mm with
# n not given.
shifted <- lm0$equation_id == "vanleeuwen2025_row010"
stopifnot(sum(shifted) == 1, isTRUE(all.equal(lm0$n[shifted], 3.33)),
          isTRUE(all.equal(lm0$length_min[shifted], 13.84)), is.na(lm0$length_max[shifted]))
lm0$length_max[shifted] <- 13.84
lm0$length_min[shifted] <- 3.33
lm0$n[shifted] <- NA
lm0$note[shifted] <- paste0(ifelse(is.na(lm0$note[shifted]), "", paste0(lm0$note[shifted], "; ")),
                            "n = 3.33, minimum 13.84 and no maximum in Table S2 read as the range 3.33-13.84 mm")
eq <- list()
eq$pkg <- data.frame(equation_id = lm0$equation_id, source_id = lm0$source,
                     study = ifelse(lm0$source == "vanleeuwen2025",
                                    sub("^van Leeuwen et al\\. 2025, Table S2: ", "", lm0$source_detail),
                                    lm0$source),
                     taxon = lm0$taxon, kind = "equation", model_domain = "invertebrate",
                     length_type = lm0$length_measure, stage = lm0$life_stage,
                     ln_a = log(lm0$a), b = lm0$b, length_min = lm0$length_min,
                     length_max = lm0$length_max, length_mean = NA_real_, mass_mean = NA_real_,
                     n = lm0$n, sigma_ln = NA_real_, var_ln_a = NA_real_, var_b = NA_real_,
                     cov_ln_a_b = NA_real_, scope = "published", note = lm0$note,
                     stringsAsFactors = FALSE)
say("## Length-mass equations")
say("")
say("- One model per domain and length type; masses in mg (invertebrates, dry mass) or g (fish and amphibians, wet mass), lengths in mm.")
say("- 0.1.0 `length_mass`: %d equations, %d of them without a length range (their range is imputed below). The row of *Ecdyonurus dispar* in Table S2 of van Leeuwen et al. (2025) has a shifted column (n = 3.33, minimum 13.84, no maximum) and is read as the range 3.33-13.84 mm.",
    nrow(lm0), sum(!(is.finite(lm0$length_min) & is.finite(lm0$length_max))))
# Equations fitted from individual data: the equations of each site where a
# taxon has them, otherwise the pooled one; not the flagged ones (r2 < 0.5 or b
# outside 1.5-4.5).
keep_src <- c("cadses_boxcorer", "cadses_leafpack", "coughlan2021", "arranz2022", "koning2024",
              "neon_fish", "santini2018", "knapp2021")
f <- lf[lf$source_id %in% keep_src, ]
native <- (f$mass_currency == "dry_mass" & !f$shell_included) | f$mass_currency == "wet_mass"
native <- native & !(f$source_id %in% c("cadses_boxcorer", "cadses_leafpack", "coughlan2021") &
                       f$mass_currency == "wet_mass")
shelled <- f$source_id %in% c("cadses_boxcorer", "cadses_leafpack") & f$mass_currency == "dry_mass" &
  f$shell_included
say("- Fitted equations: %d of %d in the native currency and %d of the CADSES with shell (converted below); the other currencies of Coughlan et al. (2021) are left out, as the same animals have shell-free dry mass; %d flagged left out.",
    sum(native), nrow(f), sum(shelled), sum((native | shelled) & !is.na(f$quality_flag)))
f <- f[(native | shelled) & is.na(f$quality_flag), ]
grp <- paste(f$source_id, f$taxon, f$length_type, f$mass_currency)
has_pop <- tapply(f$scope == "population", grp, any)
f <- f[(f$scope == "population") | !has_pop[grp], ]
stopifnot(all(f$length_unit == "mm"),
          all(f$mass_unit[f$mass_currency == "dry_mass"] == "mg"),
          all(f$mass_unit[f$mass_currency == "wet_mass"] == "g"))
eq$fit <- data.frame(equation_id = f$equation_id, source_id = f$source_id, study = f$source_id,
                     taxon = f$taxon, kind = "equation",
                     model_domain = ifelse(f$source_id == "knapp2021", "amphibian_larva",
                                           ifelse(f$mass_currency == "dry_mass", "invertebrate",
                                                  ifelse(f$source_id == "santini2018", "amphibian",
                                                         "fish"))),
                     length_type = f$length_type, stage = f$life_stage, ln_a = f$ln_a, b = f$b,
                     length_min = f$length_min, length_max = f$length_max, length_mean = NA_real_,
                     mass_mean = NA_real_, n = f$n, sigma_ln = f$sigma_ln, var_ln_a = f$var_ln_a,
                     var_b = f$var_b, cov_ln_a_b = f$cov_ln_a_b, scope = f$scope, note = NA_character_,
                     converted_from = ifelse(f$shell_included & f$mass_currency == "dry_mass",
                                             "dry_mass+shell", NA_character_),
                     stringsAsFactors = FALSE)
say("- Bagley et al. (2022) left out: masses after ethanol preservation, with no conversion to fresh mass.")
# Amphibians: the mean snout-vent length and mean mass of the adults of a
# species (Santini et al. 2018; Trochet et al. 2014) are one point each of the
# relation between the mean length and the mean mass of a population, which is
# what the model predicts. The published equations of Santini et al. (2018,
# Table 3) were fitted to the same means and are left out (they are compared
# with the model in the validation).
k <- tr$source_id %in% c("santini2018", "trochet2014") & tr$use_in_model & tr$life_stage == "adult" &
  tr$length_type %in% "snout_vent_length" & is.finite(tr$length_mean) & tr$length_mean > 0
stopifnot(all(tr$unit[k] == "g"), all(tr$length_unit[k] == "mm"))
eq$means <- data.frame(equation_id = tr$record_id[k], source_id = tr$source_id[k], study = tr$source_id[k],
                       taxon = tr$taxon[k], kind = "mean", model_domain = "amphibian",
                       length_type = "snout_vent_length", stage = "adult", ln_a = NA_real_, b = NA_real_,
                       length_min = NA_real_, length_max = NA_real_, length_mean = tr$length_mean[k],
                       mass_mean = tr$mean[k], n = tr$n[k], sigma_ln = NA_real_, var_ln_a = NA_real_,
                       var_b = NA_real_, cov_ln_a_b = NA_real_, scope = "species_mean",
                       note = tr$note[k], converted_from = NA_character_, stringsAsFactors = FALSE)
say("- Amphibians: %d adult means of snout-vent length and mass (Santini et al. 2018; Trochet et al. 2014) enter as one point each.",
    sum(k))
eq$pkg$converted_from <- NA_character_
# Published equations of zooplankton (phase 4k; data-raw/published_lwr.R):
# mass in mg from body length in mm, those marked for use (see the reasons in
# the notes of lwr_published.csv). Dry mass, except the rotifers of
# Ejsmont-Karabin (1998), in wet mass and converted below. Compilations give
# the original study of each equation.
lp <- read.csv(file.path(out_dir, "lwr_published.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
lp <- lp[lp$use_in_model, ]
wet <- lp$mass_currency == "wet_mass"
stopifnot(all(lp$mass_currency[!wet] == "dry_mass"), all(lp$source_id[wet] == "ejsmontkarabin1998"),
          all(lp$mass_unit == "mg"), all(lp$length_unit == "mm"), all(lp$length_type == "body_length"))
eq$zoo <- data.frame(equation_id = lp$equation_id, source_id = lp$source_id,
                     study = ifelse(is.na(lp$study), lp$source_id, lp$study),
                     taxon = lp$taxon, kind = "equation", model_domain = "invertebrate",
                     length_type = lp$length_type, stage = lp$life_stage, ln_a = lp$ln_a, b = lp$b,
                     length_min = lp$length_min, length_max = lp$length_max, length_mean = NA_real_,
                     mass_mean = NA_real_, n = lp$n, sigma_ln = lp$sigma_ln, var_ln_a = lp$var_ln_a,
                     var_b = lp$var_b, cov_ln_a_b = lp$cov_ln_a_b, scope = lp$scope, note = lp$note,
                     converted_from = ifelse(wet, "wet_mass (rotifer)", NA_character_),
                     stringsAsFactors = FALSE)
say("- Published equations of zooplankton: %d (%s); %d studies.", nrow(eq$zoo),
    paste(sprintf("%s %d", names(table(eq$zoo$source_id)), as.integer(table(eq$zoo$source_id))),
          collapse = ", "), length(unique(eq$zoo$study)))
eq <- do.call(rbind, eq)
rownames(eq) <- NULL
eq <- with_lineage(eq)
# The body length of copepod nauplii has relations of its own: for the same
# length they are much heavier than copepodites and adults (another shape).
naup <- eq$model_domain == "invertebrate" & eq$class %in% "Copepoda" & eq$stage %in% "nauplius" &
  eq$length_type == "body_length"
eq$length_type[naup] <- "nauplius_body_length"
say("- Copepod nauplii: %d equations with a length type of their own (nauplius_body_length).", sum(naup))
eq$model <- paste(eq$model_domain, eq$length_type, sep = "_")
# CADSES with shell: shell-free dry mass = dry mass with shell x a factor of
# the taxon. Bivalves: (ash-free dry mass with shell / dry mass with shell) of
# the taxon, from the conversion model dry_mass+shell_to_afdm+shell, which
# holds the CADSES ash measurements (transitional waters), x (shell-free dry
# mass / ash-free dry mass with shell), the mean of the model
# afdm+shell_to_dry_mass (freshwater bivalves of Coughlan et al. 2021; ln
# factor of a new species from its root and all its variance components).
# Snails (from version 0.6.0): dry mass without shell / dry mass with shell of
# the taxon, from the model dry_mass+shell_to_dry_mass, which holds the
# shares of the shell of snails (conversions.R, sections 10 and 11), by the
# lowest level of the taxon in the model (species, family, order or the
# class): the factor of the bivalves underestimates the soft tissue of snails
# (decision of 2026-11-09). That share does not follow the shell of each
# species, so it is capped by the ash of the taxon: the soft tissue cannot
# hold more organic matter than the whole snail, so its share is at most a /
# r_t, with a the ash-free dry mass with shell / dry mass with shell of the
# taxon (model dry_mass+shell_to_afdm+shell, as for the bivalves) and r_t the
# ash-free dry mass / dry mass of its soft tissue (model dry_mass_to_afdm),
# as if the shell had no organic matter (decision of 2026-11-10). The capped
# factors take the variance of both models. Barnacles are left out.
eq$conversion_ln <- NA_real_
eq$conversion_var <- NA_real_
eq$conversion_level <- NA_character_
sh <- which(eq$converted_from %in% "dry_mass+shell")
barn <- sh[eq$phylum[sh] != "Mollusca"]
if (length(barn)) say("- CADSES equations with shell left out (no shell-free factor): %s.",
                      paste(unique(paste0(eq$taxon[barn], " (", eq$class[barn], ")")), collapse = ", "))
sh <- setdiff(sh, barn)
snail <- sh[eq$class[sh] %in% "Gastropoda"]
biv <- setdiff(sh, snail)
stopifnot(all(eq$class[biv] %in% "Bivalvia"))
p1 <- convert_ln(eq[biv, ], "invertebrate dry_mass+shell_to_afdm+shell", "transitional")
t2 <- lapply(cn, function(x) x[x$model == "invertebrate afdm+shell_to_dry_mass", , drop = FALSE])
r2 <- t2$nodes[t2$nodes$level == "root", ]
v2 <- r2$path_var + sum(t2$components$value[t2$components$quantity == "variance"])
eq$conversion_ln[biv] <- p1$eta + r2$path_intercept
eq$conversion_var[biv] <- p1$v_total + v2
eq$conversion_level[biv] <- p1$matched_level
p3 <- convert_ln(eq[snail, ], "invertebrate dry_mass+shell_to_dry_mass", "transitional")
pa <- convert_ln(eq[snail, ], "invertebrate dry_mass+shell_to_afdm+shell", "transitional")
pr <- convert_ln(eq[snail, ], "invertebrate dry_mass_to_afdm", "transitional")
cap_ln <- pa$eta - pr$eta
capped <- p3$eta > cap_ln
eq$conversion_ln[snail] <- ifelse(capped, cap_ln, p3$eta)
eq$conversion_var[snail] <- ifelse(capped, pa$v_total + pr$v_total, p3$v_total)
eq$conversion_level[snail] <- ifelse(capped, paste0(pa$matched_level, " (cap by the ash)"), p3$matched_level)
eq$ln_a[sh] <- eq$ln_a[sh] + eq$conversion_ln[sh]
level_counts <- function(x) paste(sprintf("%s %d", names(table(x)), as.integer(table(x))), collapse = ", ")
say("- CADSES bivalves with shell: %d equations converted to shell-free dry mass; level of the ash factor used: %s; median factor %.3f (ash-free with shell / dry with shell, %.3f; shell-free / ash-free with shell, %.3f).",
    length(biv), level_counts(p1$matched_level), exp(stats::median(eq$conversion_ln[biv])),
    exp(stats::median(p1$eta)), exp(r2$path_intercept))
say("- CADSES snails with shell: %d equations converted to shell-free dry mass with the share of the soft tissue in the dry mass; level used: %s; share %.3f-%.3f; %d equations capped by the ash of their taxon (%s); factors %.3f-%.3f (median %.3f).",
    length(snail), level_counts(p3$matched_level), exp(min(p3$eta)), exp(max(p3$eta)), sum(capped),
    paste(unique(sprintf("*%s* %.3f instead of %.3f", eq$taxon[snail][capped], exp(cap_ln[capped]),
                         exp(p3$eta[capped]))), collapse = ", "),
    exp(min(eq$conversion_ln[snail])), exp(max(eq$conversion_ln[snail])),
    exp(stats::median(eq$conversion_ln[snail])))
eq <- eq[setdiff(seq_len(nrow(eq)), barn), ]
# Rotifers of Ejsmont-Karabin (1998): wet mass to dry mass with the convention
# of US EPA (2016) used by convert_currency() (10 %, 3.9 % in Asplanchna),
# without uncertainty.
rot <- which(eq$converted_from %in% "wet_mass (rotifer)")
stopifnot(all(eq$phylum[rot] == "Rotifera"))
eq$conversion_ln[rot] <- log(ifelse(eq$genus[rot] %in% "Asplanchna", 0.039, 0.1))
eq$conversion_var[rot] <- 0
eq$conversion_level[rot] <- "convention"
eq$ln_a[rot] <- eq$ln_a[rot] + eq$conversion_ln[rot]
say("- Rotifers: %d equations of wet mass converted to dry mass with the convention of US EPA (2016).",
    length(rot))
# Tadpoles: the domain of their taxa is "amphibian"; their equations form
# the model of their own domain, "amphibian_larva".
mismatch <- eq$domain != sub("_larva$", "", eq$model_domain)
say("- %d rows (%d without a COL taxon left out; %d whose taxon is outside the domain of the model left out%s).",
    nrow(eq) - sum(mismatch), attr(eq, "lost"), sum(mismatch),
    if (any(mismatch)) paste0(": ", paste(unique(paste(eq$source_id[mismatch], eq$taxon[mismatch])), collapse = ", ")) else "")
eq <- eq[!mismatch, ]
eq$domain <- eq$model_domain
say("- Tadpoles (Knapp et al. 2021): %d equations of total length and dry mass (by wetland where the species has them).",
    sum(eq$domain == "amphibian_larva"))

# Equations without a length range (published compilations that did not report
# it): the range is taken from the equations of the same model and genus, or
# else family, or else order, with a range (medians of the log limits); with
# one limit given, the other is set with the median log width, and a range of
# zero width (minimum = maximum) is widened around its length with that width.
# Without relatives up to the order the equation is left out.
eq$range_source <- ifelse(eq$kind == "mean", "mean",
                          ifelse(is.finite(eq$length_min) & is.finite(eq$length_max) &
                                   eq$length_max > eq$length_min, "reported", NA))
lo <- log(eq$length_min)
hi <- log(eq$length_max)
ok <- eq$range_source %in% "reported"
for (i in which(is.na(eq$range_source))) {
  rel <- FALSE
  for (r in c("genus", "family", "order")) {
    rel <- ok & eq$model == eq$model[i] & eq[[r]] == eq[[r]][i]
    if (any(rel)) break
  }
  if (!any(rel)) {
    eq$range_source[i] <- "none"
    next
  }
  w <- stats::median(hi[rel] - lo[rel])
  if (is.finite(lo[i]) && is.finite(hi[i])) {
    mid <- (lo[i] + hi[i]) / 2
    lo[i] <- mid - w / 2
    hi[i] <- mid + w / 2
    eq$range_source[i] <- paste0("width_from_", r)
  } else if (is.finite(lo[i])) {
    hi[i] <- lo[i] + w
    eq$range_source[i] <- paste0("maximum_from_", r)
  } else if (is.finite(hi[i])) {
    lo[i] <- hi[i] - w
    eq$range_source[i] <- paste0("minimum_from_", r)
  } else {
    lo[i] <- stats::median(lo[rel])
    hi[i] <- stats::median(hi[rel])
    eq$range_source[i] <- paste0("range_from_", r)
  }
}
eq$length_min <- exp(lo)
eq$length_max <- exp(hi)
say("- Length ranges of the equations: %s.",
    paste(sprintf("%s %d", names(table(eq$range_source[eq$kind == "equation"])),
                  as.integer(table(eq$range_source[eq$kind == "equation"]))), collapse = ", "))
if (any(eq$range_source == "none")) {
  say("- Left out, without a range and without relatives with one: %s.",
      paste(sprintf("%s (%s)", eq$taxon[eq$range_source == "none"], eq$source_id[eq$range_source == "none"]),
            collapse = "; "))
}
eq <- eq[eq$range_source != "none", ]
stopifnot(all(eq$kind == "mean" | eq$length_max > eq$length_min))
tab <- table(eq$model, eq$source_id)
say("")
say("| Model | %s |", paste(colnames(tab), collapse = " | "))
say("|---|%s", paste(rep("---|", ncol(tab)), collapse = ""))
for (i in seq_len(nrow(tab))) say("| %s | %s |", rownames(tab)[i], paste(tab[i, ], collapse = " | "))
say("")
for (m in sort(unique(eq$model))) {
  s <- eq[eq$model == m, ]
  say("- %s: %d rows (%d equations, %d means) of %d species nodes (%d with more than one row), %d genera, %d families, %d studies.",
      m, nrow(s), sum(s$kind == "equation"), sum(s$kind == "mean"), length(unique(s$species)),
      sum(table(s$species) > 1), length(unique(s$genus)), length(unique(s$family)),
      length(unique(s$study)))
}
eq <- eq[, c("model", "equation_id", "source_id", "study", "taxon", "col_id", "domain", "lowest_rank",
             main, "kind", "length_type", "stage", "ln_a", "b", "length_min", "length_max",
             "range_source", "length_mean", "mass_mean", "n", "sigma_ln", "var_ln_a", "var_b",
             "cov_ln_a_b", "scope", "converted_from", "conversion_ln", "conversion_var",
             "conversion_level", "note")]
write.csv(eq, file.path(out_dir, "model_lwr.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(log_lines, file.path(out_dir, "model_data_log.md"), useBytes = TRUE)
cat(log_lines, sep = "\n")
