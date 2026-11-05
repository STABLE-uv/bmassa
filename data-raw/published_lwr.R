# Published length-mass equations transcribed by hand (version 0.2.0)
#
# Builds, in data-raw/derived/:
#   lwr_published.csv       equations in the layout of lwr_fitted.csv (natural
#                           logarithms, lengths in mm)
#   lwr_published_log.md    plausibility checks
# from the transcriptions in data-raw/*_table*.csv. This is also where
# equations taken from original publications (fish, in place of FishBase) will
# be added.
#
# Run from the package root, after data-raw/typical_records.R, in a UTF-8
# locale:  Rscript data-raw/published_lwr.R
# Packages used only here: readxl.

source("data-raw/sources.R")

out_dir <- "data-raw/derived"
dir.create(out_dir, showWarnings = FALSE)
reg <- read_sources()
ln10 <- log(10)
eqs <- list()
log_lines <- c("# Published length-mass equations", "",
               paste0("Built by `data-raw/published_lwr.R` on ", format(Sys.Date()), " (",
                      R.version.string, ")."), "")
say <- function(...) log_lines <<- c(log_lines, sprintf(...))

# Santini et al. (2018), Table 3 (page 41): body mass (g) from snout-vent
# length (mm), log10 scale, phylogenetic generalized least squares on species
# means of adults. The values in parentheses are standard errors (caption:
# coefficients and their standard errors in Table S7). Transcribed from the
# page image; no covariance between intercept and slope is given.
if (use_source("santini2018", reg)) {
  t3 <- read.csv("data-raw/santini2018_table3.csv", stringsAsFactors = FALSE)
  t3$condition[!nzchar(t3$condition)] <- NA_character_
  eqs$santini <- data.frame(
    equation_id = t3$equation_id, source_id = "santini2018", study = "santini2018", taxon = t3$taxon,
    rank = t3$rank,
    condition = t3$condition, scope = "pooled", length_type = "snout_vent_length",
    length_unit = "mm", mass_currency = "wet_mass", mass_unit = "g", shell_included = FALSE,
    preservation = "unknown", life_stage = "adult", fit = "PGLS on species means",
    ln_a = t3$log10_a * ln10, b = t3$b, var_ln_a = (t3$log10_a_se * ln10)^2,
    var_b = t3$b_se^2, cov_ln_a_b = NA_real_, sigma_ln = NA_real_,
    rmse_ln = t3$rmse_log10 * ln10, r2 = t3$adj_r2, n = t3$n_species, n_unit = "species",
    length_min = NA_real_, length_max = NA_real_,
    note = paste0(ifelse(is.na(t3$note) | !nzchar(t3$note), "", paste0(t3$note, "; ")),
                  "species means of 15.99-262 mm (frogs) and 23.99-542 mm (salamanders) snout-vent length in the whole data set"),
    stringsAsFactors = FALSE)
}

# Zooplankton (phase 4k): masses in mg, lengths in mm. Dry mass of
# crustaceans (ug in the sources, ln a - ln 1000); wet mass of rotifers
# (Ejsmont-Karabin 1998), converted to dry mass in model_data.R.
# `use_in_model` is FALSE, with the reason in `note`, for equations of lengths
# the reference does not use, for those already in it through another source,
# for those pooling nauplii with later stages (other shapes), and for slopes
# outside 1.5-4.5 (the rule of the fitted equations, data-raw/fit_lwr.R). The
# equations of Santini et al. (2018) are not used in the model (it is fitted
# to the means they were fitted to).
zoo_csv <- function(f) {
  read.csv(file.path("data-raw", f), stringsAsFactors = FALSE, encoding = "UTF-8",
           na.strings = c("", "NA"))
}
rank_of <- function(taxon) {
  ifelse(grepl("^[A-Z][a-z-]+ [a-z-]+", taxon), "species",
         ifelse(grepl("(idae|ida|oda)$", taxon), "higher", "genus"))
}
add_note <- function(a, b) ifelse(is.na(a) | !nzchar(a), b, ifelse(is.na(b), a, paste0(a, "; ", b)))
var_from_ci <- function(ci, n) (ci / stats::qt(0.975, pmax(n - 2, 1)))^2
# Equations pooling nauplii with copepodites or adults.
pools_nauplii <- function(stage) grepl("nauplius", stage) & grepl("\\|", stage)

# Culver et al. (1985), Table 3 (journal p. 1384): W = a L^b, W in ug dry mass
# (4 % formalin), L in mm; fitted to ln W on ln L over weighings of 10-40
# animals of the same length (n = number of weighings); a and b checked
# against the weights of Table 2 and the values quoted in the text. Copepods:
# body length to the base of the caudal rami (nauplii without appendages);
# Daphnia: total length (anterior margin of the helmet to the base of the tail
# spine), the body length of the other sources (top of the head to the base
# of the tail spine); their standard length (from the eye, which leaves the
# helmet out) is not used; other cladocerans: head to the end of the carapace.
if (use_source("culver1985", reg)) {
  t3 <- zoo_csv("culver1985_table3.csv")
  lt <- ifelse(t3$length_measure == "total_length", "body_length", t3$length_measure)
  use <- t3$length_measure %in% c("body_length", "total_length")
  note <- add_note(paste0(t3$row_printed, " (", t3$taxon_printed, ")"), t3$taxon_note)
  note <- add_note(note, ifelse(t3$length_measure == "total_length",
                                "total length (helmet to the base of the tail spine) taken as body length", NA))
  note <- add_note(note, ifelse(use, NA, "standard length (eye to the base of the tail spine): not used"))
  note <- add_note(note, ifelse(t3$record_key == "culver1985_t3_04",
                                "printed RMS equals the standard error of b: its RMS is not used", NA))
  eqs$culver <- data.frame(
    equation_id = t3$record_key, source_id = "culver1985", study = "culver1985", taxon = t3$taxon,
    rank = rank_of(t3$taxon), condition = t3$row_printed, scope = "population", length_type = lt,
    length_unit = "mm", mass_currency = "dry_mass", mass_unit = "mg", shell_included = FALSE,
    preservation = "formalin", life_stage = t3$life_stage, fit = "OLS of ln W on ln L (weighings)",
    ln_a = log(t3$a) - log(1000), b = t3$b, var_ln_a = NA_real_,
    var_b = var_from_ci(t3$b_ci95, t3$n_weights), cov_ln_a_b = NA_real_,
    sigma_ln = ifelse(t3$record_key == "culver1985_t3_04", NA_real_, sqrt(t3$rms)),
    rmse_ln = NA_real_, r2 = t3$r2, n = t3$n_weights, n_unit = "weighings",
    length_min = t3$length_min_mm, length_max = t3$length_max_mm, use_in_model = use, note = note,
    stringsAsFactors = FALSE)
}

# Ejsmont-Karabin (1998), Table II (pp. 516-521): W = c a^3, W in mg wet
# weight (body volume in mm3, density 1), a = body length in mm (lorica or
# body, without spines, toes or foot, as given for each species), with the
# weight of the appendices (column 4, % of the body weight) added, as in the
# standard weights of the table; the mucus sheaths are left out. Range: the
# lengths of column 5. Equations of the body width (b^3, a b^2) are not used.
# Lindia truncata is printed "W = 0.030 b3", but its standard weight is that
# equation applied to its body length: taken as a^3. Geometric formulas (b =
# 3), not fitted: no errors.
if (use_source("ejsmontkarabin1998", reg)) {
  ek <- zoo_csv("ejsmontkarabin1998_table2.csv")
  lindia <- ek$record_key == "ejsmontkarabin1998_t2_145"
  stopifnot(sum(lindia) == 1, ek$taxon[lindia] == "Lindia truncata", ek$variables[lindia] == "b3")
  ek$variables[lindia] <- "a3"
  app <- ifelse(is.na(ek$appendices_pct), 0, ek$appendices_pct)
  use <- ek$variables == "a3"
  note <- paste0(ifelse(ek$state == "C", "contracted on preservation", "living"), "; a = ",
                 ek$std_length_type, ifelse(app > 0, sprintf("; appendices %s %% added", ek$appendices_printed), ""))
  note <- add_note(note, ek$taxon_note)
  note <- add_note(note, ifelse(lindia, "printed 'W=0.030 b3', taken as a^3 (its standard weight is computed with the body length)", NA))
  note <- add_note(note, ifelse(use, NA, "equation of the body width: not used"))
  eqs$ejsmontkarabin <- data.frame(
    equation_id = ek$record_key, source_id = "ejsmontkarabin1998", study = "ejsmontkarabin1998",
    taxon = ek$taxon, rank = "species", condition = ek$formula_printed, scope = "pooled",
    length_type = "body_length", length_unit = "mm", mass_currency = "wet_mass", mass_unit = "mg",
    shell_included = FALSE, preservation = ifelse(ek$state == "C", "unknown", "none"),
    life_stage = "female", fit = "geometric formula (body volume)",
    ln_a = log(ek$coef) + log1p(app / 100), b = 3, var_ln_a = NA_real_, var_b = NA_real_,
    cov_ln_a_b = NA_real_, sigma_ln = NA_real_, rmse_ln = NA_real_, r2 = NA_real_, n = NA_real_,
    n_unit = NA_character_, length_min = ek$length_min_um / 1000, length_max = ek$length_max_um / 1000,
    use_in_model = use, note = note, stringsAsFactors = FALSE)
}

# McCauley (1984), Table 7.2 (book pp. 234-237): ln W = ln a + b ln L, W in ug
# dry mass, L in mm, compiled from nine studies (`study`). Those of Dumont et
# al. (1975), Bottrell et al. (1976) and Burns (1969) are already in the
# reference (the first in the data of 0.1.0, the others transcribed from
# Bottrell et al. 1976, where all 48 coincide in ln a and b) and are not used
# again. Lengths: code 1 (total length without the caudal setae;
# for Daphnia, head to the base of the tail spine) and code 3 (copepods
# without the furcal rami) are taken as body length; codes 2 (from the eye),
# 4 (thoracic segments) and 5 (Holopedium, setae natores to the terminal
# claw) are not used. Also left out: two ranges of Persson and Ekbohm (1980)
# that cannot be in mm (Holopedium 3.01-3.37, Bosmina obtusirostris
# 2.64-3.22), the equation of Mesocyclops leuckarti of Burgis (1974), which
# repeats that of Thermocyclops hyalinus, and the equation of Daphnia parvula
# of Pace and Orcutt (1981) on preserved animals (the same animals fresh are
# used).
if (use_source("mccauley1984", reg)) {
  mc <- zoo_csv("mccauley1984_table7_2.csv")
  stopifnot(nrow(mc) == 88)
  dup <- mc$original_source %in% c("Dumont et al. 1975", "Bottrell et al. 1976", "Burns 1969")
  len_ok <- mc$length_code %in% c(1, 3)
  bad_range <- mc$record_key %in% c("mccauley1984_t7.2_33", "mccauley1984_t7.2_44")
  repeat_eq <- mc$record_key == "mccauley1984_t7.2_74"
  preserved_twin <- mc$record_key == "mccauley1984_t7.2_02"
  stopifnot(mc$taxon[repeat_eq] == "Mesocyclops leuckarti", mc$taxon[preserved_twin] == "Daphnia parvula")
  slope_ok <- mc$b >= 1.5 & mc$b <= 4.5
  pooled <- pools_nauplii(mc$life_stage)
  use <- !dup & len_ok & !bad_range & !repeat_eq & !preserved_twin & slope_ok & !pooled
  note <- add_note(mc$note, paste0("length: ", mc$length_definition))
  note <- add_note(note, ifelse(dup, paste0("from ", mc$original_source, ", already in the reference: not used"), NA))
  note <- add_note(note, ifelse(!dup & !len_ok, "length not used by the reference: not used", NA))
  note <- add_note(note, ifelse(bad_range, "range cannot be in mm: not used", NA))
  note <- add_note(note, ifelse(repeat_eq, "repeats the equation of Thermocyclops hyalinus (Burgis 1974): not used", NA))
  note <- add_note(note, ifelse(preserved_twin, "preserved animals; the equation of the same animals fresh is used: not used", NA))
  note <- add_note(note, ifelse(!dup & !slope_ok, "slope outside 1.5-4.5: not used", NA))
  note <- add_note(note, ifelse(!dup & pooled, "pools nauplii with later stages: not used", NA))
  eqs$mccauley <- data.frame(
    equation_id = mc$record_key, source_id = "mccauley1984", study = mc$original_source,
    taxon = mc$taxon, rank = rank_of(mc$taxon), condition = mc$row_printed, scope = "population",
    length_type = "body_length", length_unit = "mm", mass_currency = "dry_mass", mass_unit = "mg",
    shell_included = FALSE, preservation = ifelse(mc$preserved %in% TRUE, "preserved", "none"),
    life_stage = mc$life_stage, fit = "OLS of ln W on ln L",
    ln_a = mc$ln_a - log(1000), b = mc$b, var_ln_a = NA_real_, var_b = NA_real_,
    cov_ln_a_b = NA_real_, sigma_ln = ifelse(is.na(mc$rms), NA_real_, sqrt(mc$rms)), rmse_ln = NA_real_,
    r2 = ifelse(is.na(mc$r), NA_real_, mc$r^2), n = mc$n, n_unit = "animals or weighings",
    length_min = mc$length_min_mm, length_max = mc$length_max_mm, use_in_model = use, note = note,
    stringsAsFactors = FALSE)
}

# Bottrell et al. (1976), Tables VII, X, XI and XII (pp. 432-441): ln W = ln
# a + b ln L, W in ug dry mass, L in mm; regressions on raw data of IBP
# participants (each one a `study`) and published ones. Not used: those of
# Dumont et al. (1975), already in the reference (Table VII, the rotifers, is
# theirs); the pooled regressions of several species, genera or families,
# which reuse the raw data of the regressions of the species; Table XI, the
# regressions by date of Daphnia longispina in Lake Shirakomanoike, whose
# pooled regression is in Table X; Holopedium (setae natores to the terminal
# claw) and Arctodiaptomus spinosus (thoracic segments only), lengths the
# reference does not use; those pooling nauplii with later stages; and
# slopes outside 1.5-4.5.
if (use_source("bottrell1976", reg)) {
  bt <- zoo_csv("bottrell1976_length_mass.csv")
  contributor <- sub("^[A-Za-z]+: ", "", sub(" \\(.*$", "", bt$data_source))
  dumont <- grepl("Dumont", bt$data_source)
  pooled_sp <- grepl("^Pooled|spp|probably", bt$data_source) | grepl("\\|", bt$taxon) |
    !grepl(" ", bt$taxon)
  by_date <- bt$table == "XI"
  other_len <- bt$taxon %in% c("Holopedium gibberum", "Arctodiaptomus spinosus")
  slope_ok <- bt$b >= 1.5 & bt$b <= 4.5
  pooled <- pools_nauplii(bt$life_stage)
  use <- !dumont & !pooled_sp & !by_date & !other_len & slope_ok & !pooled
  note <- add_note(bt$note, paste0("source: ", bt$data_source))
  note <- add_note(note, ifelse(dumont, "from Dumont et al. (1975), already in the reference: not used", NA))
  note <- add_note(note, ifelse(!dumont & pooled_sp, "pooled regression of several taxa: not used", NA))
  note <- add_note(note, ifelse(!dumont & by_date, "regression of one date; the pooled one of Table X is used: not used", NA))
  note <- add_note(note, ifelse(other_len, "length not used by the reference: not used", NA))
  note <- add_note(note, ifelse(!dumont & !slope_ok, "slope outside 1.5-4.5: not used", NA))
  note <- add_note(note, ifelse(!dumont & !pooled_sp & pooled, "pools nauplii with later stages: not used", NA))
  eqs$bottrell <- data.frame(
    equation_id = bt$record_key, source_id = "bottrell1976",
    study = ifelse(grepl("^Burns", bt$data_source), "Burns 1969", paste("Bottrell et al. 1976:", contributor)),
    taxon = bt$taxon, rank = rank_of(bt$taxon), condition = bt$row_printed, scope = "population",
    length_type = "body_length", length_unit = "mm", mass_currency = "dry_mass", mass_unit = "mg",
    shell_included = FALSE, preservation = "unknown", life_stage = bt$life_stage,
    fit = "OLS of ln W on ln L", ln_a = bt$ln_a - log(1000), b = bt$b, var_ln_a = NA_real_,
    var_b = NA_real_, cov_ln_a_b = NA_real_, sigma_ln = ifelse(is.na(bt$rms), NA_real_, sqrt(bt$rms)),
    rmse_ln = NA_real_, r2 = NA_real_, n = bt$n, n_unit = "animals or weighings",
    length_min = bt$length_min_mm, length_max = bt$length_max_mm, use_in_model = use, note = note,
    stringsAsFactors = FALSE)
}

if (!is.null(eqs$santini)) {
  eqs$santini$use_in_model <- FALSE
}
pub <- do.call(rbind, eqs)
rownames(pub) <- NULL
pub$permission <- reg$permission[match(pub$source_id, reg$source_id)]
pub <- pub[, c("equation_id", "source_id", "permission", setdiff(names(pub), c("equation_id", "source_id", "permission")))]
write.csv(pub, file.path(out_dir, "lwr_published.csv"), row.names = FALSE, fileEncoding = "UTF-8")
say("`lwr_published.csv`: %d equations from %d sources.", nrow(pub), length(unique(pub$source_id)))
say("")

# Check 0: zooplankton equations used and left out, and the mass of each used
# equation at the middle of its range (geometric), by source: the median and
# the extremes, for a quick look at plausibility (crustaceans of 0.5-2 mm
# weigh about 1-50 ug of dry mass; rotifers of 0.1-0.3 mm about 0.01-1 ug of
# wet mass).
zoo <- pub[pub$source_id %in% c("culver1985", "ejsmontkarabin1998", "mccauley1984", "bottrell1976"), ]
if (nrow(zoo)) {
  say("## Zooplankton equations")
  say("")
  say("| Source | Equations | Used | Mass at the middle of the range (ug), median (min-max) |")
  say("|---|---|---|---|")
  for (id in unique(zoo$source_id)) {
    z <- zoo[zoo$source_id == id & zoo$use_in_model, ]
    mid <- sqrt(z$length_min * z$length_max)
    w <- exp(z$ln_a) * mid^z$b * 1000
    say("| %s | %d | %d | %s |", id, sum(zoo$source_id == id), nrow(z),
        if (any(is.finite(w))) sprintf("%.3g (%.3g-%.3g)", stats::median(w, na.rm = TRUE),
                                       min(w, na.rm = TRUE), max(w, na.rm = TRUE)) else "no ranges")
  }
  why <- sub("^.*; ([^;]*): not used$", "\\1", zoo$note[!zoo$use_in_model])
  say("")
  say("- Not used: %s.", paste(sprintf("%s %d", names(table(why)), as.integer(table(why))), collapse = "; "))
  say("")
}

# Check 1: mass of a 50 mm frog or salamander from each equation of Santini et
# al. (2018).
say("## Mass at 50 mm snout-vent length")
say("")
say("| Equation | Taxon | Condition | Mass (g) |")
say("|---|---|---|---|")
for (i in which(pub$source_id == "santini2018")) {
  say("| %s | %s | %s | %.2f |", pub$equation_id[i], pub$taxon[i],
      ifelse(is.na(pub$condition[i]), "", pub$condition[i]), exp(pub$ln_a[i]) * 50^pub$b[i])
}
say("")

# Check 2: adult masses and snout-vent lengths of the typical records
# (species means of Santini et al. 2018, Table S1, and of Trochet et al. 2014)
# against the equation of their family, or of their order (by locomotion in
# frogs, by paedomorphosis in salamanders) when Table 3 has none for the
# family. Order, family, locomotion and paedomorphosis from Table S1; Trochet
# species not in Table S1 are left out. The Santini records are the data the
# equations were fitted to (a check of the transcription); the Trochet
# records are an independent compilation for European species.
if (use_source("santini2018", reg)) {
  tr_file <- file.path(out_dir, "typical_records.csv")
  if (!file.exists(tr_file)) stop("Run data-raw/typical_records.R first: ", tr_file, " not found.")
  tr <- read.csv(tr_file, stringsAsFactors = FALSE, encoding = "UTF-8")
  tr <- tr[tr$source_id %in% c("santini2018", "trochet2014") & tr$use_in_model &
             tr$life_stage %in% "adult" & tr$length_type %in% "snout_vent_length" &
             !is.na(tr$mean) & !is.na(tr$length_mean) & tr$unit == "g" & tr$length_unit == "mm", ]
  s1 <- suppressMessages(as.data.frame(readxl::read_excel(source_files("santini2018", reg)[1],
                                                          sheet = "Data")))
  first <- function(v) { v <- v[!is.na(v)]; if (length(v)) names(sort(table(v), decreasing = TRUE))[1] else NA }
  sp <- do.call(rbind, lapply(split(s1, s1$Species), function(d) {
    data.frame(species = d$Species[1], order = first(d$Order), family = first(d$Family),
               locomotion = first(d$Locomotion), paedomorphic = first(as.character(d$Pedomorphic)),
               stringsAsFactors = FALSE)
  }))
  k <- match(tr$taxon, sp$species)
  tr <- cbind(tr, sp[k, c("order", "family", "locomotion", "paedomorphic")])
  missing_sp <- unique(tr$taxon[is.na(k)])
  tr <- tr[!is.na(k), ]
  anura <- c(Arboreal = "arboreal", "Semi-Arboreal" = "arboreal",
             Aquatic = "semi-aquatic and aquatic", "Semi-Aquatic" = "semi-aquatic and aquatic",
             Terrestrial = "terrestrial and fossorial")
  order_eq <- ifelse(tr$order == "Anura",
                     pub$equation_id[match(paste("Anura", anura[tr$locomotion]),
                                           paste(pub$taxon, pub$condition))],
                     ifelse(tr$order == "Urodela",
                            pub$equation_id[match(paste("Caudata", ifelse(tr$paedomorphic == "TRUE", "paedomorphic", "non-paedomorphic")),
                                                  paste(pub$taxon, pub$condition))],
                            NA))
  fam <- pub[pub$rank == "family", ]
  fam_eq <- fam$equation_id[match(tr$family, fam$taxon)]
  axo <- tr$taxon == "Ambystoma mexicanum"
  fam_eq[tr$family == "Ambystomatidae"] <- ifelse(axo[tr$family == "Ambystomatidae"], "santini2018_t3_11", "santini2018_t3_10")
  pred <- function(eq) {
    i <- match(eq, pub$equation_id)
    exp(pub$ln_a[i]) * tr$length_mean^pub$b[i]
  }
  tr$ratio_order <- tr$mean / pred(order_eq)
  tr$eq_best <- ifelse(is.na(fam_eq), order_eq, fam_eq)
  tr$ratio_best <- tr$mean / pred(tr$eq_best)
  q <- function(r) {
    r <- r[is.finite(r)]
    sprintf("%d records, median %.2f, 95%% of records %.2f-%.2f", length(r), median(r),
            quantile(r, 0.025), quantile(r, 0.975))
  }
  say("## Observed adult mass over mass predicted from snout-vent length")
  say("")
  say("- Typical records of %d species; %d Trochet species not in Table S1 left out%s.",
      length(unique(tr$taxon)), length(missing_sp),
      if (length(missing_sp)) paste0(" (", paste0("*", missing_sp, "*", collapse = ", "), ")") else "")
  for (s in c("santini2018", "trochet2014")) {
    k <- tr$source_id == s
    say("- %s, equation of the family where available: %s.", s, q(tr$ratio_best[k]))
    say("- %s, equation of the order: %s.", s, q(tr$ratio_order[k]))
  }
  say("")
  say("By equation (family where available, records of both sources):")
  say("")
  say("| Equation | Taxon | Records | Median ratio | 2.5%% | 97.5%% |")
  say("|---|---|---|---|---|---|")
  for (e in sort(unique(tr$eq_best))) {
    r <- tr$ratio_best[tr$eq_best == e & is.finite(tr$ratio_best)]
    i <- match(e, pub$equation_id)
    say("| %s | %s%s | %d | %.2f | %.2f | %.2f |", e, pub$taxon[i],
        ifelse(is.na(pub$condition[i]), "", paste0(", ", pub$condition[i])), length(r),
        median(r), quantile(r, 0.025), quantile(r, 0.975))
  }
  out_far <- tr[is.finite(tr$ratio_best) & (tr$ratio_best > 3 | tr$ratio_best < 1 / 3), ]
  say("")
  say("Records more than three times above or below the prediction: %d%s.", nrow(out_far),
      if (nrow(out_far)) paste0(" (", paste(sprintf("*%s* %.1f g at %.1f mm, ratio %.2f, %s",
                                                    out_far$taxon, out_far$mean, out_far$length_mean,
                                                    out_far$ratio_best, out_far$source_id),
                                            collapse = "; "), ")") else "")
}
writeLines(log_lines, file.path(out_dir, "lwr_published_log.md"), useBytes = TRUE)
cat(log_lines, sep = "\n")
