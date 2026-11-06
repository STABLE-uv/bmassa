# Typical masses (mean mass of the individuals of a taxon) from the sources
# added in version 0.2.0
#
# Builds, in data-raw/derived/:
#   typical_records.csv      one row per record (taxon x source x stratum)
#   typical_records_log.md   counts, decisions and checks
#
# Run from the package root:  Rscript data-raw/typical_records.R
# Originals: folder BMASSA_DATABASES (see data-raw/sources.R); MD5 checked.
# Package used only here: readxl. Run in a UTF-8 locale (LANG=C.UTF-8).
#
# Columns: `mean` and the other statistics are masses of one counting unit
# (individual, cell or natural unit) in `unit`; `n` counts `n_unit`
# (individuals, samples or unknown). `sampling`: the sampling of the record
# when a site was sampled more than once (NEON: year, season and gear).
# `use_in_model = FALSE` marks records that would count the same animals twice
# or are not individuals (see `note`).

source("data-raw/sources.R")

out_dir <- "data-raw/derived"
dir.create(out_dir, showWarnings = FALSE)
reg <- read_sources()
log_lines <- c("# Typical masses from the sources added in version 0.2.0", "",
               paste0("Built by `data-raw/typical_records.R` on ", format(Sys.Date()), " (",
                      R.version.string, ")."), "")
say <- function(...) log_lines <<- c(log_lines, sprintf(...))
num <- function(x) suppressWarnings(as.numeric(x))

cols <- c("source_id", "taxon", "name_as_given", "identified_to", "life_stage", "count_unit",
          "currency", "unit", "shell_included", "preservation", "mean", "median", "se", "sd",
          "min", "max", "n", "n_unit", "scope", "site", "sampling", "country", "length_type",
          "length_mean", "length_unit", "original_reference", "use_in_model", "note")

# A data frame with all the columns above; missing ones filled with NA.
std <- function(...) {
  x <- data.frame(..., stringsAsFactors = FALSE)
  for (k in setdiff(cols, names(x))) x[[k]] <- NA
  x[, cols]
}

records <- list()

# 1. Hrycik et al. (2024): Great Lakes benthic invertebrates -------------------
# Mean individual wet mass (g, converted to mg) of formalin-preserved animals,
# molluscs with shell; mean, SE, median, minimum and maximum over samples
# (SampleSize_Nsamples = number of samples where the taxon was present).
# "All" in Lake, DepthZone_m and Basin = pooled; otherwise one stratum.

if (use_source("hrycik2024", reg)) {
  say("## Hrycik et al. (2024)")
  say("")
  path <- source_files("hrycik2024", reg)
  a <- read.csv(path["Hrycik2024/IndividualWeights_AllData.csv"], stringsAsFactors = FALSE)
  g <- read.csv(path["Hrycik2024/IndividualWeights_MajorGroups.csv"], stringsAsFactors = FALSE)
  s <- read.csv(path["Hrycik2024/SpeciesList.csv"], stringsAsFactors = FALSE)

  # Names: "X sp."/"X spp." -> genus X; "X y group" -> genus X (species
  # group); "A/B" -> lowest common taxon; "Unid. X" -> X.
  fix <- c("Musculium/Sphaerium" = "Sphaeriidae", "Nais communis/variablis" = "Nais",
           "Potamothrix bedoti/bavaricus" = "Potamothrix",
           "Valvata sincera/piscinalis" = "Valvata",
           "immature tubificid with hairs" = "Tubificinae",
           "immature tubificid without hairs" = "Tubificinae",
           "Unid. Orthocladiinae" = "Orthocladiinae", "Unid. Chironomidae" = "Chironomidae",
           "Oligochaeta Fragments" = "Oligochaeta")
  clean_name <- function(n) {
    t <- ifelse(n %in% names(fix), unname(fix[n]), n)
    t <- sub(" (sp|spp)\\.$", "", t)
    t <- sub("^([A-Z][a-z]+) [a-z]+ group$", "\\1", t)
    t
  }
  identified <- function(n) {
    ifelse(grepl(" (sp|spp)\\.$| group$|/", n) | n %in% c("Nais communis/variablis"), "genus_or_group",
           ifelse(grepl("^[A-Z][a-z]+ [a-z]+( [a-z]+)?$", n), "species", "higher"))
  }
  grp <- s$Group[match(a$Name, s$Name)]
  grp[a$Name == "Lumbriculidae"] <- "Oligochaeta"
  shell <- grp %in% c("Dreissena", "Gastropoda", "Sphaeriidae")
  pooled <- a$Lake == "All" & a$DepthZone_m == "All" & a$Basin == "All"
  fragments <- a$Name == "Oligochaeta Fragments"
  records$hrycik <- std(
    source_id = "hrycik2024", taxon = clean_name(a$Name), name_as_given = a$Name,
    identified_to = identified(a$Name), life_stage = NA_character_, count_unit = "individual",
    currency = "wet_mass", unit = "mg", shell_included = shell, preservation = "formalin",
    mean = 1000 * a$AverageIndividualWeight_g, median = 1000 * a$MedianIndividualWeight_g,
    se = 1000 * a$SE_weight, min = 1000 * a$MinimumIndividualWeight_g,
    max = 1000 * a$MaximumIndividualWeight_g, n = a$SampleSize_Nsamples, n_unit = "samples",
    scope = ifelse(pooled, "pooled", "population"),
    site = ifelse(pooled, NA_character_, paste0("Lake ", a$Lake, "; depth ", a$DepthZone_m,
                                                " m; basin ", a$Basin)),
    country = "Canada; United States", use_in_model = !fragments,
    note = ifelse(fragments, "Fragments, not individuals",
                  ifelse(a$Name %in% names(fix), paste0("Name in the source: ", a$Name), NA_character_)))
  say("- `IndividualWeights_AllData.csv`: %d records (%d pooled over lakes, depths and basins; %d for one stratum) of %d taxa; %d with shell (molluscs). 'Oligochaeta Fragments' kept but not used (fragments, not individuals).",
      nrow(a), sum(pooled), sum(!pooled), length(unique(a$Name)), sum(shell))
  # Check: median <= maximum and minimum <= mean <= maximum.
  bad <- with(a, MinimumIndividualWeight_g > AverageIndividualWeight_g |
                AverageIndividualWeight_g > MaximumIndividualWeight_g |
                MedianIndividualWeight_g > MaximumIndividualWeight_g)
  say("- Check: minimum <= mean <= maximum and median <= maximum in %d of %d records.",
      sum(!bad, na.rm = TRUE), nrow(a))

  gshell <- g$Name %in% c("Dreissena", "Gastropoda", "Sphaeriidae")
  gpooled <- g$Lake == "All" & g$DepthZone_m == "All" & g$Basin == "All"
  records$hrycik_groups <- std(
    source_id = "hrycik2024", taxon = ifelse(g$Name == "OtherInsecta", "Insecta", g$Name),
    name_as_given = g$Name, identified_to = "higher", life_stage = NA_character_,
    count_unit = "individual", currency = "wet_mass", unit = "mg", shell_included = gshell,
    preservation = "formalin", mean = 1000 * g$AverageIndividualWeight_g,
    median = 1000 * g$MedianIndividualWeight_g, se = 1000 * g$SE_weight,
    min = 1000 * g$MinimumIndividualWeight_g, max = 1000 * g$MaximumIndividualWeight_g,
    n = g$SampleSize_Nsamples, n_unit = "samples", scope = "group",
    site = ifelse(gpooled, NA_character_, paste0("Lake ", g$Lake, "; depth ", g$DepthZone_m,
                                                 " m; basin ", g$Basin)),
    country = "Canada; United States", use_in_model = FALSE,
    note = "Mean of a major group over the same samples as the taxon records; not used, to avoid counting animals twice")
  say("- `IndividualWeights_MajorGroups.csv`: %d records of %d major groups, kept for reference and not used in the model (same samples as the taxon records).",
      nrow(g), length(unique(g$Name)))
  say("")
}

# 2. Trochet et al. (2014): European amphibians --------------------------------
# One value per species from the literature: body mass (g) of males, females,
# adults and juveniles; snout-vent length (mm) of males, females and adults.
# "DD" = data deficient. Adult records are used; sexes and juveniles are kept
# with use_in_model = FALSE (the adult value often combines them).

if (use_source("trochet2014", reg)) {
  say("## Trochet et al. (2014)")
  say("")
  path <- source_files("trochet2014", reg)
  x <- as.data.frame(suppressMessages(readxl::read_excel(path[1], sheet = 1, col_names = FALSE,
                                                         col_types = "text",
                                                         .name_repair = "unique_quiet")))
  header <- as.character(unlist(x[4, ]))
  stopifnot(header[1] == "Species")
  d <- x[-(1:4), ]
  names(d) <- make.unique(ifelse(is.na(header), "blank", header))
  d <- d[!is.na(d$Species), ]
  col <- function(name) {
    stopifnot(name %in% names(d))
    num(d[[name]])
  }
  rows <- list(
    adult = list(mass = "Adult body mass", svl = "Adult snout-to-vent length", use = TRUE),
    male = list(mass = "Body mass in males", svl = "Snout-to-vent length in males", use = FALSE),
    female = list(mass = "Body mass in females", svl = "Snout-to-vent length in females", use = FALSE),
    juvenile = list(mass = "Body mass in juveniles", svl = NA, use = FALSE))
  tr <- list()
  for (st in names(rows)) {
    m <- col(rows[[st]]$mass)
    L <- if (is.na(rows[[st]]$svl)) rep(NA_real_, nrow(d)) else col(rows[[st]]$svl)
    ok <- !is.na(m)
    tr[[st]] <- std(source_id = "trochet2014", taxon = d$Species[ok], name_as_given = d$Species[ok],
                    identified_to = "species", life_stage = st, count_unit = "individual",
                    currency = "wet_mass", unit = "g", shell_included = FALSE,
                    preservation = "unknown", mean = m[ok], n_unit = "unknown", scope = "pooled",
                    country = "Europe", length_type = ifelse(is.na(L[ok]), NA, "snout_vent_length"),
                    length_mean = L[ok], length_unit = ifelse(is.na(L[ok]), NA, "mm"),
                    use_in_model = rows[[st]]$use,
                    note = if (rows[[st]]$use) NA_character_ else "Sex or juvenile value; the adult value is used")
  }
  records$trochet <- do.call(rbind, tr)
  say("- %d species; body mass of adults in %d, males in %d, females in %d and juveniles in %d (the rest 'DD', data deficient).",
      nrow(d), nrow(tr$adult), nrow(tr$male), nrow(tr$female), nrow(tr$juvenile))
  say("")
}

# 3. Santini et al. (2018): amphibian means ------------------------------------
# Supplementary Table S1: mass (g) and snout-vent length (mm) reported
# together for the same animals, from 241 references (kept in
# `original_reference`). Here only the means (Data_type "Mean", "Mean +
# range" or "Monthly Mean") with both values, sexes pooled (M_mean, SVL_mean)
# or by sex. Individual measurements go to data-raw/fit_lwr.R.

if (use_source("santini2018", reg)) {
  say("## Santini et al. (2018), Table S1")
  say("")
  path <- source_files("santini2018", reg)
  x <- as.data.frame(suppressMessages(readxl::read_excel(path[1], sheet = "Data")))
  means <- x$Data_type %in% c("Mean", "Mean + range", "Monthly Mean")
  # A species name repeated many times in one cell ("Polypedates
  # maculatusPolypedates maculatus...") is reduced to one.
  repeated <- grepl("^(.+?)\\1+$", x$Species, perl = TRUE)
  x$Species <- sub("^(.+?)\\1+$", "\\1", x$Species, perl = TRUE)
  note <- ifelse(is.na(x$Notes), NA_character_, paste0("Note in the source: ", x$Notes))
  note[repeated] <- paste0(ifelse(is.na(note[repeated]), "", paste0(note[repeated], "; ")),
                           "species name repeated in the cell of Table S1, reduced to one")
  # Rhinella jimi (Tolledo and Toledo 2015): mass and snout-vent length are
  # swapped in Table S1, in the individuals (see fit_lwr.R) and in their means
  # (mean "SVL" 177-195 mm with a maximum of 786 mm; mean "mass" 120-134 g).
  sw <- means & x$Species == "Rhinella jimi" & grepl("^Tolledo", x$Ref)
  for (s in c("", "_m", "_f")) {
    for (stat in c("mean", "min", "max", "midrange")) {
      a <- paste0("M", s, "_", stat)
      b <- paste0("SVL", s, "_", stat)
      tmp <- x[[a]][sw]
      x[[a]][sw] <- x[[b]][sw]
      x[[b]][sw] <- tmp
    }
  }
  note[sw] <- paste0(ifelse(is.na(note[sw]), "", paste0(note[sw], "; ")),
                     "mass and snout-vent length swapped back (swapped in Table S1)")
  parts <- list(
    adult = c("M_mean", "SVL_mean", "TRUE"), male = c("M_m_mean", "SVL_m_mean", "FALSE"),
    female = c("M_f_mean", "SVL_f_mean", "FALSE"))
  sa <- list()
  for (st in names(parts)) {
    m <- num(x[[parts[[st]][1]]])
    L <- num(x[[parts[[st]][2]]])
    nn <- num(x[[switch(st, adult = "N", male = "N_m", female = "N_f")]])
    ok <- means & !is.na(m) & !is.na(L) & m > 0 & L > 0
    sa[[st]] <- std(source_id = "santini2018", taxon = x$Species[ok], name_as_given = x$Species[ok],
                    identified_to = "species", life_stage = st, count_unit = "individual",
                    currency = "wet_mass", unit = "g", shell_included = FALSE,
                    preservation = "unknown", mean = m[ok], n = nn[ok],
                    n_unit = ifelse(is.na(nn[ok]), "unknown", "individuals"), scope = "population",
                    length_type = "snout_vent_length", length_mean = L[ok], length_unit = "mm",
                    original_reference = x$Ref[ok], use_in_model = as.logical(parts[[st]][3]),
                    note = note[ok])
  }
  records$santini <- do.call(rbind, sa)
  say("- Means with mass and snout-vent length of the same animals: %d with sexes pooled (%d species), %d of males and %d of females (kept, not used: the pooled means are used).",
      nrow(sa$adult), length(unique(sa$adult$taxon)), nrow(sa$male), nrow(sa$female))
  say("- *Rhinella jimi*: mass and snout-vent length swapped back in %d means of Tolledo and Toledo (2015), as in their individuals (data-raw/fit_lwr.R).",
      sum(sw))
  say("- Species names repeated in one cell of Table S1 reduced to one: %d rows (%s).",
      sum(repeated), paste(unique(x$Species[repeated]), collapse = ", "))
  say("- Life stage: adults (the paper collected adult measurements).")
  say("")
}

# 4. Rimet and Druart (2018): phytoplankton of temperate lakes -----------------
# One value per taxon: cell biovolume (um3) and, for colonies, the cumulated
# biovolume of the cells of a colony (um3), comparable to the natural-unit
# biovolumes of Kremer et al. (2014). "Notes on biovolumes" says where the
# dimensions come from (mostly literature averages).

if (use_source("rimet2018", reg)) {
  say("## Rimet and Druart (2018)")
  say("")
  path <- source_files("rimet2018", reg)
  x <- as.data.frame(suppressMessages(readxl::read_excel(path[1], sheet = 1)))
  nm <- function(pattern) {
    k <- grep(pattern, names(x), fixed = TRUE)
    stopifnot(length(k) == 1)
    names(x)[k]
  }
  name <- x[[nm("Genus + species name")]]
  level <- x[[nm("Determination uncertainty level")]]
  cellv <- num(x[[nm("Cell biovolume")]])
  ncell <- num(x[[nm("Number of cells per colony")]])
  colv <- num(x[[nm("Cumulated biovolume of cells in a colony")]])
  notes <- x[[nm("Notes on biovolumes")]]
  ident <- ifelse(!is.na(level) & tolower(level) %in% c("genus", "family", "order", "class"),
                  tolower(level), ifelse(grepl(" ", name), "species", "genus"))
  ok_c <- !is.na(cellv) & cellv > 0
  ok_n <- !is.na(colv) & colv > 0 & !is.na(ncell) & ncell > 1
  # Entries that are not cells or colonies of one taxon are kept but not
  # used: cysts and stomatocysts, empty loricas ("vide"), mixtures of taxa
  # ("+", or "/" before a genus), size classes of centric diatoms and informal
  # groups.
  why_not <- ifelse(grepl("(^|[^a-z])(cysts?|stomatocysts?)([^a-z]|$)", tolower(name)), "cyst",
             ifelse(grepl("(^|\\s)vide$", name), "empty lorica",
             ifelse(grepl("\\+|/ ?[A-Z]", name), "mixture of taxa",
             ifelse(grepl("^Centric", name), "size class",
             ifelse(name %in% c("micro Chlorophyceae", "Chlorophyceae sp. filaments",
                                "Volvocale 4 flagella"), "informal group", NA_character_)))))
  use_r <- is.na(why_not)
  note_r <- function(k, base) {
    out <- paste0(base, ifelse(is.na(why_not[k]), "", paste0(ifelse(nzchar(base), "; ", ""),
                                                             "not used: ", why_not[k])))
    out[!nzchar(out)] <- NA_character_
    out
  }
  records$rimet <- rbind(
    std(source_id = "rimet2018", taxon = name[ok_c], name_as_given = name[ok_c],
        identified_to = ident[ok_c], count_unit = "cell", currency = "biovolume", unit = "um3",
        shell_included = FALSE, preservation = "unknown", mean = cellv[ok_c], n_unit = "unknown",
        scope = "pooled", use_in_model = use_r[ok_c],
        note = note_r(which(ok_c), ifelse(is.na(notes[ok_c]), "", notes[ok_c]))),
    std(source_id = "rimet2018", taxon = name[ok_n], name_as_given = name[ok_n],
        identified_to = ident[ok_n], count_unit = "natural_unit", currency = "biovolume",
        unit = "um3", shell_included = FALSE, preservation = "unknown", mean = colv[ok_n],
        n_unit = "unknown", scope = "pooled", use_in_model = use_r[ok_n],
        note = note_r(which(ok_n), paste0("Cumulated biovolume of the ", ncell[ok_n], " cells of a colony",
                                          ifelse(is.na(notes[ok_n]), "", paste0("; ", notes[ok_n]))))))
  # Check: cumulated colony biovolume = cell biovolume x cells per colony.
  both <- ok_c & ok_n
  rel <- abs(colv[both] - cellv[both] * ncell[both]) / colv[both]
  say("- %d taxa; cell biovolume in %d; colonies (more than one cell) in %d.",
      nrow(x), sum(ok_c), sum(ok_n))
  say("- Check: cumulated colony biovolume = cell biovolume x cells per colony within 1%% in %d of %d colonial taxa.",
      sum(rel <= 0.01), sum(both))
  say("- Origin of the dimensions: %s.", paste(sprintf("%s %d", names(head(sort(table(notes), decreasing = TRUE), 2)),
                                                       as.integer(head(sort(table(notes), decreasing = TRUE), 2))),
                                               collapse = "; "))
  say("- Entries kept but not used (not cells or colonies of one taxon): %s.",
      paste(sprintf("%s %d", names(table(why_not)), as.integer(table(why_not))), collapse = ", "))
  say("")
}

# 5. Zooplankton transcriptions (phase 4k) ---------------------------------------
# Hand transcriptions in data-raw/, values as printed (checks in each log of
# the transcription, summarized in data-raw/README.md), with the phylum and
# class of each taxon written in (hints for the taxonomy). Weights in ug,
# stored in mg.

zoo_csv <- function(f) {
  read.csv(file.path("data-raw", f), stringsAsFactors = FALSE, encoding = "UTF-8",
           na.strings = c("", "NA"))
}
identified_as <- function(taxon) {
  ifelse(grepl("^[A-Z][a-z-]+ [a-z-]+", taxon), "species",
         ifelse(grepl("^[A-Z][a-z-]+$", taxon) & !grepl("(idae|ida|oda|ini|inae)$", taxon),
                "genus_or_group", "higher"))
}
gmean2 <- function(lo, hi) sqrt(lo * hi)

# 5a. Culver et al. (1985), Table 2: lengths and dry weights of the stages of
# the copepods of Lake Erie (4 % formalin), weighed in groups of 10-40 animals.
# One record per taxon and stage group: nauplii and copepodites, the mean of
# the weights of their instars (NI-NVI, CI-CV; equal numbers in each instar);
# adults (CVI), the geometric mean of the limits of the printed range, by
# sex. Eggs are left out.
if (use_source("culver1985", reg)) {
  say("## Culver et al. (1985), Table 2")
  say("")
  t2 <- zoo_csv("culver1985_table2.csv")
  t2 <- t2[t2$life_stage != "egg", ]
  t2$w <- ifelse(is.na(t2$dw_ug), gmean2(t2$dw_min_ug, t2$dw_max_ug), t2$dw_ug)
  t2$L <- ifelse(is.na(t2$length_mm), gmean2(t2$length_min_mm, t2$length_max_mm), t2$length_mm)
  stopifnot(all(is.finite(t2$w)), all(is.finite(t2$L)))
  g <- split(t2, paste(t2$taxon, t2$life_stage), drop = TRUE)
  cu <- do.call(rbind, lapply(g, function(d) {
    instars <- d$life_stage[1] %in% c("nauplius", "copepodite")
    data.frame(taxon = d$taxon[1], name_as_given = d$taxon_printed[1], life_stage = d$life_stage[1],
               mean = mean(d$w) / 1000, length_mean = mean(d$L), n = sum(d$n),
               note = paste0(if (instars) sprintf("mean of the instars %s (%s ug)",
                                                  paste(d$stage_printed, collapse = ", "),
                                                  paste(d$w, collapse = ", "))
                             else sprintf("geometric mean of the printed range %s-%s ug at %s-%s mm",
                                          d$dw_min_ug, d$dw_max_ug, d$length_min_mm, d$length_max_mm),
                             if (!is.na(d$taxon_note[1])) paste0("; ", d$taxon_note[1]) else ""),
               stringsAsFactors = FALSE)
  }))
  records$culver <- std(source_id = "culver1985", taxon = cu$taxon, name_as_given = cu$name_as_given,
                        identified_to = identified_as(cu$taxon), life_stage = cu$life_stage,
                        count_unit = "individual", currency = "dry_mass", unit = "mg",
                        shell_included = FALSE, preservation = "formalin", mean = cu$mean, n = cu$n,
                        n_unit = "individuals", scope = "population", site = "Lake Erie",
                        country = "Canada; United States", length_type = "body_length",
                        length_mean = cu$length_mean, length_unit = "mm", use_in_model = TRUE,
                        note = cu$note)
  say("- %d printed lines without the eggs, %d records (taxon and stage group): %s.",
      nrow(t2), nrow(cu), paste(sprintf("%s %d", names(table(cu$life_stage)),
                                        as.integer(table(cu$life_stage))), collapse = ", "))
  say("")
}

# 5b. Ejsmont-Karabin (1998), Table II: standard wet weight (ug) of the
# females of 251 species of planktonic rotifers, from the equation of the
# species at its standard length (lorica or body length, without spines,
# toes or foot as stated; a body or corona width for the species whose
# equation uses it), with the appendices (spines, toes) included and the
# mucus sheaths left out. C: contracted on preservation; L: living. One record
# per species; converted to dry mass in model_data.R. Asplanchnopus multiceps
# is left out: its printed weight (3.983 ug) is a tenth of its own equation
# at its standard length (39.8 ug), probably a misplaced decimal point.
if (use_source("ejsmontkarabin1998", reg)) {
  say("## Ejsmont-Karabin (1998), Table II")
  say("")
  ek <- zoo_csv("ejsmontkarabin1998_table2.csv")
  stopifnot(nrow(ek) == 251, all(is.finite(ek$std_weight_ug_ww)))
  width <- ek$std_length_type %in% c("corona", "max. b.", "max. c.")
  bad <- ek$record_key == "ejsmontkarabin1998_t2_009"
  stopifnot(sum(bad) == 1, ek$taxon[bad] == "Asplanchnopus multiceps")
  note <- paste0(ifelse(ek$state == "C", "contracted on preservation", "living"),
                 "; standard ", ifelse(width, "width", "length"), " ", ek$std_length_printed, " um",
                 ifelse(!is.na(ek$appendices_pct), sprintf("; appendices %s %% of the weight", ek$appendices_printed), ""),
                 ifelse(!is.na(ek$taxon_note), paste0("; ", ek$taxon_note), ""),
                 ifelse(bad, "; not used: printed weight a tenth of the equation at the standard length", ""))
  records$ejsmontkarabin <- std(
    source_id = "ejsmontkarabin1998", taxon = ek$taxon, name_as_given = ek$species_printed,
    identified_to = "species", life_stage = NA_character_, count_unit = "individual",
    currency = "wet_mass", unit = "mg", shell_included = FALSE,
    preservation = ifelse(ek$state == "C", "unknown", "none"),
    mean = ek$std_weight_ug_ww / 1000, n_unit = "unknown", scope = "pooled",
    length_type = ifelse(width, NA_character_, "body_length"),
    length_mean = ifelse(width, NA_real_, ek$std_length_um / 1000), length_unit = "mm",
    use_in_model = !bad, note = note)
  say("- %d species (%d with a standard width instead of a length); %d left out (Asplanchnopus multiceps: printed weight a tenth of its equation).",
      nrow(ek), sum(width), sum(bad))
  say("- Standard wet weights from %.3f to %.1f ug (median %.3f).", min(ek$std_weight_ug_ww),
      max(ek$std_weight_ug_ww), stats::median(ek$std_weight_ug_ww))
  say("")
}

# 5c. Hébert et al. (2016): mean individual dry mass (mg) of freshwater
# crustacean zooplankton (cladocerans and copepods) from the mean body length
# of each species in the freshwater zooplankton database of the CIEE working
# group (Barnett et al. 2007; North American lake surveys) and the length-mass
# equations of McCauley (1984) and Culver et al. (1985) (references "1,18,20"
# of the data set): the typical length of a species with its equation. Their
# other dry masses (sizes of Lynch 1980 with the equations of Bottrell et al.
# 1976, and animals of physiological experiments) are left out, as are the
# marine taxa and the hybrids. Names: "sp."/"spp." and pairs of species give
# the genus; subspecies written with a hyphen are dropped ("galeata-mendotae"
# is Daphnia mendotae, "bicuspidatus-thomasi" Diacyclops thomasi); misspelt
# genera are corrected (Ceriopdaphnia, Alglaodiaptomus) and groups written as
# genera ("Bosminids sp.") take the name of their family or order.
if (use_source("hebert2016", reg)) {
  say("## Hébert et al. (2016)")
  say("")
  hb <- read.csv2(source_files("hebert2016", reg)[1], stringsAsFactors = FALSE, na.strings = "NA",
                  fileEncoding = "latin1")
  n_all <- nrow(hb)
  hb <- hb[hb$Habitat == "Freshwater" & !is.na(hb$Dry.mass) & hb$Ref.dm %in% "1,18,20", ]
  hybrid <- hb$Species %in% c("cucullataxgaleata", "hyalinaxgaleata")
  # Misspelt genera and the groups written as genera.
  genus <- hb$Genus
  genus[genus == "Ceriopdaphnia"] <- "Ceriodaphnia"
  genus[genus == "Alglaodiaptomus"] <- "Aglaodiaptomus"
  group_name <- c(Bosminids = "Bosminidae", Diaptomiids = "Diaptomidae", Harpacticoids = "Harpacticoida")
  sp <- hb$Species
  sp <- sub("^galeata-mendotae$", "mendotae", sp)
  sp <- sub("^(longimanus)-\\(total\\)$", "\\1", sp)
  sp <- sub("^([a-z]+)-[a-z]+$", "\\1", sp)
  genus[hb$Genus == "Cyclops" & hb$Species == "bicuspidatus-thomasi"] <- "Diacyclops"
  sp[hb$Genus == "Cyclops" & hb$Species == "bicuspidatus-thomasi"] <- "thomasi"
  to_genus <- sp %in% c("sp.", "spp.") | grepl("/", sp)
  taxon <- ifelse(to_genus, genus, paste(genus, sp))
  grp_row <- genus %in% names(group_name)
  stopifnot(all(to_genus[grp_row]))
  taxon[grp_row] <- unname(group_name[genus[grp_row]])
  stopifnot(!anyNA(taxon), all(is.finite(hb$Dry.mass)), all(hb$Dry.mass > 0))
  records$hebert <- std(
    source_id = "hebert2016", taxon = taxon[!hybrid], name_as_given = paste(hb$Genus, hb$Species)[!hybrid],
    identified_to = ifelse(grp_row, "higher", ifelse(to_genus, "genus_or_group", "species"))[!hybrid],
    life_stage = NA_character_,
    count_unit = "individual", currency = "dry_mass", unit = "mg", shell_included = FALSE,
    preservation = "unknown", mean = hb$Dry.mass[!hybrid], n_unit = "unknown", scope = "pooled",
    country = "Canada; United States", length_type = "body_length",
    length_mean = hb$Body.length[!hybrid], length_unit = "mm",
    original_reference = "CIEE freshwater zooplankton database; McCauley (1984); Culver et al. (1985)",
    use_in_model = TRUE,
    note = paste0(hb$Group[!hybrid], "; mean body length ", hb$Body.length[!hybrid],
                  " mm with the equations of McCauley (1984) and Culver et al. (1985)"))
  say("- %d rows in the data set; %d freshwater rows with a dry mass from the CIEE lengths (references 1,18,20): %s; %d hybrids left out.",
      n_all, nrow(hb), paste(sprintf("%s %d", names(table(hb$Group)), as.integer(table(hb$Group))), collapse = ", "),
      sum(hybrid))
  say("- Dry masses from %.2f to %.1f ug (median %.2f).", 1000 * min(hb$Dry.mass), 1000 * max(hb$Dry.mass),
      1000 * stats::median(hb$Dry.mass))
  say("")
}

# 5d. Bottrell et al. (1976): individual dry weights of rotifers (Table VIII,
# p. 433; weighings compiled from several studies) and wet weights from
# calculated volumes (Table IX, p. 434; the range of the means of 1-9 IBP
# sites; the caption says mg, but only ug fit: a Keratella would weigh
# 50-620 ug, against dry weights of 0.04-0.11 ug in Table VIII). Left out: the weights from Dumont et al. (1975), already in the
# reference, and the Protozoa (outside the domains of the package). A list or
# a range of values gives one record, at the geometric mean of its limits
# (the range of Table IX spans the means of its sites).
if (use_source("bottrell1976", reg)) {
  say("## Bottrell et al. (1976), Tables VIII and IX")
  say("")
  bt <- zoo_csv("bottrell1976_typical.csv")
  n_all <- nrow(bt)
  # Rows of Dumont et al. (1975): cited, or from their sites without the
  # citation repeated (River Sambre).
  dumont <- grepl("Dumont|Sambre", bt$locality_source)
  bt <- bt[bt$group == "Rotifera" & !dumont, ]
  stopifnot(all(bt$currency %in% c("dry_mass", "wet_mass")), all(is.finite(bt$weight_min_ug)),
            all(is.finite(bt$weight_max_ug)))
  w <- gmean2(bt$weight_min_ug, bt$weight_max_ug)
  spread <- bt$weight_max_ug > bt$weight_min_ug
  records$bottrell <- std(
    source_id = "bottrell1976", taxon = bt$taxon, name_as_given = bt$taxon_printed,
    identified_to = identified_as(bt$taxon), life_stage = bt$life_stage, count_unit = "individual",
    currency = bt$currency, unit = "mg", shell_included = FALSE, preservation = "unknown",
    mean = w / 1000, min = ifelse(spread, bt$weight_min_ug / 1000, NA), max = ifelse(spread, bt$weight_max_ug / 1000, NA),
    n = suppressWarnings(as.numeric(bt$sites)), n_unit = ifelse(bt$table == "IX", "sites", "unknown"),
    scope = "pooled", site = ifelse(bt$table == "VIII", bt$locality_source, NA),
    use_in_model = TRUE,
    note = paste0("Table ", bt$table, ": ", bt$weight_printed, " ug",
                  ifelse(bt$table == "IX", " (wet weight from calculated volumes)", " (dry weight)"),
                  ifelse(spread, "; geometric mean of the limits", "")))
  say("- %d rows; %d of rotifers not from Dumont et al. (1975): dry weights %d (Table VIII), wet weights %d (Table IX).",
      n_all, nrow(bt), sum(bt$currency == "dry_mass"), sum(bt$currency == "wet_mass"))
  say("")
}

# 5e. Fish without lengths (decision of 2026-11-05): mean individual wet mass
# (g) of each population with at least `min_fish` individuals weighed. Arranz
# et al. (2022): species and site of the electrofishing survey of Iberian
# streams (all sizes caught; masses to 0.1 g). Koning et al. (2024): species
# in the catch of the 'dai' bag-net fishery of the Tonle Sap River in
# 2020-2021, one population each (the file with total lengths; the one with
# standard lengths has the same fish). The mean of a population depends on
# the sizes caught, which the model takes as variation between populations
# and sources.
min_fish <- 5
fish_means <- function(taxon, site, mass, length) {
  ok <- is.finite(mass) & mass > 0
  key <- paste(taxon, site, sep = "\r")[ok]
  g <- split(which(ok), factor(key, levels = unique(key)))
  n <- lengths(g)
  g <- g[n >= min_fish]
  first <- vapply(g, `[`, integer(1), 1)
  data.frame(taxon = taxon[first], site = site[first], n = lengths(g),
             mean = vapply(g, function(i) mean(mass[i]), numeric(1)),
             median = vapply(g, function(i) stats::median(mass[i]), numeric(1)),
             sd = vapply(g, function(i) stats::sd(mass[i]), numeric(1)),
             min = vapply(g, function(i) min(mass[i]), numeric(1)),
             max = vapply(g, function(i) max(mass[i]), numeric(1)),
             length_mean = vapply(g, function(i) mean(length[i], na.rm = TRUE), numeric(1)),
             stringsAsFactors = FALSE, row.names = NULL)
}
if (use_source("arranz2022", reg)) {
  say("## Arranz et al. (2022): mean mass of fish populations")
  say("")
  path <- source_files("arranz2022", reg)
  x <- read.delim(path[1], fileEncoding = "latin1", stringsAsFactors = FALSE)
  m <- fish_means(x$Scientific_name, as.character(x$Code_ID), x$Mass_g, x$Length_mm)
  tl_sp <- c("Anguilla anguilla", "Gambusia holbrooki", "Salaria fluviatilis", "Silurus glanis",
             "Ameiurus melas")
  records$arranz <- std(
    source_id = "arranz2022", taxon = m$taxon, name_as_given = m$taxon,
    identified_to = "species", life_stage = NA_character_, count_unit = "individual",
    currency = "wet_mass", unit = "g", shell_included = FALSE, preservation = "fresh",
    mean = m$mean, median = m$median, sd = m$sd, min = m$min, max = m$max, n = m$n,
    n_unit = "individuals", scope = "population", site = m$site, country = "Spain",
    length_type = ifelse(m$taxon %in% tl_sp, "total_length", "fork_length"),
    length_mean = m$length_mean, length_unit = "mm", use_in_model = TRUE,
    note = "mean of the individuals of the species at the site (electrofishing)")
  say("- %d populations (species and site) with at least %d individuals, of %d species at %d sites (%d individuals); population means from %.2f to %.0f g (median %.1f).",
      nrow(m), min_fish, length(unique(m$taxon)), length(unique(m$site)), sum(m$n),
      min(m$mean), max(m$mean), stats::median(m$mean))
  say("")
}
if (use_source("koning2024", reg)) {
  say("## Koning et al. (2024): mean mass of fish in the catch")
  say("")
  path <- source_files("koning2024", reg)
  tl <- read.csv(path["Koning2024/Dai_2020_2021_TL.csv"], fileEncoding = "UTF-8-BOM",
                 stringsAsFactors = FALSE)
  sp <- gsub("_", " ", tl$Genus_species)
  m <- fish_means(sp, rep("Tonle Sap River, dai fishery", length(sp)), tl$W.g, 10 * tl$tl.cm)
  records$koning <- std(
    source_id = "koning2024", taxon = m$taxon, name_as_given = m$taxon,
    identified_to = "species", life_stage = NA_character_, count_unit = "individual",
    currency = "wet_mass", unit = "g", shell_included = FALSE, preservation = "fresh",
    mean = m$mean, median = m$median, sd = m$sd, min = m$min, max = m$max, n = m$n,
    n_unit = "individuals", scope = "population", site = m$site, country = "Cambodia",
    length_type = "total_length", length_mean = m$length_mean, length_unit = "mm",
    use_in_model = TRUE, note = "mean of the individuals of the species in the catch of 2020-2021")
  say("- %d species with at least %d individuals (%d individuals); means from %.1f to %.0f g (median %.1f).",
      nrow(m), min_fish, sum(m$n), min(m$mean), max(m$mean), stats::median(m$mean))
  say("")
}
# NEON: species in each sampling (bout: site, year and season) and gear of the
# streams and lakes of 2016-2024 (RELEASE-2026), read with data-raw/neon.R.
# The fish measured are up to 50 (30 at some sites) per species and pass,
# the rest being counted only, so the mean of those measured is that of the
# catch. Fish recorded at the minimum of the scale take the mean mass of the
# fish of their length lighter than it (neon.R); populations with such a fish
# without an equation to give it a mass are left out. A site is sampled about
# twice a year: its samplings of a species are one population
# (data-raw/model_data.R).
if (use_source("neon_fish", reg)) {
  say("## NEON: mean mass of fish populations in each sampling")
  say("")
  source("data-raw/neon.R")
  path <- source_files("neon_fish", reg)
  nf <- read_neon_fish(path[1])
  x <- nf$fish
  key <- paste(x$site, x$bout, x$gear, sep = "|")
  grp <- paste(x$taxon, key, sep = "\r")
  n_grp <- table(grp)
  no_mass <- unique(grp[is.na(x$mass)])
  big <- names(n_grp)[n_grp >= min_fish]
  keep <- !grp %in% no_mass
  m <- fish_means(x$taxon[keep], key[keep], x$mass[keep], x$length[keep])
  part <- do.call(rbind, strsplit(m$site, "|", fixed = TRUE))
  records$neon <- std(
    source_id = "neon_fish", taxon = m$taxon, name_as_given = m$taxon,
    identified_to = "species", life_stage = NA_character_, count_unit = "individual",
    currency = "wet_mass", unit = "g", shell_included = FALSE, preservation = "fresh",
    mean = m$mean, median = m$median, sd = m$sd, min = m$min, max = m$max, n = m$n,
    n_unit = "individuals", scope = "population", site = part[, 1],
    sampling = paste0(part[, 2], ", ", part[, 3]), country = "United States",
    length_type = "total_length", length_mean = m$length_mean, length_unit = "mm",
    use_in_model = TRUE, note = "mean of the individuals of the species measured in one sampling with one gear")
  imp <- tapply(x$imputed[keep], grp[keep], mean)[paste(m$taxon, m$site, sep = "\r")]
  say("- %d fish kept by data-raw/neon.R (%d with the mass of their length below the minimum of the scale: %d from the equation of their species, %d from that of their genus; %d recorded at the minimum without an equation).",
      nrow(x), sum(x$imputed), nf$counts[["imputed_species"]], nf$counts[["imputed_genus"]],
      sum(is.na(x$mass)))
  say("- %d populations (species, sampling and gear) with at least %d individuals, of %d species at %d sites (%d individuals); %d more left out for holding a fish without mass (%s). Gears: %s.",
      nrow(m), min_fish, length(unique(m$taxon)), length(unique(part[, 1])), sum(m$n),
      sum(big %in% no_mass),
      paste(names(sort(table(sub("\r.*$", "", intersect(big, no_mass))), decreasing = TRUE)),
            collapse = ", "),
      paste(sprintf("%s %d", names(table(part[, 3])), as.integer(table(part[, 3]))), collapse = ", "))
  say("- Population means from %.2f to %.0f g (median %.1f); the masses given to the fish below the minimum are more than half of the fish in %d populations.",
      min(m$mean), max(m$mean), stats::median(m$mean), sum(imp > 0.5))
  say("")
}

# Output -------------------------------------------------------------------------

rec <- do.call(rbind, records)
rownames(rec) <- NULL
rec$permission <- reg$permission[match(rec$source_id, reg$source_id)]
rec <- rec[order(rec$source_id, rec$taxon, rec$life_stage, rec$count_unit, rec$scope,
                 rec$site, rec$sampling, method = "radix"), ]
rec$record_id <- sprintf("%s_typ_%04d", rec$source_id,
                         ave(seq_len(nrow(rec)), rec$source_id, FUN = seq_along))
rec <- rec[, c("record_id", "source_id", "permission", cols[-1])]
write.csv(rec, file.path(out_dir, "typical_records.csv"), row.names = FALSE, fileEncoding = "UTF-8")

say("## Output")
say("")
say("`typical_records.csv`: %d records (%d used in the model) of %d taxa from %d sources.",
    nrow(rec), sum(rec$use_in_model), length(unique(rec$taxon)), length(unique(rec$source_id)))
tab <- table(rec$source_id, paste(rec$currency, rec$unit))
say("")
say("| Source | %s |", paste(colnames(tab), collapse = " | "))
say("|---|%s", paste(rep("---|", ncol(tab)), collapse = ""))
for (i in seq_len(nrow(tab))) say("| %s | %s |", rownames(tab)[i], paste(tab[i, ], collapse = " | "))
writeLines(log_lines, file.path(out_dir, "typical_records_log.md"), useBytes = TRUE)
cat(log_lines, sep = "\n")
