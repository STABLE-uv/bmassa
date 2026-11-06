# Length-mass relationships fitted from individual measurements
#
# Builds, in data-raw/derived/:
#   lwr_fitted.csv          one row per fitted equation ln(M) = ln(a) + b ln(L)
#   conversions_cadses.csv  dry mass -> ash-free dry mass by taxon (CADSES)
#   fit_lwr_log.md          counts, exclusions and checks against published values
#
# Run from the package root:  Rscript data-raw/fit_lwr.R
# Originals: folder BMASSA_DATABASES (see data-raw/sources.R); MD5 checked.
# Sources whose permission is pending or unverified are included unless
# BMASSA_INCLUDE_PENDING = "false" (data-raw/sources.R, use_source()).
# Package used only here: readxl. Run in a UTF-8 locale (LANG=C.UTF-8).

source("data-raw/sources.R")

# Settings ---------------------------------------------------------------------
# Provisional choices:
min_n <- 10             # individuals per equation, after exclusions
min_range_ratio <- 1.5  # largest / smallest length within an equation
max_rstudent <- 4       # gross outliers: |studentized residual| > 4, one pass
# Flags (equations are kept; the model step decides how to use them):
flag_r2 <- 0.5          # r2 below this: "low_r2"
flag_b <- c(1.5, 4.5)   # b outside this range: "b_outside_range"

out_dir <- "data-raw/derived"
dir.create(out_dir, showWarnings = FALSE)
reg <- read_sources()
log_lines <- c("# Length-mass relationships fitted from individual data", "",
               paste0("Built by `data-raw/fit_lwr.R` on ", format(Sys.Date()), " (",
                      R.version.string, ")."), "",
               sprintf("Settings: at least %d individuals and a length range of at least %.1f-fold per equation; individuals with |studentized residual| > %g removed (one pass).",
                       min_n, min_range_ratio, max_rstudent), "")
say <- function(...) log_lines <<- c(log_lines, sprintf(...))

num <- function(x) suppressWarnings(as.numeric(x))

# Fitting ----------------------------------------------------------------------

# Ordinary least squares of ln(M) on ln(L). Returns the coefficients, their
# covariance, the residual standard deviation (ln scale) and the length range,
# or the reason why no equation was fitted.
fit_one <- function(L, M) {
  keep <- is.finite(L) & is.finite(M) & L > 0 & M > 0
  x <- log(L[keep])
  y <- log(M[keep])
  check <- function(x) {
    if (length(x) < min_n) return("fewer individuals than min_n")
    if (exp(max(x) - min(x)) < min_range_ratio) return("length range narrower than min_range_ratio")
    NA_character_
  }
  why <- check(x)
  if (!is.na(why)) return(list(reason = why))
  f <- lm(y ~ x)
  rs <- rstudent(f)
  drop <- is.finite(rs) & abs(rs) > max_rstudent
  if (any(drop)) {
    x <- x[!drop]
    y <- y[!drop]
    why <- check(x)
    if (!is.na(why)) return(list(reason = why))
    f <- lm(y ~ x)
  }
  v <- vcov(f)
  s <- summary(f)
  list(reason = NA_character_, n = length(x), n_excluded = sum(drop),
       ln_a = unname(coef(f)[1]), b = unname(coef(f)[2]),
       var_ln_a = v[1, 1], var_b = v[2, 2], cov_ln_a_b = v[1, 2],
       sigma_ln = s$sigma, r2 = s$r.squared,
       length_min = exp(min(x)), length_max = exp(max(x)))
}

group_cols <- c("source_id", "taxon", "length_type", "length_unit", "mass_currency",
                "mass_unit", "shell_included", "preservation", "life_stage")

# Fits every group of a standardized data set `d` (columns of group_cols plus
# site, country, L and M). scope "pooled": all individuals of a taxon;
# "population": one site, only for taxa found at more than one site.
fit_dataset <- function(d, population = TRUE) {
  d$shell_included <- as.character(d$shell_included)
  stats_cols <- c("n", "n_excluded", "ln_a", "b", "var_ln_a", "var_b", "cov_ln_a_b",
                  "sigma_ln", "r2", "length_min", "length_max")
  one_scope <- function(d, by, scope) {
    key <- do.call(paste, c(d[by], sep = "\r"))
    idx <- split(seq_len(nrow(d)), key)
    rows <- lapply(idx, function(i) {
      r <- fit_one(d$L[i], d$M[i])
      st <- setNames(as.list(rep(NA_real_, length(stats_cols))), stats_cols)
      st[intersect(names(r), stats_cols)] <- r[intersect(names(r), stats_cols)]
      first <- d[i[1], group_cols, drop = FALSE]
      sites <- sort(unique(d$site[i]))
      countries <- sort(unique(d$country[i]))
      data.frame(first, scope = scope,
                 site = if (scope == "population") sites[1] else NA_character_,
                 n_sites = length(sites), country = paste(countries, collapse = "; "),
                 n_input = length(i), as.data.frame(st), reason = r$reason,
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }
  res <- one_scope(d, group_cols, "pooled")
  if (population) {
    n_site <- tapply(d$site, do.call(paste, c(d[group_cols], sep = "\r")),
                     function(s) length(unique(s)))
    multi <- names(n_site)[n_site > 1]
    dp <- d[do.call(paste, c(d[group_cols], sep = "\r")) %in% multi, ]
    if (nrow(dp)) res <- rbind(res, one_scope(dp, c(group_cols, "site"), "population"))
  }
  res
}

report_fits <- function(f, label) {
  ok <- is.na(f$reason)
  say("- %s: %d equations fitted (%d pooled, %d by site) for %d taxa; %d individuals excluded as gross outliers (pooled fits); groups not fitted: %s.",
      label, sum(ok), sum(ok & f$scope == "pooled"), sum(ok & f$scope == "population"),
      length(unique(f$taxon[ok])), sum(f$n_excluded[ok & f$scope == "pooled"], na.rm = TRUE),
      if (any(!ok)) paste(sprintf("%d %s", as.integer(table(f$reason[!ok])),
                                  names(table(f$reason[!ok]))), collapse = ", ") else "none")
}

all_fits <- list()

# 1. CADSES lagoons (Rosati and Basset 2020, 2021) -----------------------------
# Individual length (mm) and dry mass (mg; units assumed, see sources.csv).
# Length: total body length (worms), head to last abdominal segment (insects,
# crustaceans, polychaetes), shell height (gastropods), valve length
# (bivalves) (CASCADE D3.1.1, 2021). Dry mass of molluscs and barnacles
# includes the shell or plates.

read_cadses <- function(id) {
  path <- source_files(id, reg)
  x <- read.table(path[1], sep = ";", header = TRUE, quote = "", comment.char = "",
                  na.strings = "NA", fileEncoding = "UTF-8-BOM", colClasses = "character")
  x$source_id <- id
  x
}

cad_ids <- c("cadses_boxcorer", "cadses_leafpack")
cad_ids <- cad_ids[use_source(cad_ids, reg)]
cad <- lapply(cad_ids, read_cadses)
names(cad) <- cad_ids

if (length(cad)) {
  say("")
  say("## CADSES lagoons")
  say("")
  # Class of molluscs: from the box-corer classification (by taxon, then by
  # order), completed by hand for taxa without classification.
  all_cad <- do.call(rbind, lapply(cad, function(x) {
    data.frame(source_id = x$source_id, taxon = x$scientificName, phylum = x$phylum,
               class = if ("class" %in% names(x)) x$class else NA_character_,
               order = x$order, stringsAsFactors = FALSE)
  }))
  known <- all_cad[!is.na(all_cad$class) & all_cad$phylum %in% "Mollusca", ]
  class_by_taxon <- tapply(known$class, known$taxon, function(v) names(which.max(table(v))))
  class_by_order <- tapply(known$class, known$order, function(v) names(which.max(table(v))))
  class_by_order <- c(class_by_order, Mesogastropoda = "Gastropoda",
                      Sacoglossa = "Gastropoda", Nudibranchia = "Gastropoda")
  class_by_taxon <- c(class_by_taxon, Tellina = "Bivalvia", Mya = "Bivalvia",
                      Bivalvia = "Bivalvia", Ancylinae = "Gastropoda")
  mollusc_class <- function(taxon, order, phylum) {
    cl <- unname(class_by_taxon[taxon])
    cl[is.na(cl)] <- unname(class_by_order[order[is.na(cl)]])
    cl[!(phylum %in% "Mollusca") & !(taxon %in% names(class_by_taxon))] <- NA_character_
    cl
  }
  barnacle <- function(taxon, order) order %in% c("Sessilia", "Thoracica") | grepl("balanus", taxon, ignore.case = TRUE)

  conv_rows <- list()
  for (id in names(cad)) {
    x <- cad[[id]]
    L <- num(x$totalLength)
    bad_L <- sum(!is.na(x$totalLength) & is.na(L))
    DW <- num(x$dryWeight)
    AFDW <- num(x$ashFreeDryWeight)
    incons <- !is.na(AFDW) & !is.na(DW) & AFDW > DW
    cl <- mollusc_class(x$scientificName, x$order, x$phylum)
    hard <- !is.na(cl) | barnacle(x$scientificName, x$order)
    phylum_known <- !is.na(x$phylum) | !is.na(cl)
    shell <- ifelse(hard, TRUE, ifelse(phylum_known, FALSE, NA))
    # Insect and other taxa without phylum are not molluscs (checked by hand):
    no_shell <- c("Chironomidae", "Chironomus plumosus", "Diamesinae", "Nereis",
                  "Mysta picta", "Ceratopogonidae", "Paguroidea", "Lestes", "Trichoptera",
                  "Tipulidae", "Jaera hopeana")
    shell[is.na(shell) & x$scientificName %in% no_shell] <- FALSE
    ltype <- ifelse(cl %in% "Gastropoda", "shell_height",
                    ifelse(cl %in% "Bivalvia", "shell_length", "body_length"))
    d <- data.frame(source_id = id, taxon = x$scientificName, length_type = ltype,
                    length_unit = "mm", mass_currency = "dry_mass", mass_unit = "mg",
                    shell_included = shell,
                    preservation = if (id == "cadses_boxcorer") "formalin" else "unknown",
                    life_stage = NA_character_, site = x$locality, country = x$country,
                    L = L, M = DW, stringsAsFactors = FALSE)
    use <- !is.na(d$L) & !is.na(d$M) & d$L > 0 & d$M > 0 & !incons & !is.na(d$taxon)
    say("- %s: %d individuals; %d with length and dry mass; %d non-numeric lengths; %d with ash-free dry mass above dry mass (left out); %d used.",
        id, nrow(x), sum(!is.na(L) & !is.na(DW)), bad_L, sum(incons), sum(use))
    f <- fit_dataset(d[use, ])
    report_fits(f, id)
    all_fits[[id]] <- f

    # Dry mass -> ash-free dry mass: ash fractions come from pooled crucibles,
    # so each distinct value within taxon x lagoon x date counts once.
    # AFDW equal to DW (no ash recorded) is not a measured ash fraction.
    ok <- use & !is.na(AFDW) & AFDW > 0 & AFDW < DW
    r <- AFDW[ok] / DW[ok]
    grp <- paste(x$scientificName[ok], x$locality[ok], x$eventDate[ok], round(r, 4), sep = "\r")
    first <- !duplicated(grp)
    tx <- x$scientificName[ok][first]
    rr <- r[first]
    nind <- tapply(rep(1, sum(ok)), x$scientificName[ok], sum)
    conv_rows[[id]] <- data.frame(
      source_id = id, taxon = names(tapply(rr, tx, median)),
      from = "dry_mass", to = "afdm",
      factor_median = as.numeric(tapply(rr, tx, median)),
      factor_min = as.numeric(tapply(rr, tx, min)),
      factor_max = as.numeric(tapply(rr, tx, max)),
      n_groups = as.integer(tapply(rr, tx, length)),
      n_individuals = as.integer(nind[names(tapply(rr, tx, median))]),
      shell_included = shell[ok][first][match(names(tapply(rr, tx, median)), tx)],
      stringsAsFactors = FALSE)
  }
  conv <- do.call(rbind, conv_rows)
  rownames(conv) <- NULL
  conv$permission <- reg$permission[match(conv$source_id, reg$source_id)]
  write.csv(conv, file.path(out_dir, "conversions_cadses.csv"), row.names = FALSE,
            fileEncoding = "UTF-8")
  say("- Dry mass -> ash-free dry mass: %d taxon factors (median of the distinct ash fractions of each taxon; %d with at least 3 groups).",
      nrow(conv), sum(conv$n_groups >= 3))
}

# 2. Arranz et al. (2022): Iberian stream fish ---------------------------------
# Fork length (mm) except total length in eel, mosquitofish, blenny and catfish
# (caption of their Table 2); mass in g, weighed in the field (fresh, assumed).
# Ameiurus melas (catfish, not in their Table 2) taken as total length.

if (use_source("arranz2022", reg)) {
  say("")
  say("## Arranz et al. (2022)")
  say("")
  path <- source_files("arranz2022", reg)
  x <- read.delim(path[1], fileEncoding = "latin1", stringsAsFactors = FALSE)
  tl <- c("Anguilla anguilla", "Gambusia holbrooki", "Salaria fluviatilis", "Silurus glanis",
          "Ameiurus melas")
  d <- data.frame(source_id = "arranz2022", taxon = x$Scientific_name,
                  length_type = ifelse(x$Scientific_name %in% tl, "total_length", "fork_length"),
                  length_unit = "mm", mass_currency = "wet_mass", mass_unit = "g",
                  shell_included = FALSE, preservation = "fresh", life_stage = NA_character_,
                  site = as.character(x$Code_ID), country = "Spain",
                  L = x$Length_mm, M = x$Mass_g, stringsAsFactors = FALSE)
  say("- %d individuals of %d species from %d sites; masses recorded to 0.1 g (%d individuals below 1 g).",
      nrow(d), length(unique(d$taxon)), length(unique(d$site)), sum(d$M < 1))
  f <- fit_dataset(d)
  report_fits(f, "arranz2022")
  all_fits[["arranz2022"]] <- f

  # Check: their Table 2 (log10 M = a + b log10 L, species with > 25
  # individuals). Values read from the article page (PMC9120222); the sample
  # sizes of that reading were garbled, so only a, b and r2 are compared. The
  # paper removed 405 outliers (residuals above 2 SD); the deposited data seem
  # to be those that remained, since fits to all deposited individuals give
  # the published coefficients exactly in several species.
  pub <- data.frame(
    taxon = c("Achondrostoma arcasii", "Alburnus alburnus", "Anguilla anguilla",
              "Barbatula quignardi", "Barbus haasi", "Barbus meridionalis", "Cobitis calderoni",
              "Cyprinus carpio", "Gambusia holbrooki", "Gobio lozanoi", "Lepomis gibbosus",
              "Luciobarbus graellsii", "Parachondrostoma miegii", "Phoxinus phoxinus",
              "Salaria fluviatilis", "Salmo trutta", "Silurus glanis", "Squalius laietanus"),
    a = c(-5.494, -5.651, -6.282, -5.356, -5.141, -4.964, -5.408, -4.763, -4.586, -5.344,
          -5.099, -4.869, -5.252, -5.445, -5.344, -4.984, -5.121, -5.336),
    b = c(3.342, 3.318, 3.218, 3.140, 3.132, 3.101, 2.985, 3.022, 2.750, 3.232, 3.225,
          2.996, 3.152, 3.296, 3.198, 3.023, 2.967, 3.211),
    r2 = c(0.956, 0.961, 0.964, 0.966, 0.987, 0.993, 0.784, 0.996, 0.977, 0.949, 0.947,
           0.990, 0.981, 0.958, 0.957, 0.991, 0.999, 0.988),
    stringsAsFactors = FALSE)
  chk <- do.call(rbind, lapply(split(d, d$taxon), function(s) {
    if (nrow(s) <= 25) return(NULL)
    f2 <- lm(log10(s$M) ~ log10(s$L))
    data.frame(taxon = s$taxon[1], n = nrow(s), a = unname(coef(f2)[1]),
               b = unname(coef(f2)[2]), r2 = summary(f2)$r.squared, stringsAsFactors = FALSE)
  }))
  m <- merge(pub, chk, by = "taxon", suffixes = c("_pub", "_ours"))
  say("")
  same <- abs(round(m$a_ours, 3) - m$a_pub) < 1e-9 & abs(round(m$b_ours, 3) - m$b_pub) < 1e-9
  say("Check against their Table 2 (all deposited individuals, log10 scale; %d of %d species): %d species give the published a and b to 3 decimals; in the rest, differences up to %.3f (a) and %.3f (b), probably individuals left out of the deposited data or of their regressions.",
      nrow(m), nrow(pub), sum(same), max(abs(m$a_pub - m$a_ours)), max(abs(m$b_pub - m$b_ours)))
  say("")
  say("| Species | n deposited | a published | a ours | b published | b ours | r2 published | r2 ours |")
  say("|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(m))) {
    say("| *%s* | %d | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |", m$taxon[i], m$n[i],
        m$a_pub[i], m$a_ours[i], m$b_pub[i], m$b_ours[i], m$r2_pub[i], m$r2_ours[i])
  }
  arranz_check <- m
}

# 3. Koning et al. (2024): Mekong 'dai' fishery --------------------------------
# Total and standard length (cm, converted to mm), mass (g), fresh (assumed).

if (use_source("koning2024", reg)) {
  say("")
  say("## Koning et al. (2024)")
  say("")
  path <- source_files("koning2024", reg)
  rd <- function(p) read.csv(p, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  tl <- rd(path["Koning2024/Dai_2020_2021_TL.csv"])
  sl <- rd(path["Koning2024/Dai_2020_2021_SL.csv"])
  mk <- function(x, lcol, ltype) {
    data.frame(source_id = "koning2024", taxon = gsub("_", " ", x$Genus_species),
               length_type = ltype, length_unit = "mm", mass_currency = "wet_mass",
               mass_unit = "g", shell_included = FALSE, preservation = "fresh",
               life_stage = NA_character_, site = "Tonle Sap River, dai fishery",
               country = "Cambodia", L = 10 * x[[lcol]], M = x$W.g, stringsAsFactors = FALSE)
  }
  d <- rbind(mk(tl, "tl.cm", "total_length"), mk(sl, "sl.cm", "standard_length"))
  say("- %d (total length) and %d (standard length) individuals of %d species.",
      nrow(tl), nrow(sl), length(unique(d$taxon)))
  f <- fit_dataset(d, population = FALSE)
  report_fits(f, "koning2024")
  all_fits[["koning2024"]] <- f

  # Check: their regressions, log10 W (g) on log10 L (cm), all individuals,
  # coefficients rounded to 2 decimals.
  pubrd <- function(p, A, B) {
    x <- rd(p)
    data.frame(taxon = gsub("_", " ", x$Genus_species), n_pub = x$num_obs,
               A = x[[A]], B = x[[B]], stringsAsFactors = FALSE)
  }
  ptl <- pubrd(path["Koning2024/Dai_2020_2021_TL_reg_results.csv"], "A.tl", "B.tl")
  psl <- pubrd(path["Koning2024/Dai_2020_2021_SL_reg_results.csv"], "A.sl", "B.sl")
  ours <- function(x, lcol) {
    do.call(rbind, lapply(split(x, x$Genus_species), function(s) {
      s <- s[is.finite(s[[lcol]]) & is.finite(s$W.g) & s[[lcol]] > 0 & s$W.g > 0, ]
      f1 <- lm(log10(s$W.g) ~ log10(s[[lcol]]))
      data.frame(taxon = gsub("_", " ", s$Genus_species[1]), n = nrow(s),
                 A_ours = unname(coef(f1)[1]), B_ours = unname(coef(f1)[2]),
                 stringsAsFactors = FALSE)
    }))
  }
  for (k in list(list("total length", ptl, ours(tl, "tl.cm")),
                 list("standard length", psl, ours(sl, "sl.cm")))) {
    m <- merge(k[[2]], k[[3]], by = "taxon")
    same <- abs(round(m$A_ours, 2) - m$A) < 1e-9 & abs(round(m$B_ours, 2) - m$B) < 1e-9
    say("- Check, %s: %d of %d species give the published log10 a and b after rounding to 2 decimals; largest differences %.3f (log10 a) and %.3f (b); sample sizes equal in %d.",
        k[[1]], sum(same), nrow(m), max(abs(m$A_ours - m$A)), max(abs(m$B_ours - m$B)),
        sum(m$n == m$n_pub))
    if (any(!same)) {
      say("  Species that differ: %s.", paste(sprintf("*%s* (n %d vs %d; %.3f vs %.2f; %.3f vs %.2f)",
          m$taxon[!same], m$n[!same], m$n_pub[!same], m$A_ours[!same], m$A[!same],
          m$B_ours[!same], m$B[!same]), collapse = "; "))
    }
  }
}

# 3b. NEON: fish of North American streams and lakes --------------------------
# Total length (mm) and fresh mass (g) of the fish measured one by one in
# 2016-2024 (RELEASE-2026), read with data-raw/neon.R. Only the fish longer
# than the length at which the floor of the weights (0.3 or 0.5 g recorded
# for lighter fish) stops mattering for their species enter, so that the
# floor does not bend the equations; their length ranges start there.

if (use_source("neon_fish", reg)) {
  say("")
  say("## NEON")
  say("")
  source("data-raw/neon.R")
  path <- source_files("neon_fish", reg)
  nf <- read_neon_fish(path[1])
  cn <- nf$counts
  x <- nf$fish[nf$fish$for_lwr, ]
  d <- data.frame(source_id = "neon_fish", taxon = x$taxon, length_type = "total_length",
                  length_unit = "mm", mass_currency = "wet_mass", mass_unit = "g",
                  shell_included = FALSE, preservation = "fresh", life_stage = NA_character_,
                  site = x$site, country = "United States", L = x$length, M = x$weight,
                  stringsAsFactors = FALSE)
  say("- %d fish measured, %d of RELEASE-2026, %d of them identified to species and weighed; %d weights of %.1f g or less (%d of %.1f g), the minimum of the scale (censored: lighter than %.2f g); %d gross errors left out (data-raw/neon.R).",
      cn[["all"]], cn[["release"]], cn[["weighed"]], cn[["censored"]], neon_settings$floor_max,
      cn[["at_0.3"]], 0.3, neon_settings$limit, cn[["gross_uncensored"]] + cn[["gross_censored"]])
  fl <- nf$equations$floor_length
  say("- %d fish of %d species from %d sites, longer than the floor length of their species (the length at which it weighs %.1f g; median %.0f mm, range %.0f-%.0f mm, %d species), enter the equations.",
      nrow(d), length(unique(d$taxon)), length(unique(d$site)), 2 * neon_settings$limit,
      stats::median(fl), min(fl), max(fl), length(fl))
  f <- fit_dataset(d)
  report_fits(f, "neon_fish")
  all_fits[["neon_fish"]] <- f
}

# 4. Bagley et al. (2022): Central American fishes -----------------------------
# Standard length (mm) and mass (g) of specimens preserved in 95% ethanol;
# total length for some.

if (use_source("bagley2022", reg)) {
  say("")
  say("## Bagley et al. (2022)")
  say("")
  path <- source_files("bagley2022", reg)
  x <- as.data.frame(suppressMessages(readxl::read_excel(path[1], col_types = "text",
                                                          .name_repair = "unique_quiet")))
  names(x) <- make.names(names(x))
  # The published equations include the specimens labelled "sp." in the
  # species of their genus: their Table 2 gives n = 49 and SL 19.33-53.91 mm
  # for Poecilia gillii (47 + 2 "Poecilia sp.") and n = 14 and SL 32.82-57.43
  # mm for Atherinella hubbsi (13 + 1 "Atherinella sp."). Followed here.
  x$Species[x$Species == "Poecilia sp."] <- "Poecilia gillii"
  x$Species[x$Species == "Atherinella sp."] <- "Atherinella hubbsi"
  W <- num(x$Weight..g.)
  d <- rbind(
    data.frame(source_id = "bagley2022", taxon = x$Species, length_type = "standard_length",
               length_unit = "mm", mass_currency = "wet_mass", mass_unit = "g",
               shell_included = FALSE, preservation = "ethanol", life_stage = NA_character_,
               site = x$Locality, country = x$Country, L = num(x$SL..mm.), M = W,
               stringsAsFactors = FALSE),
    data.frame(source_id = "bagley2022", taxon = x$Species, length_type = "total_length",
               length_unit = "mm", mass_currency = "wet_mass", mass_unit = "g",
               shell_included = FALSE, preservation = "ethanol", life_stage = NA_character_,
               site = x$Locality, country = x$Country, L = num(x$TL..mm.), M = W,
               stringsAsFactors = FALSE))
  say("- %d specimens; %d with mass and standard length, %d with mass and total length.",
      nrow(x), sum(!is.na(W) & !is.na(num(x$SL..mm.))), sum(!is.na(W) & !is.na(num(x$TL..mm.))))
  f <- fit_dataset(d, population = FALSE)
  report_fits(f, "bagley2022")
  all_fits[["bagley2022"]] <- f

  # Check: their Table 2 (W = a SL^b, g and mm). Outliers were removed by eye
  # in the four poeciliids (n after removal 92, 42, 48 and 63), which cannot be
  # reproduced; the cichlids and the silverside use all specimens.
  pub <- data.frame(
    taxon = c("Alfaro cultratus", "Phallichthys amates", "Poecilia gillii",
              "Priapichthys annectens", "Parachromis dovii", "Parachromis managuensis",
              "Atherinella hubbsi"),
    a = c(2.3e-6, 5.3e-6, 1.7e-6, 1.5e-5, 1.4e-5, 1.9e-5, 9.7e-6),
    b = c(3.446, 3.439, 3.696, 3.0901, 3.1588, 3.1051, 2.9359), stringsAsFactors = FALSE)
  sl_d <- d[d$length_type == "standard_length" & is.finite(d$L) & is.finite(d$M), ]
  chk <- do.call(rbind, lapply(split(sl_d, sl_d$taxon), function(s) {
    f1 <- lm(log10(s$M) ~ log10(s$L))
    data.frame(taxon = s$taxon[1], n = nrow(s), a_ours = 10^unname(coef(f1)[1]),
               b_ours = unname(coef(f1)[2]), stringsAsFactors = FALSE)
  }))
  m <- merge(pub, chk, by = "taxon")
  say("")
  say("| Species | n | a published | a ours | b published | b ours |")
  say("|---|---|---|---|---|---|")
  for (i in seq_len(nrow(m))) {
    say("| *%s* | %d | %.2g | %.3g | %s | %.4f |", m$taxon[i], m$n[i], m$a[i], m$a_ours[i],
        format(m$b[i]), m$b_ours[i])
  }
  bagley_check <- m
}

# 5. Coughlan et al. (2021): invasive bivalves ---------------------------------
# Shell length (mm); masses (g, converted to mg): live, wet, dry (shell and
# soft tissue), shell-free dry and ash-free dry (whole specimen).

if (use_source("coughlan2021", reg)) {
  say("")
  say("## Coughlan et al. (2021)")
  say("")
  path <- source_files("coughlan2021", reg)
  species <- c("C. fluminea" = "Corbicula fluminea", "C. largillierti" = "Corbicula largillierti",
               "D. bugensis" = "Dreissena bugensis", "D. polymorpha" = "Dreissena polymorpha",
               "L. fortunei" = "Limnoperna fortunei", "S. Woodiana" = "Sinanodonta woodiana")
  currencies <- list(LW = c("live_mass", "TRUE"), WW = c("wet_mass", "TRUE"),
                     DW = c("dry_mass", "TRUE"), SFDW = c("dry_mass", "FALSE"),
                     AFDW = c("afdm", "TRUE"))
  rows <- list()
  ratio_ok <- ratio_n <- 0
  for (sh in names(species)) {
    x <- as.data.frame(readxl::read_excel(path[1], sheet = sh, col_types = "text"))
    v <- function(col) num(x[[col]])
    # Internal check: the ratio columns are 100 x the ratio of the mass columns.
    for (rc in c("DW/WW", "SFDW/DW", "AFDW/DW", "SW/DW")) {
      parts <- strsplit(rc, "/", fixed = TRUE)[[1]]
      calc <- 100 * v(parts[1]) / v(parts[2])
      given <- v(rc)
      both <- is.finite(calc) & is.finite(given)
      ratio_n <- ratio_n + sum(both)
      ratio_ok <- ratio_ok + sum(abs(calc[both] - given[both]) <= 1e-6 * pmax(1, abs(given[both])))
    }
    for (cu in names(currencies)) {
      rows[[paste(sh, cu)]] <- data.frame(
        source_id = "coughlan2021", taxon = species[[sh]], length_type = "shell_length",
        length_unit = "mm", mass_currency = currencies[[cu]][1], mass_unit = "mg",
        shell_included = currencies[[cu]][2], preservation = "fresh_or_thawed",
        life_stage = NA_character_, site = x$Site, country = x$Origin,
        L = v("Length"), M = 1000 * v(cu), stringsAsFactors = FALSE)
    }
  }
  d <- do.call(rbind, rows)
  say("- %d individuals of %d species; ratio columns equal to 100 x the ratio of the mass columns in %d of %d values.",
      sum(d$mass_currency == "live_mass"), length(species), ratio_ok, ratio_n)
  f <- fit_dataset(d)
  report_fits(f, "coughlan2021")
  all_fits[["coughlan2021"]] <- f
}

# 6. Santini et al. (2018): amphibians, individual measurements -----------------
# Supplementary Table S1, Data_type "Individual": mass (g) and snout-vent
# length (mm) of the same adult, sexes pooled or by sex (all pooled here).
# The original reference of each record is in the source (column Ref).

if (use_source("santini2018", reg)) {
  say("")
  say("## Santini et al. (2018), individual measurements")
  say("")
  path <- source_files("santini2018", reg)
  x <- as.data.frame(suppressMessages(readxl::read_excel(path[1], sheet = "Data")))
  x <- x[x$Data_type %in% "Individual", ]
  # A species name repeated many times in one cell ("Polypedates
  # maculatusPolypedates maculatus...") is reduced to one.
  x$Species <- sub("^(.+?)\\1+$", "\\1", x$Species, perl = TRUE)
  pairs <- list(c("M_mean", "SVL_mean"), c("M_m_mean", "SVL_m_mean"), c("M_f_mean", "SVL_f_mean"))
  d <- do.call(rbind, lapply(pairs, function(p) {
    m <- num(x[[p[1]]])
    L <- num(x[[p[2]]])
    ok <- !is.na(m) & !is.na(L)
    data.frame(source_id = "santini2018", taxon = x$Species[ok], length_type = "snout_vent_length",
               length_unit = "mm", mass_currency = "wet_mass", mass_unit = "g",
               shell_included = FALSE, preservation = "unknown", life_stage = "adult",
               site = x$Ref[ok], country = NA_character_, L = L[ok], M = m[ok],
               stringsAsFactors = FALSE)
  }))
  say("- %d individuals with mass and snout-vent length (sexes pooled or by sex) of %d species from %d references.",
      nrow(d), length(unique(d$taxon)), length(unique(d$site)))
  # Rhinella jimi (173 toads from Tolledo and Toledo 2015): mass and length are
  # swapped in Table S1 ("masses" of 69-177 g and "lengths" of 31-786 mm; this
  # toad reaches about 180 mm). Swapped back.
  sw <- d$taxon == "Rhinella jimi"
  if (any(sw)) {
    tmp <- d$L[sw]
    d$L[sw] <- d$M[sw]
    d$M[sw] <- tmp
    say("- *Rhinella jimi*: mass and snout-vent length swapped back in %d individuals (Table S1 gives lengths of up to 786 mm; this toad reaches about 180 mm).",
        sum(sw))
  }
  f <- fit_dataset(d, population = FALSE)
  report_fits(f, "santini2018")
  all_fits[["santini2018"]] <- f
}

# Output -------------------------------------------------------------------------

fits <- do.call(rbind, all_fits)
rownames(fits) <- NULL
fitted <- fits[is.na(fits$reason), ]
fitted$reason <- NULL
fitted <- fitted[order(fitted$source_id, fitted$taxon, fitted$length_type,
                       fitted$mass_currency, fitted$shell_included, fitted$scope,
                       fitted$site, method = "radix"), ]
fitted$equation_id <- sprintf("%s_fit_%04d", fitted$source_id,
                              ave(seq_len(nrow(fitted)), fitted$source_id, FUN = seq_along))
fitted$permission <- reg$permission[match(fitted$source_id, reg$source_id)]
fl <- cbind(ifelse(fitted$r2 < flag_r2, "low_r2", NA),
            ifelse(fitted$b < flag_b[1] | fitted$b > flag_b[2], "b_outside_range", NA))
fitted$quality_flag <- apply(fl, 1, function(v) if (all(is.na(v))) NA_character_ else paste(na.omit(v), collapse = ";"))
front <- c("equation_id", "source_id", "permission", "taxon", "scope", "site", "n_sites",
           "country")
fitted <- fitted[, c(front, setdiff(names(fitted), front))]
write.csv(fitted, file.path(out_dir, "lwr_fitted.csv"), row.names = FALSE, fileEncoding = "UTF-8")

say("")
say("## Output")
say("")
say("`lwr_fitted.csv`: %d equations (%d pooled, %d by site) for %d taxa from %d sources.",
    nrow(fitted), sum(fitted$scope == "pooled"), sum(fitted$scope == "population"),
    length(unique(fitted$taxon)), length(unique(fitted$source_id)))
pf <- fitted[fitted$scope == "pooled", ]
say("")
say("Flagged pooled equations (r2 < %.1f or b outside %.1f-%.1f): %d, mostly %s.",
    flag_r2, flag_b[1], flag_b[2], sum(!is.na(pf$quality_flag)),
    paste(names(head(sort(table(pf$source_id[!is.na(pf$quality_flag)]), decreasing = TRUE), 2)), collapse = " and "))
say("Flagged taxa: %s.", paste(sprintf("*%s* (b = %.2f, r2 = %.2f)", pf$taxon[!is.na(pf$quality_flag)],
                                        pf$b[!is.na(pf$quality_flag)], pf$r2[!is.na(pf$quality_flag)]),
                                collapse = ", "))
tab <- table(fitted$source_id, fitted$mass_currency)
say("")
say("| Source | %s |", paste(colnames(tab), collapse = " | "))
say("|---|%s", paste(rep("---|", ncol(tab)), collapse = ""))
for (i in seq_len(nrow(tab))) say("| %s | %s |", rownames(tab)[i], paste(tab[i, ], collapse = " | "))
writeLines(log_lines, file.path(out_dir, "fit_lwr_log.md"), useBytes = TRUE)
cat(log_lines, sep = "\n")
