# Conversion factors between mass currencies (version 0.2.0)
#
# Builds, in data-raw/derived/:
#   conversions.csv       one row per factor: value in `to` = factor x value in `from`
#   conversions_log.md    counts, reinterpretations and checks
# from:
#   - Hrycik et al. (2024), Table 3 (transcription in data-raw/hrycik2024_table3.csv)
#   - Coughlan et al. (2021): ratios of each bivalve, summarized by species and site
#   - CADSES lagoons: data-raw/derived/conversions_cadses.csv (run fit_lwr.R first)
#   - Brey et al. (2010) data bank, read with read_brey_databank() of
#     data-raw/brey.R; only the parameters of the models fitted to it are
#     distributed
#   - Energy densities and water contents from other articles (hand
#     transcription in data-raw/energy_literature.csv)
#   - Cummins and Wuycheck (1971), Table 3, aquatic macroconsumers (hand
#     transcription in data-raw/cummins1971_table3.csv, with the original
#     sources in data-raw/cummins1971_sources.csv)
#   - Evans-White et al. (2005), Appendix: carbon, nitrogen and phosphorus of
#     stream invertebrates (transcription in data-raw/evanswhite2005_appendix.csv)
#   - Hébert et al. (2016): carbon, nitrogen and phosphorus of crustacean
#     zooplankton (from version 0.6.0)
#   - A compilation of the carbon, nitrogen and phosphorus of freshwater
#     animals from 18 data sets, articles and theses (from version 0.6.0)
#   - A compilation of the mass ratios, carbon and energy densities of
#     freshwater molluscs from articles and theses (from version 0.6.0)
#   - A compilation of the shells and soft tissue of snails, and Kołodziejczyk
#     (1984), Table IV (hand transcription in
#     data-raw/kolodziejczyk1984_table4.csv) (from version 0.6.0)
#   - A compilation of the mass ratios, carbon and energy densities of
#     freshwater zooplankton from articles (from version 0.6.0)
#   - The Limno-STOICH database (Corman et al. 2026): carbon, nitrogen and
#     phosphorus of aquatic animals of inland waters (from version 0.6.0)
#   - A compilation of the remaining conversion factors of freshwater animals
#     (triclads, fish) from articles (from version 0.6.0)
#
# Run from the package root, in a UTF-8 locale:  Rscript data-raw/conversions.R
# Packages used only here: readxl.
#
# Currencies: live_mass (bivalves with the water of the mantle cavity),
# wet_mass, dry_mass, afdm, carbon, nitrogen, phosphorus, energy.
# `from_shell` and `to_shell` say whether that mass includes a shell or
# carbonate skeleton (FALSE also for animals without one), as
# `shell_included` in the other derived tables.
# `factor_unit` is "ratio" (mass over mass) or "J/ug".
#
# `use_in_model` = FALSE marks rows kept for the record but not to be used to
# fit the conversion model; `flag` says why (or what was reinterpreted).
# `ln_mean` and `ln_sd` are the mean and SD of the log ratio of individuals
# (only where individual data exist).

source("data-raw/sources.R")

out_dir <- "data-raw/derived"
reg <- read_sources()
log_lines <- c("# Conversion factors between mass currencies", "",
               paste0("Built by `data-raw/conversions.R` on ", format(Sys.Date()), " (",
                      R.version.string, ")."), "")
say <- function(...) log_lines <<- c(log_lines, sprintf(...))
num <- function(x) suppressWarnings(as.numeric(x))

# `sample` identifies the material behind several factors (a record of the
# data bank of Brey et al., a row of a table): factors with the same `sample`
# share their measurements, which fit_conversions.R uses to check them.
# `life_stage` is "larva" for the tadpoles (their factors have models of
# their own, apart from adult amphibians, in fit_conversions.R) and NA
# otherwise.
cols <- c("source_id", "taxon", "name_as_given", "habitat", "life_stage", "from", "from_shell", "to",
          "to_shell", "factor_unit", "statistic", "factor", "factor_sd", "factor_min",
          "factor_max", "ln_mean", "ln_sd", "n", "n_unit", "scope", "site",
          "original_reference", "sample", "use_in_model", "flag", "note")
std <- function(...) {
  x <- data.frame(..., stringsAsFactors = FALSE)
  for (k in setdiff(cols, names(x))) x[[k]] <- NA
  x[, cols]
}
conv <- list()

# 1. Hrycik et al. (2024), Table 3 ----------------------------------------------
# Mean (and SD when more than one record was available) of DW/WW ratios of
# freshwater benthic invertebrates from published and unpublished sources.
# Columns as printed: "DW/WW without shells" (caption: DW is shell-free for
# molluscs), SFDW/SFWW, and "DW/WW with shells". Reading used here:
#   - main column, molluscs: shell-free DW over wet mass WITH shell. Their wet
#     weights include shells and "the amount of soft tissue in the mollusks was
#     calculated from conversion factors based on the percentage of shell
#     weight for each genus" (Dermott, unpubl.); values 0.048-0.088.
#     Values of 0.10 or more cannot be that ratio (it is 0.029-0.066 in the six
#     bivalves of Coughlan et al. 2021) and match shell-free tissue ratios: in
#     both Dreissena they equal the SFDW/SFWW column, and Brey et al. (2010)
#     ratios refer to the shell-free body. They are read as SFDW/SFWW.
#   - "DW/WW with shells": 0.295-0.366 in gastropods and Sphaerium; 0.047 and
#     0.068 in Dreissena, where DW/WW with shell is about 0.64 (Coughlan et al.
#     2021) and SFDW/WW about 0.04: values below 0.15 are read as shell-free
#     DW over wet mass with shell.
#   - SFDW/SFWW above 0.30 (Valvata sincera, 0.372) is not plausible for soft
#     tissue (0.116-0.256 in the other rows) and is not used.
# Rows citing only Brey et al. (2010) summarize records of that data bank,
# which enters directly (section 4): kept, but not used in the model, and
# compared below with our reading of the data bank. Values repeated without SD
# by Dermott (unpubl.) for three or more taxa of the same group (e.g. 0.160 in
# 19 oligochaetes) look like group values applied to each taxon: the taxon rows
# are not used and each value is kept once at the level of the group.

if (use_source("hrycik2024", reg)) {
  t3 <- read.csv("data-raw/hrycik2024_table3.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
  t3$taxon <- sub(" (sp|spp)\\.$", "", t3$name)
  t3$taxon[t3$name == "Amnicola limnosus"] <- "Amnicola limosa"
  t3$group <- ifelse(is.na(t3$major_group) | !nzchar(t3$major_group), t3$phylum, t3$major_group)
  mollusc <- t3$phylum == "Mollusca"
  columns <- list(main = c("dw_ww_shell_free", "dw_ww_shell_free_sd", "DW/WW without shells"),
                  sf = c("sfdw_sfww", "sfdw_sfww_sd", "SFDW/SFWW"),
                  ws = c("dw_ww_with_shell", "dw_ww_with_shell_sd", "DW/WW with shells"))
  h <- do.call(rbind, lapply(names(columns), function(k) {
    v <- t3[[columns[[k]][1]]]
    i <- which(!is.na(v))
    data.frame(i = i, column = k, label = columns[[k]][3], factor = v[i],
               factor_sd = t3[[columns[[k]][2]]][i], stringsAsFactors = FALSE)
  }))
  m <- mollusc[h$i]
  h$from_shell <- m & h$column %in% c("main", "ws")
  h$to_shell <- m & h$column == "ws"
  h$use <- TRUE
  h$flag <- NA_character_
  add_flag <- function(k, text) {
    h$flag[k] <<- ifelse(is.na(h$flag[k]), text, paste(h$flag[k], text, sep = "; "))
  }
  # Reinterpretations (see above).
  k <- m & h$column == "main" & h$factor >= 0.10
  h$from_shell[k] <- FALSE
  add_flag(which(k), "read as shell-free DW over shell-free WW (value too high for DW over WW with shell)")
  sf_same_row <- t3$sfdw_sfww[h$i]
  k_dup <- which(k & !is.na(sf_same_row) & abs(h$factor - sf_same_row) < 1e-9)
  h$use[k_dup] <- FALSE
  add_flag(k_dup, "duplicate of the SFDW/SFWW value of the same row")
  k <- which(m & h$column == "ws" & h$factor < 0.15)
  h$to_shell[k] <- FALSE
  add_flag(k, "read as shell-free DW over WW with shell (DW/WW with shell is about 0.64 in Dreissena)")
  k <- which(h$column == "sf" & h$factor > 0.30)
  h$use[k] <- FALSE
  add_flag(k, "not plausible for soft tissue; probably DW/WW with shell printed in the wrong column")
  key <- paste(t3$taxon[h$i], h$column, h$factor, h$factor_sd)
  k <- which(duplicated(key) & h$column == "sf")
  h$use[k] <- FALSE
  add_flag(k, "same value as an earlier row of the same taxon")
  k <- which(!is.na(h$factor_sd) & h$factor_sd >= h$factor)
  h$use[k] <- FALSE
  add_flag(k, "SD as large as the mean")
  brey_only <- t3$authors[h$i] == "Brey et al. (2010)"
  h$use[brey_only] <- FALSE
  add_flag(which(brey_only), "summary of Brey et al. (2010) records, which enter directly")
  # Values repeated by Dermott (unpubl.) within a group.
  cand <- h$column == "main" & is.na(h$factor_sd) & t3$authors[h$i] == "Dermott (unpubl.)"
  rep_key <- paste(t3$group[h$i], h$factor)
  counts <- table(rep_key[cand])
  repeated <- names(counts)[counts >= 3]
  # Lowest taxon shared by the taxa of each cluster (checked by hand).
  cluster_taxon <- c("Diptera 0.11" = "Chironomidae", "Diptera 0.12" = "Chironomidae",
                     "Gastropoda 0.06" = "Gastropoda", "Oligochaeta 0.16" = "Oligochaeta",
                     "Unionoida 0.06" = "Pisidium", "Unionoida 0.066" = "Sphaeriidae")
  if (!setequal(repeated, names(cluster_taxon))) {
    stop("Repeated values changed; update cluster_taxon: ", paste(repeated, collapse = ", "))
  }
  k <- which(cand & rep_key %in% repeated)
  h$use[k] <- FALSE
  add_flag(k, "value repeated for several taxa of the group; kept once at the level of the group")
  q <- function(i) ifelse(is.na(t3$qualifier[i]) | !nzchar(t3$qualifier[i]), "", paste0(", ", t3$qualifier[i]))
  one <- std(source_id = "hrycik2024", taxon = t3$taxon[h$i], name_as_given = t3$name[h$i],
             habitat = "freshwater", from = "wet_mass", from_shell = h$from_shell,
             to = "dry_mass", to_shell = h$to_shell, factor_unit = "ratio",
             statistic = ifelse(is.na(h$factor_sd), "single", "mean"), factor = h$factor,
             factor_sd = h$factor_sd, n = ifelse(is.na(h$factor_sd), 1, NA), n_unit = "records",
             scope = ifelse(t3$qualifier[h$i] %in% c("Great Lakes", "superior"), "population", "pooled"),
             site = ifelse(t3$qualifier[h$i] %in% "Great Lakes", "Laurentian Great Lakes",
                           ifelse(t3$qualifier[h$i] %in% "superior", "Lake Superior", NA)),
             original_reference = t3$authors[h$i], use_in_model = h$use, flag = h$flag,
             note = paste0("Table 3, row ", t3$row[h$i], " (p. ", t3$page[h$i], "), column '",
                           h$label, "'", q(h$i)))
  groups <- lapply(repeated, function(r) {
    k <- which(cand & rep_key == r)
    std(source_id = "hrycik2024", taxon = cluster_taxon[[r]],
        name_as_given = paste(t3$name[h$i[k]], collapse = "; "), habitat = "freshwater",
        from = "wet_mass", from_shell = h$from_shell[k[1]], to = "dry_mass",
        to_shell = h$to_shell[k[1]], factor_unit = "ratio", statistic = "single",
        factor = h$factor[k[1]], n = length(k), n_unit = "taxa", scope = "pooled",
        original_reference = "Dermott (unpubl.)", use_in_model = TRUE,
        flag = "group value: same value given for each of these taxa",
        note = paste0("Table 3, rows ", paste(t3$row[h$i[k]], collapse = ", ")))
  })
  conv$hrycik <- rbind(one, do.call(rbind, groups))
  say("## Hrycik et al. (2024), Table 3")
  say("")
  say("- %d table rows (transcribed from the page images, pp. 704-708) give %d ratios: %d in the main column, %d SFDW/SFWW and %d DW/WW with shell.",
      nrow(t3), nrow(h), sum(h$column == "main"), sum(h$column == "sf"), sum(h$column == "ws"))
  say("- Molluscs, main column: %d values below 0.10 read as shell-free DW over WW with shell; %d of 0.10 or more read as shell-free DW over shell-free WW (%s).",
      sum(m & h$column == "main" & h$factor < 0.10), sum(m & h$column == "main" & h$factor >= 0.10),
      paste(sprintf("%s %.3f", t3$name[h$i][m & h$column == "main" & h$factor >= 0.10],
                    h$factor[m & h$column == "main" & h$factor >= 0.10]), collapse = ", "))
  say("- Molluscs, DW/WW with shells below 0.15 read as shell-free DW over WW with shell: %s.",
      paste(sprintf("%s %.3f", t3$name[h$i][m & h$column == "ws" & h$factor < 0.15],
                    h$factor[m & h$column == "ws" & h$factor < 0.15]), collapse = ", "))
  say("- Repeated values kept once at the level of the group: %s.",
      paste(sprintf("%s (%s, %d taxa)", names(cluster_taxon), cluster_taxon,
                    as.integer(counts[names(cluster_taxon)])), collapse = "; "))
  say("- Not used in the model: %d of %d ratios (%d Brey et al. 2010 summaries, %d repeated values, %d duplicates, %d implausible).",
      sum(!h$use), nrow(h), sum(brey_only), sum(cand & rep_key %in% repeated),
      sum(grepl("duplicate|same value as", h$flag)),
      sum(grepl("not plausible|SD as large", h$flag)))
  say("")
}

# 2. Coughlan et al. (2021): ratios of each bivalve -----------------------------
# Live mass (LW, after removing the water outside the shell and some of the
# mantle cavity), wet mass (WW, opened and blotted), dry mass of shell and soft
# tissue (DW), shell-free dry mass (SFDW) and ash-free dry mass of the whole
# specimen (AFDW = DW - ash) of the same individuals. Ratios computed from the
# masses, summarized by species and habitat (pooled, over two or more sites)
# and by species and site (population; at least 3 individuals). From version
# 0.6.0 the pooled factor of a species, habitat and conversion is used only
# when no site has 3 or more individuals, as with the length-mass equations
# fitted to individuals (those of each site, otherwise the pooled one;
# model_data.R); otherwise it repeats the individuals of the populations and
# is kept but not used. Ratios above 1 between a mass and a part of it are
# measurement errors and are left out; SFDW/AFDW can exceed 1. The sites are
# freshwater except the Szczecin Lagoon (from version 0.6.0), a brackish lagoon
# of the Baltic Sea at the mouth of the Oder, taken as transitional, as in the
# compilation of molluscs of 2026-11-09 (section 10).

if (use_source("coughlan2021", reg)) {
  path <- source_files("coughlan2021", reg)
  species <- c("C. fluminea" = "Corbicula fluminea", "C. largillierti" = "Corbicula largillierti",
               "D. bugensis" = "Dreissena bugensis", "D. polymorpha" = "Dreissena polymorpha",
               "L. fortunei" = "Limnoperna fortunei", "S. Woodiana" = "Sinanodonta woodiana")
  transitional_sites <- "Szczecin Lagoon"
  currency <- list(LW = c("live_mass", TRUE), WW = c("wet_mass", TRUE), DW = c("dry_mass", TRUE),
                   SFDW = c("dry_mass", FALSE), AFDW = c("afdm", TRUE))
  pairs <- list(c("LW", "WW"), c("LW", "DW"), c("LW", "SFDW"), c("LW", "AFDW"), c("WW", "DW"),
                c("WW", "SFDW"), c("WW", "AFDW"), c("DW", "SFDW"), c("DW", "AFDW"),
                c("AFDW", "SFDW"))
  summ <- function(v, sp, site, scope, p, habitat) {
    std(source_id = "coughlan2021", taxon = sp, name_as_given = sp, habitat = habitat,
        from = currency[[p[1]]][1], from_shell = as.logical(currency[[p[1]]][2]),
        to = currency[[p[2]]][1], to_shell = as.logical(currency[[p[2]]][2]),
        factor_unit = "ratio", statistic = "mean", factor = mean(v), factor_sd = sd(v),
        factor_min = min(v), factor_max = max(v), ln_mean = mean(log(v)), ln_sd = sd(log(v)),
        n = length(v), n_unit = "individuals", scope = scope, site = site,
        original_reference = "Coughlan et al. (2021), data", use_in_model = TRUE,
        note = paste0(p[2], " / ", p[1], " of each individual"))
  }
  rows <- list()
  dropped <- 0
  for (sh in names(species)) {
    x <- as.data.frame(readxl::read_excel(path[1], sheet = sh, col_types = "text"))
    hab <- ifelse(x$Site %in% transitional_sites, "transitional", "freshwater")
    for (p in pairs) {
      v <- num(x[[p[2]]]) / num(x[[p[1]]])
      ok <- is.finite(v) & v > 0
      bad <- ok & v > 1 & !identical(p, c("AFDW", "SFDW"))
      dropped <- dropped + sum(bad)
      ok <- ok & !bad
      for (h in unique(hab[ok])) {
        kh <- ok & hab == h
        pops <- list()
        for (s in unique(x$Site[kh])) {
          k <- kh & x$Site == s
          if (sum(k) >= 3) pops[[length(pops) + 1]] <- summ(v[k], species[[sh]], s, "population", p, h)
        }
        # A pooled factor of a single site would repeat that population; with
        # populations, it repeats their individuals and is not used.
        if (sum(kh) >= 3 && length(unique(x$Site[kh])) > 1) {
          pool <- summ(v[kh], species[[sh]], NA_character_, "pooled", p, h)
          if (length(pops)) {
            pool$use_in_model <- FALSE
            pool$flag <- "pooled over the populations of the species, which enter instead"
          }
          rows[[length(rows) + 1]] <- pool
        }
        rows <- c(rows, pops)
      }
    }
  }
  conv$coughlan <- do.call(rbind, rows)
  pooled <- conv$coughlan[conv$coughlan$scope == "pooled", ]
  say("## Coughlan et al. (2021)")
  say("")
  say("- %d factors (%d by species and habitat, %d by species and site) for %d conversions; %d ratios above 1 left out; %d pooled factors not used (the populations enter instead), %d used (no site with 3 or more individuals).",
      nrow(conv$coughlan), nrow(pooled), sum(conv$coughlan$scope == "population"),
      length(pairs), dropped, sum(!pooled$use_in_model), sum(pooled$use_in_model))
  tr <- conv$coughlan$habitat == "transitional"
  say("- Transitional (%s): %d factors of %s.", paste(transitional_sites, collapse = ", "), sum(tr),
      paste(sprintf("*%s*", unique(conv$coughlan$taxon[tr])), collapse = " and "))
  for (nt in c("DW / WW of each individual", "SFDW / WW of each individual",
               "SFDW / DW of each individual")) {
    s <- pooled[pooled$note == nt & pooled$habitat == "freshwater", ]
    say("- %s (freshwater): %s.", sub(" of each individual", "", nt),
        paste(sprintf("*%s* %.3f", s$taxon, s$factor), collapse = ", "))
  }
  say("")
}

# 3. CADSES lagoons -------------------------------------------------------------
# Dry mass to ash-free dry mass of the whole animal (shell included in
# molluscs and barnacles), from conversions_cadses.csv (fit_lwr.R): median
# over the distinct ash fractions of each taxon, which were applied to many
# individuals.

cad_file <- file.path(out_dir, "conversions_cadses.csv")
if (!file.exists(cad_file)) stop("Run data-raw/fit_lwr.R first: ", cad_file, " not found.")
cc <- read.csv(cad_file, stringsAsFactors = FALSE, encoding = "UTF-8")
cc <- cc[use_source(cc$source_id, reg), ]
conv$cadses <- std(source_id = cc$source_id, taxon = cc$taxon, name_as_given = cc$taxon,
                   habitat = "transitional", from = cc$from, from_shell = cc$shell_included,
                   to = cc$to, to_shell = cc$shell_included, factor_unit = "ratio",
                   statistic = "median", factor = cc$factor_median, factor_min = cc$factor_min,
                   factor_max = cc$factor_max, n = cc$n_groups, n_unit = "ash_fractions",
                   scope = "pooled", original_reference = "CADSES project data (LifeWatch Italy)",
                   use_in_model = TRUE,
                   note = paste0("Median over ", cc$n_groups, " distinct ash fractions of ",
                                 cc$n_individuals, " individuals"))
say("## CADSES lagoons")
say("")
say("- %d dry mass to AFDM factors (%s).", nrow(conv$cadses),
    paste(sprintf("%d %s", as.integer(table(cc$source_id)), names(table(cc$source_id))), collapse = ", "))
say("")

# 4. Brey et al. (2010) data bank -------------------------------------------------
# One row per record and ratio, in the direction given by the data bank (the
# inverted mass ratios returned by read_brey_databank() are dropped). Ratios
# refer to the shell-free body, except coralline red algae, echinoderms and
# brachiopods, where skeleton and tissue were not separated (Brey et al. 2010).

if (use_source("brey2010", reg)) {
  path <- source_files("brey2010", reg)
  tmp <- tempfile("brey")
  dir.create(tmp)
  utils::unzip(path[1], files = "Conversion03.xls", exdir = tmp)
  xls <- file.path(tmp, "Conversion03.xls")
  source("data-raw/brey.R")
  b <- read_brey_databank(xls, habitat = "both", quiet = TRUE)
  ratios <- brey_ratios
  b <- b[paste(b$from, b$to) %in% paste(ratios$from, ratios$to), ]
  # Superfamily and major taxon, not returned by read_brey_databank().
  raw <- as.data.frame(readxl::read_excel(xls, sheet = "Data", col_names = FALSE,
                                          col_types = "text", .name_repair = "minimal"))[-(1:5), ]
  b$record <- sub("^Brey et al\\. 2010, record ([^:]+).*$", "\\1", b$source_detail)
  i <- match(b$record, trimws(raw[[1]]))
  if (anyNA(i)) stop("Records of the Brey data bank not found in the raw sheet")
  b$superfamily <- trimws(raw[[6]][i])
  b$major_taxon <- trimws(raw[[9]][i])
  hard <- b$major_taxon %in% c("Echinodermata", "Brachiopoda") | b$family %in% "Corallinaceae"
  b$reference <- sub("^Brey et al\\. 2010, record [^:]+: ?", "", b$source_detail)
  conv$brey <- std(source_id = "brey2010", taxon = b$taxon, name_as_given = b$taxon,
                   habitat = b$habitat, from = b$from, from_shell = hard, to = b$to,
                   to_shell = hard, factor_unit = ifelse(b$to == "energy", "J/ug", "ratio"),
                   statistic = "single", factor = b$factor, n = 1, n_unit = "records",
                   scope = "record", original_reference = b$reference,
                   sample = paste0("brey2010:", b$record), use_in_model = TRUE,
                   note = paste0("Data bank record ", b$record,
                                 ifelse(hard, "; skeleton included", "")))
  # Energy per dry mass computed in the data bank with the carbon content of
  # the record (J/DM = J/AFDM x C/DM) instead of its AFDM/DM: left out (the
  # energy per AFDM of the same record is used). Found in the 19 records of
  # Salonen et al. (1976), whose J/DM was 0.45-0.62 times their J/AFDM.
  num_raw <- function(col) suppressWarnings(as.numeric(raw[[col]][i]))
  j_dm <- num_raw(43)
  j_af <- num_raw(45)
  af_dm <- num_raw(30)
  c_dm <- num_raw(35)
  carbon_basis <- b$from == "dry_mass" & b$to == "energy" & !is.na(j_dm + j_af + af_dm + c_dm) &
    abs(j_dm / (j_af * c_dm) - 1) < 0.01 & abs(j_dm / (j_af * af_dm) - 1) > 0.05
  carbon_basis[is.na(carbon_basis)] <- FALSE
  conv$brey$use_in_model[carbon_basis] <- FALSE
  conv$brey$flag[carbon_basis] <- "energy per dry mass of the data bank equal to J/AFDM x C/DM of the record (not x AFDM/DM)"
  say("## Brey et al. (2010)")
  say("")
  say("- %d factors (one per record and ratio) of %d taxa: %d freshwater and %d marine; %d with skeleton included (echinoderms, brachiopods, coralline algae).",
      nrow(b), length(unique(b$taxon)), sum(b$habitat == "freshwater"), sum(b$habitat == "marine"),
      sum(hard))
  say("- Energy per dry mass equal to the energy per AFDM times the carbon content of the record, not used: %d (%s).",
      sum(carbon_basis), paste(unique(b$reference[carbon_basis]), collapse = "; "))
  # Check: Hrycik et al. (2024) Table 3 rows citing only Brey et al. (2010)
  # against the mean DM/WM of the matching records of the data bank.
  if (!is.null(conv$hrycik)) {
    hb <- conv$hrycik[conv$hrycik$original_reference %in% "Brey et al. (2010)", ]
    dm <- b[b$from == "wet_mass" & b$to == "dry_mass", ]
    lin <- dm[, c("taxon", "genus", "family", "superfamily", "order", "class", "major_taxon")]
    chk <- do.call(rbind, lapply(seq_len(nrow(hb)), function(r) {
      hit <- rowSums(lin == hb$taxon[r], na.rm = TRUE) > 0
      fw <- hit & dm$habitat == "freshwater"
      data.frame(taxon = hb$taxon[r], column = sub("^.*column '([^']+)'.*$", "\\1", hb$note[r]),
                 hrycik = hb$factor[r], n_fw = sum(fw),
                 mean_fw = if (any(fw)) mean(dm$factor[fw]) else NA,
                 n_all = sum(hit), mean_all = if (any(hit)) mean(dm$factor[hit]) else NA)
    }))
    near <- function(a, b) !is.na(b) & abs(a - b) < 0.0015
    chk$match <- ifelse(near(chk$hrycik, chk$mean_fw), "freshwater",
                        ifelse(near(chk$hrycik, chk$mean_all), "all", "no"))
    say("- Check against the %d values of Hrycik et al. (2024), Table 3, that cite only Brey et al. (2010): the mean DM/WM of the matching records of the data bank gives the printed value (to 3 decimals) in %d with the freshwater records and %d more with all records. Not reproduced: %s (none of these rows is used in the model).",
        nrow(chk), sum(chk$match == "freshwater"), sum(chk$match == "all"),
        paste(sprintf("%s (%s)", chk$taxon[chk$match == "no"], chk$column[chk$match == "no"]),
              collapse = ", "))
    say("")
    say("| Taxon | Column | Hrycik | Brey freshwater (n) | Brey all (n) | Match |")
    say("|---|---|---|---|---|---|")
    for (r in seq_len(nrow(chk))) {
      say("| %s | %s | %.4f | %s (%d) | %s (%d) | %s |", chk$taxon[r], chk$column[r], chk$hrycik[r],
          ifelse(is.na(chk$mean_fw[r]), "-", sprintf("%.4f", chk$mean_fw[r])), chk$n_fw[r],
          ifelse(is.na(chk$mean_all[r]), "-", sprintf("%.4f", chk$mean_all[r])), chk$n_all[r],
          chk$match[r])
    }
  }
  say("")
}

# 5. Energy densities and mass ratios from other articles -------------------------
# Hand transcription in data-raw/energy_literature.csv: the value and unit as
# printed, converted here to J/ug (energy) or to a mass ratio: dry/wet mass
# from % water or % dry mass, AFDM/dry mass from % ash or % AFDM and
# carbon/dry mass from % carbon. Only sources registered as included; rows of
# animals prepared as prey (skin, stomach or spawn removed), of animals
# preserved in ethanol or formalin (preservation changes their wet, dry and
# ash-free dry mass by amounts that depend on the preservative, the time and
# the species; Dermott and Paterson 1974, Table 5), values already in the
# data bank of Brey et al. (2010) and values computed with the ash content of
# another study are kept but not used. Masses with shell keep `from_shell` =
# `to_shell` = TRUE (the whole animal was analysed) and are not modelled.

to_J_ug <- c("J/mg" = 1e-3, "kJ/g" = 1e-3, "J/g" = 1e-6, "cal/g" = 4.184e-6, "cal/mg" = 4.184e-3)
ratio_units <- data.frame(unit = c("% water", "% dry mass", "% ash", "% AFDM", "% carbon"),
                          from = c("wet_mass", "wet_mass", "dry_mass", "dry_mass", "dry_mass"),
                          to = c("dry_mass", "dry_mass", "afdm", "afdm", "carbon"),
                          complement = c(TRUE, FALSE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
el <- read.csv("data-raw/energy_literature.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
included <- reg$source_id[reg$status == "included"]
el <- el[el$source_id %in% included, ]
el <- el[use_source(el$source_id, reg), ]
if (nrow(el)) {
  is_energy <- el$to == "energy"
  ru <- ratio_units[match(el$unit_printed, ratio_units$unit), ]
  stopifnot(all(el$unit_printed[is_energy] %in% names(to_J_ug)),
            !anyNA(ru$unit[!is_energy]), all(el$from[!is_energy] == ru$from[!is_energy]),
            all(el$to[!is_energy] == ru$to[!is_energy]))
  fac <- ifelse(is_energy, el$value_printed * unname(to_J_ug[el$unit_printed]),
                ifelse(ru$complement, 1 - el$value_printed / 100, el$value_printed / 100))
  err <- ifelse(is.na(el$error_printed), "",
                sprintf("; %s %s %s", c(se = "SE", ci95 = "95 % CI +/-", ci99 = "99 % CI +/-",
                                        sd = "SD", pm = "+/- (type not stated)")[el$error_type],
                        el$error_printed, el$unit_printed))
  conv$literature <- std(
    source_id = el$source_id, taxon = el$taxon, name_as_given = el$name_as_given,
    habitat = el$habitat, from = el$from, from_shell = el$from_shell, to = el$to,
    to_shell = el$from_shell, factor_unit = ifelse(is_energy, "J/ug", "ratio"),
    statistic = el$statistic, factor = fac, n = el$n, n_unit = el$n_unit, scope = "pooled",
    original_reference = el$original_reference,
    sample = paste(el$source_id, el$taxon, el$life_stage, el$parts, sep = ":"),
    use_in_model = el$use_in_model, flag = ifelse(nzchar(el$flag), el$flag, NA),
    note = paste0(el$page, ": ", el$value_printed, " ", el$unit_printed, err,
                  ifelse(nzchar(el$life_stage), paste0("; ", el$life_stage), ""),
                  ifelse(nzchar(el$parts), paste0("; ", el$parts), ""),
                  ifelse(nzchar(el$note), paste0("; ", el$note), "")))
  say("## Energy densities and mass ratios from other articles")
  say("")
  say("- %d factors from %d sources (%s); %d not used (%s).", nrow(el),
      length(unique(el$source_id)),
      paste(sprintf("%s %d", names(table(el$source_id)), as.integer(table(el$source_id))), collapse = ", "),
      sum(!el$use_in_model), paste(unique(el$flag[!el$use_in_model]), collapse = "; "))
  say("")
}

# 6. Cummins and Wuycheck (1971), Table 3 ----------------------------------------
# Energy per gram of dry mass, ash-free dry mass and wet mass of the aquatic
# macroconsumers ("III. Macroconsumers, A. Aquatic", pp. 91-119 of the paper),
# transcribed from the page images into data-raw/cummins1971_table3.csv (one
# row per printed row, values as printed; `taxon` holds the name with
# abbreviations, repeat marks and misprints resolved, `habitat` the habitat of
# the row where it differs from that of its source, `curation_flag` the
# reason to leave a row out, and `edit_note` and `transcription_note` what
# was edited and the doubts of the transcription). Calories are converted with
# 4.184 J/cal. The habitat and reference of each row are otherwise those of
# its original source (data-raw/cummins1971_sources.csv). Not used in the
# model:
#   - energy densities of original sources that are in the data bank of Brey
#     et al. (2010), which enters directly (section 4);
#   - energy calculated from protein, fat and carbohydrate (method 7) or from
#     nitrogen (general note 102), not measured by combustion;
#   - eggs, embryos and yolk, exuviae and chitin, fat or fat-free material,
#     feces and parts of animals (stage codes below);
#   - terrestrial animals printed in the aquatic section;
#   - rows with a curation flag.
# Molluscs are taken with their shell unless the stage code (9b) or the notes
# (85e, 105) say that the shell was removed; echinoderms and brachiopods with
# their skeleton, as in Brey et al. (2010). The % H2O and % ash columns give
# wet to dry and dry to ash-free dry mass factors, kept for the next update of
# the mass models but not used yet.

cw_file <- "data-raw/cummins1971_table3.csv"
if (use_source("cummins1971", reg) && file.exists(cw_file)) {
  cw <- read.csv(cw_file, stringsAsFactors = FALSE, encoding = "UTF-8", colClasses = "character",
                 na.strings = "")
  cs <- read.csv("data-raw/cummins1971_sources.csv", stringsAsFactors = FALSE, encoding = "UTF-8",
                 colClasses = "character", na.strings = "")
  cw$src <- sub("[^0-9].*$", "", trimws(cw$source))
  i_src <- match(cw$src, cs$source)
  if (anyNA(i_src)) stop("Sources of Cummins and Wuycheck not in cummins1971_sources.csv: ",
                         paste(unique(cw$source[is.na(i_src)]), collapse = ", "))
  brey_energy <- cs$brey_energy[i_src] %in% "TRUE"
  brey_mass <- cs$brey_mass[i_src] %in% "TRUE"
  hab <- ifelse(is.na(cw$habitat), cs$habitat[i_src], cw$habitat)
  ref <- cs$reference[i_src]
  stage <- trimws(cw$stage)
  notes <- trimws(cw$notes)
  bad_stage <- c(paste0("16", letters[1:11]), paste0("22", letters[1:4]), "12f", "12g", "14g",
                 "18d", "18e", "18f", "10s", "10t", "10u", "24a", "24b")
  why <- rep(NA_character_, nrow(cw))
  addw <- function(k, text) why[k] <<- ifelse(is.na(why[k]), text, paste(why[k], text, sep = "; "))
  addw(which(brey_energy), "original source in the data bank of Brey et al. (2010)")
  addw(which(cw$method %in% "7" | notes %in% "102"),
       "energy calculated from organic composition, not measured")
  addw(which(stage %in% bad_stage),
       "eggs, embryos, yolk, exuviae, feces, fat or part of the animals")
  addw(which(hab %in% "terrestrial"), "terrestrial animal")
  addw(which(is.na(hab)), "habitat of the original source unknown")
  k <- which(!is.na(cw$curation_flag))
  addw(k, cw$curation_flag[k])
  mollusc <- grepl("^Mollusca", cw$higher_taxa)
  shell <- (mollusc & !(stage %in% "9b" | notes %in% c("85e", "105", "105a"))) |
    grepl("^(Echinodermata|Brachiopoda|Brachyopoda)", cw$higher_taxa)
  num <- function(x) suppressWarnings(as.numeric(x))
  # Number of samples as printed ("4*" in the rows of the Polish lakes).
  n_s <- num(sub("[^0-9.].*$", "", cw$n_samples))
  sd_of <- function(v, e, code) {
    e <- num(e)
    out <- rep(NA_real_, length(v))
    k <- code %in% "4"
    out[k] <- e[k]
    k <- code %in% "2" & !is.na(n_s)
    out[k] <- e[k] * sqrt(n_s[k])
    k <- code %in% "3"
    out[k] <- e[k] / 2
    k <- code %in% "5"
    out[k] <- v[k] * e[k] / 100
    out
  }
  where <- paste0("Table 3, p. ", cw$printed_page, " (PDF p. ", cw$pdf_page, "), row ",
                  cw$row_on_page, "; stage ", ifelse(is.na(stage), "-", stage), "; method ",
                  ifelse(is.na(cw$method), "-", cw$method), "; source ", cw$source,
                  ifelse(is.na(notes), "", paste0("; note ", notes)))
  rows <- list()
  for (col in c("dw", "afdw", "ww")) {
    v <- num(cw[[paste0("cal_", col)]])
    k <- which(!is.na(v))
    if (!length(k)) next
    from <- c(dw = "dry_mass", afdw = "afdm", ww = "wet_mass")[[col]]
    sdv <- sd_of(v, cw[[paste0("cal_", col, "_err")]], cw[[paste0("cal_", col, "_code")]])
    rows[[col]] <- std(
      source_id = "cummins1971", taxon = cw$taxon[k], name_as_given = cw$species_name[k],
      habitat = hab[k], from = from, from_shell = shell[k], to = "energy", to_shell = shell[k],
      factor_unit = "J/ug", statistic = ifelse(!is.na(n_s[k]) & n_s[k] > 1, "mean", "single"),
      factor = v[k] * 4.184e-6, factor_sd = sdv[k] * 4.184e-6, n = n_s[k], n_unit = "samples",
      scope = "pooled", original_reference = ref[k],
      sample = paste0("cummins1971:p", cw$pdf_page[k], "r", cw$row_on_page[k]),
      use_in_model = is.na(why[k]), flag = why[k],
      note = paste0(where[k], "; ", v[k], " cal/g"))
  }
  h2o <- num(cw$h2o_pct)
  ash <- num(cw$ash_pct)
  # Mass factors: not used yet; the reasons that would leave them out of the
  # mass models are kept in the flag.
  mflag <- "mass factor of Cummins and Wuycheck (1971): kept for the next update of the mass models"
  mwhy <- rep(NA_character_, nrow(cw))
  addm <- function(k, text) mwhy[k] <<- ifelse(is.na(mwhy[k]), text, paste(mwhy[k], text, sep = "; "))
  addm(which(brey_mass), "original source in the data bank of Brey et al. (2010)")
  addm(which(stage %in% bad_stage), "eggs, embryos, yolk, exuviae, feces, fat or part of the animals")
  addm(which(hab %in% "terrestrial"), "terrestrial animal")
  k <- which(!is.na(cw$curation_flag))
  addm(k, cw$curation_flag[k])
  mflag <- ifelse(is.na(mwhy), mflag, paste(mflag, mwhy, sep = "; "))
  k <- which(!is.na(h2o) & h2o > 0 & h2o < 100)
  if (length(k)) rows$h2o <- std(
    source_id = "cummins1971", taxon = cw$taxon[k], name_as_given = cw$species_name[k],
    habitat = hab[k], from = "wet_mass", from_shell = shell[k], to = "dry_mass", to_shell = shell[k],
    factor_unit = "ratio", statistic = ifelse(!is.na(n_s[k]) & n_s[k] > 1, "mean", "single"),
    factor = 1 - h2o[k] / 100, n = n_s[k], n_unit = "samples", scope = "pooled",
    original_reference = ref[k], sample = paste0("cummins1971:p", cw$pdf_page[k], "r", cw$row_on_page[k]),
    use_in_model = FALSE, flag = mflag[k],
    note = paste0(where[k], "; ", h2o[k], " % H2O"))
  k <- which(!is.na(ash) & ash >= 0 & ash < 100)
  if (length(k)) rows$ash <- std(
    source_id = "cummins1971", taxon = cw$taxon[k], name_as_given = cw$species_name[k],
    habitat = hab[k], from = "dry_mass", from_shell = shell[k], to = "afdm", to_shell = shell[k],
    factor_unit = "ratio", statistic = ifelse(!is.na(n_s[k]) & n_s[k] > 1, "mean", "single"),
    factor = 1 - ash[k] / 100, n = n_s[k], n_unit = "samples", scope = "pooled",
    original_reference = ref[k], sample = paste0("cummins1971:p", cw$pdf_page[k], "r", cw$row_on_page[k]),
    use_in_model = FALSE, flag = mflag[k],
    note = paste0(where[k], "; ", ash[k], " % ash"))
  conv$cummins <- do.call(rbind, rows)
  e <- conv$cummins[conv$cummins$to == "energy", ]
  say("## Cummins and Wuycheck (1971), Table 3")
  say("")
  say("- %d printed rows of aquatic macroconsumers (PDF pp. %s-%s) give %d energy densities (%d per dry mass, %d per ash-free dry mass, %d per wet mass) of %d names, from %d original sources.",
      nrow(cw), min(as.integer(cw$pdf_page)), max(as.integer(cw$pdf_page)), nrow(e),
      sum(e$from == "dry_mass"), sum(e$from == "afdm"), sum(e$from == "wet_mass"),
      length(unique(cw$taxon)), length(unique(cw$src)))
  for (w in c("original source in the data bank of Brey et al. (2010)",
              "energy calculated from organic composition, not measured",
              "eggs, embryos, yolk, exuviae, feces, fat or part of the animals",
              "terrestrial animal", "habitat of the original source unknown", "row inconsistent")) {
    say("- Not used, %s: %d.", w, sum(grepl(w, e$flag, fixed = TRUE)))
  }
  say("- Used: %d energy densities (%d with shell or skeleton, not modelled); habitats %s.",
      sum(e$use_in_model), sum(e$use_in_model & e$from_shell),
      paste(sprintf("%s %d", names(table(e$habitat[e$use_in_model])),
                    as.integer(table(e$habitat[e$use_in_model]))), collapse = ", "))
  say("- Kept for later (not used): %d wet to dry mass factors (%% H2O) and %d dry to ash-free dry mass factors (%% ash).",
      sum(conv$cummins$to == "dry_mass"), sum(conv$cummins$to == "afdm"))
  say("")
}

# 7. Evans-White et al. (2005), Appendix --------------------------------------------
# Total carbon and nitrogen (CHN analyser) and total phosphorus (combustion at
# 500 C, HCl and the ascorbic acid method) as % of dry mass of 46 taxa of
# benthic invertebrates of streams of Indiana-Michigan and Wisconsin: mean, SD
# and range among streams (n = number of streams; in each stream one to three
# pooled samples), transcribed from the Appendix (pp. 1798-1799) into
# data-raw/evanswhite2005_appendix.csv (columns c_*, n_* and p_*; nitrogen and
# phosphorus from version 0.6.0), with `taxon` correcting two misprints
# (`edit_note`). The animals were kept overnight at 4 C to clear their guts,
# frozen, and freeze-dried or oven-dried at 45 C; insects and crustaceans were
# analysed whole and molluscs without shell (pp. 1788-1789).

ew_file <- "data-raw/evanswhite2005_appendix.csv"
if (use_source("evanswhite2005", reg) && file.exists(ew_file)) {
  ew <- read.csv(ew_file, stringsAsFactors = FALSE, encoding = "UTF-8", colClasses = "character",
                 na.strings = "")
  n_ew <- as.integer(ew$n)
  stopifnot(!anyNA(n_ew))
  sites <- c(IM = "Indiana-Michigan", W = "Wisconsin", "IM, W" = "Indiana-Michigan and Wisconsin")
  stopifnot(all(ew$region %in% names(sites)))
  elements <- c(carbon = "c", nitrogen = "n", phosphorus = "p")
  symbol <- c(carbon = "C", nitrogen = "N", phosphorus = "P")
  conv$evanswhite <- do.call(rbind, lapply(names(elements), function(el) {
    col <- function(k) ew[[paste0(elements[[el]], "_", k)]]
    v <- num(col("mean"))
    stopifnot(!anyNA(v), all(v > 0 & v < 100))
    std(source_id = "evanswhite2005", taxon = ew$taxon, name_as_given = ew$taxon_printed,
        habitat = "freshwater", from = "dry_mass", from_shell = FALSE, to = el, to_shell = FALSE,
        factor_unit = "ratio", statistic = ifelse(n_ew > 1, "mean", "single"), factor = v / 100,
        factor_sd = num(col("sd")) / 100, factor_min = num(col("min")) / 100,
        factor_max = num(col("max")) / 100, n = n_ew, n_unit = "streams", scope = "pooled",
        site = unname(sites[ew$region]), original_reference = "Evans-White et al. (2005)",
        sample = paste0("evanswhite2005:p", ew$page, "r", ew$row_on_page), use_in_model = TRUE,
        note = paste0("Appendix, p. ", ew$page, ", row ", ew$row_on_page, "; ", col("mean"), " % ",
                      symbol[[el]],
                      ifelse(is.na(col("sd")), "", paste0(" (SD ", col("sd"), "; range ", col("min"),
                                                          "-", col("max"), ")")),
                      "; ", ifelse(ew$group == "mollusc", "soft tissue, shell removed", "whole animals"),
                      "; total ", el))
  }))
  say("## Evans-White et al. (2005), Appendix")
  say("")
  say("- %d taxa (%s), carbon, nitrogen and phosphorus as %% of dry mass; %d with one stream, the others means of 2-%d streams.",
      nrow(ew), paste(sprintf("%s %d", names(table(ew$group)), as.integer(table(ew$group))), collapse = ", "),
      sum(n_ew == 1), max(n_ew))
  say("")
}

# Range of carbon, nitrogen and phosphorus in whole animals, in % of dry mass,
# used to screen the values of sections 8 and 9: no animal is more than 65 %
# carbon (lipid, the richest component, is 77 %) or 16 % nitrogen (protein),
# nor more than 8 % phosphorus (fish, with bone); below 20 % carbon, 2 %
# nitrogen or 0.05 % phosphorus the sample is mostly mineral (sediment in the
# gut, or a shell). The values of animals with calcified cuticles (ostracods,
# amphipods, soldier fly larvae) lie inside it.
elemental_bounds <- list(carbon = c(20, 65), nitrogen = c(2, 16), phosphorus = c(0.05, 8))
bounds_text <- "C 20-65 %, N 2-16 %, P 0.05-8 % of dry mass"
# First and last of the places a factor was read from (sections 9 to 11).
first_last <- function(x) {
  x <- unique(x[!is.na(x)])
  if (!length(x)) "" else if (length(x) == 1L) x else paste0(x[1], " ... ", x[length(x)], " (", length(x), ")")
}
# Factors from a table with the columns of the compilations of molluscs and
# zooplankton of 2026-11-09 and 2026-11-10 (sections 10 to 12), to which the
# caller has added `source_id`, `taxon`, the conversion (`conv_from_shell`,
# `conv_to`, `conv_to_shell`), the factor (`f`, `f_sd`, `f_min`, `f_max`) and
# why it is not to be used (`why`, NA otherwise). The values of `pool` (by
# default the single values) are averaged by source, study, taxon, name as
# given, locality, life stage and conversion, with the mean and SD of their
# logs; the other values (means) enter as given. `tag` starts the `sample` of
# each factor.
compilation_factors <- function(mo, tag, pool = mo$statistic == "single") {
  where <- paste0(ifelse(is.na(mo$page), "", paste0("p. ", mo$page, ", ")), mo$table_or_figure)
  given_text <- paste0(mo$value_as_given, " ", mo$unit_as_given,
                       ifelse(is.na(mo$se), "", paste0(" (SE ", mo$se, ")")),
                       ifelse(mo$measured_or_derived == "derived", "; derived from published means", ""))
  stage_text <- ifelse(is.na(mo$life_stage), "", paste0("; ", mo$life_stage))
  key <- paste(mo$source_id, mo$original_reference, mo$taxon, mo$name_as_given, mo$locality,
               mo$life_stage, mo$habitat, mo$from, mo$conv_from_shell, mo$conv_to, mo$conv_to_shell,
               mo$why, sep = "\u001f")
  agg <- lapply(split(seq_len(nrow(mo))[pool], key[pool], drop = TRUE), function(k) {
    r <- mo[k[1], ]
    f <- mo$f[k]
    one <- length(k) == 1L
    units <- unique(mo$n_unit[k])
    std(source_id = r$source_id, taxon = r$taxon, name_as_given = r$name_as_given,
        habitat = r$habitat, from = r$from, from_shell = r$conv_from_shell, to = r$conv_to,
        to_shell = r$conv_to_shell, factor_unit = if (r$conv_to == "energy") "J/ug" else "ratio",
        statistic = if (one) r$statistic else "mean", factor = mean(f),
        factor_sd = if (one) NA else stats::sd(f),
        factor_min = if (one && !is.na(r$f_min)) r$f_min else min(f),
        factor_max = if (one && !is.na(r$f_max)) r$f_max else max(f),
        ln_mean = if (one) NA else mean(log(f)), ln_sd = if (length(k) > 2) stats::sd(log(f)) else NA,
        n = if (one) r$n else length(k),
        n_unit = if (one) r$n_unit else if (length(units) == 1L && !is.na(units) && all(mo$n[k] %in% 1)) units else "values",
        scope = if (is.na(r$locality)) "pooled" else "population", site = r$locality,
        original_reference = r$original_reference,
        sample = paste(tag, r$source_id, r$taxon, r$name_as_given, r$locality, r$life_stage, sep = ":"),
        use_in_model = is.na(r$why), flag = r$why,
        note = if (one) paste0(where[k], "; ", given_text[k], stage_text[k]) else
          paste0(first_last(where[k]), stage_text[k[1]], "; ", length(k), " values"))
  })
  means <- which(!pool)
  given <- if (length(means)) {
    std(source_id = mo$source_id[means], taxon = mo$taxon[means],
        name_as_given = mo$name_as_given[means], habitat = mo$habitat[means],
        from = mo$from[means], from_shell = mo$conv_from_shell[means], to = mo$conv_to[means],
        to_shell = mo$conv_to_shell[means],
        factor_unit = ifelse(mo$conv_to[means] == "energy", "J/ug", "ratio"),
        statistic = "mean", factor = mo$f[means], factor_sd = mo$f_sd[means],
        factor_min = mo$f_min[means], factor_max = mo$f_max[means], n = mo$n[means],
        n_unit = mo$n_unit[means],
        scope = ifelse(is.na(mo$locality[means]), "pooled", "population"),
        site = mo$locality[means], original_reference = mo$original_reference[means],
        sample = paste(tag, mo$source_id[means], mo$table_or_figure[means],
                       mo$name_as_given[means], mo$locality[means], sep = ":"),
        use_in_model = is.na(mo$why[means]), flag = mo$why[means],
        note = paste0(where[means], "; ", given_text[means], stage_text[means]))
  }
  rbind(do.call(rbind, agg), given)
}
pair_of <- function(x) paste0(x$from, ifelse(x$from_shell, " (shell)", ""), " -> ", x$to,
                              ifelse(x$to_shell, " (shell)", ""))

# 8. Hébert et al. (2016): carbon, nitrogen and phosphorus of crustacean zooplankton --
# Carbon, nitrogen and phosphorus as % of dry mass (columns Cct.prop, Nct.prop
# and Pct.prop of zooplankton_traits.csv, from version 0.6.0), with the
# reference of each value (Ref.C, Ref.N, Ref.P: codes of references.csv). The
# rows of a species, habitat and reference (replicates in the compilation)
# are averaged. Values of studies whose first author and year are among the
# references of the carbon, nitrogen or phosphorus records of the data bank of
# Brey et al. (2010) are kept but not used (they enter through the data bank),
# as are means outside the range of whole animals (above).

if (use_source("hebert2016", reg)) {
  hp <- source_files("hebert2016", reg)
  ht <- read.csv2(hp[1], stringsAsFactors = FALSE, na.strings = c("", "NA"))
  hr <- read.csv2(hp[2], stringsAsFactors = FALSE, fileEncoding = "latin1")
  first_author <- function(x) tolower(iconv(sub("[ ,].*$", "", trimws(x)), "UTF-8", "ASCII//TRANSLIT"))
  # First author and year of the references of the data bank (as "Andersen &
  # Hessen '91 LimnOc 36: 807-814").
  brey_keys <- character(0)
  if (use_source("brey2010", reg)) {
    tmp <- tempfile("brey")
    dir.create(tmp)
    utils::unzip(source_files("brey2010", reg)[1], files = "Conversion03.xls", exdir = tmp)
    raw_b <- as.data.frame(readxl::read_excel(file.path(tmp, "Conversion03.xls"), sheet = "Data",
                                              col_names = FALSE, col_types = "text",
                                              .name_repair = "minimal"))[-(1:5), ]
    elem <- !is.na(num(raw_b[[35]])) | !is.na(num(raw_b[[36]])) | !is.na(num(raw_b[[37]]))
    rb <- trimws(raw_b[[51]][elem])
    yy <- num(sub("^[^']*'([0-9]{2}).*$", "\\1", rb))
    brey_keys <- unique(paste(first_author(rb), ifelse(yy > 20, 1900 + yy, 2000 + yy)))
  }
  hr$key <- paste(first_author(hr$Authors), hr$Year)
  elements <- c(carbon = "Cct.prop", nitrogen = "Nct.prop", phosphorus = "Pct.prop")
  refs <- c(carbon = "Ref.C", nitrogen = "Ref.N", phosphorus = "Ref.P")
  rows <- list()
  for (el in names(elements)) {
    v <- num(ht[[elements[[el]]]])
    k <- which(!is.na(v) & v > 0)
    sp <- ifelse(is.na(ht$Species[k]) | ht$Species[k] %in% c("sp.", "spp."), "", ht$Species[k])
    x <- data.frame(taxon = trimws(paste(ht$Genus[k], sp)), group = ht$Group[k],
                    habitat = tolower(ht$Habitat[k]), value = v[k],
                    ref = gsub("\\s+", "", ht[[refs[[el]]]][k]), stringsAsFactors = FALSE)
    # Rows without habitat or reference cannot be placed or checked: left out.
    x <- x[!is.na(x$habitat) & !is.na(x$ref), ]
    for (g in split(x, paste(x$taxon, x$habitat, x$ref), drop = TRUE)) {
      codes <- strsplit(g$ref[1], ",", fixed = TRUE)[[1]]
      rr <- hr[match(codes, hr$Ref.code), ]
      cite <- paste(sprintf("%s (%s)", sub("[,.] .*$", "", trimws(rr$Authors)), rr$Year),
                    collapse = "; ")
      in_brey <- any(rr$key %in% brey_keys)
      vv <- g$value / 100
      b <- elemental_bounds[[el]]
      out_b <- mean(g$value) < b[1] || mean(g$value) > b[2]
      why_h <- c(if (in_brey) "original source in the data bank of Brey et al. (2010)",
                 if (out_b) paste0("outside the range of whole animals (", bounds_text, ")"))
      rows[[length(rows) + 1L]] <- std(
        source_id = "hebert2016", taxon = g$taxon[1], name_as_given = g$taxon[1],
        habitat = g$habitat[1], from = "dry_mass", from_shell = FALSE, to = el, to_shell = FALSE,
        factor_unit = "ratio", statistic = if (nrow(g) > 1) "mean" else "single",
        factor = mean(vv), factor_sd = if (nrow(g) > 1) stats::sd(vv) else NA,
        factor_min = min(vv), factor_max = max(vv), n = nrow(g), n_unit = "rows of the compilation",
        scope = "pooled", original_reference = cite,
        sample = paste("hebert2016", g$taxon[1], g$habitat[1], g$ref[1], sep = ":"),
        use_in_model = !length(why_h),
        flag = if (length(why_h)) paste(why_h, collapse = "; ") else NA,
        note = sprintf("%s; %s; references %s; %d row(s), mean %.3g %% of dry mass", g$group[1],
                       elements[[el]], g$ref[1], nrow(g), mean(g$value)))
    }
  }
  conv$hebert <- do.call(rbind, rows)
  hb <- conv$hebert
  say("## Hébert et al. (2016), carbon, nitrogen and phosphorus")
  say("")
  say("- %d factors (species, habitat and reference) of %d taxa: %s; %d freshwater.",
      nrow(hb), length(unique(hb$taxon)),
      paste(sprintf("%s %d", names(table(hb$to)), as.integer(table(hb$to))), collapse = ", "),
      sum(hb$habitat == "freshwater"))
  in_b <- grepl("Brey", hb$flag, fixed = TRUE)
  say("- Not used, original source in the data bank of Brey et al. (2010): %d (%s).",
      sum(in_b), paste(sort(unique(hb$original_reference[in_b])), collapse = "; "))
  out_h <- grepl("outside the range", hb$flag, fixed = TRUE)
  say("- Not used, outside the range of whole animals: %d (%s).", sum(out_h),
      paste(sprintf("%s %s %.3g", hb$taxon[out_h], hb$to[out_h], 100 * hb$factor[out_h]), collapse = "; "))
  say("")
}

# 9. Elemental composition of freshwater animals (compilation of 2026-11-09) ---------
# Carbon, nitrogen and phosphorus as % of dry mass of freshwater animals, read
# from 18 sources (data sets of Beck et al. 2022, González et al. 2018 and
# 2025, Knapp et al. 2021, Moody et al. 2018, May et al. 2026, Evangelista et
# al. 2024 and Musin et al. 2023; tables of articles and theses) into one table
# (Stoichiometry2026/ under BMASSA_DATABASES), with the file and line, or the
# page and table, of each value and its text as given, the name as given and
# the accepted name in the COL (26.10 XR), life stage, locality, preservation,
# gut clearance and methods. The sources that were already used here (Brey et
# al. 2010, Hébert et al. 2016, Evans-White et al. 2005 and the other sources
# of the conversion factors) were left out of the search. Each value is
# credited to the source it was read from (`source_id`), with the study that
# measured it in `original_reference` (Benstead et al. 2010 and Vanni et al.
# 2002 through the compilation of González et al. 2025, Singer and Battin 2007
# through that of Beck et al. 2022). The tadpoles (Knapp et al. 2021; Vanni
# et al. 2002 through González et al. 2025) are used, with `life_stage`
# "larva", for the models of tadpoles (all their stages, including the
# metamorphic climax). Kept but not used:
#   - adult amphibians (Xenopus laevis of Brenes-Soto et al. 2019, two
#     values): too few for a model;
#   - shells of molluscs measured alone (not the animal);
#   - studies that are in the data bank of Brey et al. (2010): Walve and
#     Larsson (1999), Liess and Hillebrand (2005), and the thesis of Ventura
#     (2005), whose Lake Redon zooplankton is that of Ventura and Catalan (2005);
#   - first-instar larvae reared from eggs, and snails processed for food;
#   - Zhang et al. (2017), until its values are checked against the article
#     (it could not be opened);
#   - values outside the range of whole animals (above; e.g. oligochaetes
#     with 11-18 % C, whose guts held sediment), and the other elements of
#     the same sample.
# Molluscs are taken without shell (soft tissue) and crustaceans whole, as in
# the other sources. Brackish values (a tidal river) are "transitional". The
# single values (individuals or pooled samples) are averaged by source, study,
# taxon, name as given, locality and life stage (as for the bivalves of
# Coughlan et al. 2021): the factors are means of populations; the means given
# by the sources enter as given (SD from the SE when n is known).

blocka_ids <- c(
  "10.5061/dryad.g4f4qrfr5" = "beck2022", "10.5061/dryad.3cd2611" = "gonzalez2018",
  "10.5061/dryad.3tx95x6r2" = "gonzalez2025", "10.6084/m9.figshare.13683910.v1" = "knapp2021",
  "10.5061/dryad.7f753f0" = "moody2018", "10.5061/dryad.6m905qg99" = "may2026",
  "10.6084/m9.figshare.23822535.v1" = "evangelista2024",
  "10.1371/journal.pone.0262972" = "musin2023", "10.1016/j.jenvman.2020.111889" = "goedkoop2021",
  "10.1093/plankt/21.12.2309" = "walve1999", "10.18307/2017.0617" = "zhang2017",
  "10.32792/utq/utjsci/v6i2.83" = "abuturab2017", "10.1073/pnas.2008223118" = "li2021",
  "10.1899/12-181.1" = "backking2013", "10.1016/j.anres.2017.12.007" = "ghosh2017",
  "10.7717/peerj.7365" = "brenessoto2019")
# Names without an accepted name in the COL: the family given by the source
# (corrected spelling), or the higher taxon of the morphospecies.
blocka_names <- c(
  "Anisoptera_sp1" = "Anisoptera", "Brachycera_sp1" = "Brachycera", "Eriopterini" = "Limoniidae",
  "Hexatomini" = "Limoniidae", "Limoniini" = "Limoniidae", "Pediciini" = "Limoniidae",
  "Hirudinea.sp3" = "Hirudinidae", "Hirudinea_sp1" = "Hirudinidae", "Hirudinea_sp2" = "Hirudinidae",
  "Pomacea canaliculata" = "Pomacea canaliculata", "Silo-Lithax" = "Goeridae",
  "Stenophylacini-Chaetopterygini" = "Limnephilidae", "Tentrepohlia_sp1" = "Tipulidae",
  "Tentrepohlia_sp2" = "Tipulidae")
blocka_sources <- unique(unname(c(blocka_ids, "ventura2005", "back2013")))
if (any(use_source(blocka_sources, reg))) {
  # One file for all the sources (registered with each of them).
  ba_id <- blocka_sources[use_source(blocka_sources, reg)][1]
  ba <- read.csv(source_files(ba_id, reg)[1], stringsAsFactors = FALSE, check.names = FALSE,
                 encoding = "UTF-8", na.strings = "")
  stopifnot(all(ba$from == "dry_mass"), all(ba$value_unit == "percent"),
            all(ba$to %in% names(elemental_bounds)), all(ba$measured_or_derived == "measured"),
            all(ba$digitized == "no"))
  ba$source_id <- unname(blocka_ids[ba$doi])
  ba$source_id[is.na(ba$doi) & startsWith(ba$citation, "Ventura (2005)")] <- "ventura2005"
  ba$source_id[is.na(ba$doi) & startsWith(ba$citation, "Back (2013)")] <- "back2013"
  if (anyNA(ba$source_id)) stop("Rows of an unknown source: ", paste(unique(ba$citation[is.na(ba$source_id)]), collapse = "; "))
  ba <- ba[use_source(ba$source_id, reg), ]
  unnamed <- is.na(ba$accepted_name)
  if (!all(ba$name_as_given[unnamed] %in% names(blocka_names))) {
    stop("Names without accepted name not in blocka_names: ",
         paste(unique(setdiff(ba$name_as_given[unnamed], names(blocka_names))), collapse = ", "))
  }
  ba$taxon <- ifelse(unnamed, unname(blocka_names[ba$name_as_given]), ba$accepted_name)
  ba$habitat <- ifelse(ba$habitat == "brackish", "transitional", ba$habitat)
  stopifnot(all(ba$habitat %in% c("freshwater", "transitional")))
  mollusc <- ba$group %in% c("Gastropods", "Bivalves")
  shell_only <- mollusc & ba$from_shell == "included"
  # The sample of each value: the line or table row it was read from.
  ba$sample_key <- paste(ba$source_id, ba$table_or_figure, ba$page, ba$verbatim, ba$name_as_given,
                         ba$life_stage, sep = "\u001f")
  lo <- vapply(ba$to, function(e) elemental_bounds[[e]][1], numeric(1))
  hi <- vapply(ba$to, function(e) elemental_bounds[[e]][2], numeric(1))
  out_range <- !shell_only & (ba$value < lo | ba$value > hi)
  same_sample <- !out_range & ba$sample_key %in% ba$sample_key[out_range]
  ref <- ba$original_reference
  why <- rep(NA_character_, nrow(ba))
  addw <- function(k, text) why[k] <<- ifelse(is.na(why[k]), text, paste(why[k], text, sep = "; "))
  larva <- ba$group == "Amphibians" & startsWith(ba$life_stage, "tadpole")
  addw(which(ba$group == "Amphibians" & !larva),
       "adult amphibians: two values of one species, too few for a model")
  addw(which(shell_only), "shell measured alone, not the animal")
  addw(which(startsWith(ref, "Walve & Larsson") | startsWith(ref, "Liess & Hillebrand") |
               ba$source_id == "ventura2005"),
       "original source in the data bank of Brey et al. (2010)")
  addw(which(grepl("first-instar larva reared from egg", ba$life_stage, fixed = TRUE)),
       "first-instar larvae reared from eggs")
  addw(which(ba$source_id == "ghosh2017"), "snails processed for food")
  addw(which(ba$source_id == "zhang2017"), "values not checked against the source")
  addw(which(out_range), paste0("outside the range of whole animals (", bounds_text, ")"))
  addw(which(same_sample), "another element of the same sample outside the range of whole animals")
  ba$why <- why
  # Single values averaged by source, study, taxon, name, locality and stage.
  single <- ba$statistic == "single"
  key <- paste(ba$source_id, ref, ba$taxon, ba$name_as_given, ba$locality, ba$life_stage, ba$habitat,
               ba$to, ba$why, sep = "\u001f")
  agg <- lapply(split(seq_len(nrow(ba))[single], key[single], drop = TRUE), function(k) {
    v <- ba$value[k] / 100
    r <- ba[k[1], ]
    units <- unique(ba$n_unit[k])
    std(source_id = r$source_id, taxon = r$taxon, name_as_given = r$name_as_given,
        habitat = r$habitat, life_stage = if (larva[k[1]]) "larva" else NA_character_,
        from = "dry_mass", from_shell = FALSE, to = r$to, to_shell = FALSE,
        factor_unit = "ratio", statistic = if (length(k) > 1) "mean" else "single",
        factor = mean(v), factor_sd = if (length(k) > 1) stats::sd(v) else NA,
        factor_min = min(v), factor_max = max(v), ln_mean = if (length(k) > 1) mean(log(v)) else NA,
        ln_sd = if (length(k) > 2) stats::sd(log(v)) else NA, n = length(k),
        n_unit = if (length(units) == 1L && !is.na(units)) units else "samples",
        scope = if (is.na(r$locality)) "pooled" else "population", site = r$locality,
        original_reference = r$original_reference,
        sample = paste("stoich2026", r$source_id, r$taxon, r$name_as_given, r$locality, r$life_stage,
                       sep = ":"),
        use_in_model = is.na(r$why), flag = r$why,
        note = paste0(first_last(paste0(ifelse(is.na(ba$page[k]), "", paste0("p. ", ba$page[k], ", ")),
                                        ba$table_or_figure[k])),
                      ifelse(is.na(r$life_stage), "", paste0("; ", r$life_stage)),
                      ifelse(is.na(r$preservation), "", paste0("; preservation: ", r$preservation)),
                      ifelse(is.na(r$gut), "", paste0("; gut: ", r$gut)),
                      "; ", length(k), " value(s)"))
  })
  means <- which(!single)
  sd_m <- ifelse(!is.na(ba$sd[means]), ba$sd[means],
                 ifelse(!is.na(ba$se[means]) & !is.na(ba$n[means]), ba$se[means] * sqrt(ba$n[means]), NA))
  given <- std(source_id = ba$source_id[means], taxon = ba$taxon[means],
               name_as_given = ba$name_as_given[means], habitat = ba$habitat[means],
               life_stage = ifelse(larva[means], "larva", NA_character_), from = "dry_mass",
               from_shell = FALSE, to = ba$to[means], to_shell = FALSE, factor_unit = "ratio",
               statistic = "mean", factor = ba$value[means] / 100, factor_sd = sd_m / 100,
               factor_min = ba$min[means] / 100, factor_max = ba$max[means] / 100, n = ba$n[means],
               n_unit = ba$n_unit[means], scope = ifelse(is.na(ba$locality[means]), "pooled", "population"),
               site = ba$locality[means], original_reference = ba$original_reference[means],
               sample = paste("stoich2026", ba$source_id[means], ba$table_or_figure[means],
                              ba$name_as_given[means], ba$life_stage[means], sep = ":"),
               use_in_model = is.na(ba$why[means]), flag = ba$why[means],
               note = paste0(ifelse(is.na(ba$page[means]), "", paste0("p. ", ba$page[means], ", ")),
                             ba$table_or_figure[means], "; ", ba$value_as_given[means], " ",
                             ba$unit_as_given[means],
                             ifelse(is.na(ba$se[means]), "", paste0(" (SE ", ba$se[means], ")")),
                             ifelse(is.na(ba$life_stage[means]), "", paste0("; ", ba$life_stage[means]))))
  conv$stoich2026 <- rbind(do.call(rbind, agg), given)
  st <- conv$stoich2026
  say("## Elemental composition of freshwater animals (compilation of 2026-11-09)")
  say("")
  say("- %d values (%s) of %d sources; %d single values averaged into %d factors, %d means as given.",
      nrow(ba), paste(sprintf("%s %d", names(table(ba$to)), as.integer(table(ba$to))), collapse = ", "),
      length(unique(ba$source_id)), sum(single), length(agg), length(means))
  for (w in c("adult amphibians: two values of one species, too few for a model",
              "shell measured alone, not the animal",
              "original source in the data bank of Brey et al. (2010)",
              "first-instar larvae reared from eggs", "snails processed for food",
              "values not checked against the source")) {
    say("- Values not used, %s: %d.", w, sum(grepl(w, ba$why, fixed = TRUE)))
  }
  say("- Values not used, outside the range of whole animals: %d; other elements of the same samples: %d.",
      sum(out_range), sum(same_sample))
  say("- Factors used: %d (%s); %d freshwater, %d transitional.", sum(st$use_in_model),
      paste(sprintf("%s %d", names(table(st$to[st$use_in_model])),
                    as.integer(table(st$to[st$use_in_model]))), collapse = ", "),
      sum(st$use_in_model & st$habitat == "freshwater"),
      sum(st$use_in_model & st$habitat == "transitional"))
  tl <- st$use_in_model & st$life_stage %in% "larva"
  say("- Of them, tadpoles: %d (%s) of %d species.", sum(tl),
      paste(sprintf("%s %d", names(table(st$to[tl])), as.integer(table(st$to[tl]))), collapse = ", "),
      length(unique(st$taxon[tl])))
  say("- Values outside the range of whole animals: %s.",
      paste(sprintf("%s %s %.1f %%", ba$name_as_given[out_range], ba$to[out_range], ba$value[out_range]),
            collapse = "; "))
  say("")
}

# 10. Freshwater molluscs (compilation of 2026-11-09) ------------------------------
# Mass ratios of the soft tissue (dry over wet mass, AFDM over dry mass), its
# carbon and its energy per dry mass and AFDM, and the share of the shell in
# the whole dry or wet mass, of freshwater bivalves and snails, read from 15
# sources into one table (Molluscs2026/ under BMASSA_DATABASES), with the page
# and table of each value and its text as given, the name as given and the
# accepted name in the COL (26.10 XR), locality, preservation, drying, ashing
# and method, and whether each ratio was measured or derived from published
# means. Rows that repeat values read elsewhere are not read: those of the
# workbook of Coughlan et al. (2021) (section 2), the carbon of Goedkoop et
# al. (2021) (section 9) and the bivalves of Alimov and Shadrin (1977), read
# from Alimov (1981), whose means are in the data bank of Brey et al. (2010)
# (section 4). A share s of the shell in the dry (wet) mass gives the dry
# (wet) mass without shell over that with shell, 1 - s. Kept but not used:
#   - the edible flesh of snails (Debnath et al. 2016, Baby et al. 2010), not
#     the whole soft body (the edible flesh of bivalves is their soft body);
#   - animals exposed to diuron or to high calcium in experiments (the
#     controls are used) or starved for six months (Doherty 1986);
#   - the share of the valves in the wet mass of the quagga mussels of Baer et
#     al. (2022), at odds with the smaller shares that follow from the masses
#     in their section 2.4 (and without a model of that conversion);
#   - carbon outside the range of whole animals (section 8).
# Single values are averaged as in section 9 (by source, study, taxon, name
# as given, locality, life stage and conversion); means enter as given (SD
# from the SE when n is known).

molluscs_file <- "Molluscs2026/bmassa_block_b_mollusc_factors.csv"
molluscs_ids <- c(
  "10.4081/jlimnol.2007.15" = "ravera2007", "10.1111/jfb.15043" = "baer2022",
  "10.1016/j.jenvman.2020.111889" = "goedkoop2021",
  "10.17582/journal.pjz/2018.50.6.2153.2158" = "sereflisan2018",
  "10.23818/limn.33.16" = "lopezdoval2014", "10.2478/cjf-2025-0014" = "mazuran2025",
  "10.56093/ft.v53i4.65545" = "debnath2016", "10.3329/jsr.v2i2.3362" = "baby2010",
  "10.2307/1468447" = "benke1999", "10.5061/dryad.j3tx95xfg" = "moreau2021")
# Sources without a DOI, by the start of their citation.
molluscs_cited <- c("Begam & Khan (2002)" = "begam2002", "Tanjung (2015)" = "tanjung2015",
                    "Doherty (1986)" = "doherty1986")
molluscs_sources <- unique(unname(c(molluscs_ids, molluscs_cited)))
if (any(use_source(molluscs_sources, reg))) {
  mo_id <- molluscs_sources[use_source(molluscs_sources, reg)][1]
  mo <- read.csv(source_files(mo_id, reg)[[molluscs_file]], stringsAsFactors = FALSE,
                 check.names = FALSE, encoding = "UTF-8", na.strings = "")
  stopifnot(all(mo$digitized == "no"), all(mo$habitat %in% c("freshwater", "brackish")),
            all(mo$value_unit %in% c("fraction", "percent", "J_per_mg")),
            all(mo$from_shell %in% c("included", "removed")),
            all(mo$group %in% c("Bivalves", "Gastropods")), !anyNA(mo$accepted_name))
  n_rows <- nrow(mo)
  not_read <- list(
    coughlan = mo$doi %in% "10.5061/dryad.tdz08kq0c",
    alimov = is.na(mo$doi) & startsWith(mo$citation, "Alimov (1981)"),
    goedkoop_carbon = mo$doi %in% "10.1016/j.jenvman.2020.111889" & mo$to == "carbon")
  mo <- mo[!Reduce(`|`, not_read), ]
  mo$source_id <- unname(molluscs_ids[mo$doi])
  for (k in names(molluscs_cited)) {
    mo$source_id[is.na(mo$doi) & startsWith(mo$citation, k)] <- molluscs_cited[[k]]
  }
  if (anyNA(mo$source_id)) {
    stop("Rows of an unknown source: ", paste(unique(mo$citation[is.na(mo$source_id)]), collapse = "; "))
  }
  mo <- mo[use_source(mo$source_id, reg), ]
  mo$taxon <- mo$accepted_name
  mo$habitat <- ifelse(mo$habitat == "brackish", "transitional", mo$habitat)
  # Conversion of each value: shares of the shell as mass without shell over
  # mass with shell.
  share <- mo$to %in% c("shell_dry_mass", "shell_wet_mass")
  stopifnot(all(mo$from_shell[share] == "included"), all(mo$from_shell[!share] == "removed"),
            all(mo$from[share] == sub("^shell_", "", mo$to[share])),
            all(mo$value_unit[share] == "fraction"))
  mo$conv_to <- ifelse(share, mo$from, mo$to)
  mo$conv_from_shell <- mo$from_shell == "included"
  scale <- unname(c(fraction = 1, percent = 0.01, J_per_mg = 1e-3)[mo$value_unit])
  mo$f <- ifelse(share, 1 - mo$value * scale, mo$value * scale)
  mo$f_sd <- scale * ifelse(!is.na(mo$sd), mo$sd,
                            ifelse(!is.na(mo$se) & !is.na(mo$n), mo$se * sqrt(mo$n), NA))
  mo$f_min <- ifelse(share, 1 - mo$max * scale, mo$min * scale)
  mo$f_max <- ifelse(share, 1 - mo$min * scale, mo$max * scale)
  is_ratio <- mo$conv_to != "energy"
  stopifnot(all(mo$f[is_ratio] > 0 & mo$f[is_ratio] <= 1), all(mo$f[!is_ratio] > 0))
  why <- rep(NA_character_, nrow(mo))
  addw <- function(k, text) why[k] <<- ifelse(is.na(why[k]), text, paste(why[k], text, sep = "; "))
  edible <- mo$group == "Gastropods" & grepl("edible", mo$notes, ignore.case = TRUE)
  addw(which(edible), "edible flesh of snails, not the whole soft body")
  addw(which(mo$source_id == "lopezdoval2014" & !grepl("Measured diuron 0 ", mo$notes, fixed = TRUE)),
       "exposed to diuron in an experiment (the control is used)")
  addw(which(mo$source_id == "mazuran2025" & grepl("Ca exposure", mo$notes, fixed = TRUE)),
       "exposed to high calcium in an experiment (the controls are used)")
  addw(which(mo$source_id == "doherty1986" & grepl("without food", mo$table_or_figure, fixed = TRUE)),
       "starved for six months in the laboratory")
  addw(which(mo$source_id == "baer2022" & mo$to == "shell_wet_mass"),
       "share of the valves in the wet mass at odds with the masses in section 2.4 of the paper")
  elem <- mo$conv_to %in% names(elemental_bounds)
  lo <- vapply(mo$conv_to, function(e) if (e %in% names(elemental_bounds)) elemental_bounds[[e]][1] else NA_real_, 0)
  hi <- vapply(mo$conv_to, function(e) if (e %in% names(elemental_bounds)) elemental_bounds[[e]][2] else NA_real_, 0)
  out_range <- elem & (100 * mo$f < lo | 100 * mo$f > hi)
  addw(which(out_range), paste0("outside the range of whole animals (", bounds_text, ")"))
  mo$why <- why
  mo$conv_to_shell <- FALSE
  single <- mo$statistic == "single"
  conv$molluscs2026 <- compilation_factors(mo, "molluscs2026")
  mc <- conv$molluscs2026
  say("## Freshwater molluscs (compilation of 2026-11-09)")
  say("")
  say("- %d rows; not read: %d of the workbook of Coughlan et al. (2021) (read from the workbook, section 2), %d of Alimov (1981) (Alimov and Shadrin 1977, whose means are in the data bank of Brey et al. 2010) and %d carbon values of Goedkoop et al. (2021) (section 9).",
      n_rows, sum(not_read$coughlan), sum(not_read$alimov), sum(not_read$goedkoop_carbon))
  say("- %d values of %d sources (%s); %d single values averaged into %d factors, %d means as given.",
      nrow(mo), length(unique(mo$source_id)),
      paste(sprintf("%s %d", names(table(mo$source_id)), as.integer(table(mo$source_id))), collapse = ", "),
      sum(single), nrow(mc) - sum(!single), sum(!single))
  for (w in sort(unique(unlist(strsplit(mo$why[!is.na(mo$why)], "; ", fixed = TRUE))))) {
    say("- Values not used, %s: %d.", w, sum(grepl(w, mo$why, fixed = TRUE)))
  }
  u_mc <- mc[mc$use_in_model, ]
  say("- Factors used: %d of %d taxa (%s).", nrow(u_mc), length(unique(u_mc$taxon)),
      paste(sprintf("%s %d", names(table(pair_of(u_mc))), as.integer(table(pair_of(u_mc)))), collapse = ", "))
  say("")
}

# 11. Shells and soft tissue of snails (compilation of 2026-11-09) --------------------
# The share of the shell in the dry mass of snails, the dry mass of the whole
# snail (shell included) over its wet mass and the dry mass of the soft tissue
# over the wet mass of the whole snail, from a second compilation of
# freshwater and brackish molluscs (Molluscs2026/ under BMASSA_DATABASES; the
# columns of section 10, with whether each mass includes the shell), and the
# two snails of Kołodziejczyk (1984, Table IV: mean dry mass of one individual
# with and without shell), transcribed with the same columns in
# data-raw/kolodziejczyk1984_table4.csv. Not read:
#   - the organic content of empty shells (Strayer and Malcom 2007) and the dry
#     shell over the live mass of Bithynia tentaculata (Kamler 1969), which are
#     not conversions between masses of the animal;
#   - the share of the shell in the wet mass of Lymnaea stagnalis (Grosell and
#     Brix 2009), a conversion without a model;
#   - the ratios of Sphaerium lacustre of Yan and Liang (2003), whose
#     calibration specimens and their handling are not described (those of the
#     survey were kept in formalin);
#   - the shares of the shell of Theodoxus macri and Melanopsis praemorsa of
#     Elkarmi and Ismail (2006a, b), from means of the dry shell and of the dry
#     soft body that disagree with the regressions printed in the same papers
#     and leave 39-45 % of the dry mass to the soft body (12-25 % in the other
#     snails of this section).
# A share s of the shell gives the dry mass without shell over that with
# shell, 1 - s. Single values are averaged and means enter as given, as in
# section 10: the 60 snails weighed one by one by Crichton (2003) give one
# factor with the mean and SD of their log ratios. Brackish sites are
# transitional.

shells_file <- "Molluscs2026/bmassa_block_b2_mollusc_shells.csv"
shells_ids <- c("10.1016/S0025-326X(03)00241-8" = "sokolowski2003",
                "10.1080/00288330.1998.9516826" = "james1998")
shells_cited <- c("Crichton (2003)" = "crichton2003", "Kamler (1969)" = "kamler1969",
                  "Kołodziejczyk (1984)" = "kolodziejczyk1984")
shells_sources <- unique(unname(c(shells_ids, shells_cited)))
if (any(use_source(shells_sources, reg))) {
  read_tab <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                                      encoding = "UTF-8", na.strings = "")
  sh2 <- read_tab("data-raw/kolodziejczyk1984_table4.csv")
  n_kolo <- nrow(sh2)
  in_file <- setdiff(shells_sources, "kolodziejczyk1984")
  in_file <- in_file[use_source(in_file, reg)]
  if (length(in_file)) {
    b2 <- read_tab(source_files(in_file[1], reg)[[shells_file]])
    stopifnot(identical(names(b2), names(sh2)))
    sh2 <- rbind(b2, sh2)
  }
  stopifnot(all(sh2$digitized == "no"), all(sh2$habitat %in% c("freshwater", "brackish")),
            all(sh2$value_unit == "fraction"), all(sh2$group %in% c("Bivalves", "Gastropods")),
            !anyNA(sh2$accepted_name))
  n_rows <- nrow(sh2)
  shell_alone <- c("shell_wet_mass", "shell_dry_mass", "shell_afdm")
  not_read <- list(
    shell_masses = sh2$from %in% shell_alone |
      (sh2$to %in% shell_alone & sh2$from != sub("^shell_", "", sh2$to)),
    wet_share = sh2$from == "wet_mass" & sh2$to == "shell_wet_mass",
    yan = startsWith(sh2$citation, "Yan & Liang (2003)"),
    elkarmi = startsWith(sh2$citation, "Elkarmi & Ismail (2006)"))
  stopifnot(all(Reduce(`+`, not_read) <= 1))
  sh2 <- sh2[!Reduce(`|`, not_read), ]
  sh2$source_id <- unname(shells_ids[sh2$doi])
  for (k in names(shells_cited)) {
    sh2$source_id[is.na(sh2$doi) & startsWith(sh2$citation, k)] <- shells_cited[[k]]
  }
  if (anyNA(sh2$source_id)) {
    stop("Rows of an unknown source: ", paste(unique(sh2$citation[is.na(sh2$source_id)]), collapse = "; "))
  }
  sh2 <- sh2[use_source(sh2$source_id, reg), ]
  sh2$taxon <- sh2$accepted_name
  sh2$habitat <- ifelse(sh2$habitat == "brackish", "transitional", sh2$habitat)
  share <- sh2$to == "shell_dry_mass"
  stopifnot(all(sh2$from_shell == "included"), all(sh2$from[share] == "dry_mass"),
            all(sh2$to_shell[!share] %in% c("included", "removed")))
  sh2$conv_to <- ifelse(share, "dry_mass", sh2$to)
  sh2$conv_from_shell <- TRUE
  sh2$conv_to_shell <- !share & sh2$to_shell == "included"
  sh2$f <- ifelse(share, 1 - sh2$value, sh2$value)
  sh2$f_sd <- ifelse(!is.na(sh2$sd), sh2$sd, ifelse(!is.na(sh2$se) & !is.na(sh2$n), sh2$se * sqrt(sh2$n), NA))
  sh2$f_min <- ifelse(share, 1 - sh2$max, sh2$min)
  sh2$f_max <- ifelse(share, 1 - sh2$min, sh2$max)
  stopifnot(all(sh2$f > 0 & sh2$f < 1))
  sh2$why <- NA_character_
  conv$shells2026 <- compilation_factors(sh2, "shells2026")
  sc <- conv$shells2026
  say("## Shells and soft tissue of snails (compilation of 2026-11-09)")
  say("")
  say("- %d rows (%d of Kołodziejczyk 1984, Table IV); not read: %d of masses of the shell alone, %d shares of the shell in the wet mass, %d of Yan and Liang (2003) and %d of Elkarmi and Ismail (2006a, b).",
      n_rows, n_kolo, sum(not_read$shell_masses), sum(not_read$wet_share), sum(not_read$yan),
      sum(not_read$elkarmi))
  say("- %d values of %d sources (%s) give %d factors: %s.", nrow(sh2), length(unique(sh2$source_id)),
      paste(sprintf("%s %d", names(table(sh2$source_id)), as.integer(table(sh2$source_id))), collapse = ", "),
      nrow(sc), paste(sprintf("*%s* %s %.3f", sc$taxon, pair_of(sc), sc$factor), collapse = ", "))
  say("")
}

# 12. Freshwater zooplankton (compilation of 2026-11-10) -----------------------------
# Dry over wet mass, AFDM over dry mass, carbon and energy per dry mass of
# freshwater rotifers, copepods, cladocerans and ostracods, read from 11
# sources into one table (Zooplankton2026/ under BMASSA_DATABASES), with the
# columns of section 10 and how the wet mass was weighed. Not read:
#   - values of studies in the data bank of Brey et al. (2010) (section 4):
#     Salonen et al. (1976), Schindler et al. (1971) and Comita and Schindler
#     (1963), the last two through Moreau et al. (2021), and the thesis of
#     Ventura (2005), whose Lake Redon zooplankton is that of Ventura and
#     Catalan (2005);
#   - the carbon of González et al. (2018), read with the elemental
#     composition (section 9);
#   - Işık et al. (1999), excluded in version 0.2.0: their rotifers were
#     weighed with the algae filtered from the cultures;
#   - energy densities through Moreau et al. (2021) of Sitaramaiah (1967),
#     calculated from nitrogen (and in Cummins and Wuycheck 1971, section 6),
#     and of Anderson and Smith (1998), whose ostracods may or may not include
#     their valves;
#   - mixed zooplankton without a taxon (Karpowicz et al. 2020).
# Kept but not used:
#   - the energy densities of Riccardi and Mangoni (1999), calculated from the
#     elemental composition (equations of Gnaiger and Bitterlich 1984), not
#     measured, as in section 6;
#   - the moisture of Ovie and Ovie (2006), of zooplankton air-dried on filter
#     paper and stored before it was weighed, not of live animals.
# The five feeds of the cultures of Moina macrocopa of Hasan et al. (2023), a
# culture without a control feed, are averaged into one factor per
# conversion. The other values are combined as in section 10 (single values
# averaged by source, taxon, locality, life stage and conversion; means as
# given). Subgenera are dropped from the names.

zoo_file <- "Zooplankton2026/bmassa_block_c_zooplankton_factors.csv"
zoo_ids <- c("10.4081/jlimnol.1999.58" = "riccardi1999", "10.5962/bhl.title.31524" = "birge1922",
             "10.2331/suisan.44.979" = "watanabe1978", "10.22271/fish.2023.v11.i2b.2798" = "hasan2023",
             "10.46989/001c.20428" = "ovie2006")
zoo_sources <- unname(zoo_ids)
if (any(use_source(zoo_sources, reg))) {
  zo_id <- zoo_sources[use_source(zoo_sources, reg)][1]
  zo <- read.csv(source_files(zo_id, reg)[[zoo_file]], stringsAsFactors = FALSE, check.names = FALSE,
                 encoding = "UTF-8", na.strings = "")
  stopifnot(all(zo$digitized == "no"), all(zo$habitat == "freshwater"),
            all(zo$value_unit %in% c("fraction", "percent", "J_per_mg")))
  n_rows <- nrow(zo)
  from_moreau <- zo$doi %in% "10.5061/dryad.j3tx95xfg"
  ref_starts <- function(x) !is.na(zo$original_reference) & startsWith(zo$original_reference, x)
  not_read <- list(
    brey = zo$doi %in% "10.4319/lo.1976.21.5.0724" | startsWith(zo$citation, "Ventura (2005)") |
      (from_moreau & (ref_starts("Schindler") | ref_starts("Comita"))),
    gonzalez = zo$doi %in% "10.5061/dryad.3cd2611",
    isik = zo$doi %in% "10.1016/S0044-8486(99)00013-7",
    energy = from_moreau & (ref_starts("Sitaramaiah") | ref_starts("Anderson")),
    mixed = zo$doi %in% "10.1007/s10750-020-04271-5")
  stopifnot(all(Reduce(`+`, not_read) <= 1), all(!from_moreau | Reduce(`|`, not_read)))
  zo <- zo[!Reduce(`|`, not_read), ]
  zo$source_id <- unname(zoo_ids[zo$doi])
  if (anyNA(zo$source_id)) {
    stop("Rows of an unknown source: ", paste(unique(zo$citation[is.na(zo$source_id)]), collapse = "; "))
  }
  zo <- zo[use_source(zo$source_id, reg), ]
  stopifnot(!anyNA(zo$accepted_name), all(zo$from_shell %in% c("not_applicable", NA)))
  zo$taxon <- gsub(" \\([A-Z][a-z]+\\)", "", zo$accepted_name)
  zo$conv_to <- zo$to
  zo$conv_from_shell <- FALSE
  zo$conv_to_shell <- FALSE
  scale <- unname(c(fraction = 1, percent = 0.01, J_per_mg = 1e-3)[zo$value_unit])
  zo$f <- zo$value * scale
  zo$f_sd <- scale * ifelse(!is.na(zo$sd), zo$sd, ifelse(!is.na(zo$se) & !is.na(zo$n), zo$se * sqrt(zo$n), NA))
  zo$f_min <- zo$min * scale
  zo$f_max <- zo$max * scale
  is_ratio <- zo$conv_to != "energy"
  stopifnot(all(zo$f[is_ratio] > 0 & zo$f[is_ratio] <= 1), all(zo$f[!is_ratio] > 0))
  why <- rep(NA_character_, nrow(zo))
  why[zo$source_id == "riccardi1999" & zo$conv_to == "energy"] <-
    "energy calculated from organic composition, not measured"
  why[zo$source_id == "ovie2006"] <- "wet mass after air-drying on filter paper and storage, not of live animals"
  zo$why <- why
  feeds <- zo$source_id == "hasan2023"
  conv$zooplankton2026 <- compilation_factors(zo, "zooplankton2026", pool = zo$statistic == "single" | feeds)
  zc <- conv$zooplankton2026
  say("## Freshwater zooplankton (compilation of 2026-11-10)")
  say("")
  say("- %d rows; not read: %d of studies in the data bank of Brey et al. (2010), %d carbon values of González et al. (2018) (section 9), %d of Işık et al. (1999) (excluded in version 0.2.0), %d energy densities of Sitaramaiah (1967) and Anderson and Smith (1998) through Moreau et al. (2021) and %d of mixed zooplankton.",
      n_rows, sum(not_read$brey), sum(not_read$gonzalez), sum(not_read$isik), sum(not_read$energy),
      sum(not_read$mixed))
  say("- %d values of %d sources (%s) give %d factors (the five feeds of Hasan et al. 2023 averaged into one per conversion).",
      nrow(zo), length(unique(zo$source_id)),
      paste(sprintf("%s %d", names(table(zo$source_id)), as.integer(table(zo$source_id))), collapse = ", "),
      nrow(zc))
  for (w in sort(unique(zo$why[!is.na(zo$why)]))) {
    say("- Values not used, %s: %d.", w, sum(zo$why %in% w))
  }
  u_zc <- zc[zc$use_in_model, ]
  say("- Factors used: %d of %d taxa (%s).", nrow(u_zc), length(unique(u_zc$taxon)),
      paste(sprintf("%s %d", names(table(pair_of(u_zc))), as.integer(table(pair_of(u_zc)))), collapse = ", "))
  say("")
}

# 13. Limno-STOICH (Corman et al. 2026) ------------------------------------------
# Carbon, nitrogen and phosphorus as % of dry mass of the animals of the
# Limno-STOICH database (release of 2025-12-11, package edi.2231.2 of the
# Environmental Data Initiative; LimnoSTOICH2026/ under BMASSA_DATABASES):
# organisms of inland waters entered from more than 190 data sources (data
# sets of their authors, some unpublished, and articles) and checked by the
# project, which entered only whole organisms (molluscs without shell count as
# whole) and replaced by 999999999 the carbon below 20 %, the nitrogen below
# 4 % and the phosphorus below detection (Corman et al. 2026; here NA). The
# values are credited to the database (`source_id` corman2026), with their
# data source in `original_reference`: first author, year and DOI (or title),
# or the contributor of unpublished data. Decisions of the maintainer of
# 2026-11-10. Not read:
#   - data sources used already: Beck et al. (2022), Moody et al. (2018) and
#     Knapp et al. (2021) (section 9), and González et al. (2017), whose
#     invertebrates of bromeliads are those of González et al. (2018);
#   - values that are not % of dry mass (masses per individual or per litre);
#   - parts of animals, also when entered as whole organisms (the feet of
#     snails and mussels, the tail of crayfish, fins, liver, eggs of birds);
#   - other vertebrates (a grebe, a crocodile) and records without a genus,
#     family or order;
#   - eggs, embryos and parasites (metacercariae of trematodes).
# Kept but not used:
#   - fish analysed as muscle: those of Kiffney et al. (2002; dorsal muscle,
#     p. 14 of their report), the trout of Volk and Kiffney (2012), with 1.0 %
#     of phosphorus as in muscle (the whole fish of the other sources have 2.9
#     to 4.5 %), and those of the unpublished data of C. Volk, without
#     phosphorus, of the same species and year, with the same carbon and
#     nitrogen;
#   - amphibians other than the tadpoles of an unpublished thesis on ponds of
#     Alberta: juveniles emerging after metamorphosis (Fritz and Whiles 2018),
#     adult salamanders and juveniles of tailed frogs (C. Volk), juvenile frogs
#     of Trinidad without genus, and juveniles of the thesis without their
#     Gosner stage;
#   - animals that are not aquatic: samples of the riparian zone, spiders,
#     grasshoppers, ants and wasps, and the adults of mayflies, stoneflies,
#     caddisflies, alderflies, dragonflies and flies, whose larvae are the
#     aquatic stage;
#   - molluscs with less than 33 % carbon and 6 % nitrogen (or without
#     nitrogen): their soft tissue has 7 to 12 % nitrogen and 38 to 47 % carbon
#     in most data sources, and so much mineral matter points to the shell (the
#     snails of Watson et al. 2011, entered with and without shell, have 30 %
#     and 47 % carbon) or to sediment in the gut;
#   - values outside the range of whole animals (section 9), and the other
#     elements of the same record.
# The tadpoles of the thesis (Gosner 25-42, entered as juveniles and
# metamorphs) have `life_stage` "larva": models of tadpoles, all their stages
# as in section 9. Estuaries, the Odense Fjord, the lagoons of the Beaufort Sea
# and a salt marsh of Virginia are "transitional"; the floodplain lagoons of the
# Mitchell River (Jardine et al. 2011) and the rivers of the West Coast of New
# Zealand (B. Crichton), entered as coastal, are freshwater. The single records
# are averaged by data source, taxon, name, site, life stage (the Gosner stage
# of tadpoles) and habitat, as in section 9, and so are the means of several
# samples entered as aggregates, which can repeat a population (sampling
# dates): a factor is the composition of a population.

ls_file <- "LimnoSTOICH2026/STOICH_Release_Full_Join_20251211_shortened_notes.csv"
if (use_source("corman2026", reg)) {
  oc <- function(x) paste0(x, ".OrganismStoichiometry")
  lim <- read.csv(source_files("corman2026", reg)[[ls_file]], stringsAsFactors = FALSE, check.names = FALSE,
                 encoding = "UTF-8", na.strings = c("", "NA"), colClasses = "character")
  animal_types <- c("Insecta", "Crustacean", "Mollusca", "Annelida", "Zooplankton", "Invertebrates (Other)",
                    "Other", "Fish", "Amphibian", "Vertebrates (Other)")
  lim <- lim[lim[[oc("Type")]] %in% animal_types, ]
  stopifnot(!anyDuplicated(lim[[oc("Id")]]))
  n_rows <- nrow(lim)
  n_src <- length(unique(lim$Id.Source))
  type <- lim[[oc("Type")]]
  src <- lim$Id.Source
  lower <- function(x) tolower(trimws(x))
  doi_of <- function(x) {
    m <- regexpr("10\\.[0-9]{4,9}/[^ ]+", x)
    out <- rep(NA_character_, length(x))
    out[!is.na(m) & m > 0] <- regmatches(x, m)
    out
  }
  doi <- doi_of(lim$DOI.Source)
  # Identities of the data sources named in the rules (checked, as their
  # numbers could change in another release).
  who <- function(id) unique(c(lim$FirstAuthor.Source[src == id], lim$LastName.Contact[src == id]))
  stopifnot("Kiffney, Peter" %in% who("125"), "Volk, Carol" %in% who("139"), "Volk" %in% who("140"),
            "Whiting" %in% who("74"), "Jardine, Timothy" %in% who("174"), "Crichton, Ben" %in% who("188"),
            "Dunton, Kenneth" %in% who("150"), "Johnson, David" %in% who("94"))
  # Values: % of dry mass only; 999999999 (left out by the project) and values
  # not above 0 are missing.
  pct <- function(e, what = "Mean") {
    v <- num(lim[[oc(paste0(e, what))]])
    v[!(lim[[oc(paste0(e, "Units"))]] %in% "percent_dm") | !is.finite(v) | v <= 0 | abs(v) >= 999999999] <- NA
    v
  }
  elements <- c(carbon = "Carbon", nitrogen = "Nitrogen", phosphorus = "Phosphorus")
  vals <- sapply(elements, pct)
  sds <- sapply(elements, pct, what = "Stdev")
  # Taxon: genus and species, or the genus, family or order (classes of
  # molluscs entered as orders); names as given are in lower case.
  low <- function(x) {
    x <- lower(lim[[oc(x)]])
    x[!is.na(x) & !nzchar(x)] <- NA
    x
  }
  gen <- low("TaxonomicGenus")
  spp <- low("TaxonomicSpecies")
  fam <- low("TaxonomicFamily")
  ord <- low("TaxonomicOrder")
  spp[spp %in% c("sp", "sp.", "spp", "spp.")] <- NA
  mollusc_class <- c(class_bivalva = "Bivalvia", class_gastropoda = "Gastropoda")
  stopifnot(all(ord[startsWith(ord, "class_") %in% TRUE] %in% names(mollusc_class)))
  cls <- unname(mollusc_class[ord])
  ord[startsWith(ord, "class_") %in% TRUE] <- NA
  up <- function(x) ifelse(is.na(x), NA, paste0(toupper(substr(x, 1, 1)), substring(x, 2)))
  lim$taxon <- ifelse(!is.na(gen) & !is.na(spp), paste(up(gen), spp), ifelse(!is.na(gen), up(gen),
                     ifelse(!is.na(fam), up(fam), ifelse(!is.na(ord), up(ord), cls))))
  lim$name_as_given <- ifelse(!is.na(gen) & !is.na(spp), paste(gen, spp),
                             ifelse(!is.na(gen), gen, ifelse(!is.na(fam), fam, low("TaxonomicOrder"))))
  stage <- lower(lim[[oc("DevelopmentStage")]])
  notes <- lower(lim[[oc("Notes")]])
  tissue <- lower(lim[[oc("TissueSampled")]])
  sampled <- lim[[oc("OrganismSampled")]]
  order_l <- low("TaxonomicOrder")
  habitat_ev <- lim$Habitat.SampleEvent
  # Not read.
  dup_dois <- c("10.5061/dryad.7f753f0", "10.1002/ecs2.3466", "10.1111/geb.13413", "10.3389/fevo.2017.00110")
  soft <- type == "Mollusca" & grepl("shell removed|removed shell|without shell|no shell", tissue)
  not_read <- list(
    dup = doi %in% dup_dois,
    unit = rowSums(!is.na(vals)) == 0,
    part = !((sampled %in% "whole_organism") | soft) |
      grepl("foot|muscle|tail|meat|fin|liver|albumen|yolk|needle", tissue),
    vertebrate = type == "Vertebrates (Other)",
    no_taxon = is.na(lim$taxon),
    egg = stage %in% c("egg", "eggs", "embryos", "metacercariae"))
  read <- !Reduce(`|`, not_read)
  # Kept but not used.
  why <- rep(NA_character_, nrow(lim))
  addw <- function(k, text) why[k] <<- ifelse(is.na(why[k]), text, paste(why[k], text, sep = "; "))
  gosner <- as.integer(sub("^.*gosner ([0-9]+).*$", "\\1", ifelse(grepl("gosner [0-9]+", notes), notes, NA)))
  # Tadpoles: records of the thesis with their Gosner stage (those without it,
  # entered as juveniles, are of an uncertain stage).
  tadpole <- type == "Amphibian" & src == "74" & stage %in% c("juvenile", "metamorph") & !is.na(gosner)
  stopifnot(all(gosner[tadpole] >= 25))
  addw(which(type == "Fish" & src %in% c("125", "139", "140")), "fish analysed as muscle")
  addw(which(type == "Amphibian" & !tadpole),
       "amphibians other than tadpoles, or of an uncertain stage")
  not_aquatic <- habitat_ev %in% "Riparian" | order_l %in% c("araneae", "orthoptera", "hymenoptera") |
    (stage %in% c("adult", "adults") &
       order_l %in% c("ephemeroptera", "plecoptera", "trichoptera", "megaloptera", "odonata", "diptera"))
  addw(which(not_aquatic), "not aquatic (riparian samples, terrestrial orders, emerged adults of insects)")
  mineral <- type == "Mollusca" & !is.na(vals[, "carbon"]) & vals[, "carbon"] < 33 &
    (is.na(vals[, "nitrogen"]) | vals[, "nitrogen"] < 6)
  addw(which(mineral), "molluscs with the carbon and nitrogen of a shell or of other mineral matter")
  lo <- sapply(names(elements), function(e) elemental_bounds[[e]][1])
  hi <- sapply(names(elements), function(e) elemental_bounds[[e]][2])
  out_el <- sweep(vals, 2, lo, "<") | sweep(vals, 2, hi, ">")
  out_el[is.na(out_el)] <- FALSE
  out_rec <- rowSums(out_el) > 0
  # Habitat.
  coastal <- lim$EcosystemType.Site %in% "Coastal"
  lim$habitat <- ifelse((coastal & !src %in% c("174", "188")) | src %in% c("150", "94"),
                       "transitional", "freshwater")
  # Data source of each record.
  surname <- ifelse(grepl(",", lim$FirstAuthor.Source), sub(",.*$", "", lim$FirstAuthor.Source),
                    ifelse(grepl("_", lim$FirstAuthor.Source), sub("_.*$", "", lim$FirstAuthor.Source),
                           sub("^.* ", "", lim$FirstAuthor.Source)))
  surname[lim$FirstAuthor.Source %in% "Rana, El-Sabaawi"] <- "El-Sabaawi"
  surname[lim$FirstAuthor.Source %in% "Dörte, Köster"] <- "Köster"
  year <- sub("\\.0$", "", lim$YearPublished.Source)
  title <- up(lim$Title.Source)
  unpublished <- lim$Type.Source %in% "raw_unpublished" | is.na(lim$FirstAuthor.Source)
  contact <- trimws(paste(ifelse(is.na(lim$FirstName.Contact), "", lim$FirstName.Contact), lim$LastName.Contact))
  lim$original_reference <- ifelse(unpublished,
                                  paste0("Unpublished ", ifelse(grepl("thesis", lower(lim$Notes.Source)), "thesis", "data"),
                                         " of ", contact, " (Limno-STOICH source ", src, ")"),
                                  paste0(surname, " (", year, ")",
                                         ifelse(!is.na(doi), paste0(", ", doi), ifelse(is.na(title), "", paste0(", ", title))),
                                         " (Limno-STOICH source ", src, ")"))
  site_of <- function(...) {
    p <- cbind(...)
    apply(p, 1, function(z) paste(z[!is.na(z) & nzchar(z)], collapse = ", "))
  }
  lim$locality <- site_of(lim$Name.Site, lim$State.Site, lim$Country.Site)
  lim$stage_text <- ifelse(tadpole, paste0(ifelse(stage == "metamorph", "metamorph", "tadpole"), ", Gosner ", gosner),
                          stage)
  lim$why <- why
  # One row per record and element.
  long <- do.call(rbind, lapply(names(elements), function(e) {
    k <- which(read & !is.na(vals[, e]))
    data.frame(i = k, to = e, value = vals[k, e], sd = sds[k, e],
               why = ifelse(out_el[k, e], paste0("outside the range of whole animals (", bounds_text, ")"),
                            ifelse(out_rec[k], "another element of the same record outside the range of whole animals",
                                   NA_character_)),
               stringsAsFactors = FALSE)
  }))
  long$why <- ifelse(is.na(lim$why[long$i]), long$why,
                     ifelse(is.na(long$why), lim$why[long$i], paste(lim$why[long$i], long$why, sep = "; ")))
  r <- lim[long$i, ]
  aggregate <- r[[oc("StatisticalPooling")]] %in% "Aggregate"
  key <- paste(r$Id.Source, r$taxon, r$name_as_given, r$Id.Site, r$stage_text, r$habitat, long$to, long$why,
               sep = "\u001f")
  larva_of <- function(k) ifelse(tadpole[long$i[k]], "larva", NA_character_)
  gut <- r[[oc("GutClearance")]]
  agg <- lapply(split(seq_len(nrow(r)), key, drop = TRUE), function(k) {
    v <- long$value[k] / 100
    one <- length(k) == 1L
    std(source_id = "corman2026", taxon = r$taxon[k[1]], name_as_given = r$name_as_given[k[1]],
        habitat = r$habitat[k[1]], life_stage = larva_of(k[1]), from = "dry_mass", from_shell = FALSE,
        to = long$to[k[1]], to_shell = FALSE, factor_unit = "ratio",
        statistic = if (one && !aggregate[k]) "single" else "mean",
        factor = mean(v), factor_sd = if (one) long$sd[k] / 100 else stats::sd(v),
        factor_min = min(v), factor_max = max(v),
        ln_mean = if (one) NA else mean(log(v)), ln_sd = if (length(k) > 2) stats::sd(log(v)) else NA,
        n = if (one && aggregate[k]) num(r[[oc("SampleSize")]][k]) else length(k),
        n_unit = if (one && aggregate[k]) "samples" else "records", scope = "population",
        site = r$locality[k[1]], original_reference = r$original_reference[k[1]],
        sample = paste("limnostoich2026", r$Id.Source[k[1]], r$taxon[k[1]], r$name_as_given[k[1]],
                       r$Id.Site[k[1]], r$stage_text[k[1]], sep = ":"),
        use_in_model = is.na(long$why[k[1]]), flag = long$why[k[1]],
        note = paste0(if (one) "record " else "records ", first_last(r[[oc("Id")]][k]),
                      ifelse(is.na(r$stage_text[k[1]]), "", paste0("; ", r$stage_text[k[1]])),
                      ifelse(is.na(gut[k[1]]), "", paste0("; gut: ", gut[k[1]])), "; ",
                      if (one) "" else paste0(length(k), " records, "),
                      if (any(aggregate[k])) paste0(sum(aggregate[k]), " of them aggregates of samples")
                      else "single samples"))
  })
  conv$limnostoich2026 <- do.call(rbind, agg)
  lc <- conv$limnostoich2026
  say("## Limno-STOICH (Corman et al. 2026)")
  say("")
  say("- %d records of animals of %d data sources; not read: %d of data sources used already, %d without values in %% of dry mass, %d parts of animals, %d other vertebrates, %d without a genus, family or order, %d eggs, embryos and parasites (a record can fall in several).",
      n_rows, n_src, sum(not_read$dup), sum(not_read$unit), sum(not_read$part), sum(not_read$vertebrate),
      sum(not_read$no_taxon), sum(not_read$egg))
  say("- %d records read, of %d data sources: %d values (%s), %d of them means of aggregated samples; averaged by population into %d factors.",
      sum(read), length(unique(src[read])), nrow(long),
      paste(sprintf("%s %d", names(table(long$to)), as.integer(table(long$to))), collapse = ", "),
      sum(aggregate), length(agg))
  for (w in c("fish analysed as muscle", "amphibians other than tadpoles, or of an uncertain stage",
              "not aquatic (riparian samples, terrestrial orders, emerged adults of insects)",
              "molluscs with the carbon and nitrogen of a shell or of other mineral matter")) {
    say("- Records not used, %s: %d.", w, sum(grepl(w, lim$why[read], fixed = TRUE)))
  }
  say("- Values not used, outside the range of whole animals: %d; other elements of the same records: %d.",
      sum(out_el[read, ] & !is.na(vals[read, ])), sum(grepl("another element", long$why, fixed = TRUE)))
  say("- Factors used: %d (%s); %d freshwater, %d transitional; of %d data sources.", sum(lc$use_in_model),
      paste(sprintf("%s %d", names(table(lc$to[lc$use_in_model])), as.integer(table(lc$to[lc$use_in_model]))),
            collapse = ", "),
      sum(lc$use_in_model & lc$habitat == "freshwater"), sum(lc$use_in_model & lc$habitat == "transitional"),
      length(unique(lc$original_reference[lc$use_in_model])))
  tl <- lc$use_in_model & lc$life_stage %in% "larva"
  say("- Of them, tadpoles: %d (%s) of %d species.", sum(tl),
      paste(sprintf("%s %d", names(table(lc$to[tl])), as.integer(table(lc$to[tl]))), collapse = ", "),
      length(unique(lc$taxon[tl])))
  say("")
}

# 14. Remaining conversion factors of freshwater animals (compilation of 2026-11-10)
# Mass ratios, carbon and energy densities of groups with few factors (water
# mites, triclads, leeches, whole fish, tadpoles), read from seven articles
# into one table (Remaining2026/ under BMASSA_DATABASES), with the columns of
# section 10. Not read, as sources excluded in the registry:
#   - Morante et al. (2012), excluded since version 0.2.0 for the quality of
#     its masses (here: those of water mites near the error of the balance,
#     and two families of triclads with the same masses);
#   - Salonen et al. (1976), in the data bank of Brey et al. (2010);
#   - farmed animals fed formulated diets: the bass of Wang et al. (2021),
#     fattier than wild fish (48-54 % C of the dry mass, against 42-48 % in
#     the wild fish of Penczak et al. 1985), and the bullfrog tadpoles of Gao
#     et al. (2024) and Seixas Filho et al. (2022).
# Kept but not used: the triclads starved for two weeks (Calow and Jennings
# 1974) and the rainbow trout farmed in cages (Penczak et al. 1985). The
# energy per dry mass of Calow and Jennings (1974) is their energy per AFDM
# (micro-bomb calorimetry) times their AFDM/DM. The values are combined as in
# section 10 (single values averaged by source, taxon, locality, life stage
# and conversion; means as given).

rem_file <- "Remaining2026/bmassa_block_e_remaining_factors.csv"
rem_ids <- c("10.23818/limn.31.08" = "morante2012", "10.2307/1540570" = "calow1974",
             "10.4319/lo.1976.21.5.0724" = "salonen1976", "10.1155/2024/5513496" = "gao2024",
             "10.1016/j.aqrep.2021.100757" = "wang2021", "10.1590/S1519-99402122312022" = "seixas2022")
rem_sources <- c("calow1974", "penczak1985")
if (any(use_source(rem_sources, reg))) {
  re_id <- rem_sources[use_source(rem_sources, reg)][1]
  re <- read.csv(source_files(re_id, reg)[[rem_file]], stringsAsFactors = FALSE, check.names = FALSE,
                 encoding = "UTF-8", na.strings = "")
  stopifnot(all(re$digitized == "no"), all(re$habitat == "freshwater"),
            all(re$value_unit %in% c("fraction", "percent", "J_per_mg")),
            all(re$from_shell %in% "not_applicable"))
  n_rows <- nrow(re)
  re$source_id <- unname(rem_ids[re$doi])
  re$source_id[is.na(re$doi) & startsWith(re$citation, "Penczak et al. (1985)")] <- "penczak1985"
  if (anyNA(re$source_id)) {
    stop("Rows of an unknown source: ", paste(unique(re$citation[is.na(re$source_id)]), collapse = "; "))
  }
  read_it <- use_source(re$source_id, reg)
  not_read <- table(re$source_id[!read_it])
  re <- re[read_it, ]
  re$taxon <- re$accepted_name
  re$conv_to <- re$to
  re$conv_from_shell <- FALSE
  re$conv_to_shell <- FALSE
  scale <- unname(c(fraction = 1, percent = 0.01, J_per_mg = 1e-3)[re$value_unit])
  re$f <- re$value * scale
  re$f_sd <- scale * ifelse(!is.na(re$sd), re$sd, ifelse(!is.na(re$se) & !is.na(re$n), re$se * sqrt(re$n), NA))
  re$f_min <- re$min * scale
  re$f_max <- re$max * scale
  stopifnot(all(re$f[re$conv_to != "energy"] > 0 & re$f[re$conv_to != "energy"] <= 1))
  why <- rep(NA_character_, nrow(re))
  why[re$source_id == "calow1974" & grepl("starved", re$treatment)] <- "starved for two weeks before the analysis"
  why[re$source_id == "penczak1985" & re$treatment %in% "caged"] <- "farmed in cages, not a wild population"
  re$why <- why
  conv$remaining2026 <- compilation_factors(re, "remaining2026")
  rc <- conv$remaining2026
  say("## Remaining conversion factors of freshwater animals (compilation of 2026-11-10)")
  say("")
  say("- %d rows; not read (sources excluded in the registry): %s.", n_rows,
      paste(sprintf("%s %d", names(not_read), as.integer(not_read)), collapse = ", "))
  say("- %d values of %d sources (%s) give %d factors.", nrow(re), length(unique(re$source_id)),
      paste(sprintf("%s %d", names(table(re$source_id)), as.integer(table(re$source_id))), collapse = ", "),
      nrow(rc))
  for (w in sort(unique(re$why[!is.na(re$why)]))) {
    say("- Values not used, %s: %d.", w, sum(re$why %in% w))
  }
  u_rc <- rc[rc$use_in_model, ]
  say("- Factors used: %d of %d taxa (%s).", nrow(u_rc), length(unique(u_rc$taxon)),
      paste(sprintf("%s %d", names(table(pair_of(u_rc))), as.integer(table(pair_of(u_rc)))), collapse = ", "))
  say("")
}

# Output -------------------------------------------------------------------------

out <- do.call(rbind, conv)
rownames(out) <- NULL
out$permission <- reg$permission[match(out$source_id, reg$source_id)]
out <- out[order(out$source_id, out$taxon, out$from, out$to, out$scope, out$site,
                 method = "radix"), ]
out$conversion_id <- sprintf("%s_conv_%05d", out$source_id,
                             ave(seq_len(nrow(out)), out$source_id, FUN = seq_along))
out <- out[, c("conversion_id", "source_id", "permission", cols[-1])]
write.csv(out, file.path(out_dir, "conversions.csv"), row.names = FALSE, fileEncoding = "UTF-8")

say("## Output")
say("")
say("`conversions.csv`: %d factors of %d taxa from %d sources; %d to be used in the model.",
    nrow(out), length(unique(out$taxon)), length(unique(out$source_id)), sum(out$use_in_model))
lab <- function(cur, shell) paste0(cur, ifelse(shell, " (shell)", ""))
tab <- table(paste(lab(out$from, out$from_shell), "->", lab(out$to, out$to_shell))[out$use_in_model],
             out$source_id[out$use_in_model])
say("")
say("Factors used in the model, by conversion and source:")
say("")
say("| Conversion | %s |", paste(colnames(tab), collapse = " | "))
say("|---|%s", paste(rep("---|", ncol(tab)), collapse = ""))
for (i in seq_len(nrow(tab))) say("| %s | %s |", rownames(tab)[i], paste(tab[i, ], collapse = " | "))
writeLines(log_lines, file.path(out_dir, "conversions_log.md"), useBytes = TRUE)
cat(log_lines, sep = "\n")
