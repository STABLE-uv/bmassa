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
#   - Evans-White et al. (2005), Appendix: carbon content of stream
#     invertebrates (transcription in data-raw/evanswhite2005_appendix.csv)
#
# Run from the package root, in a UTF-8 locale:  Rscript data-raw/conversions.R
# Packages used only here: readxl.
#
# Currencies: live_mass (bivalves with the water of the mantle cavity),
# wet_mass, dry_mass, afdm, carbon, energy. `from_shell` and `to_shell` say
# whether that mass includes a shell or carbonate skeleton (FALSE also for
# animals without one), as `shell_included` in the other derived tables.
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
cols <- c("source_id", "taxon", "name_as_given", "habitat", "from", "from_shell", "to",
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
# masses, summarized by species (pooled) and by species and site (population;
# at least 3 individuals). Ratios above 1 between a mass and a part of it are
# measurement errors and are left out; SFDW/AFDW can exceed 1.

if (use_source("coughlan2021", reg)) {
  path <- source_files("coughlan2021", reg)
  species <- c("C. fluminea" = "Corbicula fluminea", "C. largillierti" = "Corbicula largillierti",
               "D. bugensis" = "Dreissena bugensis", "D. polymorpha" = "Dreissena polymorpha",
               "L. fortunei" = "Limnoperna fortunei", "S. Woodiana" = "Sinanodonta woodiana")
  currency <- list(LW = c("live_mass", TRUE), WW = c("wet_mass", TRUE), DW = c("dry_mass", TRUE),
                   SFDW = c("dry_mass", FALSE), AFDW = c("afdm", TRUE))
  pairs <- list(c("LW", "WW"), c("LW", "DW"), c("LW", "SFDW"), c("LW", "AFDW"), c("WW", "DW"),
                c("WW", "SFDW"), c("WW", "AFDW"), c("DW", "SFDW"), c("DW", "AFDW"),
                c("AFDW", "SFDW"))
  summ <- function(v, sp, site, scope, p) {
    std(source_id = "coughlan2021", taxon = sp, name_as_given = sp, habitat = "freshwater",
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
    for (p in pairs) {
      v <- num(x[[p[2]]]) / num(x[[p[1]]])
      ok <- is.finite(v) & v > 0
      bad <- ok & v > 1 & !identical(p, c("AFDW", "SFDW"))
      dropped <- dropped + sum(bad)
      ok <- ok & !bad
      if (sum(ok) < 3) next
      rows[[length(rows) + 1]] <- summ(v[ok], species[[sh]], NA_character_, "pooled", p)
      for (s in unique(x$Site[ok])) {
        k <- ok & x$Site == s
        if (sum(k) >= 3) rows[[length(rows) + 1]] <- summ(v[k], species[[sh]], s, "population", p)
      }
    }
  }
  conv$coughlan <- do.call(rbind, rows)
  pooled <- conv$coughlan[conv$coughlan$scope == "pooled", ]
  say("## Coughlan et al. (2021)")
  say("")
  say("- %d factors (%d by species, %d by species and site) for %d conversions; %d ratios above 1 left out.",
      nrow(conv$coughlan), nrow(pooled), sum(conv$coughlan$scope == "population"),
      length(pairs), dropped)
  for (nt in c("DW / WW of each individual", "SFDW / WW of each individual",
               "SFDW / DW of each individual")) {
    s <- pooled[pooled$note == nt, ]
    say("- %s: %s.", sub(" of each individual", "", nt),
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
# Total carbon (CHN analyser) as % of dry mass of 46 taxa of benthic
# invertebrates of streams of Indiana-Michigan and Wisconsin: mean, SD and
# range among streams (n = number of streams; in each stream one to three
# pooled samples), transcribed from the Appendix (pp. 1798-1799) into
# data-raw/evanswhite2005_appendix.csv, with `taxon` correcting two misprints
# (`edit_note`). The animals were kept overnight at 4 C to clear their guts,
# frozen, and freeze-dried or oven-dried at 45 C; insects and crustaceans were
# analysed whole and molluscs without shell (pp. 1788-1789).

ew_file <- "data-raw/evanswhite2005_appendix.csv"
if (use_source("evanswhite2005", reg) && file.exists(ew_file)) {
  ew <- read.csv(ew_file, stringsAsFactors = FALSE, encoding = "UTF-8", colClasses = "character",
                 na.strings = "")
  c_mean <- num(ew$c_mean)
  n_ew <- as.integer(ew$n)
  stopifnot(!anyNA(c_mean), !anyNA(n_ew), all(c_mean > 0 & c_mean < 100))
  sites <- c(IM = "Indiana-Michigan", W = "Wisconsin", "IM, W" = "Indiana-Michigan and Wisconsin")
  stopifnot(all(ew$region %in% names(sites)))
  conv$evanswhite <- std(
    source_id = "evanswhite2005", taxon = ew$taxon, name_as_given = ew$taxon_printed,
    habitat = "freshwater", from = "dry_mass", from_shell = FALSE, to = "carbon", to_shell = FALSE,
    factor_unit = "ratio", statistic = ifelse(n_ew > 1, "mean", "single"), factor = c_mean / 100,
    factor_sd = num(ew$c_sd) / 100, factor_min = num(ew$c_min) / 100, factor_max = num(ew$c_max) / 100,
    n = n_ew, n_unit = "streams", scope = "pooled", site = unname(sites[ew$region]),
    original_reference = "Evans-White et al. (2005)",
    sample = paste0("evanswhite2005:p", ew$page, "r", ew$row_on_page), use_in_model = TRUE,
    note = paste0("Appendix, p. ", ew$page, ", row ", ew$row_on_page, "; ", ew$c_mean, " % C",
                  ifelse(is.na(ew$c_sd), "", paste0(" (SD ", ew$c_sd, "; range ", ew$c_min, "-", ew$c_max, ")")),
                  "; ", ifelse(ew$group == "mollusc", "soft tissue, shell removed", "whole animals"),
                  "; total carbon"))
  say("## Evans-White et al. (2005), Appendix")
  say("")
  say("- %d taxa (%s), carbon as %% of dry mass; %d with one stream, the others means of 2-%d streams.",
      nrow(ew), paste(sprintf("%s %d", names(table(ew$group)), as.integer(table(ew$group))), collapse = ", "),
      sum(n_ew == 1), max(n_ew))
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
