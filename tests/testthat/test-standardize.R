# standardize_taxa(), with recorded answers of the COL XR ----------------------------

skip_if_not_installed("jsonlite")
skip_if_not_installed("curl")

test_that("each kind of name gets its COL XR taxon (recorded answers)", {
  fx <- local_colxr_answers()
  expect_warning(s <- standardize_taxa(fx$names, quiet = TRUE), "2 of 16 names were not found",
                 class = "bmassa_warning_not_estimable")
  row <- function(name) s[s$scientific_name == name, ]
  # Exact species.
  b <- row("Baetis rhodani")
  expect_equal(b$match_type, "exact")
  expect_equal(unlist(b[c("family", "genus", "species", "col_id")], use.names = FALSE),
               c("Baetidae", "Baetis", "Baetis rhodani", "KD2J"))
  # A synonym gets the accepted name and classification.
  e <- row("Ephemerella ignita")
  expect_equal(c(e$col_name, e$col_status, e$genus), c("Serratella ignita", "synonym", "Serratella"))
  expect_match(e$note, "synonym of Serratella ignita")
  # An epithet with another ending is the same name.
  v <- row("Chironomus plumosa")
  expect_equal(c(v$match_type, v$species, v$col_id), c("variant", "Chironomus plumosus", "8BGDB"))
  # A misspelt epithet is not corrected: the taxon goes to its genus, as
  # placed in the COL XR (Daphnidae, not a homonymous record of Daphnia).
  d <- row("Daphnia magnaa")
  expect_equal(c(d$match_type, d$genus, d$family, d$col_id), c("higher_rank", "Daphnia", "Daphnidae", "3ZT6"))
  expect_true(is.na(d$species))
  expect_match(d$note, "closest name: Daphnia magna")
  # "sp." and "cf." go to the genus, not to the homonymous subgenus.
  for (nm in c("Chironomus sp.", "Chironomus cf. plumosus")) {
    g <- row(nm)
    expect_equal(c(g$col_rank, g$col_id, g$genus), c("genus", "BXWQR", "Chironomus"))
    expect_true(is.na(g$species))
  }
  expect_match(row("Chironomus cf. plumosus")$note, "cf./aff.")
  # A one-word name of a rank that is not a main rank, with a homonymous
  # placeholder genus in the COL XR.
  o <- row("Oligochaeta")
  expect_equal(c(o$col_rank, o$class), c("subclass", "Clitellata"))
  expect_true(is.na(o$order))
  expect_equal(row("Chironomidae")$col_rank, "family")
  expect_match(row("Gammarus pulex larvae")$note, "life-stage words removed")
  expect_equal(row("Gammarus pulex larvae")$species, "Gammarus pulex")
  # Phytoplankton outside Animalia (a cyanobacterium in Bacillati).
  expect_equal(row("Microcystis aeruginosa")$kingdom, "Bacillati")
  expect_equal(row("Cyclotella meneghiniana")$class, "Bacillariophyceae")
  # A rank given that differs from the COL XR is reported.
  expect_match(row("Squalius pyrenaicus")$note, "family given Leuciscidae, COL XR Cyprinidae")
  expect_equal(row("Pelophylax perezi")$class, "Amphibia")
  # Homonyms in other groups are reported.
  cl <- row("Cladocera")
  expect_equal(cl$col_name, "Cladoceromorpha")
  expect_match(cl$note, "other taxa with this name in the COL XR: Cladocera \\(genus, Insecta")
  # Not found: classification as given, and the closest name if any.
  h <- row("Hirudinea")
  expect_equal(h$match_type, "none")
  expect_true(all(is.na(unlist(h[c("kingdom", "col_id", "col_release")]))))
  expect_match(h$note, "closest name: Hirundinea \\(genus, not used\\)")
  expect_equal(row("Xyzzyqwerty foo")$match_type, "none")
  # Output contract.
  expect_equal(s$taxon_id, fx$names$taxon_id)
  expect_equal(s$scientific_name, fx$names$scientific_name)
  expect_true(all(s$col_release[s$match_type != "none"] == "COL26.6 XR"))
  expect_true(all(c("col_id", "col_name", "col_rank", "col_status", "match_type",
                    "match_confidence", "col_release", "note", "kingdom", "species") %in% names(s)))
  expect_identical(s, fx$expected)
})

test_that("its result goes into check_coverage() and estimate_biomass()", {
  fx <- local_colxr_answers()
  s <- suppressWarnings(standardize_taxa(fx$names, quiet = TRUE))
  ok <- s[s$match_type != "none", ]
  ok$length_mean <- ifelse(ok$taxon_id %in% c("x01", "x02", "x03"),
                           c(x01 = 6, x02 = 5, x03 = 12)[ok$taxon_id], NA)
  ok$length_unit <- ifelse(is.na(ok$length_mean), NA, "mm")
  ok$length_type <- ifelse(is.na(ok$length_mean), NA, "body_length")
  cov <- as.data.frame(check_coverage(ok))
  expect_true(all(cov$taxon_id %in% ok$taxon_id))
  by_length <- cov[cov$method == "length", ]
  # The synonym Ephemerella ignita is found as Serratella ignita.
  expect_equal(by_length$match_name, c("Baetis rhodani", "Serratella ignita", "Chironomus plumosus"))
  expect_true(all(by_length$match_level == "species"))
})

test_that("col_id finds the node of the reference even under another name", {
  tx <- data.frame(taxon_id = "a", broad_group = "macroinvertebrate", kingdom = "Animalia",
                   phylum = "Arthropoda", class = "Insecta", order = "Diptera",
                   family = "Chironomidae", genus = "Chironomus",
                   species = "Chironomus plumosa", col_id = "8BGDB", length_mean = 12,
                   length_unit = "mm", length_type = "body_length")
  cov <- as.data.frame(check_coverage(tx))
  cov <- cov[cov$method == "length", ]
  expect_equal(cov$match_level, "species")
  expect_equal(cov$match_name, "Chironomus plumosus")
  expect_match(cov$note, "matched by col_id")
  # Without col_id the misspelt species is not found, and the genus is used.
  tx$col_id <- NA_character_
  cov <- as.data.frame(check_coverage(tx))
  expect_equal(cov$match_level[cov$method == "length"], "genus")
})

test_that("names are parsed into the name sent and the rank it implies", {
  p <- parse_taxon_names(c("Daphnia (Ctenodaphnia) magna", "Asellus aquaticus (Linnaeus, 1758)",
                           "Baetis sp.", "Baetis spp", "Cyclops vicinus vicinus",
                           "Navicula radiosa var. tenella", "Tanypodinae", "Nitzschia"))
  expect_equal(p$query, c("Daphnia magna", "Asellus aquaticus", "Baetis", "Baetis", "Cyclops vicinus",
                          "Navicula radiosa", "Tanypodinae", "Nitzschia"))
  expect_equal(p$rank, c("SPECIES", "SPECIES", "GENUS", "GENUS", "SPECIES", "SPECIES", "SUBFAMILY",
                         "GENUS"))
  expect_match(p$note[2], "authorship removed")
  expect_match(p$note[5], "infraspecific")
  expect_true(same_taxon_name("Chironomus plumosa", "Chironomus plumosus"))
  expect_false(same_taxon_name("Daphnia magnaa", "Daphnia magna"))
  expect_false(same_taxon_name("Hirudinea", "Hirundinea"))
})

test_that("zero rows, input errors, missing packages and no connection", {
  z <- standardize_taxa(data.frame(taxon_id = character(), broad_group = character(),
                                   scientific_name = character()))
  expect_equal(nrow(z), 0L)
  expect_true(all(c("col_id", "match_type", "species") %in% names(z)))
  cls <- "bmassa_error_input"
  bad <- data.frame(taxon_id = c("a", "a"), broad_group = "zooplankton",
                    scientific_name = c("Daphnia magna", "Daphnia pulex"))
  expect_error(standardize_taxa(bad), class = cls)
  bad$taxon_id <- c("a", "b")
  bad$broad_group <- c("zooplankton", "plankton")
  expect_error(standardize_taxa(bad), class = cls)
  bad$broad_group <- "zooplankton"
  bad$scientific_name[2] <- NA
  expect_error(standardize_taxa(bad), class = cls)
  expect_error(standardize_taxa(bad[, 1:2]), class = cls)
  ok <- data.frame(taxon_id = "a", broad_group = "zooplankton", scientific_name = "Daphnia magna")
  local_mocked_bindings(is_installed = function(pkg) FALSE)
  expect_error(standardize_taxa(ok), class = "bmassa_error_missing_package")
})

test_that("a service that cannot be reached gives a network error", {
  expect_error(colxr_fetch("https://name-matching.invalid/", tries = 1L, wait = 0),
               class = "bmassa_error_network")
})

test_that("an answer that is not JSON is a network error", {
  skip_if_not_installed("jsonlite")
  expect_error(colxr_json("<html>Service unavailable</html>"), class = "bmassa_error_network")
})
