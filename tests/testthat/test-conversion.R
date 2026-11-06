# convert_currency() for results of estimate_biomass() -------------------------------

ex <- example_survey
est <- suppressWarnings(estimate_biomass(ex$abundance, ex$taxa, ex$lengths, bias_correction = "lognormal"))
# A taxon unit without its length; the diatom t01 counted in cells and in
# natural units (two units).
without_length <- function(tx) {
  tx[c("length_mean", "length_unit", "length_type", "length_sd", "length_cv")] <- NA
  tx
}
cells_and_units <- rbind(ex$taxa[ex$taxa$taxon_id == "t01", ],
                         transform(ex$taxa[ex$taxa$taxon_id == "t01", ], taxon_id = "t01n",
                                   count_unit = "natural_unit"))
ref <- get_reference()
e0 <- as.data.frame(est)
row_of <- function(id, site = "L1") which(e0$taxon_id == id & e0$site == site)
tau2_of <- function(e, i) 2 * log(e$ind_mean[i] / e$ind_median[i])
comp <- function(m, l) {
  cc <- ref$conversion_components
  v <- cc$value[cc$model == m & cc$level == l & cc$parameter == "intercept"]
  if (length(v)) v else 0
}
# Deepest node of a conversion model for a taxon of `taxa`, by name from the
# species upwards (a simple lookup, independent of the matching code).
conv_node <- function(m, id) {
  tx <- ex$taxa[ex$taxa$taxon_id == id, ]
  n <- ref$conversion_nodes
  for (l in c("species", "genus", "family", "order", "class", "phylum")) {
    k <- which(n$model == m & n$level == l & n$node == tx[[l]])
    if (length(k) == 1L) return(k)
  }
  NA_integer_
}
missing_var <- function(m, level) {
  lv <- strsplit(ref$conversion_models$levels[ref$conversion_models$model == m], ",")[[1]]
  sum(vapply(lv[seq_along(lv) > match(level, lv)], function(l) comp(m, l), numeric(1)))
}

test_that("the conversion models have no effect of source", {
  expect_false(any(ref$conversion_components$level == "source"))
  expect_false(any(grepl("source", ref$conversion_models$context, fixed = TRUE)))
  expect_true(any(ref$conversion_models$context == "habitat"))
})

test_that("a factor conversion adds the path and the variances of its model (hand calculation)", {
  cv <- suppressWarnings(convert_currency(est, to = "carbon"))
  e1 <- as.data.frame(cv)
  m <- "invertebrate dry_mass_to_carbon"
  for (id in c("t10", "t13", "t06")) {
    i <- row_of(id)
    j <- conv_node(m, id)
    n <- ref$conversion_nodes[j, ]
    # Without habitat: its variance is added (the conversion models have no
    # effect of source).
    v_f <- n$path_var + missing_var(m, n$level) + comp(m, "residual") + comp(m, "habitat")
    expect_equal(log(e1$ind_median[i]) - log(e0$ind_median[i]), n$path_intercept, tolerance = 1e-12)
    expect_equal(tau2_of(e1, i) - tau2_of(e0, i), v_f, tolerance = 1e-10)
    expect_equal(e1$unit[i], "mg")
    expect_equal(e1$converted_from[i], "dry_mass")
    expect_match(e1$uncertainty[i], "conversion$")
  }
})

test_that("an inverse conversion changes the sign of the factor, not its variance", {
  wm <- suppressWarnings(convert_currency(est, to = "wet_mass"))
  e1 <- as.data.frame(wm)
  m <- "invertebrate wet_mass_to_dry_mass"
  i <- row_of("t10")
  n <- ref$conversion_nodes[conv_node(m, "t10"), ]
  v_f <- n$path_var + missing_var(m, n$level) + comp(m, "residual") + comp(m, "habitat")
  expect_equal(log(e1$ind_median[i]) - log(e0$ind_median[i]), -n$path_intercept, tolerance = 1e-12)
  expect_equal(tau2_of(e1, i) - tau2_of(e0, i), v_f, tolerance = 1e-10)
  expect_match(e1$conversion[i], "^dry_mass->wet_mass: inverse of invertebrate wet_mass_to_dry_mass")
})

test_that("fish go to carbon through dry mass, and the two factors add up", {
  cv <- suppressWarnings(convert_currency(est, to = "carbon", habitat = "freshwater"))
  e1 <- as.data.frame(cv)
  i <- row_of("t16")
  eta <- 0
  v <- 0
  for (m in c("fish wet_mass_to_dry_mass", "fish dry_mass_to_carbon")) {
    n <- ref$conversion_nodes[conv_node(m, "t16"), ]
    h <- ref$conversion_nodes[ref$conversion_nodes$model == m &
                                ref$conversion_nodes$level == "habitat" &
                                ref$conversion_nodes$node == "freshwater", ]
    eta <- eta + n$path_intercept + h$effect
    v <- v + n$path_var + missing_var(m, n$level) + comp(m, "residual") + h$effect_var
  }
  expect_equal(log(e1$ind_median[i]) - log(e0$ind_median[i]), eta, tolerance = 1e-12)
  expect_equal(tau2_of(e1, i) - tau2_of(e0, i), v, tolerance = 1e-10)
  expect_match(e1$conversion[i], "wet_mass->dry_mass.*; dry_mass->carbon")
})

test_that("phytoplankton cells follow Menden-Deuer and Lessard (2000), Table 4", {
  cv <- suppressWarnings(convert_currency(est, to = "carbon"))
  e1 <- as.data.frame(cv)
  # Cyclotella meneghiniana, a diatom: log10 C (pg) = -0.541 + 0.811 log10 V.
  i <- row_of("t01")
  v <- e0$ind_median[i]
  expect_equal(e1$ind_median[i], 10^(-0.541 + 0.811 * log10(v)), tolerance = 1e-12)
  expect_equal(e1$unit[i], "pg")
  # The width of the interval shrinks by b on the log scale.
  expect_equal(log(e1$ind_upper[i] / e1$ind_median[i]),
               0.811 * log(e0$ind_upper[i] / e0$ind_median[i]), tolerance = 1e-12)
  expect_match(e1$note[i], "not quantified")
  expect_false(grepl("conversion", e1$uncertainty[i]))
  # Hand calculation for 1000 um3: 77.98 pg (diatom), 141.9 pg (other protists).
  expect_equal(10^(-0.541 + 3 * 0.811), 77.98, tolerance = 1e-4)
  # Microcystis (cyanobacterium): protist equation, marked as an extrapolation.
  k <- row_of("t02")
  expect_equal(e1$ind_median[k], 10^(-0.665 + 0.939 * log10(e0$ind_median[k])), tolerance = 1e-12)
  expect_match(e1$conversion[k], "extrapolated: prokaryote")
})

test_that("phytoplankton wet mass is the biovolume at 1 g/cm3, for cells and natural units", {
  wm <- suppressWarnings(convert_currency(est, to = "wet_mass"))
  e1 <- as.data.frame(wm)
  phy <- which(e0$broad_group == "phytoplankton" & !is.na(e0$ind_median))
  # 1 um3 = 1 pg at the density of water (EN 16695, section 6.6): the same
  # numbers in pg, and the same interval.
  for (col in c("ind_median", "ind_mean", "ind_lower", "ind_upper", "biomass_median")) {
    expect_equal(e1[[col]][phy], e0[[col]][phy])
  }
  expect_true(all(e1$currency[phy] == "wet_mass"))
  expect_true(all(e1$unit[phy] == "pg"))
  expect_true(all(e1$converted_from[phy] == "biovolume"))
  expect_true(all(e1$conversion_sources[phy] == "cen2015"))
  expect_match(e1$conversion[phy[1]], "density 1 g/cm3")
  expect_match(e1$note[phy[1]], "not quantified")
  expect_false(any(grepl("conversion", e1$uncertainty[phy])))
  # 1 mm3/L of biovolume is 1 mg/L of wet mass: 1 um3 = 1e-9 mm3, 1 pg = 1e-9 mg.
  mg <- suppressWarnings(convert_currency(est, to = "wet_mass", unit = "mg"))
  e2 <- as.data.frame(mg)
  expect_equal(e2$ind_median[phy], 1e-9 * e0$ind_median[phy])
  expect_equal(e2$biomass_upper[phy], 1e-9 * e0$biomass_upper[phy])
  expect_equal(e2$biomass_unit[phy], paste0("mg/", e0$per[phy]))
  # The density has no error: the draws are those of the biovolumes, and the
  # mean of a total is the sum of its rows.
  ld1 <- withr::with_seed(5, draw_log_mass(mg, phy, 500, "shared"))
  ld0 <- withr::with_seed(5, draw_log_mass(est, phy, 500, "shared"))
  expect_equal(ld1, ld0 + log(1e-9), tolerance = 1e-12)
  s <- summarize_biomass(mg, by = c("site", "broad_group"), na_rm = TRUE, seed = 4)
  s0 <- summarize_biomass(est, by = c("site", "broad_group"), na_rm = TRUE, seed = 4)
  k <- s$broad_group == "phytoplankton"
  expect_equal(s$biomass_mean[k], 1e-9 * s0$biomass_mean[s0$broad_group == "phytoplankton"])
  expect_true(all(s$biomass_unit[k] == "mg/mL"))
  # Natural units are converted too (the conversion is linear).
  ab <- data.frame(site = "a", taxon_id = c("t01", "t01n"), abundance = 10, per = "mL")
  e <- estimate_biomass(ab, cells_and_units)
  w <- as.data.frame(convert_currency(e, to = "wet_mass"))
  expect_equal(w$ind_median, as.data.frame(e)$ind_median)
  expect_true(all(w$unit == "pg"))
  expect_true("cen2015" %in% cite_sources(wm)$source_id)
})

test_that("rotifers go to wet mass with dry mass 10 % of it, 3.9 % in Asplanchna", {
  wm <- suppressWarnings(convert_currency(est, to = "wet_mass"))
  e1 <- as.data.frame(wm)
  rot <- c(row_of("t06"), row_of("t07"))
  # US EPA (2016), section 8.2.2: wet mass = dry mass / 0.1, with the same
  # interval on the log scale.
  for (col in c("ind_median", "ind_mean", "ind_lower", "ind_upper", "biomass_median")) {
    expect_equal(e1[[col]][rot], e0[[col]][rot] / 0.1)
  }
  expect_true(all(e1$currency[rot] == "wet_mass"))
  expect_true(all(e1$unit[rot] == "mg"))
  expect_true(all(e1$converted_from[rot] == "dry_mass"))
  expect_true(all(e1$conversion_sources[rot] == "usepa2016"))
  expect_match(e1$conversion[rot[1]], "rotifers, dry mass 10 % of wet mass")
  expect_match(e1$note[rot[1]], "not quantified")
  expect_false(any(grepl("conversion", e1$uncertainty[rot])))
  ld1 <- withr::with_seed(6, draw_log_mass(wm, rot, 300, "shared"))
  ld0 <- withr::with_seed(6, draw_log_mass(est, rot, 300, "shared"))
  expect_equal(ld1, ld0 - log(0.1), tolerance = 1e-12)
  # Asplanchna: 3.9 % (Dumont et al. 1975, p. 94).
  tx <- ex$taxa[ex$taxa$taxon_id == "t07", ]
  tx$taxon_id <- "asp"
  tx$scientific_name <- "Asplanchna priodonta"
  tx$family <- "Asplanchnidae"
  tx$genus <- "Asplanchna"
  tx$species <- "Asplanchna priodonta"
  tx$col_id <- NA_character_
  ab <- data.frame(site = "a", taxon_id = "asp", abundance = 5, per = "L")
  e <- suppressWarnings(estimate_biomass(ab, tx))
  w <- as.data.frame(convert_currency(e, to = "wet_mass"))
  expect_equal(w$ind_median, as.data.frame(e)$ind_median / 0.039)
  expect_match(w$conversion, "3.9 %.*Asplanchna")
  expect_equal(w$conversion_sources, "dumont1975;usepa2016")
  # Their AFDM, carbon and energy still come from the conversion models.
  cv <- as.data.frame(suppressWarnings(convert_currency(est, to = "carbon")))
  expect_match(cv$conversion[rot[1]], "invertebrate dry_mass_to_carbon")
})

test_that("amphibians go to dry mass and energy with the values of frogs and toads", {
  amph <- which(e0$broad_group == "amphibian" & !is.na(e0$ind_median))
  expect_length(amph, 2L)
  # US EPA (1993, Table 4-1, p. 4-13), frogs and toads, mean (SD) of 3 studies:
  # 85 % water (4.7) and 4.6 kcal/g of dry mass (0.45). Lognormal with that
  # mean: variance between populations s2 = ln(1 + CV^2), plus s2/3 for the
  # mean of the 3 studies; ln-median ln(mean) - (s2 + s2/3)/2.
  s2 <- log(1 + (0.047 / 0.15)^2)
  v <- s2 * (1 + 1 / 3)
  expect_equal(v, 0.1249, tolerance = 1e-3)
  expect_equal(0.15 * exp(-v / 2), 0.1409, tolerance = 1e-3)
  s2e <- log(1 + (0.45 / 4.6)^2)
  ve <- s2e * (1 + 1 / 3)
  dm <- suppressWarnings(convert_currency(est, to = "dry_mass"))
  e1 <- as.data.frame(dm)
  expect_equal(e1$ind_mean[amph], 0.15 * e0$ind_mean[amph])
  expect_equal(e1$biomass_mean[amph], 0.15 * e0$biomass_mean[amph])
  expect_equal(e1$ind_median[amph], 0.15 * exp(-v / 2) * e0$ind_median[amph])
  expect_equal(tau2_of(e1, amph) - tau2_of(e0, amph), rep(v, 2L))
  expect_true(all(e1$currency[amph] == "dry_mass"))
  expect_true(all(e1$unit[amph] == "g"))
  expect_true(all(e1$converted_from[amph] == "wet_mass"))
  expect_true(all(e1$conversion_sources[amph] == "usepa1993"))
  expect_equal(e1$conversion[amph[1]], paste(
    "wet_mass->dry_mass: amphibians, dry mass 15 % (SD 4.7) of wet mass (frogs and toads,",
    "US EPA 1993)"))
  expect_true(all(grepl("; conversion$", e1$uncertainty[amph])))
  expect_true("usepa1993" %in% cite_sources(dm)$source_id)
  # Energy: wet mass (g) * 1e6 ug/g * 0.15 * 4.6 kcal/g * 4.184e-3 J/ug per kcal/g,
  # 70.1 kJ for a mean wet mass of 24.3 g.
  en <- suppressWarnings(convert_currency(est, to = "energy"))
  j <- as.data.frame(en)
  expect_equal(j$ind_mean[amph], e0$ind_mean[amph] * 1e6 * 0.15 * 4.6 * 4.184e-3)
  expect_equal(j$ind_mean[amph[1]] / 1000, 70.1, tolerance = 1e-3)
  expect_equal(tau2_of(j, amph) - tau2_of(e0, amph), rep(v + ve, 2L))
  expect_true(all(j$unit[amph] == "J"))
  expect_match(j$conversion[amph[1]], paste0(
    "^wet_mass->dry_mass: amphibians.*; dry_mass->energy: amphibians, 4.6 kcal/g \\(SD 0.45\\) ",
    "of dry mass \\(frogs and toads, US EPA 1993\\)$"))
  # Draws: the errors added to the log mass have the mean and variance of the
  # convention. With population_effect = "shared" the two rows (one taxon at
  # two sites) share all of it; with "independent" only the error of the mean
  # of the studies (s2/3), which all the amphibian rows share.
  d <- withr::with_seed(8, draw_log_mass(dm, amph, 20000, "shared")) -
    withr::with_seed(8, draw_log_mass(est, amph, 20000, "shared"))
  expect_equal(colMeans(d), rep(log(0.15) - v / 2, 2L), tolerance = 0.01)
  expect_equal(stats::var(d[, 1]), v, tolerance = 0.05)
  expect_equal(d[, 1], d[, 2])
  di <- withr::with_seed(8, draw_log_mass(dm, amph, 20000, "independent")) -
    withr::with_seed(8, draw_log_mass(est, amph, 20000, "independent"))
  expect_equal(stats::var(di[, 2]), v, tolerance = 0.05)
  expect_equal(stats::cov(di[, 1], di[, 2]), s2 / 3, tolerance = 0.1)
  de <- withr::with_seed(9, draw_log_mass(en, amph, 20000, "shared")) -
    withr::with_seed(9, draw_log_mass(est, amph, 20000, "shared"))
  expect_equal(mean(de[, 1]), log(0.15) - v / 2 + log(1e6 * 4.6 * 4.184e-3) - ve / 2,
               tolerance = 0.01)
  expect_equal(stats::var(de[, 1]), v + ve, tolerance = 0.05)
  # AFDM and carbon: not converted.
  cv <- as.data.frame(suppressWarnings(convert_currency(est, to = "carbon")))
  expect_true(all(is.na(cv$ind_median[amph])))
  expect_match(cv$reason[amph[1]], "too few conversion factors for amphibians to carbon")
  # From dry to wet mass (an estimate in dry mass): the opposite shift and the
  # same variance, so the mean wet mass is the mean dry mass / 0.15 * exp(v).
  x <- est
  x$estimates$currency[amph] <- "dry_mass"
  w <- as.data.frame(suppressWarnings(convert_currency(x, to = "wet_mass")))
  expect_equal(w$ind_median[amph], e0$ind_median[amph] / (0.15 * exp(-v / 2)))
  expect_equal(w$ind_mean[amph], e0$ind_mean[amph] / 0.15 * exp(v))
  expect_match(w$conversion[amph[1]], "^dry_mass->wet_mass: amphibians")
  # Other amphibians get the values of frogs and toads, as an extrapolation.
  tx <- without_length(ex$taxa[ex$taxa$taxon_id == "t19", ])
  tx[c("taxon_id", "scientific_name", "order", "family", "genus", "species", "col_id")] <-
    list("sal", "Triturus marmoratus", "Caudata", "Salamandridae", "Triturus",
         "Triturus marmoratus", NA_character_)
  ab <- data.frame(site = "a", taxon_id = "sal", abundance = 3, per = "100 m2")
  e <- estimate_biomass(ab, tx)
  s <- as.data.frame(convert_currency(e, to = "dry_mass"))
  expect_equal(s$ind_mean, 0.15 * as.data.frame(e)$ind_mean)
  expect_match(s$conversion, "US EPA 1993; extrapolated: not a frog or toad\\)$")
  expect_match(s$note, "values of frogs and toads, extrapolated to other amphibians")
})

test_that("rows that cannot be converted get NA and a reason, with one warning", {
  expect_warning(cv <- convert_currency(est, to = "dry_mass"), "could not be converted",
                 class = "bmassa_warning_not_estimable")
  e1 <- as.data.frame(cv)
  phy <- e0$broad_group == "phytoplankton" & !is.na(e0$ind_median)
  expect_true(all(is.na(e1$ind_median[phy])))
  expect_true(all(grepl("no conversion from biovolume", e1$reason[phy])))
  amph <- which(e0$broad_group == "amphibian" & !is.na(e0$ind_median))
  # Amphibians go only between wet and dry mass: to AFDM they are not converted.
  af <- as.data.frame(suppressWarnings(convert_currency(est, to = "afdm")))
  expect_true(all(is.na(af$ind_median[amph])))
  expect_true(all(af$status[amph] == "not_estimated"))
  expect_match(af$reason[amph[1]], "too few conversion factors for amphibians to afdm")
  # Wet mass: amphibians and fish stay as they are.
  wm <- suppressWarnings(convert_currency(est, to = "wet_mass"))
  e2 <- as.data.frame(wm)
  vert <- e0$broad_group %in% c("fish", "amphibian")
  expect_equal(e2$ind_median[vert], e0$ind_median[vert])
  expect_true(all(is.na(e2$converted_from[vert])))
})

test_that("natural units of phytoplankton go to carbon only with cells_per_unit", {
  ab <- data.frame(site = "a", taxon_id = c("t01", "t01n"), abundance = 10, per = "mL")
  e <- estimate_biomass(ab, cells_and_units)
  cv <- suppressWarnings(convert_currency(e, to = "carbon"))
  out <- as.data.frame(cv)
  expect_false(is.na(out$ind_median[1]))
  expect_true(is.na(out$ind_median[2]))
  expect_match(out$reason[2], "single cells")
  # With 8 cells per natural unit: the biovolume of 8 cells, from the model of
  # cells, with its interval; and the carbon of 8 cells (hand calculation:
  # ln C = ln10 log10(a) + (1 - b) ln 8 + b ln(8 V) = ln 8 + ln C of one cell).
  tx8 <- cells_and_units
  tx8$cells_per_unit <- c(NA, 8)
  e8 <- estimate_biomass(ab, tx8)
  d8 <- as.data.frame(e8)
  expect_equal(d8$ind_median[2], 8 * d8$ind_median[1])
  expect_equal(d8$ind_upper[2] / d8$ind_median[2], d8$ind_upper[1] / d8$ind_median[1])
  expect_equal(d8$count_unit[2], "natural_unit")
  expect_match(d8$note[2], "natural unit of 8 cells")
  c8 <- as.data.frame(suppressWarnings(convert_currency(e8, to = "carbon")))
  expect_equal(c8$ind_median[2], 8 * c8$ind_median[1], tolerance = 1e-12)
  # Colonial genera counted in natural units without cells_per_unit: a note.
  mic <- transform(ex$taxa[ex$taxa$taxon_id == "t02", ], count_unit = "natural_unit")
  em <- as.data.frame(estimate_biomass(data.frame(site = "a", taxon_id = "t02", abundance = 1,
                                                  per = "mL"), mic))
  expect_match(em$note, "not whole colonies")
  # cells_per_unit only for phytoplankton in natural units, and at least 1.
  bad <- cells_and_units
  bad$cells_per_unit <- c(8, NA)
  expect_error(estimate_biomass(ab, bad), class = "bmassa_error_input")
  bad$cells_per_unit <- c(NA, 0.5)
  expect_error(estimate_biomass(ab, bad), class = "bmassa_error_input")
})

test_that("phytoplankton cells go to energy through their carbon (hand calculation)", {
  # Platt and Irwin (1973, p. 309 and Table 3): 11.40 cal/mg C, SD 1.108 among
  # 10 samples. Lognormal with that mean, as the conventions of amphibians:
  # variance between populations s2 = ln(1 + CV^2), plus s2/10 for the mean.
  s2 <- log(1 + (1.108 / 11.40)^2)
  v <- s2 * (1 + 1 / 10)
  f <- 11.40 * 4.184e-9 # J per pg of carbon
  c1 <- as.data.frame(suppressWarnings(convert_currency(est, to = "carbon")))
  en <- suppressWarnings(convert_currency(est, to = "energy"))
  j <- as.data.frame(en)
  i <- row_of("t01")
  expect_equal(j$ind_mean[i], c1$ind_mean[i] * f)
  expect_equal(j$ind_median[i], c1$ind_median[i] * f * exp(-v / 2))
  expect_equal(tau2_of(j, i) - tau2_of(c1, i), v)
  expect_equal(j$unit[i], "J")
  expect_match(j$conversion[i], "Platt and Irwin 1973")
  expect_true(grepl("conversion", j$uncertainty[i]))
  expect_true("platt1973" %in% cite_sources(en)$source_id)
  # 1000 um3 of a diatom: 77.98 pg of carbon, 3.72 uJ.
  expect_equal(77.98 * f * 1e6, 3.720, tolerance = 1e-3)
})

test_that("units, abundances and repeated conversions", {
  mg <- suppressWarnings(convert_currency(est, to = "carbon"))
  ug <- suppressWarnings(convert_currency(est, to = "carbon", unit = "ug"))
  a <- as.data.frame(mg)
  b <- as.data.frame(ug)
  inv <- which(a$unit %in% "mg")
  expect_equal(b$ind_median[inv], 1000 * a$ind_median[inv])
  expect_equal(b$biomass_upper[inv], 1000 * a$biomass_upper[inv])
  pg <- which(a$unit %in% "pg")
  expect_equal(b$ind_median[pg], 1e-6 * a$ind_median[pg])
  expect_true(all(b$unit[!is.na(b$unit)] == "ug"))
  expect_equal(b$biomass_unit[inv], paste0("ug/", a$per[inv]))
  expect_equal(a$biomass_mean, a$abundance * a$ind_mean)
  # A conversion always starts from the currency of the estimate.
  twice <- suppressWarnings(convert_currency(convert_currency(est, to = "wet_mass"), to = "carbon"))
  expect_equal(as.data.frame(twice), a)
  # Same rows, same order.
  expect_equal(a[, c("site", "taxon_id")], e0[, c("site", "taxon_id")])
})

test_that("the draws of converted rows reproduce their distribution, and totals their means", {
  cv <- suppressWarnings(convert_currency(est, to = "carbon", unit = "ug"))
  e1 <- as.data.frame(cv)
  ok <- which(!is.na(e1$ind_median))
  ld <- withr::with_seed(21, draw_log_mass(cv, ok, 20000, "shared"))
  expect_equal(apply(ld, 2, stats::var), 2 * log(e1$ind_mean[ok] / e1$ind_median[ok]),
               tolerance = 0.05)
  expect_equal(colMeans(ld), log(e1$ind_median[ok]), tolerance = 0.05)
  s <- summarize_biomass(cv, by = c("site", "broad_group"), na_rm = TRUE, seed = 3)
  kb <- !is.na(e1$biomass_mean)
  hand <- tapply(e1$biomass_mean[kb], paste(e1$site, e1$broad_group)[kb], sum)
  expect_equal(s$biomass_mean[!is.na(s$currency)],
               as.vector(hand[paste(s$site, s$broad_group)[!is.na(s$currency)]]))
  expect_true(all(s$currency[!is.na(s$currency)] == "carbon"))
})

test_that("rows sharing a conversion node have correlated factors", {
  # Two species absent from the reference, of the same genus: both use the
  # node of Daphnia in the carbon model, and share its error.
  tx <- without_length(ex$taxa[ex$taxa$taxon_id == "t08", ])
  tx <- rbind(tx, tx)
  tx$taxon_id <- c("d1", "d2")
  tx$species <- c("Daphnia alpha", "Daphnia beta")
  tx$col_id <- NA_character_
  ab <- data.frame(site = "a", taxon_id = c("d1", "d2"), abundance = 1, per = "L")
  e <- estimate_biomass(ab, tx)
  cv <- suppressWarnings(convert_currency(e, to = "carbon"))
  st <- cv$simulation$conversion$steps
  expect_equal(nrow(st), 2L)
  expect_equal(st$node[1], st$node[2])
  expect_equal(ref$conversion_nodes$node[st$node[1]], "Daphnia")
  ld0 <- withr::with_seed(5, draw_log_mass(e, 1:2, 20000, "shared"))
  ld1 <- withr::with_seed(5, draw_log_mass(cv, 1:2, 20000, "shared"))
  # The added covariance is the variance of the shared parts of the factor:
  # the path of the node and the habitat.
  m <- "invertebrate dry_mass_to_carbon"
  shared <- ref$conversion_nodes$path_var[st$node[1]] + comp(m, "habitat")
  extra <- stats::cov(ld1)[1, 2] - stats::cov(ld0)[1, 2]
  expect_equal(extra, shared, tolerance = 0.25)
})

test_that("zero rows, and arguments are checked", {
  z <- suppressWarnings(estimate_biomass(ex$abundance[0, ], ex$taxa))
  cz <- convert_currency(z, to = "carbon")
  est <- suppressWarnings(estimate_biomass(ex$abundance[ex$abundance$taxon_id == "t10", ], ex$taxa,
                                           bias_correction = "lognormal"))
  expect_equal(nrow(as.data.frame(cz)), 0L)
  expect_true(all(c("converted_from", "conversion", "conversion_sources") %in%
                    names(as.data.frame(cz))))
  cls <- "bmassa_error_input"
  expect_error(convert_currency(est, to = "biovolume"), class = cls)
  expect_error(convert_currency(est, to = "carbon", habitat = "lake"), class = cls)
  expect_error(convert_currency(est, to = "carbon", unit = "kg"), class = cls)
  expect_error(convert_currency(est, to = "carbon", unit = "J"), class = cls)
  expect_error(convert_currency(est, to = "energy", unit = "mg"), class = cls)
  expect_error(convert_currency(est, to = "carbon", colour = "red"), class = cls)
})

# Energy --------------------------------------------------------------------------------

test_that("energy is the mass in ug times the energy density of the taxon (hand calculation)", {
  en <- suppressWarnings(convert_currency(est, to = "energy"))
  e1 <- as.data.frame(en)
  # Invertebrates: dry mass in mg, energy per dry mass in J/ug.
  m <- "invertebrate dry_mass_to_energy"
  for (id in c("t10", "t13", "t08")) {
    i <- row_of(id)
    n <- ref$conversion_nodes[conv_node(m, id), ]
    v_f <- n$path_var + missing_var(m, n$level) + comp(m, "residual") + comp(m, "habitat")
    # ln E (J) = ln M (mg) + ln(1000 ug/mg) + ln f (J/ug)
    expect_equal(log(e1$ind_median[i]), log(e0$ind_median[i]) + log(1000) + n$path_intercept,
                 tolerance = 1e-12)
    expect_equal(tau2_of(e1, i) - tau2_of(e0, i), v_f, tolerance = 1e-10)
    expect_equal(e1$unit[i], "J")
    expect_equal(e1$biomass_unit[i], paste0("J/", e0$per[i]))
    expect_match(e1$conversion[i], "^dry_mass->energy: invertebrate dry_mass_to_energy")
  }
  # Fish: wet mass in g, energy per wet mass.
  m <- "fish wet_mass_to_energy"
  i <- row_of("t16")
  n <- ref$conversion_nodes[conv_node(m, "t16"), ]
  expect_equal(log(e1$ind_median[i]), log(e0$ind_median[i]) + log(1e6) + n$path_intercept,
               tolerance = 1e-12)
  # A dry mass of 1 mg of an invertebrate holds about 17 J (median of the
  # energy densities per dry mass in the fit, conversion_fit_log.md): the
  # model's mean at its root is of that order.
  root <- ref$conversion_nodes[ref$conversion_nodes$model == "invertebrate dry_mass_to_energy" &
                                 ref$conversion_nodes$level == "root", ]
  expect_true(exp(root$path_intercept) * 1000 > 10 && exp(root$path_intercept) * 1000 < 25)
})

test_that("energy units rescale every converted row", {
  j <- as.data.frame(suppressWarnings(convert_currency(est, to = "energy")))
  kj <- as.data.frame(suppressWarnings(convert_currency(est, to = "energy", unit = "kJ")))
  mj <- as.data.frame(suppressWarnings(convert_currency(est, to = "energy", unit = "mJ")))
  ok <- which(!is.na(j$ind_median))
  expect_equal(kj$ind_median[ok], j$ind_median[ok] / 1000)
  expect_equal(mj$biomass_upper[ok], j$biomass_upper[ok] * 1000)
  expect_true(all(kj$unit[ok] == "kJ"))
  phy <- which(e0$broad_group == "phytoplankton" & !is.na(e0$ind_median))
  expect_true(all(kj$unit[phy] == "kJ"))
  expect_equal(kj$ind_mean[phy], j$ind_mean[phy] / 1000)
  amph <- which(e0$broad_group == "amphibian" & !is.na(e0$ind_median))
  expect_true(all(kj$unit[amph] == "kJ"))
  expect_equal(kj$ind_mean[amph], j$ind_mean[amph] / 1000)
  # Totals in energy add up the row means, and keep the units apart.
  en <- suppressWarnings(convert_currency(est, to = "energy"))
  s <- summarize_biomass(en, by = c("site", "broad_group"), na_rm = TRUE, seed = 7)
  kb <- !is.na(j$biomass_mean)
  hand <- tapply(j$biomass_mean[kb], paste(j$site, j$broad_group)[kb], sum)
  has <- !is.na(s$currency)
  expect_equal(s$biomass_mean[has], as.vector(hand[paste(s$site, s$broad_group)[has]]))
  expect_true(all(s$currency[has] == "energy"))
  expect_true(all(grepl("^J/", s$biomass_unit[has])))
})

test_that("the energy models come from the screened data, without carbon", {
  cm <- ref$conversion_models
  expect_true(all(c("invertebrate dry_mass_to_energy", "invertebrate afdm_to_energy",
                    "invertebrate wet_mass_to_energy", "fish wet_mass_to_energy",
                    "fish dry_mass_to_energy", "fish afdm_to_energy") %in% cm$model))
  expect_false(any(cm$from == "carbon" & cm$to == "energy"))
  # Energy densities per ash-free dry mass of the model roots lie within the
  # bounds of the screening (15-45 J/mg).
  r <- ref$conversion_nodes[ref$conversion_nodes$level == "root" &
                              grepl("afdm_to_energy$", ref$conversion_nodes$model), ]
  expect_true(all(exp(r$path_intercept) * 1000 > 15 & exp(r$path_intercept) * 1000 < 45))
})

test_that("print() reports the conversion", {
  skip_on_cran()
  cv <- suppressWarnings(convert_currency(est, to = "carbon", unit = "ug"))
  expect_snapshot(print(cv))
  en <- suppressWarnings(convert_currency(est, to = "energy"))
  expect_snapshot(print(en))
})

test_that("convert_currency() needs the result of estimate_biomass()", {
  expect_error(convert_currency(data.frame(taxon = "Daphnia"), to = "carbon"),
               class = "bmassa_error_input")
})
