# estimate_biomass() ---------------------------------------------------------------

ref <- get_reference()
taxon <- function(id, group, lineage) {
  x <- as.list(c(taxon_id = id, broad_group = group))
  lin <- stats::setNames(as.list(rep(NA_character_, 7)), c("kingdom", "phylum", "class", "order",
                                                            "family", "genus", "species"))
  lin[names(lineage)] <- lineage
  as.data.frame(c(x, lin), stringsAsFactors = FALSE)
}
limno <- taxon("w1", "macroinvertebrate",
               c(kingdom = "Animalia", phylum = "Annelida", class = "Clitellata", order = "Tubificida",
                 family = "Naididae", genus = "Limnodrilus", species = "Limnodrilus hoffmeisteri"))
chiro <- taxon("c1", "macroinvertebrate",
               c(kingdom = "Animalia", phylum = "Arthropoda", class = "Insecta", order = "Diptera",
                 family = "Chironomidae", genus = "Chironomus", species = "Chironomus fictus"))
trout <- taxon("f1", "fish",
               c(kingdom = "Animalia", phylum = "Chordata", class = "Teleostei",
                 order = "Salmoniformes", family = "Salmonidae", genus = "Salmo", species = "Salmo trutta"))
# The mean length of a taxon unit in all its samples, as a table of lengths.
len <- function(tx, mean, unit, type, ...) {
  data.frame(taxon_id = tx$taxon_id, length_mean = mean, length_unit = unit, length_type = type,
             ..., stringsAsFactors = FALSE)
}
comp <- function(m, l, p) {
  v <- ref$components$value[ref$components$model == m & ref$components$level == l &
                              ref$components$parameter == p]
  if (length(v)) v else 0
}

test_that("a species matched in the reference is predicted from its path (hand calculation)", {
  ab <- data.frame(site = "s", taxon_id = "w1", abundance = 10, per = "m2")
  e <- as.data.frame(estimate_biomass(ab, limno, currency = "dry_mass"))
  m <- "invertebrate_individual"
  nd <- ref$nodes[ref$nodes$model == m & ref$nodes$level == "species" &
                    ref$nodes$node == "Limnodrilus hoffmeisteri", ]
  st <- ref$nodes[ref$nodes$model == m & ref$nodes$level == "stage" & ref$nodes$node == "non_insect", ]
  eta <- nd$path_intercept + st$effect
  tau2 <- nd$path_var + st$effect_var + comp(m, "residual", "intercept") + comp(m, "source", "intercept")
  expect_equal(e$match_level, "species")
  expect_equal(e$ind_median, exp(eta))
  expect_equal(e$ind_mean, exp(eta + tau2 / 2))
  expect_equal(e$ind_lower, exp(eta - stats::qnorm(0.975) * sqrt(tau2)))
  expect_equal(e$biomass_mean, 10 * e$ind_mean)
  expect_equal(e$biomass_unit, "mg/m2")
  expect_equal(e$stage_used, "non_insect")
  # The model of invertebrates without lengths has no effect of source.
  expect_equal(e$uncertainty, "parameter; between_population")
})

test_that("a species missing from the reference adds the variance of its level at its length", {
  ab <- data.frame(site = "s", taxon_id = "c1", abundance = 1, per = "m2")
  e <- as.data.frame(estimate_biomass(ab, chiro, len(chiro, 10, "mm", "body_length"),
                                      currency = "dry_mass", bias_correction = "none"))
  m <- "invertebrate_lwr"
  nd <- ref$nodes[ref$nodes$model == m & ref$nodes$level == "genus" & ref$nodes$node == "Chironomus", ]
  lt <- ref$lwr_types[ref$lwr_types$model == m & ref$lwr_types$length_type == "body_length", ]
  xc <- log(10) - lt$centre
  eta <- nd$path_intercept + nd$path_slope * xc
  v <- nd$path_var + 2 * xc * nd$path_cov + xc^2 * nd$path_var_slope +
    comp(m, "species", "intercept") + 2 * xc * comp(m, "species", "intercept_slope") +
    xc^2 * comp(m, "species", "slope") + comp(m, "residual", "intercept") +
    comp(m, "equation", "intercept") + comp(m, "study", "intercept")
  expect_equal(e$match_level, "genus")
  expect_equal(e$ind_median, exp(eta))
  expect_equal(e$ind_mean, exp(eta + v / 2))
  expect_match(e$uncertainty, "taxonomic_imputation")
  expect_equal(e$length_status, "within")
  expect_match(e$note, "without its spread")
})

test_that("a published equation lies inside the interval (Salmo trutta, Arranz et al. pooled)", {
  # ln W = -11.46273 + 3.021065 ln FL (g, mm; data-raw/derived/lwr_fitted.csv,
  # arranz2022_fit_0267, 1461 fish): 94.1 g at 200 mm.
  ab <- data.frame(site = "s", taxon_id = "f1", abundance = 1, per = "100 m2")
  e <- as.data.frame(estimate_biomass(ab, trout, len(trout, 200, "mm", "fork_length"),
                                      unit = "g", bias_correction = "none"))
  w <- exp(-11.46273 + 3.021065 * log(200))
  expect_equal(w, 94.1, tolerance = 0.001)
  expect_true(e$ind_lower < w && w < e$ind_upper)
  expect_equal(e$unit, "g")
})

test_that("the corrections for mean lengths follow their formulas", {
  ab <- data.frame(site = "s", taxon_id = "c1", abundance = 1, per = "m2")
  ln0 <- len(chiro, 10, "mm", "body_length")
  ln1 <- transform(ln0, length_cv = 0.3)
  e0 <- as.data.frame(estimate_biomass(ab, chiro, ln0, currency = "dry_mass", bias_correction = "none"))
  e1 <- as.data.frame(estimate_biomass(ab, chiro, ln1, currency = "dry_mass", bias_correction = "none"))
  e2 <- as.data.frame(estimate_biomass(ab, chiro, ln1, currency = "dry_mass",
                                       bias_correction = "lognormal"))
  nd <- ref$nodes[ref$nodes$model == "invertebrate_lwr" & ref$nodes$node == "Chironomus", ]
  b <- nd$path_slope
  # Lognormal lengths: E[L^b] / E[L]^b = (1 + CV^2)^(b (b - 1) / 2); with b = 3
  # and CV = 0.3 this is 1.09^3 = 1.295.
  expect_equal(e1$ind_median / e0$ind_median, (1 + 0.3^2)^(b * (b - 1) / 2))
  expect_equal((1 + 0.3^2)^3, 1.295029, tolerance = 1e-6)
  s <- ref$lwr_types$sigma_ind[ref$lwr_types$model == "invertebrate_lwr" &
                                 ref$lwr_types$length_type == "body_length"]
  expect_equal(e2$ind_median / e1$ind_median, exp(s^2 / 2))
  expect_equal(e1$corrections, "jensen")
  expect_equal(e2$corrections, "jensen; lognormal")
  expect_equal(transform(ln1, length_cv = NULL, length_sd = 3)$length_sd / 10, 0.3)
  e3 <- as.data.frame(estimate_biomass(ab, chiro, transform(ln0, length_sd = 3), currency = "dry_mass",
                                       bias_correction = "none"))
  expect_equal(e3$ind_median, e1$ind_median)
})

test_that("units: the same length in mm and cm, and abundances scale the biomass", {
  ab <- data.frame(site = "s", taxon_id = "f1", abundance = c(1), per = "m2")
  ln_mm <- len(trout, 150, "mm", "fork_length")
  ln_cm <- len(trout, 15, "cm", "fork_length")
  a <- as.data.frame(estimate_biomass(ab, trout, ln_mm, bias_correction = "none"))
  b <- as.data.frame(estimate_biomass(ab, trout, ln_cm, bias_correction = "none"))
  expect_equal(a$ind_median, b$ind_median)
  ab10 <- transform(ab, abundance = 10)
  c10 <- as.data.frame(estimate_biomass(ab10, trout, ln_mm, bias_correction = "none"))
  expect_equal(c10$biomass_mean, 10 * a$biomass_mean)
})

test_that("zeros and missing abundances", {
  ab <- data.frame(site = c("a", "b"), taxon_id = "w1", abundance = c(0, NA), per = "m2")
  e <- suppressWarnings(as.data.frame(estimate_biomass(ab, limno)))
  expect_equal(e$biomass_mean[1], 0)
  expect_equal(e$status, c("estimated", "not_estimated"))
  expect_true(is.na(e$biomass_mean[2]))
  expect_equal(e$reason[2], "abundance is missing")
  expect_equal(e$ind_median[1], e$ind_median[2])
})

test_that("lengths outside the reference move up the tree, or are not estimated", {
  ab <- data.frame(site = "s", taxon_id = "c1", abundance = 1, per = "m2")
  tiny <- len(chiro, 0.02, "mm", "body_length")
  e <- suppressWarnings(as.data.frame(estimate_biomass(ab, chiro, tiny, bias_correction = "none")))
  expect_equal(e$status, "not_estimated")
  expect_match(e$reason, "outside the lengths")
  e2 <- as.data.frame(estimate_biomass(ab, chiro, tiny, bias_correction = "none", extrapolate = TRUE))
  expect_equal(e2$length_status, "outside")
  expect_equal(e2$match_level, "genus")
  ex <- as.data.frame(estimate_biomass(example_survey$abundance[12, ], example_survey$taxa,
                                       example_survey$lengths, bias_correction = "none"))
  expect_equal(ex$match_name, "Coenagrionidae")
  expect_match(ex$note, "estimated from Coenagrionidae")
})

test_that("fish without length: mean mass of a population (hand calculation)", {
  ab <- data.frame(site = "s", taxon_id = c("f1", "f2", "f3"), abundance = 2, per = "m2")
  smelt <- taxon("f2", "fish",
                 c(kingdom = "Animalia", phylum = "Chordata", class = "Teleostei",
                   order = "Osmeriformes", family = "Osmeridae", genus = "Osmerus",
                   species = "Osmerus eperlanus"))
  stingray <- taxon("f3", "fish",
                    c(kingdom = "Animalia", phylum = "Chordata", class = "Elasmobranchii",
                      order = "Myliobatiformes", family = "Potamotrygonidae",
                      genus = "Potamotrygon", species = "Potamotrygon motoro"))
  e <- as.data.frame(estimate_biomass(ab, rbind(trout, smelt, stingray), unit = "g"))
  expect_equal(e$status, rep("estimated", 3))
  expect_equal(e$method, rep("typical", 3))
  expect_equal(e$unit, rep("g", 3))
  expect_match(e$note[1], "fish without length")
  # The trout from its node of the fish model: path, its variance, the
  # variation between populations (populations and samplings of a
  # population) and that of a new source.
  m <- "fish_individual"
  pred_tau2 <- function(k) {
    nd <- ref$nodes[ref$nodes$model == m & ref$nodes$node == e$match_name[k] &
                      ref$nodes$level == e$match_level[k], ]
    tau2 <- nd$path_var + comp(m, "residual", "intercept") + comp(m, "population", "intercept") +
      comp(m, "source", "intercept")
    lv <- strsplit(ref$models$levels[ref$models$model == m], ",")[[1]]
    miss <- lv[seq_along(lv) > match(e$match_level[k], c("root", lv)) - 1L]
    for (l in miss) tau2 <- tau2 + comp(m, l, "intercept")
    list(eta = nd$path_intercept, tau2 = tau2)
  }
  expect_gt(comp(m, "population", "intercept"), 0)
  p1 <- pred_tau2(1)
  expect_equal(e$ind_median[1], exp(p1$eta))
  expect_equal(e$ind_upper[1], exp(p1$eta + stats::qnorm(0.975) * sqrt(p1$tau2)))
  expect_equal(e$biomass_median[1], 2 * e$ind_median[1])
  # An order absent from the reference: from the fish of its class, with the
  # variances of the levels below.
  expect_equal(e$match_level[2], "class")
  expect_equal(e$match_name[2], "Teleostei")
  p2 <- pred_tau2(2)
  expect_equal(e$ind_upper[2], exp(p2$eta + stats::qnorm(0.975) * sqrt(p2$tau2)))
  expect_gt(e$ind_upper[2] / e$ind_median[2], e$ind_upper[1] / e$ind_median[1])
  # A class absent too: from all the fish, with the variances of all the
  # levels.
  expect_equal(e$match_level[3], "root")
  expect_match(e$note[3], "class not in the reference")
  expect_gte(e$ind_upper[3] / e$ind_median[3], e$ind_upper[2] / e$ind_median[2])
})

test_that("rows that cannot be estimated say why, with one warning", {
  ab <- data.frame(site = "s", taxon_id = c("f1", "w1"), abundance = 1, per = "m2")
  tx <- rbind(transform(trout, life_stage = "juvenile"), transform(limno, life_stage = NA))
  expect_warning(est <- estimate_biomass(ab, tx), class = "bmassa_warning_not_estimable")
  e <- as.data.frame(est)
  expect_equal(e$status, c("not_estimated", "estimated"))
  expect_match(e$reason[1], "whole populations")
  # Amphibians without a stage are larvae: newts, without masses of their
  # larvae, are not estimated; adults are.
  newt <- taxon("a1", "amphibian", c(kingdom = "Animalia", phylum = "Chordata", class = "Amphibia",
                                     order = "Caudata", family = "Salamandridae", genus = "Triturus",
                                     species = "Triturus marmoratus"))
  newts <- rbind(transform(newt, life_stage = NA), transform(newt, taxon_id = "a2", life_stage = "adult"))
  ab2 <- data.frame(site = "s", taxon_id = c("a1", "a2"), abundance = 1, per = "m2")
  e2 <- suppressWarnings(as.data.frame(estimate_biomass(ab2, newts)))
  expect_equal(e2$status, c("not_estimated", "estimated"))
  expect_match(e2$reason[1], "taken as larvae")
  odd <- taxon("x1", "macroinvertebrate", c(kingdom = "Animalia", phylum = "Bryozoa"))
  e3 <- suppressWarnings(as.data.frame(estimate_biomass(data.frame(site = "s", taxon_id = "x1",
                                                                   abundance = 1, per = "L"), odd)))
  expect_match(e3$reason, "no taxon of the reference model")
})

test_that("insects without a stage are taken as larvae", {
  units <- rbind(transform(chiro, life_stage = NA), transform(chiro, taxon_id = "c2", life_stage = "larva"),
                 transform(chiro, taxon_id = "c3", life_stage = "adult"))
  ab <- data.frame(site = "s", taxon_id = c("c1", "c2", "c3"), abundance = 1, per = "m2")
  e <- as.data.frame(estimate_biomass(ab, units))
  expect_equal(e$stage_used, c("larva (assumed)", "larva", "adult"))
  expect_equal(e$ind_median[1], e$ind_median[2])
  expect_match(e$note[1], "taken as larva")
})

test_that("copepods are taken as adults without a stage, with contrasts for copepodites and nauplii", {
  cyc <- taxon("k1", "zooplankton",
               c(kingdom = "Animalia", phylum = "Arthropoda", class = "Copepoda", order = "Cyclopoida",
                 family = "Cyclopidae", genus = "Cyclops", species = "Cyclops abyssorum"))
  units <- rbind(transform(cyc, life_stage = NA), transform(cyc, taxon_id = "k2", life_stage = "adult"),
                 transform(cyc, taxon_id = "k3", life_stage = "copepodite"),
                 transform(cyc, taxon_id = "k4", life_stage = "nauplius"),
                 transform(cyc, taxon_id = "k5", life_stage = "juvenile"))
  ab <- data.frame(site = "s", taxon_id = paste0("k", 1:5), abundance = 1, per = "L")
  e <- suppressWarnings(as.data.frame(estimate_biomass(ab, units, currency = "dry_mass")))
  expect_equal(e$stage_used, c("copepod_adult (assumed)", "copepod_adult", "copepod_copepodite",
                               "copepod_nauplius", NA))
  expect_equal(e$ind_median[1], e$ind_median[2])
  expect_match(e$note[1], "taken as adult")
  expect_true(e$ind_median[4] < e$ind_median[3] && e$ind_median[3] < e$ind_median[2])
  expect_equal(e$status[5], "not_estimated")
  expect_match(e$reason[5], "copepodite")
  # Hand calculation: path of the node, effect of the group of the other
  # invertebrates and fixed contrast, with the variance of the contrast and
  # twice its covariance with the path.
  m <- "invertebrate_individual"
  nd <- ref$nodes[ref$nodes$model == m & ref$nodes$level == e$match_level[4] &
                    ref$nodes$node == e$match_name[4], ]
  st <- ref$nodes[ref$nodes$model == m & ref$nodes$level == "stage" & ref$nodes$node == "non_insect", ]
  for (k in c("copepodite", "nauplius")) {
    i <- if (k == "copepodite") 3 else 4
    fx <- ref$fixed[ref$fixed$model == m & ref$fixed$term == paste0("c_copepod_", k), ]
    expect_equal(nrow(fx), 1L)
    expect_equal(e$ind_median[i], exp(nd$path_intercept + st$effect + fx$estimate))
    lv <- strsplit(ref$models$levels[ref$models$model == m], ",")[[1]]
    below <- lv[seq_along(lv) > match(e$match_level[i], lv)]
    v <- nd$path_var + st$effect_var + fx$se^2 + 2 * nd[[paste0("contrast_cov_copepod_", k)]] +
      sum(vapply(below, function(l) comp(m, l, "intercept"), numeric(1))) +
      comp(m, "residual", "intercept") + comp(m, "source", "intercept")
    expect_equal(e$ind_upper[i], exp(log(e$ind_median[i]) + stats::qnorm(0.975) * sqrt(v)))
  }
  # The adults have no contrast: path and group only.
  expect_equal(e$ind_median[2], exp(nd$path_intercept + st$effect))
  # Older classifications without the class Copepoda: the order tells.
  old <- transform(cyc, taxon_id = "k4", class = "Maxillopoda", life_stage = "nauplius")
  e2 <- as.data.frame(estimate_biomass(ab[4, ], old))
  expect_equal(e2$stage_used, "copepod_nauplius")
})

test_that("draws of copepod nauplii have the variance of their prediction", {
  nau <- taxon("n1", "zooplankton",
               c(kingdom = "Animalia", phylum = "Arthropoda", class = "Copepoda", order = "Cyclopoida"))
  nau$life_stage <- "nauplius"
  ab <- data.frame(site = c("a", "b"), taxon_id = "n1", abundance = 1, per = "L")
  e <- estimate_biomass(ab, nau)
  s <- simulate(e, nsim = 20000, seed = 3)
  ln <- log(s$ind_mass[1, ])
  es <- as.data.frame(e)
  sd_pred <- log(es$ind_upper[1] / es$ind_median[1]) / stats::qnorm(0.975)
  expect_equal(stats::sd(ln), sd_pred, tolerance = 0.03)
  expect_lt(abs(stats::median(ln) - log(es$ind_median[1])), 0.05 * sd_pred)
})

test_that("the body lengths of copepod nauplii use relations of their own", {
  cal <- taxon("n1", "zooplankton",
               c(kingdom = "Animalia", phylum = "Arthropoda", class = "Copepoda", order = "Calanoida"))
  units <- rbind(transform(cal, life_stage = "nauplius"),
                 transform(cal, taxon_id = "n2", life_stage = "copepodite"))
  ln <- data.frame(taxon_id = c("n1", "n2"), length_mean = 0.3, length_unit = "mm",
                   length_type = "body_length")
  ab <- data.frame(site = "s", taxon_id = c("n1", "n2"), abundance = 1, per = "L")
  e <- suppressWarnings(as.data.frame(estimate_biomass(ab, units, ln, currency = "dry_mass",
                                                       bias_correction = "none", extrapolate = TRUE)))
  expect_equal(e$length_type, c("body_length", "body_length"))
  # Hand calculation from the node of the nauplii.
  m <- "invertebrate_lwr"
  nd <- ref$nodes[ref$nodes$model == m & ref$nodes$level == e$match_level[1] &
                    ref$nodes$node == e$match_name[1] &
                    ref$nodes$length_type %in% "nauplius_body_length", ]
  lt <- ref$lwr_types[ref$lwr_types$model == m & ref$lwr_types$length_type == "nauplius_body_length", ]
  expect_equal(nrow(nd), 1L)
  expect_equal(e$ind_median[1], exp(nd$path_intercept + nd$path_slope * (log(0.3) - lt$centre)))
  expect_false(isTRUE(all.equal(e$ind_median[1], e$ind_median[2])))
})

test_that("the most specific lengths apply: sample, site, all the samples", {
  ex <- example_survey
  e <- suppressWarnings(as.data.frame(estimate_biomass(ex$abundance, ex$taxa, ex$lengths,
                                                       bias_correction = "none")))
  sq <- e[e$taxon_id == "t17", ]
  expect_equal(sq$length_mm, c(110, 85))
  # The frog as adults (t19, with a length) and as larvae (t20, without: from
  # the typical length of tadpoles).
  fr <- e[e$taxon_id %in% c("t19", "t20"), ]
  expect_equal(fr$method, c("length", "typical_length", "length", "typical_length"))
  expect_equal(fr$life_stage, c("adult", "larva", "adult", "larva"))
  # Gammarus measured in L1 only: from its individuals there, without length in L2.
  gm <- e[e$taxon_id == "t13", ]
  expect_equal(gm$method, c("length", "typical"))
  expect_equal(gm$length_basis, c("individuals", NA))
  expect_equal(gm$length_n, c(12, NA))
  # A length for all the samples, replaced in one site and on one date.
  ab <- data.frame(site = c("a", "a", "b", "b"), date = c(1, 2, 1, 2), taxon_id = "c1",
                   abundance = 1, per = "m2")
  ln <- data.frame(taxon_id = "c1", site = c(NA, "a", "b"), date = c(NA, NA, 2),
                   length_mean = c(10, 8, 12), length_unit = "mm", length_type = "body_length")
  e2 <- as.data.frame(estimate_biomass(ab, chiro, ln, bias_correction = "none"))
  expect_equal(e2$length_mm, c(8, 8, 10, 12))
  expect_equal(e2$length_basis, rep("mean", 4))
  expect_identical(e2$date, ab$date)
})

test_that("individual lengths give the mean mass of the individuals (hand calculation)", {
  ab <- data.frame(site = c("a", "b"), taxon_id = "c1", abundance = 5, per = "m2")
  L <- c(6, 8, 10, 12, 14)
  n <- c(1, 3, 4, 3, 1)
  ln <- rbind(data.frame(taxon_id = "c1", site = "a", length = L, n = n, length_unit = "mm",
                         length_type = "body_length"),
              # One row per individual gives the same as their counts.
              data.frame(taxon_id = "c1", site = "b", length = rep(L, n), n = 1, length_unit = "mm",
                         length_type = "body_length"))
  e <- as.data.frame(estimate_biomass(ab, chiro, ln, currency = "dry_mass", bias_correction = "none"))
  expect_equal(e$length_basis, c("individuals", "individuals"))
  expect_equal(e$length_n, c(12, 12))
  expect_equal(e$ind_median[1], e$ind_median[2])
  m <- "invertebrate_lwr"
  nd <- ref$nodes[ref$nodes$model == m & ref$nodes$level == e$match_level[1] &
                    ref$nodes$node == e$match_name[1] & ref$nodes$length_type %in% "body_length", ]
  lt <- ref$lwr_types[ref$lwr_types$model == m & ref$lwr_types$length_type == "body_length", ]
  b <- nd$path_slope
  mass <- function(x) exp(nd$path_intercept + b * (log(x) - lt$centre))
  # The mean of the masses of the individuals is the mass at the length of
  # their mean mass, which the estimate uses.
  expect_equal(e$ind_median[1], sum(n * mass(L)) / sum(n))
  expect_equal(e$length_mm[1], (sum(n * L^b) / sum(n))^(1 / b))
  xc <- log(e$length_mm[1]) - lt$centre
  v <- nd$path_var + 2 * xc * nd$path_cov + xc^2 * nd$path_var_slope +
    comp(m, "species", "intercept") + 2 * xc * comp(m, "species", "intercept_slope") +
    xc^2 * comp(m, "species", "slope") + comp(m, "residual", "intercept") +
    comp(m, "equation", "intercept") + comp(m, "study", "intercept")
  expect_equal(e$ind_upper[1], exp(log(e$ind_median[1]) + stats::qnorm(0.975) * sqrt(v)))
  expect_true(is.na(e$corrections[1]))
  expect_match(e$note[1], "mean mass of 12 individuals measured \\(6-14 mm\\)")
  # Individuals beyond the range of the node move it up, as a mean length.
  wide <- transform(ln[ln$site == "a", ], length = c(0.02, 8, 10, 12, 14))
  e2 <- suppressWarnings(as.data.frame(estimate_biomass(ab[1, ], chiro, wide, bias_correction = "none")))
  expect_equal(e2$status, "not_estimated")
  expect_match(e2$reason, "lengths of individuals 0.02-14 mm outside")
})

test_that("lengths in taxa are deprecated but still used, after those of lengths", {
  ab <- data.frame(site = c("a", "b"), taxon_id = "c1", abundance = 1, per = "m2")
  old <- transform(chiro, length_mean = 10, length_unit = "mm", length_type = "body_length")
  expect_warning(e_old <- estimate_biomass(ab, old, bias_correction = "none"),
                 class = "bmassa_warning_deprecated")
  e_new <- estimate_biomass(ab, chiro, len(chiro, 10, "mm", "body_length"), bias_correction = "none")
  expect_equal(as.data.frame(e_old)$ind_median, as.data.frame(e_new)$ind_median)
  ln <- data.frame(taxon_id = "c1", site = "b", length_mean = 8, length_unit = "mm",
                   length_type = "body_length")
  e2 <- suppressWarnings(as.data.frame(estimate_biomass(ab, old, ln, bias_correction = "none")))
  expect_equal(e2$length_mm, c(10, 8))
})

test_that("all the rows come out in one currency and unit", {
  ex <- example_survey
  e <- suppressWarnings(as.data.frame(estimate_biomass(ex$abundance, ex$taxa, ex$lengths,
                                                       bias_correction = "lognormal")))
  ok <- e$status == "estimated"
  expect_equal(unique(e$currency[ok]), "wet_mass")
  expect_equal(unique(e$unit[ok]), "mg")
  expect_true(all(is.na(e$converted_from[ok & e$broad_group == "fish"])))
  expect_true(all(e$converted_from[ok & e$broad_group == "phytoplankton"] == "biovolume"))
  # Phytoplankton biovolume at the density of water: 1 um3 = 1 pg = 1e-9 mg.
  pg <- suppressWarnings(as.data.frame(estimate_biomass(ex$abundance, ex$taxa, ex$lengths,
                                                        unit = "pg", bias_correction = "lognormal")))
  phy <- e$broad_group == "phytoplankton" & ok
  expect_equal(pg$ind_median[phy], 1e9 * e$ind_median[phy])
  # Dry mass: phytoplankton cannot be converted, and says why.
  expect_warning(d <- as.data.frame(estimate_biomass(ex$abundance, ex$taxa, ex$lengths,
                                                     currency = "dry_mass", bias_correction = "lognormal")),
                 class = "bmassa_warning_not_estimable")
  expect_true(all(d$status[d$broad_group == "phytoplankton"] == "not_estimated"))
  expect_match(d$reason[d$broad_group == "phytoplankton"][1], "no conversion from biovolume to dry mass")
  expect_equal(unique(d$unit[d$status == "estimated"]), "mg")
  # convert_currency() starts again from the currency of the models.
  est <- suppressWarnings(estimate_biomass(ex$abundance, ex$taxa, ex$lengths, bias_correction = "none"))
  back <- suppressWarnings(as.data.frame(convert_currency(convert_currency(est, to = "carbon"),
                                                          to = "wet_mass")))
  expect_equal(back$ind_median, as.data.frame(est)$ind_median)
  expect_error(estimate_biomass(ex$abundance, ex$taxa, currency = "biovolume"), class = "bmassa_error_input")
  expect_error(estimate_biomass(ex$abundance, ex$taxa, unit = "kJ"), class = "bmassa_error_input")
})

test_that("a site sampled on several dates has one row per date", {
  ab <- data.frame(site = "s", date = as.Date(c("2022-06-01", "2022-06-15")), taxon_id = "w1",
                   abundance = c(10, 20), per = "m2")
  e <- as.data.frame(estimate_biomass(ab, limno))
  expect_equal(nrow(e), 2L)
  expect_equal(e$biomass_mean[2], 2 * e$biomass_mean[1])
  s <- summarize_biomass(estimate_biomass(ab, limno), by = c("site", "date"), seed = 1)
  expect_equal(nrow(s), 2L)
  expect_error(estimate_biomass(rbind(ab, ab), limno), class = "bmassa_error_input")
})

test_that("the sex and the stage of a unit are given in taxa", {
  units <- rbind(transform(chiro, life_stage = "adult", sex = "female"),
                 transform(chiro, taxon_id = "c2", life_stage = NA, sex = "male"))
  ab <- data.frame(site = "s", taxon_id = c("c1", "c2"), abundance = 1, per = "m2")
  e <- as.data.frame(estimate_biomass(ab, units))
  expect_equal(e$stage_used, c("adult", "adult"))
  expect_equal(e$sex, c("female", "male"))
  cls <- "bmassa_error_input"
  expect_error(estimate_biomass(ab[1, ], transform(chiro, life_stage = "female")), "sex", class = cls)
  expect_error(estimate_biomass(ab[1, ], transform(chiro, sex = "hermaphrodite")), class = cls)
  expect_error(estimate_biomass(transform(ab[1, ], life_stage = "adult"), chiro), "taxa",
               class = cls)
  expect_error(estimate_biomass(ab[1, ], chiro, data.frame(taxon_id = "c1", site = "s",
                                                           life_stage = "larva", length_mean = 5,
                                                           length_unit = "mm",
                                                           length_type = "body_length"),
                                bias_correction = "none"), "life_stage", class = cls)
})

test_that("the result keeps the rows and their order; zero rows give zero rows", {
  ex <- example_survey
  e <- suppressWarnings(as.data.frame(estimate_biomass(ex$abundance, ex$taxa, ex$lengths,
                                                       bias_correction = "none")))
  expect_equal(e$site, ex$abundance$site)
  expect_equal(e$taxon_id, ex$abundance$taxon_id)
  o <- rev(seq_len(nrow(ex$abundance)))
  e2 <- suppressWarnings(as.data.frame(estimate_biomass(ex$abundance[o, ], ex$taxa, ex$lengths,
                                                        bias_correction = "none")))
  expect_equal(e2$ind_median, e$ind_median[o])
  z <- as.data.frame(estimate_biomass(ex$abundance[0, ], ex$taxa, ex$lengths, bias_correction = "none"))
  expect_equal(nrow(z), 0L)
  expect_identical(names(z), names(e))
})

test_that("factors and tibbles give the same result", {
  ex <- example_survey
  base <- suppressWarnings(as.data.frame(estimate_biomass(ex$abundance, ex$taxa, ex$lengths,
                                                          bias_correction = "none")))
  tx_f <- ex$taxa
  chr <- vapply(tx_f, is.character, logical(1))
  tx_f[chr] <- lapply(tx_f[chr], factor)
  ab_f <- ex$abundance
  ab_f$site <- factor(ab_f$site)
  f <- suppressWarnings(as.data.frame(estimate_biomass(ab_f, tx_f, ex$lengths, bias_correction = "none")))
  expect_equal(f$ind_median, base$ind_median)
  # A tibble (without depending on tibble): a data frame with its classes.
  as_tbl <- function(d) structure(d, class = c("tbl_df", "tbl", "data.frame"))
  t <- suppressWarnings(as.data.frame(estimate_biomass(as_tbl(ex$abundance), as_tbl(ex$taxa),
                                                       as_tbl(ex$lengths), bias_correction = "none")))
  expect_equal(t$ind_median, base$ind_median)
})

test_that("columns of abundance that the result would overwrite are an error", {
  # Regression: the list of output columns of 0.1.0 (same object name)
  # replaced that of estimate_biomass(), so these went unnoticed.
  ab <- data.frame(site = "s", taxon_id = "w1", abundance = 1, per = "m2")
  for (col in c("biomass_mean", "unit", "match_level", "conversion")) {
    bad <- ab
    bad[[col]] <- 1
    expect_error(estimate_biomass(bad, limno), col, class = "bmassa_error_input")
  }
})

test_that("input errors are classed and informative", {
  ab <- data.frame(site = "s", taxon_id = "w1", abundance = 1, per = "m2")
  cls <- "bmassa_error_input"
  expect_error(estimate_biomass(ab[, -4], limno), class = cls)
  expect_error(estimate_biomass(ab, rbind(limno, limno)), class = cls)
  expect_error(estimate_biomass(transform(ab, taxon_id = "zz"), limno), class = cls)
  expect_error(estimate_biomass(ab, transform(limno, broad_group = "benthos")), class = cls)
  expect_error(estimate_biomass(transform(ab, abundance = -1), limno), class = cls)
  expect_error(estimate_biomass(transform(ab, count_unit = "cell"), limno), class = cls)
  expect_error(estimate_biomass(ab, transform(limno, count_unit = "cell")), class = cls)
  expect_error(estimate_biomass(rbind(ab, ab), limno), class = cls)
  expect_error(estimate_biomass(ab, limno, level = 1.5), class = cls)
  phy <- taxon("p1", "phytoplankton", c(kingdom = "Chromista", genus = "Cyclotella"))
  expect_error(estimate_biomass(transform(ab, taxon_id = "p1"), phy), class = cls)
  ln5 <- len(limno, 5, "mm", "body_length")
  expect_error(estimate_biomass(ab, limno, ln5), class = cls)
  expect_error(estimate_biomass(ab, limno, transform(ln5, length_unit = "inch"), bias_correction = "none"),
               class = cls)
  expect_error(estimate_biomass(ab, limno, transform(ln5, length_sd = 1, length_cv = 0.2),
                                bias_correction = "none"), class = cls)
  expect_error(estimate_biomass(ab, limno, transform(ln5, length_unit = NA), bias_correction = "none"),
               "length_unit", class = cls)
  expect_error(estimate_biomass(ab, transform(limno, length_unit = "mm"), bias_correction = "none"),
               "length_mean", class = cls)
  ln <- data.frame(taxon_id = "w1", site = "s", length_mean = 5, length_unit = "mm",
                   length_type = "body_length")
  expect_error(estimate_biomass(ab, limno, ln[, c("taxon_id", "site", "length_unit", "length_type")],
                                bias_correction = "none"), "length_mean", class = cls)
  expect_error(estimate_biomass(ab, limno, transform(ln, date = 1), bias_correction = "none"),
               "date", class = cls)
  expect_error(estimate_biomass(ab, limno, rbind(ln, ln), bias_correction = "none"), class = cls)
  # Individual lengths: not with a mean, no spread, n only for them, one kind
  # and one type per sample, and a date needs a site.
  ind <- data.frame(taxon_id = "w1", site = "s", length = c(4, 6), length_unit = "mm",
                    length_type = "body_length")
  expect_error(estimate_biomass(ab, limno, transform(ind, length_mean = 5), bias_correction = "none"),
               "both", class = cls)
  expect_error(estimate_biomass(ab, limno, transform(ind, length_sd = 1), bias_correction = "none"),
               "spread", class = cls)
  expect_error(estimate_biomass(ab, limno, transform(ln, n = 3), bias_correction = "none"),
               "`n`", class = cls)
  expect_error(estimate_biomass(ab, limno, transform(ind, n = 0), bias_correction = "none"),
               class = cls)
  mixed <- rbind(transform(ind, length_mean = NA), transform(ln, length = NA)[, c(names(ind), "length_mean")])
  expect_error(estimate_biomass(ab, limno, mixed, bias_correction = "none"), "one or the other",
               class = cls)
  expect_error(estimate_biomass(ab, limno, transform(ind, length_type = c("body_length", "shell_length")),
                                bias_correction = "none"), "length type", class = cls)
  ab_d <- transform(ab, date = 1)
  expect_error(estimate_biomass(ab_d, limno, transform(ind, site = NA, date = 1), bias_correction = "none"),
               "without `site`", class = cls)
  expect_error(estimate_biomass(ab, limno, ind[, -3], bias_correction = "none"), "length_mean",
               class = cls)
  grouped <- structure(ab, class = c("grouped_df", "tbl_df", "tbl", "data.frame"))
  expect_error(estimate_biomass(grouped, limno), class = cls)
  expect_snapshot(estimate_biomass(ab, limno, ln5), error = TRUE)
})

test_that("homonyms are told apart by the higher ranks", {
  cand <- data.frame(level = c("genus", "genus"), node = c("Xus", "Xus"),
                     path = c("A|P1|C1|O1|F1|Xus", "A|P2|C2|O2|F2|Xus"), stringsAsFactors = FALSE)
  lin <- c(kingdom = "A", phylum = "P2", class = "C2", order = NA, family = "F2", genus = "Xus",
           species = NA)
  m <- match_lineage(lin, cand, c("phylum", "class", "order", "family", "genus", "species"))
  expect_equal(m$row, 2L)
  lin2 <- c(kingdom = "A", phylum = NA, class = NA, order = NA, family = NA, genus = "Xus", species = NA)
  m2 <- match_lineage(lin2, cand, c("genus", "species"))
  expect_null(m2$row)
  expect_match(m2$note, "ambiguous")
})

test_that("print() summarises the estimate", {
  ex <- example_survey
  est <- suppressWarnings(estimate_biomass(ex$abundance, ex$taxa, ex$lengths, bias_correction = "lognormal"))
  expect_snapshot(print(est))
})
