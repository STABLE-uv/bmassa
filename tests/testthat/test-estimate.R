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
# A taxon unit with its mean length.
with_length <- function(tx, mean, unit, type, ...) {
  tx$length_mean <- mean
  tx$length_unit <- unit
  tx$length_type <- type
  extra <- list(...)
  for (k in names(extra)) tx[[k]] <- extra[[k]]
  tx
}
comp <- function(m, l, p) {
  v <- ref$components$value[ref$components$model == m & ref$components$level == l &
                              ref$components$parameter == p]
  if (length(v)) v else 0
}

test_that("a species matched in the reference is predicted from its path (hand calculation)", {
  ab <- data.frame(site = "s", taxon_id = "w1", abundance = 10, per = "m2")
  e <- as.data.frame(estimate_biomass(ab, limno))
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
  e <- as.data.frame(estimate_biomass(ab, with_length(chiro, 10, "mm", "body_length"),
                                      bias_correction = "none"))
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
  e <- as.data.frame(estimate_biomass(ab, with_length(trout, 200, "mm", "fork_length"),
                                      bias_correction = "none"))
  w <- exp(-11.46273 + 3.021065 * log(200))
  expect_equal(w, 94.1, tolerance = 0.001)
  expect_true(e$ind_lower < w && w < e$ind_upper)
  expect_equal(e$unit, "g")
})

test_that("the corrections for mean lengths follow their formulas", {
  ab <- data.frame(site = "s", taxon_id = "c1", abundance = 1, per = "m2")
  tx0 <- with_length(chiro, 10, "mm", "body_length")
  tx1 <- transform(tx0, length_cv = 0.3)
  e0 <- as.data.frame(estimate_biomass(ab, tx0, bias_correction = "none"))
  e1 <- as.data.frame(estimate_biomass(ab, tx1, bias_correction = "none"))
  e2 <- as.data.frame(estimate_biomass(ab, tx1, bias_correction = "lognormal"))
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
  expect_equal(transform(tx1, length_cv = NULL, length_sd = 3)$length_sd / 10, 0.3)
  e3 <- as.data.frame(estimate_biomass(ab, transform(tx0, length_sd = 3), bias_correction = "none"))
  expect_equal(e3$ind_median, e1$ind_median)
})

test_that("units: the same length in mm and cm, and abundances scale the biomass", {
  ab <- data.frame(site = "s", taxon_id = "f1", abundance = c(1), per = "m2")
  tx_mm <- with_length(trout, 150, "mm", "fork_length")
  tx_cm <- with_length(trout, 15, "cm", "fork_length")
  a <- as.data.frame(estimate_biomass(ab, tx_mm, bias_correction = "none"))
  b <- as.data.frame(estimate_biomass(ab, tx_cm, bias_correction = "none"))
  expect_equal(a$ind_median, b$ind_median)
  ab10 <- transform(ab, abundance = 10)
  c10 <- as.data.frame(estimate_biomass(ab10, tx_mm, bias_correction = "none"))
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
  tiny <- with_length(chiro, 0.02, "mm", "body_length")
  e <- suppressWarnings(as.data.frame(estimate_biomass(ab, tiny, bias_correction = "none")))
  expect_equal(e$status, "not_estimated")
  expect_match(e$reason, "outside the lengths")
  e2 <- as.data.frame(estimate_biomass(ab, tiny, bias_correction = "none", extrapolate = TRUE))
  expect_equal(e2$length_status, "outside")
  expect_equal(e2$match_level, "genus")
  ex <- as.data.frame(estimate_biomass(example_survey$abundance[12, ], example_survey$taxa,
                                       example_survey$lengths, bias_correction = "none"))
  expect_equal(ex$match_name, "Coenagrionidae")
  expect_match(ex$note, "estimated from Coenagrionidae")
})

test_that("fish without length: mean mass of a population (hand calculation)", {
  ab <- data.frame(site = "s", taxon_id = c("f1", "f2"), abundance = 2, per = "m2")
  mudminnow <- taxon("f2", "fish",
                     c(kingdom = "Animalia", phylum = "Chordata", class = "Teleostei",
                       order = "Esociformes", family = "Umbridae", genus = "Umbra",
                       species = "Umbra limi"))
  e <- as.data.frame(estimate_biomass(ab, rbind(trout, mudminnow)))
  expect_equal(e$status, c("estimated", "estimated"))
  expect_equal(e$method, c("typical", "typical"))
  expect_equal(e$unit, c("g", "g"))
  expect_match(e$note[1], "fish without length")
  # The trout from its node of the fish model: path, its variance, the
  # variation between populations and that of a new source.
  m <- "fish_individual"
  nd <- ref$nodes[ref$nodes$model == m & ref$nodes$node == e$match_name[1] &
                    ref$nodes$level == e$match_level[1], ]
  tau2 <- nd$path_var + comp(m, "residual", "intercept") + comp(m, "source", "intercept")
  lv <- strsplit(ref$models$levels[ref$models$model == m], ",")[[1]]
  miss <- lv[seq_along(lv) > match(e$match_level[1], lv)]
  for (l in miss) tau2 <- tau2 + comp(m, l, "intercept")
  expect_equal(e$ind_median[1], exp(nd$path_intercept))
  expect_equal(e$ind_upper[1], exp(nd$path_intercept + stats::qnorm(0.975) * sqrt(tau2)))
  expect_equal(e$biomass_median[1], 2 * e$ind_median[1])
  # An order absent from the reference: from all the fish, with the variances
  # of all the levels.
  expect_equal(e$match_level[2], "root")
  expect_match(e$note[2], "order not in the reference")
  expect_gt(e$ind_upper[2] / e$ind_median[2], e$ind_upper[1] / e$ind_median[1])
})

test_that("rows that cannot be estimated say why, with one warning", {
  ab <- data.frame(site = "s", taxon_id = c("f1", "w1"), abundance = 1, per = "m2")
  tx <- rbind(transform(trout, life_stage = "juvenile"), transform(limno, life_stage = NA))
  expect_warning(est <- estimate_biomass(ab, tx), class = "bmassa_warning_not_estimable")
  e <- as.data.frame(est)
  expect_equal(e$status, c("not_estimated", "estimated"))
  expect_match(e$reason[1], "whole populations")
  frog <- taxon("a1", "amphibian", c(kingdom = "Animalia", phylum = "Chordata", class = "Amphibia",
                                     order = "Anura", family = "Ranidae", genus = "Rana",
                                     species = "Rana temporaria"))
  frogs <- rbind(transform(frog, life_stage = NA), transform(frog, taxon_id = "a2", life_stage = "adult"))
  ab2 <- data.frame(site = "s", taxon_id = c("a1", "a2"), abundance = 1, per = "m2")
  e2 <- suppressWarnings(as.data.frame(estimate_biomass(ab2, frogs)))
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
  e <- suppressWarnings(as.data.frame(estimate_biomass(ab, units)))
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
  units <- rbind(with_length(transform(cal, life_stage = "nauplius"), 0.3, "mm", "body_length"),
                 with_length(transform(cal, taxon_id = "n2", life_stage = "copepodite"), 0.3, "mm",
                             "body_length"))
  ab <- data.frame(site = "s", taxon_id = c("n1", "n2"), abundance = 1, per = "L")
  e <- suppressWarnings(as.data.frame(estimate_biomass(ab, units, bias_correction = "none",
                                                       extrapolate = TRUE)))
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

test_that("the lengths of a sample take precedence over those of the unit", {
  ex <- example_survey
  e <- suppressWarnings(as.data.frame(estimate_biomass(ex$abundance, ex$taxa, ex$lengths,
                                                       bias_correction = "none")))
  sq <- e[e$taxon_id == "t17", ]
  expect_equal(sq$length_mm, c(110, 85))
  # The frog as adults (t19, with a length) and as larvae (t20, without).
  fr <- e[e$taxon_id %in% c("t19", "t20"), ]
  expect_equal(fr$method, c("length", "typical", "length", "typical"))
  expect_equal(fr$life_stage, c("adult", "larva", "adult", "larva"))
  # A length in the unit, replaced in one site and on one date.
  tx <- with_length(chiro, 10, "mm", "body_length")
  ab <- data.frame(site = c("a", "a", "b", "b"), date = c(1, 2, 1, 2), taxon_id = "c1",
                   abundance = 1, per = "m2")
  ln <- data.frame(taxon_id = "c1", site = c("a", "b"), date = c(NA, 2), length_mean = c(8, 12),
                   length_unit = "mm", length_type = "body_length")
  e2 <- as.data.frame(estimate_biomass(ab, tx, ln, bias_correction = "none"))
  expect_equal(e2$length_mm, c(8, 8, 10, 12))
  expect_identical(e2$date, ab$date)
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
  limno_l <- with_length(limno, 5, "mm", "body_length")
  expect_error(estimate_biomass(ab, limno_l), class = cls)
  expect_error(estimate_biomass(ab, transform(limno_l, length_unit = "inch"), bias_correction = "none"),
               class = cls)
  expect_error(estimate_biomass(ab, transform(limno_l, length_sd = 1, length_cv = 0.2),
                                bias_correction = "none"), class = cls)
  expect_error(estimate_biomass(ab, transform(limno_l, length_unit = NA), bias_correction = "none"),
               "length_unit", class = cls)
  expect_error(estimate_biomass(ab, transform(limno, length_unit = "mm"), bias_correction = "none"),
               "length_mean", class = cls)
  ln <- data.frame(taxon_id = "w1", site = "s", length_mean = 5, length_unit = "mm",
                   length_type = "body_length")
  expect_error(estimate_biomass(ab, limno, ln[, -2], bias_correction = "none"), class = cls)
  expect_error(estimate_biomass(ab, limno, transform(ln, date = 1), bias_correction = "none"),
               "date", class = cls)
  expect_error(estimate_biomass(ab, limno, rbind(ln, ln), bias_correction = "none"), class = cls)
  grouped <- structure(ab, class = c("grouped_df", "tbl_df", "tbl", "data.frame"))
  expect_error(estimate_biomass(grouped, limno), class = cls)
  expect_snapshot(estimate_biomass(ab, limno_l), error = TRUE)
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
