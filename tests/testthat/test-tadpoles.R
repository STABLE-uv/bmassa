# Tadpoles (version 0.6.0): length-mass relation of total length to dry mass,
# typical lengths for tadpoles without lengths, and their conversions --------

ref <- get_reference()
taxon <- function(id, lineage, life_stage = "larva") {
  lin <- stats::setNames(as.list(rep(NA_character_, 7)), c("kingdom", "phylum", "class", "order",
                                                            "family", "genus", "species"))
  lin[names(lineage)] <- lineage
  as.data.frame(c(list(taxon_id = id, broad_group = "amphibian"), lin, list(life_stage = life_stage)),
                stringsAsFactors = FALSE)
}
anura <- c(kingdom = "Animalia", phylum = "Chordata", class = "Amphibia", order = "Anura")
leopard <- taxon("t1", c(anura, family = "Ranidae", genus = "Lithobates",
                         species = "Lithobates sphenocephalus"))
toad <- taxon("t2", c(anura, family = "Bufonidae", genus = "Bufo", species = "Bufo spinosus"))
newt <- taxon("t3", c(kingdom = "Animalia", phylum = "Chordata", class = "Amphibia",
                      order = "Caudata", family = "Salamandridae", genus = "Triturus",
                      species = "Triturus marmoratus"))
comp <- function(m, l) {
  v <- ref$components$value[ref$components$model == m & ref$components$level == l &
                              ref$components$parameter == "intercept"]
  if (length(v)) v else 0
}
node_of <- function(m, level, name) {
  k <- which(ref$nodes$model == m & ref$nodes$level == level & ref$nodes$node == name)
  stopifnot(length(k) == 1L)
  ref$nodes[k, ]
}
tau2_of <- function(e) 2 * log(e$ind_mean / e$ind_median)
lwr_m <- "amphibian_larva_lwr"
len_m <- "amphibian_larva_length"
lt <- ref$lwr_types[ref$lwr_types$model == lwr_m & ref$lwr_types$length_type == "total_length", ]
ab1 <- function(id) data.frame(site = "s", taxon_id = id, abundance = 10, per = "m2")

test_that("a tadpole with its total length gets the relation of its species (hand calculation)", {
  ln <- data.frame(taxon_id = "t1", length_mean = 40, length_unit = "mm", length_type = "total_length")
  e <- as.data.frame(estimate_biomass(ab1("t1"), leopard, ln, currency = "dry_mass",
                                      bias_correction = "lognormal"))
  nd <- node_of(lwr_m, "species", "Lithobates sphenocephalus")
  xc <- log(40) - lt$centre
  eta <- nd$path_intercept + nd$path_slope * xc + lt$sigma_ind^2 / 2
  tau2 <- nd$path_var + 2 * xc * nd$path_cov + xc^2 * nd$path_var_slope +
    comp(lwr_m, "residual") + comp(lwr_m, "equation")
  expect_equal(e$model, lwr_m)
  expect_equal(e$method, "length")
  expect_equal(e$stage_used, "larva")
  expect_equal(e$match_level, "species")
  expect_equal(c(e$currency, e$unit), c("dry_mass", "mg"))
  expect_equal(log(e$ind_median), eta, tolerance = 1e-12)
  expect_equal(tau2_of(e), tau2, tolerance = 1e-10)
  expect_equal(e$biomass_mean, 10 * e$ind_mean)
  # The relation is for tadpoles before the metamorphic climax: a 40 mm
  # leopard frog tadpole weighs some tens of mg of dry mass.
  expect_true(e$ind_median > 10 && e$ind_median < 100)
})

test_that("a tadpole without lengths gets the mass at the typical length of its taxon (hand calculation)", {
  e <- as.data.frame(estimate_biomass(ab1("t1"), leopard, currency = "dry_mass",
                                      bias_correction = "lognormal"))
  tm <- ref$models[ref$models$model == len_m, ]
  tn <- node_of(len_m, "species", "Lithobates sphenocephalus")
  L <- exp(tn$path_intercept)
  nd <- node_of(lwr_m, "species", "Lithobates sphenocephalus")
  expect_true(L >= nd$length_min && L <= nd$length_max)
  b <- nd$path_slope
  xc <- log(L) - lt$centre
  eta <- nd$path_intercept + b * xc + b * (b - 1) / 2 * log1p(tm$length_cv^2) + lt$sigma_ind^2 / 2
  tau2 <- nd$path_var + 2 * xc * nd$path_cov + xc^2 * nd$path_var_slope +
    comp(lwr_m, "residual") + comp(lwr_m, "equation") +
    b^2 * (tn$path_var + comp(len_m, "residual") + comp(len_m, "population"))
  expect_equal(e$method, "typical_length")
  expect_equal(e$model, lwr_m)
  expect_equal(e$length_mm, L)
  expect_equal(e$length_type, "total_length")
  expect_equal(e$length_basis, "typical")
  expect_equal(log(e$ind_median), eta, tolerance = 1e-12)
  expect_equal(tau2_of(e), tau2, tolerance = 1e-10)
  expect_match(e$corrections, "jensen; lognormal")
  expect_match(e$uncertainty, "typical_length")
  expect_match(e$note, "typical length of the tadpoles of species Lithobates sphenocephalus")
  # Without a stage, an anuran is a tadpole too (assumed), with the same mass.
  e0 <- as.data.frame(estimate_biomass(ab1("t1"), transform(leopard, life_stage = NA),
                                       currency = "dry_mass", bias_correction = "lognormal"))
  expect_equal(e0$stage_used, "larva (assumed)")
  expect_equal(e0$ind_median, e$ind_median)
  expect_match(e0$note, "taken as larva")
  # The typical length needs the back-transformation to be chosen.
  expect_error(estimate_biomass(ab1("t1"), leopard), "typical lengths",
               class = "bmassa_error_input")
})

test_that("tadpoles of families absent from the reference come from all its tadpoles", {
  e <- as.data.frame(estimate_biomass(ab1("t2"), toad, currency = "dry_mass",
                                      bias_correction = "none"))
  tm <- ref$models[ref$models$model == len_m, ]
  root_l <- ref$nodes[ref$nodes$model == len_m & ref$nodes$level == "root", ]
  root_m <- ref$nodes[ref$nodes$model == lwr_m & ref$nodes$level == "root", ]
  L <- exp(root_l$path_intercept)
  b <- root_m$path_slope
  xc <- log(L) - lt$centre
  lv_l <- strsplit(tm$levels, ",")[[1]]
  lv_m <- strsplit(ref$models$levels[ref$models$model == lwr_m], ",")[[1]]
  miss_m <- sum(vapply(lv_m, function(l) {
    v <- ref$components$value[ref$components$model == lwr_m & ref$components$level == l]
    p <- ref$components$parameter[ref$components$model == lwr_m & ref$components$level == l]
    sum(v[p == "intercept"]) + 2 * xc * sum(v[p == "intercept_slope"]) + xc^2 * sum(v[p == "slope"])
  }, numeric(1)))
  miss_l <- sum(vapply(lv_l, function(l) comp(len_m, l), numeric(1)))
  eta <- root_m$path_intercept + b * xc + b * (b - 1) / 2 * log1p(tm$length_cv^2)
  tau2 <- root_m$path_var + 2 * xc * root_m$path_cov + xc^2 * root_m$path_var_slope + miss_m +
    comp(lwr_m, "residual") + comp(lwr_m, "equation") +
    b^2 * (root_l$path_var + miss_l + comp(len_m, "residual") + comp(len_m, "population"))
  expect_equal(e$match_level, "root")
  expect_equal(e$length_mm, L)
  expect_equal(log(e$ind_median), eta, tolerance = 1e-12)
  expect_equal(tau2_of(e), tau2, tolerance = 1e-10)
  expect_match(e$note, "all the tadpoles of the reference")
})

test_that("larvae of newts, juvenile frogs and tadpoles measured otherwise are not estimated", {
  tx <- rbind(newt, transform(leopard, taxon_id = "j1", life_stage = "juvenile"),
              transform(leopard, taxon_id = "s1"))
  ab <- data.frame(site = "s", taxon_id = c("t3", "j1", "s1"), abundance = 1, per = "m2")
  ln <- data.frame(taxon_id = "s1", length_mean = 20, length_unit = "mm",
                   length_type = "snout_vent_length")
  expect_warning(e <- as.data.frame(estimate_biomass(ab, tx, ln, bias_correction = "none")),
                 class = "bmassa_warning_not_estimable")
  expect_equal(e$status, rep("not_estimated", 3))
  expect_match(e$reason[1], "tadpoles\\) only, not of salamanders")
  expect_match(e$reason[2], "given as \"juvenile\" have no reference masses")
  expect_match(e$reason[3], "estimated from their total length")
})

test_that("tadpoles go to carbon, nitrogen and phosphorus, not to wet mass or energy", {
  est <- estimate_biomass(ab1("t1"), leopard, currency = "dry_mass", bias_correction = "lognormal")
  e0 <- as.data.frame(est)
  for (el in c("carbon", "nitrogen", "phosphorus")) {
    e1 <- as.data.frame(convert_currency(est, to = el, unit = "ug"))
    m <- paste0("amphibian_larva dry_mass_to_", el)
    k <- which(ref$conversion_nodes$model == m & ref$conversion_nodes$level == "species" &
                 ref$conversion_nodes$node == "Lithobates sphenocephalus")
    expect_length(k, 1L)
    n <- ref$conversion_nodes[k, ]
    cc <- ref$conversion_components
    v_f <- n$path_var + sum(cc$value[cc$model == m & cc$level == "residual"])
    # Dry mass in mg to the element in ug: x 1000.
    expect_equal(log(e1$ind_median) - log(e0$ind_median * 1000), n$path_intercept, tolerance = 1e-12)
    expect_equal(tau2_of(e1) - tau2_of(e0), v_f, tolerance = 1e-10)
    expect_match(e1$conversion, paste0("dry_mass->", el, ": ", m, ", species"))
    expect_match(e1$conversion_sources, "knapp2021")
  }
  # Carbon, nitrogen and phosphorus of the tadpoles of each family: 40-50 %,
  # 7-12 % and 0.8-1.4 % of the dry mass, P below that of fish (with bones).
  share <- function(el, level = "family") {
    k <- ref$conversion_nodes$model == paste0("amphibian_larva dry_mass_to_", el) &
      ref$conversion_nodes$level == level
    exp(ref$conversion_nodes$path_intercept[k])
  }
  expect_true(all(share("carbon") > 0.40 & share("carbon") < 0.50))
  expect_true(all(share("nitrogen") > 0.07 & share("nitrogen") < 0.12))
  expect_true(all(share("phosphorus") > 0.008 & share("phosphorus") < 0.014))
  fish_p <- ref$conversion_nodes$model == "fish dry_mass_to_phosphorus" &
    ref$conversion_nodes$level == "root"
  expect_lt(share("phosphorus", "root"), exp(ref$conversion_nodes$path_intercept[fish_p]))
  for (to in c("wet_mass", "energy", "afdm")) {
    expect_warning(cv <- convert_currency(est, to = to), class = "bmassa_warning_not_estimable")
    e1 <- as.data.frame(cv)
    expect_true(is.na(e1$ind_median))
    expect_match(e1$reason, paste("no conversion from dry mass to", gsub("_", " ", to), "for tadpoles"))
  }
  # In the default currency (wet mass) they are not estimated, with the reason.
  expect_warning(w <- as.data.frame(estimate_biomass(ab1("t1"), leopard, bias_correction = "none")),
                 class = "bmassa_warning_not_estimable")
  expect_match(w$reason, "no conversion from dry mass to wet mass for tadpoles")
})

test_that("draws of tadpoles without lengths have the mean and variance of their estimate", {
  tx <- rbind(leopard, toad)
  ab <- data.frame(site = rep(c("a", "b"), each = 2), taxon_id = c("t1", "t2"), abundance = 1,
                   per = "m2")
  est <- estimate_biomass(ab, tx, currency = "dry_mass", bias_correction = "lognormal")
  e <- as.data.frame(est)
  d <- withr::with_seed(3, draw_log_mass(est, seq_len(nrow(e)), 40000, "shared"))
  expect_equal(colMeans(d), log(e$ind_median), tolerance = 0.02)
  expect_equal(apply(d, 2, stats::var), tau2_of(e), tolerance = 0.05)
  # The two sites of a taxon share its population (and typical length).
  expect_equal(d[, 1], d[, 3])
  s <- summarize_biomass(est, by = "site", seed = 4, n_draws = 4000)
  expect_equal(s$biomass_mean, as.vector(tapply(e$biomass_mean, e$site, sum)))
  # check_coverage() gives the width of the same interval.
  cov <- as.data.frame(check_coverage(tx))
  expect_equal(cov$method, c("typical_length", "typical_length"))
  expect_equal(cov$factor, e$ind_upper[1:2] / e$ind_median[1:2], tolerance = 1e-10)
})
