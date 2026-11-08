# check_coverage() -------------------------------------------------------------------

ex <- example_survey

test_that("one row per taxon unit, plus one per lengths of a site or sample", {
  cov <- check_coverage(ex$taxa, ex$lengths)
  d <- as.data.frame(cov)
  # The chub in two sites and the Gammarus measured in one (its individuals
  # are one row).
  expect_equal(nrow(d), nrow(ex$taxa) + 3L)
  expect_s3_class(cov, "bmassa_coverage")
  expect_false(inherits(d, "bmassa_coverage"))
  # A fish without length: mean mass of a population of all sizes.
  expect_equal(d$status[d$taxon_id == "t18" & d$method == "typical"], "estimable")
  expect_true(all(d$factor[d$status == "estimable"] > 1))
  # Units with a length are checked with it; the lengths of a sample say where.
  expect_equal(d$method[d$taxon_id == "t08"], "length")
  expect_equal(d$site[d$taxon_id == "t17"], c(NA, "L1", "L2"))
  expect_equal(d$count_unit[d$broad_group == "phytoplankton"], rep("cell", 5))
  gm <- d[d$taxon_id == "t13", ]
  expect_equal(gm$method, c("typical", "length"))
  expect_equal(gm$length_basis, c(NA, "individuals"))
  expect_equal(gm$length_n, c(NA, 12))
  # The frog: adults with a length, larvae not estimable without one.
  expect_equal(d$status[d$taxon_id %in% c("t19", "t20")], c("estimable", "not_estimable"))
})

test_that("the factor is the width of the interval of estimate_biomass()", {
  cov <- as.data.frame(check_coverage(ex$taxa))
  ab <- data.frame(site = "s", taxon_id = "t15", abundance = 1, per = "m2")
  # In the currency of the model (the conversion to another one widens it).
  e <- as.data.frame(estimate_biomass(ab, ex$taxa, currency = "dry_mass"))
  expect_equal(cov$factor[cov$taxon_id == "t15"], e$ind_upper / e$ind_median)
  cov90 <- as.data.frame(check_coverage(ex$taxa, level = 0.9))
  expect_lt(cov90$factor[cov90$taxon_id == "t15"], cov$factor[cov$taxon_id == "t15"])
})

test_that("print() summarises by broad group", {
  expect_snapshot(print(check_coverage(ex$taxa, ex$lengths)))
})
