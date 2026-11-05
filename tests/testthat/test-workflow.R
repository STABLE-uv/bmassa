# The flow of dev/ejemplos_uso.R (also in the README and the vignette).

test_that("taxa and lengths -> coverage -> estimates -> carbon -> totals -> sources", {
  ab <- example_survey$abundance
  tx <- example_survey$taxa
  ln <- example_survey$lengths
  cov <- check_coverage(tx, ln)
  expect_s3_class(cov, "bmassa_coverage")
  expect_warning(est <- estimate_biomass(ab, tx, ln, bias_correction = "lognormal"),
                 class = "bmassa_warning_not_estimable")
  e <- as.data.frame(est)
  expect_equal(nrow(e), nrow(ab))
  expect_equal(sum(e$status == "not_estimated"), 5L)
  # Fish without lengths, amphibian larvae and a missing abundance.
  expect_setequal(e$taxon_id[e$status == "not_estimated"], c("t18", "t20", "t08"))
  # Carbon in a common unit; amphibians are not converted.
  expect_warning(est_c <- convert_currency(est, to = "carbon", unit = "ug"),
                 class = "bmassa_warning_not_estimable")
  ec <- as.data.frame(est_c)
  expect_true(all(ec$unit[!is.na(ec$unit)] == "ug"))
  expect_true(all(is.na(ec$ind_median[ec$broad_group == "amphibian"])))
  tot <- summarize_biomass(est_c, by = c("site", "broad_group"), na_rm = TRUE, seed = 1)
  expect_equal(nrow(tot), 10L)
  ok <- !is.na(tot$currency)
  expect_true(all(tot$biomass_lower[ok] < tot$biomass_upper[ok]))
  expect_true(all(tot$currency[ok] == "carbon"))
  # Joint draws for other packages, all per m2, with one population per habitat.
  sites <- data.frame(site = c("L1", "L2"), depth = c(1.5, 0.8), depth_unit = "m")
  ab_h <- convert_abundance(ab, to = "m2", sites = sites)
  expect_true(all(ab_h$per == "m2"))
  ab_h$habitat <- ifelse(ab_h$site == "L1", "lake", "spring")
  est_w <- suppressWarnings(convert_currency(
    suppressWarnings(estimate_biomass(ab_h, tx, ln, bias_correction = "lognormal")),
    to = "wet_mass", unit = "mg"))
  sim <- simulate(est_w, nsim = 100, seed = 1, population_by = "habitat")
  d <- as.data.frame(sim, draws = 1:2)
  expect_equal(nrow(d), 2L * nrow(ab))
  expect_true(all(d$unit[!is.na(d$unit)] == "mg"))
  expect_true(all(d$biomass_unit[!is.na(d$biomass_unit)] == "mg/m2"))
  expect_equal(d$biomass, d$abundance * d$ind_mass)
  src <- cite_sources(est_c)
  expect_true(all(c("kremer2014", "vanleeuwen2025", "arranz2022", "brey2010", "mdl2000") %in%
                    src$source_id))
  out <- tempfile(fileext = ".csv")
  utils::write.csv(ec, out, row.names = FALSE)
  expect_equal(nrow(utils::read.csv(out)), nrow(ab))
})
