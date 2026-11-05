# summarize_biomass() for estimates: totals by simulation ----------------------------

ex <- example_survey
est <- suppressWarnings(estimate_biomass(ex$abundance, ex$taxa, ex$lengths, bias_correction = "lognormal"))

test_that("the draws of each row reproduce its analytic distribution", {
  e <- as.data.frame(est)
  ok <- which(!is.na(e$ind_median))
  ld <- withr::with_seed(11, draw_log_mass(est, ok, 20000, "shared"))
  tau2 <- 2 * log(e$ind_mean[ok] / e$ind_median[ok])
  # Monte Carlo error of a variance with 20000 draws: about 1 %.
  expect_equal(apply(ld, 2, stats::var), tau2, tolerance = 0.05)
  expect_equal(colMeans(ld), log(e$ind_median[ok]), tolerance = 0.05)
})

test_that("the total of one row has the mean of its biomass", {
  one <- est
  one$estimates <- est$estimates[est$estimates$taxon_id == "t16" & est$estimates$site == "L1", ]
  # The simulation reads the rows by position: keep the matching row.
  k <- which(est$estimates$taxon_id == "t16" & est$estimates$site == "L1")
  one$simulation$rows <- est$simulation$rows[k, ]
  one$simulation$missing_keys <- est$simulation$missing_keys[k]
  s <- summarize_biomass(one, by = "site", n_draws = 20000, seed = 3)
  expect_identical(s$biomass_mean, one$estimates$biomass_mean)
  expect_equal(s$biomass_median, one$estimates$biomass_median, tolerance = 0.03)
  expect_equal(s$biomass_lower, one$estimates$biomass_lower, tolerance = 0.03)
  expect_equal(s$biomass_upper, one$estimates$biomass_upper, tolerance = 0.03)
})

test_that("the mean of a total is the sum of the means of its rows, as that of the draws", {
  e <- as.data.frame(est)
  s <- summarize_biomass(est, by = c("site", "broad_group"), na_rm = TRUE, seed = 2)
  # Sums by hand, with the same keys.
  ok <- !is.na(e$biomass_mean)
  hand <- aggregate(biomass_mean ~ site + broad_group + currency, data = e[ok, ], FUN = sum)
  m <- merge(s, hand, by = c("site", "broad_group", "currency"), suffixes = c("", "_hand"))
  expect_equal(nrow(m), nrow(s))
  expect_equal(m$biomass_mean, m$biomass_mean_hand)
  # The mean of the simulated totals converges to it (fish and amphibians:
  # narrow intervals, small Monte Carlo error).
  narrow <- e$broad_group %in% c("fish", "amphibian") & ok
  k <- which(narrow)
  cid <- match(paste(e$site[k], e$broad_group[k]), unique(paste(e$site[k], e$broad_group[k])))
  d <- withr::with_seed(4, draw_totals(est, k, cid, max(cid), 20000, "shared"))
  sims <- colMeans(d)
  exact <- vapply(split(e$biomass_mean[k], cid), sum, numeric(1))
  expect_equal(unname(sims), unname(exact), tolerance = 0.01)
})

test_that("totals never add currencies, and missing rows make a total incomplete", {
  s <- summarize_biomass(est, by = "site", seed = 1)
  expect_setequal(unique(s$currency), c("biovolume", "dry_mass", "wet_mass"))
  expect_true(all(is.na(s$biomass_mean[s$status == "incomplete"])))
  s2 <- summarize_biomass(est, by = c("site", "broad_group"), na_rm = TRUE, seed = 1)
  expect_true(all(!is.na(s2$biomass_mean)))
  expect_equal(sum(s2$n_rows), nrow(as.data.frame(est)))
  expect_equal(s2$n_missing[s2$site == "L2" & s2$broad_group == "zooplankton"], 1L)
  expect_true(all(s2$biomass_lower <= s2$biomass_median & s2$biomass_median <= s2$biomass_upper))
})

test_that("a shared population effect widens sums over sites", {
  ab <- data.frame(site = c("a", "b"), taxon_id = "t15", abundance = 100, per = "m2", group = "all")
  e <- estimate_biomass(ab, ex$taxa)
  sh <- summarize_biomass(e, by = "group", n_draws = 4000, seed = 5)
  ind <- summarize_biomass(e, by = "group", n_draws = 4000, population_effect = "independent", seed = 5)
  expect_gt(log(sh$biomass_upper / sh$biomass_lower), log(ind$biomass_upper / ind$biomass_lower))
  one_sh <- summarize_biomass(e, by = "site", n_draws = 4000, seed = 5)
  one_ind <- summarize_biomass(e, by = "site", n_draws = 4000, population_effect = "independent", seed = 5)
  expect_equal(one_sh$biomass_upper / one_sh$biomass_lower,
               one_ind$biomass_upper / one_ind$biomass_lower, tolerance = 0.15)
})

test_that("a seed gives the same draws and restores the random state", {
  set.seed(42)
  before <- .Random.seed
  a <- summarize_biomass(est, by = "site", seed = 9)
  expect_identical(.Random.seed, before)
  b <- summarize_biomass(est, by = "site", seed = 9)
  expect_identical(a, b)
})

test_that("zero rows give zero rows, and arguments are checked", {
  z <- est
  z$estimates <- est$estimates[0, ]
  z$simulation$rows <- est$simulation$rows[0, ]
  z$simulation$missing_keys <- list()
  s <- summarize_biomass(z)
  expect_equal(nrow(s), 0L)
  cls <- "bmassa_error_input"
  expect_error(summarize_biomass(est, by = 1), class = cls)
  expect_error(summarize_biomass(est, by = "nowhere"), class = cls)
  expect_error(summarize_biomass(est, n_draws = 10), class = cls)
  expect_error(summarize_biomass(est, population_effect = "mixed"), class = cls)
  expect_error(summarize_biomass(est, by = "currency"), class = cls)
  expect_error(summarize_biomass(est, nadraws = 100), class = cls)
})

# simulate() for estimates: joint draws of every row --------------------------------

ab_h <- ex$abundance
ab_h$habitat <- ifelse(ab_h$site == "L1", "lake", "spring")
est_h <- suppressWarnings(estimate_biomass(ab_h, ex$taxa, ex$lengths, bias_correction = "lognormal"))

test_that("simulate() keeps the rows in order, with NA where there is no estimate", {
  sim <- simulate(est_h, nsim = 50, seed = 1)
  e <- as.data.frame(est_h)
  expect_s3_class(sim, "bmassa_simulation")
  expect_equal(dim(sim$ind_mass), c(nrow(e), 50L))
  expect_identical(sim$rows$site, e$site)
  expect_identical(sim$rows$taxon_id, e$taxon_id)
  expect_identical(sim$rows$habitat, e$habitat)
  expect_identical(is.na(sim$ind_mass[, 1]), is.na(e$ind_mean))
  expect_identical(sim$settings$nsim, 50L)
  expect_output(print(sim), "50 draws of the mean individual mass of 40 rows")
  d <- as.data.frame(sim)
  expect_equal(nrow(d), 50L * nrow(e))
  expect_identical(d$draw, rep(1:50, each = nrow(e)))
  expect_identical(d$ind_mass, as.vector(sim$ind_mass))
  expect_equal(d$biomass, d$abundance * d$ind_mass)
  expect_true(all(c("habitat", "broad_group", "currency", "unit", "biomass_unit", "status") %in%
                    names(d)))
  expect_false(any(c("ind_mean", "biomass_mean", "reason") %in% names(d)))
  # A missing abundance: the mass is drawn and the biomass is NA.
  k <- which(is.na(e$abundance) & !is.na(e$ind_mean))
  expect_length(k, 1L)
  dk <- d[d$site == e$site[k] & d$taxon_id == e$taxon_id[k], ]
  expect_true(all(!is.na(dk$ind_mass)) && all(is.na(dk$biomass)))
  d2 <- as.data.frame(sim, draws = c(3, 1))
  expect_identical(unique(d2$draw), c(3L, 1L))
  expect_identical(d2$ind_mass[d2$draw == 3L], unname(sim$ind_mass[, 3]))
})

test_that("the draws of each row reproduce its mean, median and interval", {
  sim <- simulate(est_h, nsim = 20000, seed = 11)
  e <- as.data.frame(est_h)
  ok <- which(!is.na(e$ind_mean))
  m <- sim$ind_mass[ok, , drop = FALSE]
  tau2 <- 2 * log(e$ind_mean[ok] / e$ind_median[ok])
  expect_equal(apply(log(m), 1, stats::var), tau2, tolerance = 0.05)
  # Monte Carlo error with 20000 draws: a few per cent of each statistic.
  expect_true(all(abs(log(rowMeans(m) / e$ind_mean[ok])) < 0.15))
  expect_true(all(abs(log(apply(m, 1, stats::median) / e$ind_median[ok])) < 0.05))
  lo <- apply(m, 1, stats::quantile, probs = 0.025, names = FALSE)
  hi <- apply(m, 1, stats::quantile, probs = 0.975, names = FALSE)
  expect_true(all(abs(log(lo / e$ind_lower[ok])) < 0.15))
  expect_true(all(abs(log(hi / e$ind_upper[ok])) < 0.15))
})

test_that("summarize_biomass() gives quantiles of the sums of the simulated biomasses", {
  for (pb in list(NULL, "habitat")) {
    s <- summarize_biomass(est_h, by = "site", n_draws = 300, seed = 4, population_by = pb,
                           na_rm = TRUE)
    d <- as.data.frame(simulate(est_h, nsim = 300, seed = 4, population_by = pb))
    d <- d[!is.na(d$biomass), ]
    tot <- aggregate(biomass ~ draw + site + currency + biomass_unit, data = d, FUN = sum)
    for (i in seq_len(nrow(s))) {
      v <- tot$biomass[tot$site == s$site[i] & tot$currency == s$currency[i] &
                         tot$biomass_unit == s$biomass_unit[i]]
      expect_length(v, 300L)
      expect_equal(s$biomass_median[i], stats::median(v), tolerance = 1e-10)
      expect_equal(s$biomass_lower[i], stats::quantile(v, 0.025, names = FALSE), tolerance = 1e-10)
      expect_equal(s$biomass_upper[i], stats::quantile(v, 0.975, names = FALSE), tolerance = 1e-10)
    }
  }
})

# One taxon without lengths in three sites, two of them in the same habitat.
ab3 <- data.frame(site = c("a", "b", "c"), taxon_id = "t15", abundance = 100, per = "m2",
                  habitat = c("x", "x", "y"))
est3 <- estimate_biomass(ab3, ex$taxa)

test_that("population_by shares a population within a group and not between groups", {
  vp <- est3$simulation$rows$v_population[1]
  expect_gt(vp, 0)
  sh <- simulate(est3, nsim = 20000, seed = 2)
  # Shared by all the rows: the same draws in the three sites.
  expect_identical(sh$ind_mass[1, ], sh$ind_mass[3, ])
  hb <- simulate(est3, nsim = 20000, seed = 2, population_by = "habitat")
  expect_identical(hb$ind_mass[1, ], hb$ind_mass[2, ])
  # Between habitats only the population differs: Var(ln a - ln c) = 2 v_population.
  expect_equal(stats::var(log(hb$ind_mass[1, ]) - log(hb$ind_mass[3, ])), 2 * vp,
               tolerance = 0.05)
  # Each row keeps its distribution.
  expect_equal(stats::var(log(hb$ind_mass[3, ])), stats::var(log(sh$ind_mass[3, ])),
               tolerance = 0.05)
  # After a conversion, the population of each conversion factor too.
  wm <- convert_currency(est3, to = "wet_mass")
  st <- wm$simulation$conversion$steps
  vc <- sum(st$v_population[st$row == 1L & st$type == "factor"])
  expect_gt(vc, 0)
  hw <- simulate(wm, nsim = 20000, seed = 3, population_by = "habitat")
  expect_identical(hw$ind_mass[1, ], hw$ind_mass[2, ])
  expect_equal(stats::var(log(hw$ind_mass[1, ]) - log(hw$ind_mass[3, ])), 2 * (vp + vc),
               tolerance = 0.05)
})

test_that("population_by narrows sums over groups and is reported", {
  s_sh <- summarize_biomass(est3, by = "per", n_draws = 4000, seed = 5)
  s_hb <- summarize_biomass(est3, by = "per", n_draws = 4000, seed = 5, population_by = "habitat")
  expect_lt(s_hb$biomass_upper / s_hb$biomass_lower, s_sh$biomass_upper / s_sh$biomass_lower)
  expect_identical(s_hb$biomass_mean, s_sh$biomass_mean)
  expect_identical(s_hb$population_by, "habitat")
  expect_true(is.na(s_sh$population_by))
  # One site alone does not change.
  one_sh <- summarize_biomass(est3, by = "site", n_draws = 4000, seed = 5)
  one_hb <- summarize_biomass(est3, by = "site", n_draws = 4000, seed = 5, population_by = "habitat")
  expect_equal(one_hb$biomass_upper / one_hb$biomass_lower,
               one_sh$biomass_upper / one_sh$biomass_lower, tolerance = 0.15)
})

test_that("simulate() is reproducible, restores the random state and handles zero rows", {
  set.seed(42)
  before <- .Random.seed
  a <- simulate(est_h, nsim = 20, seed = 9)
  expect_identical(.Random.seed, before)
  b <- simulate(est_h, nsim = 20, seed = 9)
  expect_identical(a$ind_mass, b$ind_mass)
  z <- est_h
  z$estimates <- est_h$estimates[0, ]
  z$simulation$rows <- est_h$simulation$rows[0, ]
  z$simulation$missing_keys <- list()
  sz <- simulate(z, nsim = 5, seed = 1)
  expect_equal(dim(sz$ind_mass), c(0L, 5L))
  dz <- as.data.frame(sz)
  expect_equal(nrow(dz), 0L)
  expect_true(all(c("draw", "site", "taxon_id", "ind_mass", "biomass") %in% names(dz)))
  expect_equal(nrow(summarize_biomass(z, population_by = "habitat")), 0L)
})

test_that("simulate() and summarize_biomass() check their arguments", {
  cls <- "bmassa_error_input"
  expect_error(simulate(est_h, nsim = 0), class = cls)
  expect_error(simulate(est_h, nsim = 2.5), class = cls)
  expect_error(simulate(est_h, nsim = 10, seed = "a"), class = cls)
  expect_error(simulate(est_h, nsim = 10, population_effect = "mixed"), class = cls)
  expect_error(simulate(est_h, nsim = 10, nonsense = 1), class = cls)
  expect_error(simulate(est_h, nsim = 10, population_by = "nowhere"), class = cls)
  expect_error(simulate(est_h, nsim = 10, population_by = "currency"), class = cls)
  expect_error(simulate(est_h, nsim = 10, population_by = "abundance"), class = cls)
  expect_error(simulate(est_h, nsim = 10, population_by = c("habitat", "habitat")), class = cls)
  expect_error(simulate(est_h, nsim = 10, population_by = "habitat",
                        population_effect = "independent"), class = cls)
  expect_error(summarize_biomass(est_h, population_by = 1), class = cls)
  expect_error(summarize_biomass(est_h, by = "population_by"), class = cls)
  ab_na <- ab_h
  ab_na$habitat[3] <- NA
  est_na <- suppressWarnings(estimate_biomass(ab_na, ex$taxa, ex$lengths,
                                              bias_correction = "lognormal"))
  expect_error(simulate(est_na, nsim = 10, population_by = "habitat"), "row 3", class = cls)
  ab_b <- ab3
  ab_b$biomass <- 1
  sb <- simulate(estimate_biomass(ab_b, ex$taxa), nsim = 2, seed = 1)
  expect_error(as.data.frame(sb), "biomass", class = cls)
  sim <- simulate(est_h, nsim = 5, seed = 1)
  expect_error(as.data.frame(sim, draws = 6), class = cls)
  expect_error(as.data.frame(sim, draws = 1.5), class = cls)
  expect_error(as.data.frame(sim, nonsense = 1), class = cls)
})

test_that("summarize_biomass() needs the result of estimate_biomass()", {
  expect_error(summarize_biomass(data.frame(site = "s1", biomass = 1)),
               class = "bmassa_error_input")
})
