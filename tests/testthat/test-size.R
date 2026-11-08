# summarize_size() -----------------------------------------------------------------

size_fixture <- function() {
  tx <- example_survey$taxa
  tx <- tx[tx$taxon_id %in% c("t10", "t11", "t13", "t14", "t15", "t19", "t20"), ]
  ab <- data.frame(
    site = c("S1", "S1", "S1", "S1", "S2", "S3", "S4", "S4", "S4", "S5", "S5", "S6", "S6"),
    taxon_id = c("t13", "t10", "t11", "t14", "t13", "t11", "t14", "t15", "t13", "t19", "t20",
                 "t13", "t14"),
    abundance = c(10, 30, 20, 40, 5, 8, 40, 25, 10, 4, 30, NA, 10),
    per = "m2", stringsAsFactors = FALSE)
  ln <- data.frame(
    taxon_id = c("t13", "t13", "t13", "t13", "t10", "t11", "t19"),
    site = c("S1", "S1", "S2", "S2", NA, NA, NA),
    length_mean = c(NA, NA, NA, NA, 12, 6, 58), length = c(4, 8, 4, 8, NA, NA, NA),
    n = c(1, 3, 1, 3, NA, NA, NA), length_unit = "mm",
    length_type = c(rep("body_length", 6), "snout_vent_length"),
    length_cv = c(NA, NA, NA, NA, NA, 0.2, 0.15), stringsAsFactors = FALSE)
  est <- suppressWarnings(estimate_biomass(ab, tx, ln, bias_correction = "lognormal"))
  list(est = est, ab = ab, tx = tx, ln = ln)
}

row_of <- function(s, site) s$summary[s$summary$site == site, ]

test_that("lengths: shares, mean and median by hand", {
  f <- size_fixture()
  s <- summarize_size(f$est, size = "length")
  expect_s3_class(s, "bmassa_size")
  expect_named(s$summary, c("site", "per", "size", "size_type", "size_unit", "abundance",
                            "share_excluded", "share_measured", "share_lognormal",
                            "share_taxon_mean", "mean_size", "mean_size_lower",
                            "mean_size_upper", "median_size", "median_size_lower",
                            "median_size_upper", "approximation", "level", "interval",
                            "n_rows", "n_sized", "n_missing", "status", "n_draws"))
  r1 <- row_of(s, "S1")
  # 10 Gammarus measured (1 of 4 mm, 3 of 8 mm), 30 Chironomus with a mean of
  # 12 mm, 20 Cloeon with a mean of 6 mm and CV 0.2, 40 Physella without length.
  expect_equal(r1$abundance, 100)
  expect_equal(r1$share_excluded, 0.4)
  expect_equal(c(r1$share_measured, r1$share_lognormal, r1$share_taxon_mean),
               c(10, 20, 30) / 60)
  expect_equal(r1$mean_size, (10 * 7 + 30 * 12 + 20 * 6) / 60)
  # The cumulative share of the individuals reaches one half at 12 mm (the
  # Chironomus): below it, 10 Gammarus and at most 20 Cloeon.
  expect_equal(r1$median_size, 12)
  expect_equal(c(r1$median_size_lower, r1$median_size_upper), c(12, 12))
  expect_true(r1$approximation)
  expect_true(is.na(r1$n_draws))
  expect_identical(r1$size_type, "body_length")
  expect_identical(r1$size_unit, "mm")
  # Measured lengths only: the median is a measured length.
  r2 <- row_of(s, "S2")
  expect_equal(c(r2$mean_size, r2$median_size), c(7, 8))
  expect_false(r2$approximation)
  # A mean length with its spread: lognormal, median below the mean.
  r3 <- row_of(s, "S3")
  expect_equal(r3$mean_size, 6)
  expect_equal(r3$median_size, 6 / sqrt(1 + 0.2^2), tolerance = 1e-9)
  expect_equal(r3$share_lognormal, 1)
  # Groups without lengths have no size.
  r4 <- row_of(s, "S4")
  expect_equal(r4$share_excluded, 1)
  expect_true(is.na(r4$mean_size))
})

test_that("length classes: shares of individuals and of biomass by hand", {
  f <- size_fixture()
  s <- summarize_size(f$est, size = "length", breaks = c(5, 10), n_draws = 200, seed = 1)
  cl <- s$classes
  expect_identical(unique(cl$class), c("[0, 5)", "[5, 10)", "[10, Inf)"))
  expect_equal(unique(cl$lower), c(0, 5, 10))
  shares <- function(site, col) cl[[col]][cl$site == site]
  e <- f$est$estimates
  sim <- f$est$simulation$rows
  slope <- function(i) f$est$reference$nodes$path_slope[sim$node[i]]
  # Measured lengths (S2): a quarter of the individuals in [0, 5); the shares
  # of biomass are those of the masses at each length, the same in every draw.
  expect_equal(shares("S2", "share_individuals"), c(0.25, 0.75, 0))
  b <- slope(which(e$site == "S2"))
  bio <- c(4^b, 3 * 8^b, 0) / (4^b + 3 * 8^b)
  expect_equal(shares("S2", "share_biomass"), bio)
  expect_equal(shares("S2", "share_biomass_lower"), bio)
  expect_equal(shares("S2", "share_biomass_upper"), bio)
  # Mean length with its spread (S3): lognormal lengths, and the biomass of
  # each class by numerical integration of the masses.
  sl <- sqrt(log(1 + 0.2^2))
  mu <- log(6) - sl^2 / 2
  ind <- diff(stats::plnorm(c(0, 5, 10, Inf), mu, sl))
  expect_equal(shares("S3", "share_individuals"), ind, tolerance = 1e-9)
  b <- slope(which(e$site == "S3"))
  part <- function(lo, hi) {
    stats::integrate(function(x) stats::dlnorm(x, mu, sl) * x^b, lo, hi)$value
  }
  tot <- exp(b * mu + b^2 * sl^2 / 2)
  expect_equal(shares("S3", "share_biomass"), c(part(0, 5), part(5, 10), part(10, Inf)) / tot,
               tolerance = 1e-6)
  # Several rows (S1): individuals by hand.
  z <- stats::pnorm((log(c(5, 10)) - mu) / sl)
  expect_equal(shares("S1", "share_individuals"),
               c(2.5 + 20 * z[1], 7.5 + 20 * (z[2] - z[1]), 30 + 20 * (1 - z[2])) / 60,
               tolerance = 1e-9)
  # The shares add up to 1 in every group with sizes.
  ok <- !is.na(cl$share_individuals)
  expect_equal(as.vector(tapply(cl$share_individuals[ok], cl$site[ok], sum)), rep(1, 4))
  expect_equal(as.vector(tapply(cl$share_biomass[ok], cl$site[ok], sum)), rep(1, 4),
               tolerance = 1e-9)
  expect_identical(s$summary$n_draws[1], 200L)
})

test_that("masses: indicators from the draws of simulate()", {
  f <- size_fixture()
  s <- summarize_size(f$est, n_draws = 200, seed = 3, breaks = c(0.1, 1))
  sim <- simulate(f$est, nsim = 200, seed = 3)
  e <- f$est$estimates
  q <- function(x, centre = stats::median) {
    c(centre(x), stats::quantile(x, c(0.025, 0.975), names = FALSE))
  }
  # The mean mass of the individuals: the biomass of the rows over their
  # individuals, whatever the spread of sizes within each row.
  r <- which(e$site == "S1")
  m_d <- colSums(e$abundance[r] * sim$ind_mass[r, ]) / sum(e$abundance[r])
  r1 <- row_of(s, "S1")
  expect_equal(c(r1$mean_size, r1$mean_size_lower, r1$mean_size_upper), q(m_d))
  # Masses without lengths (S4): every individual at the mean mass of its row,
  # so the median and the classes follow from the draws directly.
  r <- which(e$site == "S4")
  a <- e$abundance[r]
  med <- apply(sim$ind_mass[r, ], 2, function(m) {
    o <- order(m)
    m[o][which(cumsum(a[o]) / sum(a) >= 0.5)[1]]
  })
  r4 <- row_of(s, "S4")
  expect_equal(c(r4$median_size, r4$median_size_lower, r4$median_size_upper), q(med))
  expect_equal(r4$share_taxon_mean, 1)
  expect_true(r4$approximation)
  in_class <- function(m, lo, hi) m >= lo & m < hi
  cl <- s$classes[s$classes$site == "S4", ]
  for (j in seq_len(nrow(cl))) {
    ind <- apply(sim$ind_mass[r, ], 2, function(m) sum(a[in_class(m, cl$lower[j], cl$upper[j])]) / sum(a))
    bio <- apply(sim$ind_mass[r, ], 2, function(m) {
      sum((a * m)[in_class(m, cl$lower[j], cl$upper[j])]) / sum(a * m)
    })
    # Shares: the mean of the draws.
    expect_equal(c(cl$share_individuals[j], cl$share_individuals_lower[j],
                   cl$share_individuals_upper[j]), q(ind, mean))
    expect_equal(c(cl$share_biomass[j], cl$share_biomass_lower[j], cl$share_biomass_upper[j]),
                 q(bio, mean))
  }
  expect_identical(r4$size_type, "wet_mass")
  expect_identical(r4$size_unit, "mg")
  # One row: its mean mass, with the spread of its prediction interval.
  big <- summarize_size(f$est, n_draws = 20000, seed = 1)
  r2 <- row_of(big, "S2")
  i <- which(e$site == "S2")
  expect_equal(c(r2$mean_size, r2$mean_size_lower, r2$mean_size_upper),
               c(e$ind_median[i], e$ind_lower[i], e$ind_upper[i]), tolerance = 0.05)
})

test_that("masses of measured individuals and of lognormal lengths", {
  f <- size_fixture()
  e <- f$est$estimates
  sim <- f$est$simulation$rows
  slope <- function(i) f$est$reference$nodes$path_slope[sim$node[i]]
  s <- summarize_size(f$est, n_draws = 20000, seed = 2)
  # Measured individuals (S2): masses M (L / L*)^b, with L* the length of their
  # mean mass; the median is the mass of the 8 mm individuals.
  i <- which(e$site == "S2")
  b <- slope(i)
  lstar <- (0.25 * 4^b + 0.75 * 8^b)^(1 / b)
  expect_equal(sim$length_used[i], lstar)
  expect_equal(row_of(s, "S2")$median_size, e$ind_median[i] * (8 / lstar)^b, tolerance = 0.05)
  # Lognormal lengths (S3): lognormal masses with log-SD b times that of the
  # lengths, whose mean is the mass of the row.
  i <- which(e$site == "S3")
  sm <- slope(i) * sqrt(log(1 + 0.2^2))
  expect_equal(row_of(s, "S3")$median_size, e$ind_median[i] * exp(-sm^2 / 2), tolerance = 0.05)
  expect_equal(row_of(s, "S3")$share_lognormal, 1)
})

test_that("rows without a size are counted, not used", {
  f <- size_fixture()
  # S5: 4 adult frogs with a mean length and 30 larvae, not estimated.
  m <- row_of(summarize_size(f$est, n_draws = 100, seed = 1), "S5")
  expect_equal(m$share_excluded, 30 / 34)
  expect_equal(m$n_rows, 2L)
  expect_equal(m$n_sized, 1L)
  l <- row_of(summarize_size(f$est, size = "length"), "S5")
  expect_equal(l$share_excluded, 30 / 34)
  expect_equal(l$mean_size, 58)
  expect_identical(l$size_type, "snout_vent_length")
})

test_that("missing abundances make a group incomplete; na_rm leaves them out", {
  f <- size_fixture()
  s <- summarize_size(f$est, n_draws = 100, seed = 1)
  r6 <- row_of(s, "S6")
  expect_identical(r6$status, "incomplete")
  expect_equal(r6$n_missing, 1L)
  expect_true(is.na(r6$mean_size))
  expect_true(is.na(r6$abundance))
  expect_true(is.na(r6$share_excluded))
  s2 <- summarize_size(f$est, n_draws = 100, seed = 1, na_rm = TRUE)
  r6 <- row_of(s2, "S6")
  expect_identical(r6$status, "incomplete")
  expect_equal(r6$abundance, 10)
  expect_equal(r6$share_excluded, 0)
  expect_false(is.na(r6$mean_size))
  expect_identical(row_of(s2, "S1")$status, "complete")
})

test_that("groups are sorted, and units of abundance are kept apart", {
  f <- size_fixture()
  ab <- f$ab[f$ab$site %in% c("S1", "S2"), ]
  ab$site[ab$site == "S1"] <- "Z"
  ab$per[ab$site == "Z" & ab$taxon_id == "t14"] <- "L"
  ln <- f$ln
  ln$site[ln$site %in% "S1"] <- "Z"
  est <- suppressWarnings(estimate_biomass(ab, f$tx, ln, bias_correction = "lognormal"))
  s <- summarize_size(est, n_draws = 100, seed = 1)
  expect_identical(s$summary$site, c("S2", "Z", "Z"))
  expect_identical(s$summary$per, c("m2", "L", "m2"))
  expect_equal(s$summary$abundance, c(5, 40, 60))
  # By another column of the estimates.
  b <- summarize_size(f$est, by = c("site", "taxon_id"), size = "length")
  expect_identical(b$summary$site[1:4], rep("S1", 4))
  expect_identical(b$summary$taxon_id[1:4], c("t10", "t11", "t13", "t14"))
})

test_that("zero rows give zero rows", {
  f <- size_fixture()
  est0 <- estimate_biomass(f$ab[0, ], f$tx)
  s <- summarize_size(est0, breaks = 1)
  expect_identical(nrow(s$summary), 0L)
  expect_identical(nrow(s$classes), 0L)
  expect_identical(names(s$summary), names(summarize_size(size_fixture()$est, n_draws = 100,
                                                          seed = 1)$summary))
  expect_output(print(s), "0 groups")
})

test_that("input errors", {
  f <- size_fixture()
  est <- f$est
  expect_error(summarize_size(data.frame()), class = "bmassa_error_input")
  expect_error(summarize_size(est, by = 1), class = "bmassa_error_input")
  expect_error(summarize_size(est, by = "nowhere"), class = "bmassa_error_input")
  expect_error(summarize_size(est, by = "share_excluded"), class = "bmassa_error_input")
  expect_error(summarize_size(est, size = "weight"), class = "bmassa_error_input")
  expect_error(summarize_size(est, breaks = c(2, 1)), class = "bmassa_error_input")
  expect_error(summarize_size(est, breaks = -1), class = "bmassa_error_input")
  expect_error(summarize_size(est, breaks = c(1, NA)), class = "bmassa_error_input")
  expect_error(summarize_size(est, breaks = "1"), class = "bmassa_error_input")
  expect_error(summarize_size(est, n_draws = 10), class = "bmassa_error_input")
  expect_error(summarize_size(est, na_rm = NA), class = "bmassa_error_input")
  expect_error(summarize_size(est, extra = 1), class = "bmassa_error_input")
  # Lengths of two types in one group.
  ab <- f$ab
  ab$site[ab$site == "S5"] <- "S1"
  est2 <- suppressWarnings(estimate_biomass(ab, f$tx, f$ln, bias_correction = "lognormal"))
  expect_error(summarize_size(est2, size = "length"), "more than one type",
               class = "bmassa_error_input")
  expect_no_error(summarize_size(est2, by = c("site", "length_type"), size = "length"))
})

test_that("the seed restores the random state; lengths without classes draw nothing", {
  f <- size_fixture()
  set.seed(10)
  r0 <- .Random.seed
  a <- summarize_size(f$est, n_draws = 100, seed = 5)
  expect_identical(.Random.seed, r0)
  b <- summarize_size(f$est, n_draws = 100, seed = 5)
  expect_identical(a, b)
  s <- summarize_size(f$est, size = "length")
  expect_identical(.Random.seed, r0)
  expect_null(s$classes)
})

test_that("other currencies and units", {
  f <- size_fixture()
  dm <- suppressWarnings(convert_currency(f$est, to = "dry_mass", unit = "ug"))
  s <- summarize_size(dm, n_draws = 100, seed = 1)
  expect_identical(unique(s$summary$size_type[!is.na(s$summary$mean_size)]), "dry_mass")
  expect_identical(unique(s$summary$size_unit), "ug")
})

test_that("print", {
  f <- size_fixture()
  expect_snapshot(print(summarize_size(f$est, breaks = c(1, 10), n_draws = 200, seed = 1)))
  expect_snapshot(print(summarize_size(f$est, size = "length")))
})
