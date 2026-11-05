# convert_abundance(): abundances on a common basis ----------------------------------

ab <- data.frame(site = c("a", "a", "a", "b", "b"), taxon_id = c("t01", "t06", "t10", "t06", "t10"),
                 abundance = c(2000, 50, 120, NA, 0), per = c("mL", "L", "0.25 m2", "L", "m2"),
                 habitat = "x", stringsAsFactors = FALSE)
tanks <- data.frame(site = c("a", "b"), volume = c(1.23, 1.23), volume_unit = "m3",
                    area = c(1.54, 1.54), area_unit = "m2")

test_that("the sizes in `per` and `to` are read as volumes, areas or the site", {
  p <- parse_basis(c("L", "5 L", "0.25 m2", "100 m2", "mL", "m3", "ha", "site", "1e3 L", "ml",
                     ".5L", "100m2"))
  expect_identical(p$dim, c("volume", "volume", "area", "area", "volume", "volume", "area", "site",
                            "volume", "volume", "volume", "area"))
  expect_equal(p$size, c(1, 5, 0.25, 100, 0.001, 1000, 10000, 1, 1000, 0.001, 0.5, 100))
  expect_identical(p$label[c(1, 2, 4, 10, 12)], c("L", "5 L", "100 m2", "mL", "100 m2"))
  bad <- parse_basis(c("2 site", "0,5 L", "sample", "0 L", "-1 L", "L2", ""))
  expect_true(all(is.na(bad$dim)))
})

test_that("within one dimension the factor is the ratio of the sizes", {
  x <- ab[ab$per %in% c("mL", "L"), ]
  out <- convert_abundance(x, to = "L")
  expect_equal(out$per_factor, c(1000, 1, 1))
  expect_equal(out$abundance, c(2e6, 50, NA))
  expect_identical(out$per, rep("L", 3))
  expect_true(all(is.na(out$per_assumption)))
  out2 <- convert_abundance(data.frame(site = "a", abundance = 7, per = c("m2", "100 m2", "m3")[1:2]),
                            to = "100 m2")
  expect_equal(out2$per_factor, c(100, 1))
  expect_identical(out2$per, rep("100 m2", 2))
  out3 <- convert_abundance(data.frame(site = "a", abundance = 7, per = "m3"), to = "L")
  expect_equal(out3$abundance, 7 / 1000)
})

test_that("between volume and area the depth of the site is used", {
  d <- data.frame(site = c("a", "b"), depth = c(80, 0.5), depth_unit = c("cm", "m"))
  out <- convert_abundance(ab, to = "m2", sites = d)
  # per mL: x 1000 (mL/L) x 1000 (L/m3) x 0.8 m; per L: x 1000 x 0.8; 0.25 m2: x 4.
  expect_equal(out$per_factor, c(800000, 800, 4, 500, 1))
  expect_equal(out$abundance, c(1.6e9, 40000, 480, NA, 0))
  expect_match(out$per_assumption[1], "depth 0.8 m")
  expect_true(is.na(out$per_assumption[3]))
  # Without a depth: the volume over the area.
  out2 <- convert_abundance(ab[2, ], to = "m2", sites = tanks)
  expect_equal(out2$per_factor, 1000 * (1230 / 1000 / 1.54))
  # And back: per m2 to per L.
  out3 <- convert_abundance(ab[5, ], to = "L", sites = d)
  expect_equal(out3$per_factor, 1 / 500)
})

test_that("to and from the whole site the volume or the area of the site is used", {
  out <- convert_abundance(ab, to = "site", sites = tanks)
  expect_equal(out$per_factor, c(1.23e6, 1230, 1.54 / 0.25, 1230, 1.54))
  expect_identical(out$per, rep("site", 5))
  expect_identical(out$abundance_original, ab$abundance)
  expect_identical(out$per_original, ab$per)
  expect_identical(out[, c("site", "taxon_id", "habitat")], ab[, c("site", "taxon_id", "habitat")])
  back <- convert_abundance(data.frame(site = c("a", "a"), abundance = 1230, per = "site"),
                            to = "L", sites = tanks)
  expect_equal(back$abundance, c(1, 1))
  back2 <- convert_abundance(data.frame(site = "a", abundance = 1.54, per = "site"), to = "m2",
                             sites = tanks)
  expect_equal(back2$abundance, 1)
  same <- convert_abundance(data.frame(site = "a", abundance = 3, per = "site"), to = "site")
  expect_equal(same$per_factor, 1)
})

test_that("the biomasses of the converted table share one unit", {
  ex <- example_survey
  d <- data.frame(site = c("L1", "L2"), depth = c(1.5, 0.8), depth_unit = "m")
  a2 <- convert_abundance(ex$abundance, to = "m2", sites = d)
  est <- suppressWarnings(estimate_biomass(a2, ex$taxa, ex$lengths, bias_correction = "lognormal"))
  wm <- suppressWarnings(convert_currency(est, to = "wet_mass", unit = "mg"))
  u <- unique(stats::na.omit(as.data.frame(wm)$biomass_unit))
  expect_identical(u, "mg/m2")
  s <- summarize_biomass(wm, by = "site", na_rm = TRUE, seed = 1)
  expect_equal(nrow(s), 2L)
  # A biomass per m2 is the biomass per L times 1000 times the depth.
  e0 <- as.data.frame(suppressWarnings(estimate_biomass(ex$abundance, ex$taxa, ex$lengths,
                                                         bias_correction = "lognormal")))
  k <- which(e0$site == "L1" & ex$abundance$per == "L" & !is.na(e0$biomass_mean))[1]
  expect_equal(as.data.frame(est)$biomass_mean[k], e0$biomass_mean[k] * 1500)
})

test_that("zero rows give zero rows with the new columns", {
  out <- convert_abundance(ab[0, ], to = "site", sites = tanks)
  expect_equal(nrow(out), 0L)
  expect_true(all(c(names(ab), "abundance_original", "per_original", "per_factor",
                    "per_assumption") %in% names(out)))
})

test_that("convert_abundance() checks its input", {
  cls <- "bmassa_error_input"
  expect_error(convert_abundance(ab, to = "m2"), "needs `sites`", class = cls)
  expect_error(convert_abundance(ab, to = "sample"), class = cls)
  expect_error(convert_abundance(ab, to = c("L", "m2")), class = cls)
  bad <- ab
  bad$per[2] <- "sample"
  expect_error(convert_abundance(bad, to = "L"), "row 2", class = cls)
  bad$per[2] <- "0,5 L"
  expect_error(convert_abundance(bad, to = "L"), class = cls)
  out <- convert_abundance(ab, to = "site", sites = tanks)
  expect_error(convert_abundance(out, to = "L", sites = tanks), "converted", class = cls)
  expect_error(convert_abundance(ab, to = "site", sites = tanks[1, ]), "'b'", class = cls)
  expect_error(convert_abundance(ab, to = "site", sites = tanks[, c("site", "area", "area_unit")]),
               "volume", class = cls)
  expect_error(convert_abundance(ab, to = "m2", sites = data.frame(site = c("a", "b"))),
               "depth", class = cls)
  expect_error(convert_abundance(ab, to = "site", sites = rbind(tanks, tanks)), class = cls)
  neg <- tanks
  neg$volume[1] <- -1
  expect_error(convert_abundance(ab, to = "site", sites = neg), class = cls)
  wrong <- tanks
  wrong$volume_unit <- "gallon"
  expect_error(convert_abundance(ab, to = "site", sites = wrong), class = cls)
  nounit <- tanks[, c("site", "volume", "area", "area_unit")]
  expect_error(convert_abundance(ab, to = "site", sites = nounit), "volume_unit", class = cls)
  expect_error(convert_abundance(ab[, c("site", "abundance")], to = "L"), class = cls)
  neg_ab <- ab
  neg_ab$abundance[1] <- -3
  expect_error(convert_abundance(neg_ab, to = "L", sites = tanks), class = cls)
})
