# cite_sources() -----------------------------------------------------------------

test_that("all the sources, or those of an estimate", {
  all <- cite_sources()
  expect_named(all, c("source_id", "citation", "doi", "licence", "n_rows"))
  expect_false(is.unsorted(all$source_id))
  expect_true(all(is.na(all$n_rows)))
  ex <- example_survey
  est <- suppressWarnings(estimate_biomass(ex$abundance, ex$taxa, ex$lengths, bias_correction = "none"))
  s <- cite_sources(est)
  e <- as.data.frame(est)
  # The sources of the models and of the conversion to wet mass.
  ok <- !is.na(e$ind_median)
  used <- unique(unlist(strsplit(c(e$sources[ok], e$conversion_sources[ok & !is.na(e$conversion_sources)]),
                                 ";")))
  expect_setequal(s$source_id, used)
  expect_true(all(s$n_rows >= 1L))
  expect_true(all(s$source_id %in% all$source_id))
})

test_that("the sources of conversions are cited too", {
  ex <- example_survey
  est <- suppressWarnings(estimate_biomass(ex$abundance, ex$taxa, ex$lengths, bias_correction = "none"))
  cv <- suppressWarnings(convert_currency(est, to = "carbon"))
  s <- cite_sources(cv)
  expect_true(all(c("brey2010", "mdl2000") %in% s$source_id))
  e <- as.data.frame(cv)
  expect_equal(s$n_rows[s$source_id == "mdl2000"],
               sum(grepl("mdl2000", e$conversion_sources) & !is.na(e$ind_median)))
})

test_that("x must be an estimate", {
  expect_error(cite_sources(data.frame()), class = "bmassa_error_input")
})
