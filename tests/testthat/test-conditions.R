# Tests for the condition helpers.

test_that("abort_input() signals a classed error that names the caller", {
  f <- function(x) abort_input("`x` must be positive.", argument = "x")
  err <- expect_error(f(-1), class = "bmassa_error_input")
  expect_s3_class(err, "bmassa_error")
  expect_identical(conditionCall(err), quote(f(-1)))
  expect_identical(conditionMessage(err), "`x` must be positive.")
  expect_identical(err$argument, "x")
})

test_that("abort_internal() is a different class and asks for a report", {
  err <- expect_error(abort_internal("rows were lost"), class = "bmassa_error_internal")
  expect_false(inherits(err, "bmassa_error_input"))
  expect_match(conditionMessage(err), "bug in bmassa", fixed = TRUE)
})

test_that("warn_not_estimable() signals a classed warning", {
  expect_warning(
    warn_not_estimable("2 samples could not be estimated; see `status`."),
    class = "bmassa_warning_not_estimable"
  )
})

test_that("check_installed() explains how to install a missing package", {
  expect_true(check_installed("stats", "fit models"))
  err <- expect_error(
    check_installed("bmassanotapackage", "draw the figure"),
    class = "bmassa_error_missing_package"
  )
  expect_match(conditionMessage(err), "required to draw the figure", fixed = TRUE)
  expect_match(conditionMessage(err), 'install.packages("bmassanotapackage")', fixed = TRUE)
  expect_identical(err$package, "bmassanotapackage")
})

test_that("format_items() and format_rows() abbreviate long lists", {
  expect_identical(format_items(character()), "")
  expect_identical(format_items(letters[1:3]), "a, b, c")
  expect_identical(format_items(letters[1:7], max = 5), "a, b, c, d, e and 2 more")
  expect_identical(format_rows(3L), "row 3")
  expect_identical(format_rows(c(3L, 8L)), "rows 3, 8")
})
