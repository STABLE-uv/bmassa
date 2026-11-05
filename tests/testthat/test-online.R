# Tests that use the internet (GBIF API): not on CRAN, and only when online -------

test_that("standardize_taxa() queries the COL XR", {
  skip_on_cran()
  skip_if_not_installed("curl")
  skip_if_not_installed("jsonlite")
  skip_if_offline("api.gbif.org")
  rm(list = ls(colxr_memory, all.names = TRUE), envir = colxr_memory)
  taxa <- data.frame(taxon_id = c("a", "b"), broad_group = "macroinvertebrate",
                     scientific_name = c("Baetis rhodani", "Ephemerella ignita"))
  s <- standardize_taxa(taxa, quiet = TRUE)
  expect_equal(s$match_type, c("exact", "exact"))
  expect_equal(s$family[1], "Baetidae")
  expect_equal(s$col_status[2], "synonym")
  expect_match(s$col_release[1], "XR$")
})
