# Answers of the GBIF name-matching API (COL XR) recorded on 2026-11-02 for the
# names of fixtures/colxr_answers.rds, reduced to the fields the package reads.
# local_colxr_answers() makes standardize_taxa() use them instead of the
# internet, and forgets the answers of the session before and after.
local_colxr_answers <- function(env = parent.frame()) {
  fx <- readRDS(test_path("fixtures", "colxr_answers.rds"))
  forget <- function() rm(list = ls(colxr_memory, all.names = TRUE), envir = colxr_memory)
  forget()
  withr::defer(forget(), envir = env)
  testthat::local_mocked_bindings(
    colxr_fetch = function(urls, ...) {
      miss <- setdiff(urls, names(fx$answers))
      if (length(miss)) stop("query not recorded: ", miss[1])
      unname(fx$answers[urls])
    },
    .env = env
  )
  fx
}
