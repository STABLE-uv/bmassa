# get_reference() and the structure of the reference tables ------------------

ref <- get_reference()

test_that("the reference has its tables and versions", {
  expect_s3_class(ref, "bmassa_reference")
  expect_named(ref, c("models", "lwr_types", "nodes", "components", "fixed", "conversion_models",
                      "conversion_nodes", "conversion_components", "conversion_fixed",
                      "carbon_equations", "sources", "info"))
  expect_setequal(ref$conversion_models$domain, c("invertebrate", "fish"))
  expect_true(all(c(ref$conversion_models$from, ref$conversion_models$to) %in%
                    c("wet_mass", "dry_mass", "afdm", "carbon", "energy")))
  # Energy only from masses (no model of energy per carbon).
  to_energy <- ref$conversion_models$to == "energy"
  expect_true(any(to_energy))
  expect_true(all(ref$conversion_models$from[to_energy] %in% c("wet_mass", "dry_mass", "afdm")))
  expect_equal(ref$carbon_equations$log10_a, c(-0.541, -0.665))
  expect_equal(ref$carbon_equations$b, c(0.811, 0.939))
  expect_equal(ref$info$col_release, "COL26.6 XR")
  expect_setequal(ref$models$domain, c("phytoplankton", "invertebrate", "fish", "amphibian"))
  expect_output(print(ref), "bmassa_reference")
})

test_that("every node has its parent, and conditional variances are valid", {
  n <- ref$nodes
  tax <- which(n$level %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
  parents <- vapply(tax, parent_row, integer(1), nodes = n)
  expect_false(anyNA(parents))
  expect_true(all(n$path_var[n$level != "stage"] > 0))
  lwr <- which(!is.na(n$path_var_slope))
  det <- n$path_var[lwr] * n$path_var_slope[lwr] - n$path_cov[lwr]^2
  expect_true(all(det > -1e-12))
  # The conditional variance of a node given its parent is not negative.
  for (k in tax[seq(1, length(tax), length.out = 50)]) {
    V <- node_cov(n, k)
    p <- parent_row(k, n)
    K <- parent_cov_matrix(n, k, ncol(V))
    C <- V - K %*% pinv_sym(node_cov(n, p)) %*% t(K)
    expect_true(all(eigen(C, symmetric = TRUE)$values > -1e-8))
  }
})

test_that("the conversion nodes have their parents, and col_id follows the COL", {
  n <- ref$conversion_nodes
  tax <- which(n$level %in% c("kingdom", "phylum", "class", "order", "family", "genus", "species"))
  some <- tax[seq(1, length(tax), length.out = 200)]
  expect_false(anyNA(vapply(some, parent_row, integer(1), nodes = n)))
  expect_true(all(n$path_var[tax] > 0))
  expect_true(all(n$effect_var[n$level == "habitat"] >= 0))
  # Identifiers of the COL26.6 XR (answers of the GBIF API v2 on 2026-11-02).
  m <- ref$nodes
  expect_equal(m$col_id[m$node == "Daphnia magna"][1], "6CCSV")
  expect_equal(m$col_id[m$node == "Chironomus" & m$level == "genus"][1], "BXWQR")
  expect_equal(n$col_id[n$node == "Chironomus plumosus"][1], "8BGDB")
  expect_gt(mean(!is.na(m$col_id[m$level %in% c("genus", "species")])), 0.95)
})

test_that("rmvnorm0() draws with the given covariance", {
  S <- matrix(c(1, 0.6, 0.6, 2), 2)
  d <- withr::with_seed(1, rmvnorm0(20000, S))
  expect_equal(stats::cov(d), S, tolerance = 0.05)
  expect_equal(pinv_sym(S) %*% S, diag(2), tolerance = 1e-10)
})

test_that("check_reference() rejects other objects", {
  expect_error(check_reference(list()), class = "bmassa_error_input")
})
