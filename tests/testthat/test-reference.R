# get_reference() and the structure of the reference tables ------------------

ref <- get_reference()

test_that("the reference has its tables and versions", {
  expect_s3_class(ref, "bmassa_reference")
  expect_named(ref, c("models", "lwr_types", "nodes", "components", "fixed", "conversion_models",
                      "conversion_nodes", "conversion_components", "conversion_fixed",
                      "carbon_equations", "colonial_genera", "sources", "info"))
  expect_true(all(c("Microcystis", "Dinobryon") %in% ref$colonial_genera))
  expect_setequal(ref$conversion_models$domain, c("invertebrate", "fish", "amphibian_larva"))
  expect_true(all(c(ref$conversion_models$from, ref$conversion_models$to) %in%
                    c("wet_mass", "dry_mass", "afdm", "carbon", "nitrogen", "phosphorus", "energy")))
  # Nitrogen and phosphorus only from dry mass, for invertebrates, fish and
  # tadpoles; tadpoles only to carbon, nitrogen and phosphorus.
  elem <- ref$conversion_models$to %in% c("nitrogen", "phosphorus")
  expect_equal(sum(elem), 6L)
  expect_true(all(ref$conversion_models$from[elem] == "dry_mass"))
  tad <- ref$conversion_models$domain == "amphibian_larva"
  expect_setequal(ref$conversion_models$to[tad], c("carbon", "nitrogen", "phosphorus"))
  expect_true(all(ref$conversion_models$from[tad] == "dry_mass"))
  # The model of typical lengths of tadpoles goes to their length-mass model.
  tl <- ref$models[ref$models$kind == "typical_length", ]
  expect_equal(tl$model, "amphibian_larva_length")
  expect_equal(tl$lwr_model, "amphibian_larva_lwr")
  expect_equal(tl$currency, "length")
  expect_equal(tl$unit, "mm")
  expect_equal(tl$length_type, "total_length")
  expect_true(tl$length_cv > 0.05 && tl$length_cv < 0.5)
  tlwr <- ref$models[ref$models$model == "amphibian_larva_lwr", ]
  expect_equal(c(tlwr$currency, tlwr$unit, tlwr$life_stage), c("dry_mass", "mg", "larva"))
  # Energy only from masses (no model of energy per carbon).
  to_energy <- ref$conversion_models$to == "energy"
  expect_true(any(to_energy))
  expect_true(all(ref$conversion_models$from[to_energy] %in% c("wet_mass", "dry_mass", "afdm")))
  expect_equal(ref$carbon_equations$log10_a, c(-0.541, -0.665))
  expect_equal(ref$carbon_equations$b, c(0.811, 0.939))
  expect_equal(ref$info$col_release, "COL26.10 XR")
  expect_setequal(ref$models$domain, c("phytoplankton", "invertebrate", "fish", "amphibian",
                                       "amphibian_larva"))
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
  # Identifiers of the COL26.10 XR (answers of the GBIF API v2 on 2026-11-06; the
  # same in the COL26.6 XR).
  m <- ref$nodes
  expect_equal(m$col_id[m$node == "Daphnia magna"][1], "6CCSV")
  expect_equal(m$col_id[m$node == "Chironomus" & m$level == "genus"][1], "BXWQR")
  expect_equal(n$col_id[n$node == "Chironomus plumosus"][1], "8BGDB")
  expect_gt(mean(!is.na(m$col_id[m$level %in% c("genus", "species")])), 0.95)
})

test_that("the conversion models hold the soft tissue of freshwater molluscs", {
  n <- ref$conversion_nodes
  has <- function(m, sp) any(n$model == m & n$level == "species" & n$node == sp)
  # From the compilation of molluscs of 2026-11-09 (data-raw/conversions.R,
  # section 10): species without factors before.
  expect_true(has("invertebrate wet_mass_to_dry_mass", "Unio mancus"))
  expect_true(has("invertebrate dry_mass_to_afdm", "Planorbarius corneus"))
  expect_true(has("invertebrate dry_mass_to_afdm", "Elimia fascinans"))
  expect_true(has("invertebrate dry_mass_to_energy", "Parreysia corrugata"))
  expect_true(has("invertebrate afdm_to_energy", "Lymnaea stagnalis"))
  # Their sources are cited; the edible flesh of snails (Debnath et al. 2016)
  # is not used.
  expect_true(all(c("ravera2007", "benke1999", "mazuran2025", "doherty1986") %in% ref$sources$source_id))
  expect_false("debnath2016" %in% ref$sources$source_id)
})

test_that("the conversion models hold the ash and water of freshwater zooplankton", {
  n <- ref$conversion_nodes
  has <- function(m, sp) any(n$model == m & n$level == "species" & n$node == sp)
  # From the compilation of zooplankton of 2026-11-10 (data-raw/conversions.R,
  # section 12): species without factors before.
  expect_true(has("invertebrate dry_mass_to_afdm", "Bythotrephes longimanus"))
  expect_true(has("invertebrate dry_mass_to_afdm", "Eudiaptomus padanus"))
  expect_true(has("invertebrate wet_mass_to_dry_mass", "Moina macrocopa"))
  expect_true(all(c("riccardi1999", "birge1922", "watanabe1978", "hasan2023") %in% ref$sources$source_id))
  # The moisture of zooplankton dried on filter paper (Ovie and Ovie 2006) is
  # not used.
  expect_false("ovie2006" %in% ref$sources$source_id)
})

test_that("the models of carbon, nitrogen and phosphorus hold the animals of Limno-STOICH", {
  n <- ref$conversion_nodes
  has <- function(m, sp) any(n$model == m & n$level == "species" & n$node == sp)
  # From the Limno-STOICH database (Corman et al. 2026; data-raw/conversions.R,
  # section 13): species without factors before.
  expect_true(has("invertebrate dry_mass_to_carbon", "Aquarius remigis"))
  expect_true(has("invertebrate dry_mass_to_phosphorus", "Fusconaia cerina"))
  expect_true(has("fish dry_mass_to_nitrogen", "Anablepsoides hartii"))
  expect_true(has("fish dry_mass_to_phosphorus", "Anablepsoides hartii"))
  # The wood frog (Lithobates sylvaticus of the source) in the COL XR.
  expect_true(has("amphibian_larva dry_mass_to_carbon", "Boreorana sylvatica"))
  expect_true(has("amphibian_larva dry_mass_to_nitrogen", "Pseudacris maculata"))
  # The juvenile tailed frogs, of an uncertain stage, are not used.
  expect_false(any(n$node == "Ascaphus truei"))
  expect_true("corman2026" %in% ref$sources$source_id)
  expect_true("corman2026" %in% cite_sources()$source_id)
})

test_that("the conversion models hold the triclads and fish of the remaining factors", {
  n <- ref$conversion_nodes
  has <- function(m, sp) any(n$model == m & n$level == "species" & n$node == sp)
  # From the compilation of the remaining conversion factors of 2026-11-10
  # (data-raw/conversions.R, section 14): species without these factors before.
  expect_true(has("invertebrate dry_mass_to_afdm", "Schmidtea lugubris"))
  expect_true(has("invertebrate dry_mass_to_energy", "Planaria torva"))
  expect_true(has("fish dry_mass_to_carbon", "Coregonus albula"))
  expect_true(has("fish dry_mass_to_carbon", "Leuciscus idus"))
  expect_true(all(c("calow1974", "penczak1985") %in% ref$sources$source_id))
  # Farmed animals fed formulated diets are not used.
  expect_false(any(c("wang2021", "gao2024", "seixas2022") %in% ref$sources$source_id))
})

test_that("the sources of the shares of the shell of snails are cited", {
  # They convert the dry masses with shell of the CADSES snails (and the wet
  # masses with shell of the molluscs of Hrycik et al. 2024) before the fits
  # (data-raw/conversions.R, section 11; model_data.R).
  shells <- c("sokolowski2003", "james1998", "kamler1969", "kolodziejczyk1984", "crichton2003")
  expect_true(all(shells %in% ref$sources$source_id))
  expect_true(all(shells %in% cite_sources()$source_id))
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
