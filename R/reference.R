# Reference tables of the hierarchical models -----------------------------------

#' Reference tables of the hierarchical models
#'
#' `get_reference()` returns the tables from which [estimate_biomass()]
#' predicts and [convert_currency()] converts: the hierarchical models fitted
#' to the reference data of the package, one per domain and kind of data,
#' the models of the conversion factors between currencies, and the sources of
#' those data. Use it to inspect what a prediction rests on.
#'
#' @return A list of class `bmassa_reference` with these data frames (and a
#'   list):
#'   * `models`: one row per model, with `model`, `domain` (`"phytoplankton"`,
#'     `"invertebrate"`, `"fish"`, `"amphibian"`), `kind` (`"typical"`, masses
#'     without lengths, or `"lwr"`, length-mass relations), `count_unit`,
#'     `life_stage` (the stage of the data, if only one), `currency`, `unit`
#'     (`"um3"`, `"mg"` or `"g"`) and `levels` (taxonomic levels of the model,
#'     comma-separated, from the highest).
#'   * `lwr_types`: one row per length-mass model and length type, with
#'     `centre` (mean ln length, mm, of the data), `length_min` and
#'     `length_max` (mm), and `sigma_ind`, the median residual standard
#'     deviation (ln) of the length-mass equations fitted from individual
#'     data, used by `bias_correction = "lognormal"` (`n_sigma` equations).
#'   * `nodes`: one row per node of each model: `level` (`"root"`, a
#'     taxonomic level, or `"stage"`), `length_type` (length-mass models),
#'     `node` (name), `col_id` (identifier of the name at that rank in the
#'     Catalogue of Life eXtended Release), `path` and `parent` (names from
#'     the kingdom down, separated by `|`), `n_units` (records or equations
#'     below the node),
#'     `sources`, `length_min` and `length_max` (mm; range of lengths of the
#'     equations below the node), the estimate of the path from the root to
#'     the node (`path_intercept`; `path_slope` in length-mass models, the
#'     slope of ln mass on ln length), its conditional variances and
#'     covariance (`path_var`, `path_var_slope`, `path_cov`), the covariances
#'     with the path of the parent (`parent_cov`, `parent_cov_is`,
#'     `parent_cov_si`, `parent_cov_ss`), the effect of the node itself
#'     (`effect`, `effect_slope`; `effect_var` for stage effects), and, in
#'     the model of invertebrates without lengths, the covariance of the path
#'     with the fixed contrasts of copepodites and nauplii
#'     (`contrast_cov_copepod_copepodite`, `contrast_cov_copepod_nauplius`).
#'   * `components`: variance components (`value`, natural logs) by `model`,
#'     `level` (taxonomic levels, `"equation"`, `"study"`, `"source"`,
#'     `"stage"`, `"residual"`) and `parameter` (`"intercept"`, `"slope"`,
#'     `"intercept_slope"`).
#'   * `fixed`: fixed effects by model, with standard errors and covariances
#'     (intercepts and slopes by length type; in the model of invertebrates
#'     without lengths, the contrasts `c_copepod_copepodite` and
#'     `c_copepod_nauplius` against adult copepods).
#'   * `conversion_models`: one row per model of the conversion factors (ln
#'     of the ratio between two currencies of the same organisms), with
#'     `model`, `domain` (`"invertebrate"` or `"fish"`), `from` and `to`
#'     (currencies), `levels` and `context` (crossed effects: `habitat`,
#'     `source`).
#'   * `conversion_nodes`, `conversion_components`, `conversion_fixed`: as
#'     `nodes`, `components` and `fixed`, for the conversion models; the
#'     nodes of level `"habitat"` hold the habitat effects.
#'   * `carbon_equations`: the carbon-volume equations of Menden-Deuer and
#'     Lessard (2000, Table 4) for diatoms and the other protists, with the
#'     95 % confidence intervals of their coefficients.
#'   * `sources`: `source_id`, `citation`, `doi`, `licence`.
#'   * `info`: versions (`data_version`, `col_release` of the Catalogue of
#'     Life eXtended Release used for the taxonomy, `fitted_on`,
#'     `conversions_fitted_on`, `lme4`).
#'
#' @details
#' The models were fitted outside the package (`data-raw/fit_models.R` and
#' `data-raw/fit_conversions.R` in the source repository) with lme4 by
#' REML: masses (or biovolumes) without lengths, and conversion factors, with
#' nested random intercepts by taxonomic level, and length-mass relations with
#' random intercepts and slopes by level within length type, from three points
#' of each published or fitted equation inside its range of lengths. The estimate and variance of each path come from the joint
#' conditional distribution of the fixed and random effects, so the effects of
#' a node and of its ancestors are not treated as independent. Records
#' identified above species have nodes of their own in the fit, which are not
#' kept here because no user taxon can match them.
#'
#' @export
#' @examples
#' ref <- get_reference()
#' ref
#' ref$models[, c("model", "kind", "currency", "unit", "levels")]
#' head(ref$nodes[ref$nodes$level == "genus", c("model", "node", "n_units", "sources")])
get_reference <- function() {
  structure(ref_tables, class = "bmassa_reference")
}

#' @export
print.bmassa_reference <- function(x, ...) {
  cat("<bmassa_reference> ", x$info$data_version, ": ", nrow(x$models), " models fitted on ",
      x$info$fitted_on, " (lme4 ", x$info$lme4, "); taxonomy of the ", x$info$col_release,
      "\n", sep = "")
  n <- x$nodes[x$nodes$level %in% rank_levels, ]
  for (i in seq_len(nrow(x$models))) {
    m <- x$models[i, ]
    nm <- n[n$model == m$model, ]
    sp <- nm[nm$level == "species", ]
    cat(sprintf("  %-27s %-7s %-9s %-3s  %4d species, %4d genera; %s\n", m$model, m$kind, m$currency,
                m$unit, nrow(sp), sum(nm$level == "genus"),
                paste(sort(unique(unlist(strsplit(nm$sources, ";")))), collapse = ", ")))
  }
  cm <- x$conversion_models
  if (!is.null(cm) && nrow(cm)) {
    cat("Conversion models (fitted on ", x$info$conversions_fitted_on, "):\n", sep = "")
    cn <- x$conversion_nodes
    for (i in seq_len(nrow(cm))) {
      nm <- cn[cn$model == cm$model[i] & cn$level %in% rank_levels, ]
      cat(sprintf("  %-34s %4d species; %s\n", cm$model[i], sum(nm$level == "species"),
                  paste(sort(unique(unlist(strsplit(nm$sources, ";")))), collapse = ", ")))
    }
    cat("Carbon of phytoplankton cells: Menden-Deuer and Lessard (2000), diatoms and other protists\n")
  }
  invisible(x)
}

# Reference given to a function: the default one or an object of the same
# class.
check_reference <- function(reference, arg = "reference", call = sys.call(-1)) {
  if (!inherits(reference, "bmassa_reference")) {
    abort_input(sprintf(
      "`%s` must be the result of get_reference(), not %s.", arg, describe_value(reference)
    ), call = call)
  }
  invisible(reference)
}
