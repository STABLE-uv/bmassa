# Sources of the reference data -----------------------------------------------

#' Sources to cite for an estimate
#'
#' `cite_sources()` lists the publications and data sets behind the
#' reference data: those under the nodes used by the rows of an estimate (and
#' by its currency conversions), or all of them. Cite them together with the
#' package.
#'
#' @param x `NULL` (default: all the sources of the reference) or the result of
#'   [estimate_biomass()] or [convert_currency()].
#' @param reference Reference tables, from [get_reference()] (default);
#'   used when `x` is `NULL`.
#'
#' @return A data frame with one row per source, sorted by `source_id`:
#'   `source_id`, `citation`, `doi`, `licence` (of the data of the source)
#'   and `n_rows` (rows of `x` whose estimate or conversion rests on data of
#'   the source; `NA` when `x` is `NULL`).
#'
#' @details
#' An estimate rests on the data below the node it comes from; the variance
#' components of the model use all of its data. Some sources are compilations:
#' van Leeuwen et al. (2025) and Moreau et al. (2021) give the original
#' reference of each value. *bmassa* holds only the parameters of the models
#' fitted to these data (and the equations published in the sources), not the
#' data themselves.
#'
#' @export
#' @examples
#' cite_sources()[, c("source_id", "licence")]
#' est <- estimate_biomass(example_survey$abundance, example_survey$taxa,
#'                         example_survey$lengths, bias_correction = "lognormal")
#' cite_sources(est)[, c("source_id", "n_rows")]
cite_sources <- function(x = NULL, reference = get_reference()) {
  if (is.null(x)) {
    check_reference(reference)
    out <- reference$sources
    out$n_rows <- rep(NA_integer_, nrow(out))
  } else {
    if (!inherits(x, "bmassa_estimate")) {
      abort_input(sprintf("`x` must be NULL or the result of estimate_biomass(), not %s.",
                          describe_value(x)))
    }
    e <- x$estimates
    ok <- !is.na(e$ind_median)
    per_row <- lapply(which(ok), function(i) {
      s <- c(if (!is.na(e$sources[i])) strsplit(e$sources[i], ";", fixed = TRUE)[[1]],
             if (!is.null(e$conversion_sources) && !is.na(e$conversion_sources[i]))
               strsplit(e$conversion_sources[i], ";", fixed = TRUE)[[1]])
      unique(s)
    })
    n <- table(unlist(per_row))
    out <- x$reference$sources[x$reference$sources$source_id %in% names(n), ]
    out$n_rows <- as.integer(n[out$source_id])
  }
  out <- out[order(out$source_id, method = "radix"), ]
  rownames(out) <- NULL
  out
}
