# Totals by group ------------------------------------------------------------------

#' Total biomass by site or any other grouping
#'
#' `summarize_biomass()` adds up the biomass of the rows of each group, with
#' one total per currency and unit of biomass ([estimate_biomass()] gives all
#' its rows in one currency and unit; abundances per volume and per area stay
#' apart unless put on one basis with [convert_abundance()]). It gives the
#' mean, median and prediction interval of each total by simulation, keeping
#' the correlations between rows.
#'
#' @param data Result of [estimate_biomass()] or [convert_currency()].
#' @param by Character vector with the names of the grouping columns, e.g.
#'   `"site"` (default) or `c("site", "broad_group")`.
#' @param na_rm If `FALSE` (default), a group with rows whose biomass is not
#'   known gets `NA` totals and `status = "incomplete"`. If `TRUE`, those rows
#'   are left out: the total is then the biomass of the estimated rows only (a
#'   lower bound) and `n_missing` says how many were left out.
#' @param level Level of the prediction intervals (default 0.95).
#' @param n_draws Number of Monte Carlo draws (default 1000).
#' @param population_effect How the variation between populations is shared
#'   between the sites: `"shared"` (default; one draw per taxon and life stage
#'   for all its rows, as the sites of a study usually share populations and
#'   season) or `"independent"` (one draw per row). The total of one site does
#'   not change; sums over several sites and differences between sites do.
#' @param population_by `NULL` (default) or the names of columns of the
#'   abundance table (they are kept in the estimates) whose values separate
#'   populations, such as a habitat or a region: with
#'   `population_effect = "shared"`, the rows of a taxon and stage with the
#'   same values share their population, and those with different values have
#'   independent ones. Not allowed with `"independent"`; the columns must not
#'   have missing values.
#' @param seed `NULL` (default: the current random state is used) or an
#'   integer: the draws are then made with that seed and the random state of
#'   the session is restored afterwards (`withr::with_seed()`).
#' @param ... Not used; an error if anything is given.
#'
#' @return A data frame with one row per combination of the `by` columns,
#'   currency and unit, sorted by those columns (C locale order), with
#'   columns `by...`, `currency`, `biomass_unit`,
#'   `biomass_mean` (expectation of the total), `biomass_median`,
#'   `biomass_lower`, `biomass_upper` (prediction interval), `level`,
#'   `interval` (`"prediction"`), `n_rows`, `n_estimated`, `n_missing`,
#'   `status` (`"complete"` or `"incomplete"`), `n_draws`,
#'   `population_effect` and `population_by` (the columns, separated by
#'   commas, or `NA`). A group without any estimated row appears once
#'   with `currency = NA`. A zero-row input gives a zero-row result.
#'
#' @details
#' **Simulation.** Each draw takes, once, the error of every node of the
#' reference used by the rows, along the tree from the root (the error of a
#' node is drawn given that of its parent, from their joint covariance); the
#' deviation of every taxonomic level missing from the reference, once per
#' missing node (two species of a genus absent from the reference share the
#' deviation of the genus); the variation between populations, once per taxon
#' and stage (`population_effect = "shared"`; within each group of
#' `population_by`, if given) or per row; the error of the effect of a stage
#' group, once per group, and that of the contrast of copepodites or nauplii,
#' once per contrast and correlated with the error of the node of each row;
#' and the deviation of a new source or study, once per model, since all the
#' rows come from the same study. Two taxa
#' estimated from the same family therefore have correlated errors, and the
#' totals reflect it. Each draw of the mean individual mass is multiplied by
#' the abundance and the draws are added up by group; `biomass_median` and the
#' interval are quantiles of those draws. The draws are those of
#' [simulate()][simulate.bmassa_estimate()], which returns them row by row.
#' `biomass_mean` is not simulated: it is the sum of the `biomass_mean` of the
#' rows, which is exact because the expectation of a sum does not depend on
#' the correlations.
#'
#' **Monte Carlo error.** The median and the limits of the interval change
#' slightly with the seed. In the example data, with 1000 draws, their
#' standard deviation between seeds was up to 5 % of the median and 1-12 % of
#' each limit (the wider the interval, as for taxa without lengths, the
#' larger), and with 10000 draws about a third of that. Use more draws (e.g. `n_draws = 10000`) for final
#' results; time and memory grow with `n_draws` times the number of rows.
#'
#' @section Missing values and zeros:
#' Rows without biomass (taxa not estimated or abundances missing) cannot be
#' assigned to a currency, so they make every total of their group
#' incomplete. Include in `by` a column that separates organism groups (e.g.
#' `broad_group`) if a missing fish should not affect the phytoplankton total.
#' Rows with a zero abundance and known size count as zero.
#'
#' @section Reproducibility:
#' With the same `seed`, the draws are the same on the same platform; other
#' versions of R or other BLAS and LAPACK libraries may give slightly
#' different values.
#'
#' @export
#' @examples
#' est <- estimate_biomass(example_survey$abundance, example_survey$taxa,
#'                         example_survey$lengths, bias_correction = "lognormal")
#' summarize_biomass(est, by = c("site", "broad_group"), na_rm = TRUE, seed = 1)
summarize_biomass <- function(data, ...) {
  UseMethod("summarize_biomass")
}

#' @rdname summarize_biomass
#' @export
summarize_biomass.default <- function(data, ...) {
  abort_input(sprintf(
    "`data` must be the result of estimate_biomass() or convert_currency(), not %s.",
    describe_value(data)
  ))
}
