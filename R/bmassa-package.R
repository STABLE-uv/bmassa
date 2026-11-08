#' *bmassa*: Biomass of Aquatic Taxa from Abundances and Body Sizes
#'
#' Estimates the mean individual mass and the biomass of aquatic taxa
#' (phytoplankton, zooplankton, macroinvertebrates, fish and amphibians) from
#' abundance data, with prediction intervals, from hierarchical models fitted
#' to the reference data of the package.
#'
#' Workflow:
#' * [standardize_taxa()]: classification and identifier of each taxon in the
#'   Catalogue of Life (needs internet).
#' * [check_coverage()]: from which taxonomic level each taxon will be
#'   estimated, and how wide its interval will be.
#' * [convert_abundance()]: the abundances on one basis (per volume, per area
#'   or per site), so that the biomasses of all the groups can be added up.
#' * [estimate_biomass()]: mean individual mass of a population of each taxon
#'   and biomass of each row, from the taxon units (with their stages, sexes
#'   and counting units), their abundances by site and date and their lengths,
#'   in one currency (wet mass by default).
#' * [convert_currency()]: the estimates in wet mass, dry mass, ash-free dry
#'   mass, carbon or energy, with the uncertainty of the conversion factors.
#' * [summarize_biomass()]: totals by site or group with prediction intervals.
#' * [summarize_size()]: size structure by site or group: mean and median
#'   size of the individuals and shares of individuals and biomass by size
#'   class.
#' * [simulate()][simulate.bmassa_estimate()]: joint draws of the mass and the
#'   biomass of every row, for models that use them (such as the energy fluxes
#'   of a food web).
#' * [cite_sources()] and [get_reference()]: the sources to cite and the
#'   reference tables.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
