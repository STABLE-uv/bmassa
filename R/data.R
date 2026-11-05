# Data sets -----------------------------------------------------------------------

#' Example survey in the input format of estimate_biomass()
#'
#' Constructed data (not field data) for two sites of a shallow lake, in the
#' input format of [estimate_biomass()] and [check_coverage()]: phytoplankton
#' counted in cells, zooplankton, macroinvertebrates, fish and amphibians,
#' with mean lengths for some of the animals.
#'
#' @format A list of three data frames:
#' \describe{
#'   \item{taxa}{20 taxon units (a taxon, or a taxon in one life stage):
#'     `taxon_id`, `broad_group`, `scientific_name`, `rank`, the
#'     classification `kingdom` to `species` (Catalogue of Life, COL26.6 XR,
#'     with the orders it leaves empty taken from other classifications, and
#'     *Squalius* given in Leuciscidae), `col_id` (identifier of the taxon in
#'     the COL26.6 XR), `life_stage` (`NA` when not stated; the frog
#'     *Pelophylax perezi* is twice, as adults and as larvae), `sex` (not
#'     stated), `count_unit` (`"cell"` for phytoplankton, `"individual"`) and
#'     the mean length of seven units in both sites: `length_mean`,
#'     `length_unit` (`"mm"`, `"cm"`), `length_type`, `length_sd` and
#'     `length_cv`.}
#'   \item{abundance}{40 rows: `site` (`"L1"`, `"L2"`), `taxon_id`,
#'     `abundance` (counting units per `per`; 0 for a taxon looked for and not
#'     seen, `NA` for one not counted) and `per` (`"mL"`, `"L"`, `"m2"`,
#'     `"100 m2"`).}
#'   \item{lengths}{2 rows of mean lengths measured in each site (the chub
#'     *Squalius pyrenaicus*): `taxon_id`, `site`, `length_mean`,
#'     `length_unit`, `length_type` and `length_sd`.}
#' }
#' @source Constructed for this package (`data-raw/examples.R`).
#' @examples
#' str(example_survey, max.level = 1)
#' check_coverage(example_survey$taxa, example_survey$lengths)
"example_survey"
