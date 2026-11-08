# Data sets -----------------------------------------------------------------------

#' Example survey in the input format of estimate_biomass()
#'
#' Constructed data (not field data) for two sites of a shallow lake, in the
#' input format of [estimate_biomass()] and [check_coverage()]: phytoplankton
#' counted in cells, zooplankton, macroinvertebrates, fish and amphibians,
#' with lengths for some of the animals.
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
#'     stated) and `count_unit` (`"cell"` for phytoplankton, `"individual"`).}
#'   \item{abundance}{40 rows: `site` (`"L1"`, `"L2"`), `taxon_id`,
#'     `abundance` (counting units per `per`; 0 for a taxon looked for and not
#'     seen, `NA` for one not counted) and `per` (`"mL"`, `"L"`, `"m2"`,
#'     `"100 m2"`).}
#'   \item{lengths}{13 rows: `taxon_id`, `site`, `length_mean`, `length`,
#'     `n`, `length_unit` (`"mm"`, `"cm"`), `length_type`, `length_sd` and
#'     `length_cv`. Six mean lengths without site (they apply to both sites),
#'     the mean length of the chub *Squalius pyrenaicus* in each site, and the
#'     lengths of the *Gammarus* measured in L1, by length class (`n`
#'     individuals of each `length`).}
#' }
#' @source Constructed for this package (`data-raw/examples.R`).
#' @examples
#' str(example_survey, max.level = 1)
#' check_coverage(example_survey$taxa, example_survey$lengths)
"example_survey"

#' Raw data of a mesocosm experiment with the fungicide azoxystrobin
#'
#' Files with the counts of phytoplankton, zooplankton and macroinvertebrates
#' of the mesocosm experiment of Amador et al. (2024), as they were entered
#' (one row per sample and one column per taxon), and its design. They are the
#' starting point of `vignette("azoxystrobin", package = "bmassa")`, which
#' puts them in the tables of [estimate_biomass()].
#'
#' The experiment: 24 outdoor mesocosms (140 cm in diameter, 80 cm of water,
#' 1230 L) with the communities of two habitats of the Albufera Natural Park
#' (Valencia, Spain), a eutrophic lake (`W` in the files) and a water spring
#' (`S`), with controls and 2, 20 and 200 µg/L of azoxystrobin (applied twice,
#' 14 days apart), three mesocosms each, sampled 7 days before the first
#' application and 7, 21, 35 and 56 days after it. The files have no
#' identifier of the mesocosm: within a day, a treatment and a concentration,
#' the three rows are the replicates 1 to 3 of `design.csv`, in the same order
#' on every date and in the three files.
#'
#' @format Comma-separated files (UTF-8) in
#'   `system.file("extdata", "azoxystrobin", package = "bmassa")`:
#' \describe{
#'   \item{design.csv}{24 mesocosms: `tank`, `treatment` (`"W"`, `"S"`),
#'     `habitat` (`"eutrophic_lake"`, `"water_spring"`), `concentration_ug_L`
#'     (0, 2, 20, 200), `replicate` (1 to 3), `volume_L`, `depth_m` and
#'     `area_m2`.}
#'   \item{phytoplankton.csv}{120 samples: `Day`, `Treatment`,
#'     `Concentration` and the cells per litre of 45 taxa (units of the
#'     metadata of the original workbook).}
#'   \item{zooplankton.csv}{The same samples, individuals per litre of 34
#'     taxa.}
#'   \item{macroinvertebrates.csv}{The same samples, individuals of 29 taxa
#'     caught in the two traps (30 x 10 cm each) of each mesocosm and date,
#'     together.}
#'   \item{taxa.csv}{The classification of the taxon units of the vignette in
#'     the Catalogue of Life (COL26.10 XR), the result of [standardize_taxa()]
#'     of 2026-11-07, kept so that the vignette does not need internet.}
#' }
#' The names of the taxa are those of the original workbook (misspellings
#' included), except the shrimp *Dugastella valentina*, labelled
#' *Palaemonetes zariquieyi* in the workbook, and the label "Dytiscidae (small
#' larva)", in Spanish there; an unlabelled column of the macroinvertebrates
#' (a single individual) is left out.
#' @source Amador P, Vega C, Navarro Pacheco NI, Moratalla-López J, Palacios J,
#'   Crettaz Minaglia MC, López I, Díaz M, Rico A (2024). Effects of the
#'   fungicide azoxystrobin in two habitats representative of Mediterranean
#'   coastal wetlands: a mesocosm experiment. *Aquatic Toxicology* 267,
#'   106828. \doi{10.1016/j.aquatox.2023.106828}. Licence of the data: CC BY
#'   4.0; cite the article. Prepared with `data-raw/azoxystrobin.R`.
#' @name azoxystrobin
#' @keywords datasets
#' @examples
#' dir <- system.file("extdata", "azoxystrobin", package = "bmassa")
#' list.files(dir)
#' read.csv(file.path(dir, "design.csv"))[1:4, ]
#' read.csv(file.path(dir, "zooplankton.csv"), check.names = FALSE)[1:3, 1:6]
NULL

#' Body masses and abundances of the pelagic community of Tuesday Lake
#'
#' The communities of Tuesday Lake (Michigan, USA) sampled in 1984 and 1986,
#' with the mean body mass and the abundance of each taxon, as distributed in
#' the R package *cheddar* (Hudson et al. 2013), communities TL84 and TL86.
#' They are the test of
#' `vignette("tuesday-lake", package = "bmassa")`: none of these masses is in
#' the reference data of the package.
#'
#' @format Comma-separated files (UTF-8) in
#'   `system.file("extdata", "tuesday_lake", package = "bmassa")`:
#' \describe{
#'   \item{nodes.csv}{113 rows, one per year and taxon: `year`, `node` (name
#'     in *cheddar*), `category` (`"producer"`, `"invertebrate"`,
#'     `"vert.ecto"`), `family`, `M_kg` (mean wet mass of an individual, kg;
#'     in the phytoplankton, of an algal unit, a cell or a colony, from its
#'     volume with a density of 1 kg/L) and `N_m3` (individuals per m3).}
#'   \item{taxa.csv}{The classification of the taxa in the Catalogue of Life
#'     (COL26.10 XR), the result of [standardize_taxa()] of 2026-11-07, with
#'     the names cleaned as in the vignette.}
#' }
#' @source Jonsson T, Cohen JE, Carpenter SR (2005). Food webs, body size, and
#'   species abundance in ecological community description. *Advances in
#'   Ecological Research* 36, 1-84; Cohen JE, Jonsson T, Carpenter SR (2003).
#'   Ecological community description using the food web, species abundance,
#'   and body size. *PNAS* 100, 1781-1786. Through the R package *cheddar*
#'   (Hudson et al. 2013), BSD 2-Clause licence, Copyright (c) 2014-2016
#'   Lawrence Hudson, Dan Reuman and Rob Emerson (see the file `COPYRIGHTS` of
#'   the package). Prepared with `data-raw/tuesday_lake.R`.
#' @references Hudson LN, Emerson R, Jenkins GB, Layer K, Ledger ME, Pichler
#'   DE, Thompson MSA, O'Gorman EJ, Woodward G, Reuman DC (2013). *Cheddar*:
#'   analysis and visualisation of ecological communities in R. *Methods in
#'   Ecology and Evolution* 4(1), 99-104. \doi{10.1111/2041-210X.12005}
#' @name tuesday_lake
#' @keywords datasets
#' @examples
#' nodes <- read.csv(system.file("extdata", "tuesday_lake", "nodes.csv", package = "bmassa"))
#' table(nodes$year, nodes$category)
NULL
