# Biomass from the hierarchical reference models (version 0.2.0) ---------------

broad_group_levels <- c("phytoplankton", "zooplankton", "macroinvertebrate", "fish", "amphibian")
broad_group_domain <- c(phytoplankton = "phytoplankton", zooplankton = "invertebrate",
                        macroinvertebrate = "invertebrate", fish = "fish", amphibian = "amphibian")
main_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
length_unit_mm <- c(um = 0.001, mm = 1, cm = 10)
length_type_levels <- c("body_length", "shell_length", "shell_height", "fork_length",
                        "total_length", "standard_length", "snout_vent_length")
bias_correction_levels <- c("none", "lognormal")

# Life stages and sexes of the taxon units (version 0.2.0). Internally, an adult
# (or unstaged) unit with a sex is coded by its sex, as in the reference data.
unit_stage_levels <- c("adult", "juvenile", "copepodite", "nauplius", "larva", "nymph", "pupa")
sex_levels <- c("female", "male")
length_input_columns <- c("length_mean", "length_unit", "length_type", "length_sd", "length_cv")

# Stage groups of the insects in the model of masses without lengths; other
# invertebrates form one group. An insect without a stage is taken as a larva
# (decision of 2026-11-02: aquatic samples hold larvae; adults must
# be declared).
insect_stage_group <- c(adult = "adult", female = "adult", female_ovigerous = "adult",
                        male = "adult", juvenile = "larva", copepodite = "larva",
                        nauplius = "larva", larva = "larva", nymph = "larva", pupa = "pupa")
# Stages of the copepods in the same model: the group of the other
# invertebrates ("non_insect") and, for copepodites and nauplii, a fixed
# contrast against the adults (`copepod_contrast`). A copepod without a stage
# is taken as an adult (counts of copepods without a stage are usually of
# adults; copepodites and nauplii are counted apart); "juvenile" is either
# copepodites or nauplii, and is not estimated.
copepod_stage_group <- c(adult = "copepod_adult", female = "copepod_adult",
                         female_ovigerous = "copepod_adult", male = "copepod_adult",
                         copepodite = "copepod_copepodite", nauplius = "copepod_nauplius")
copepod_contrast <- c(copepod_copepodite = "copepod_copepodite", copepod_nauplius = "copepod_nauplius")
# Orders of copepods, for classifications that do not give the class Copepoda
# (Maxillopoda, Hexanauplia).
copepod_orders <- c("Calanoida", "Cyclopoida", "Harpacticoida", "Siphonostomatoida",
                    "Poecilostomatoida", "Misophrioida", "Monstrilloida", "Mormonilloida",
                    "Platycopioida", "Gelyelloida", "Canuelloida")

# Columns added by estimate_biomass().
estimate_biomass_columns <- c(
  "broad_group", "life_stage", "sex", "count_unit", "domain", "method", "model", "stage_used",
  "length_mm", "length_type", "length_basis", "length_n", "length_status", "match_level", "match_name",
  "n_reference", "sources", "currency", "unit",
  "ind_mean", "ind_median", "ind_lower", "ind_upper", "biomass_mean", "biomass_median",
  "biomass_lower", "biomass_upper", "biomass_unit", "level", "interval", "uncertainty",
  "corrections", "status", "reason", "note"
)

#' Individual mass and biomass of each taxon, with prediction intervals
#'
#' `estimate_biomass()` gives, for each row of an abundance table, the mean
#' individual mass of a population of its taxon (mean, median and prediction
#' interval) and the biomass, from hierarchical models fitted to the
#' reference data of the package: phytoplankton biovolume, dry mass of
#' invertebrates (without shell in molluscs) and wet mass of fish and
#' amphibians. With lengths of the taxon it uses the length-mass models;
#' otherwise, the models of mass without lengths. A taxon absent from the
#' reference is predicted from its closest relatives, and the interval widens
#' accordingly. All the rows are given in one currency and unit (`currency`,
#' wet mass in mg by default), converted with the uncertainty of the factors
#' as in [convert_currency()].
#'
#' @param abundance A data frame with one row per sample and taxon unit, with
#'   columns:
#'   * `site` (character): site or sample.
#'   * `date` (optional; any type, such as dates or days of an experiment):
#'     with it, a site is sampled on several dates, and each site and date is
#'     a sample.
#'   * `taxon_id`: key of the taxon unit in `taxa`.
#'   * `abundance` (numeric, >= 0 or `NA`): counting units per `per`.
#'   * `per` (character): what the abundance refers to, the size of the
#'     sample: a volume or an area, with or without a number before the unit
#'     (`"mL"`, `"L"`, `"m3"`, `"cm2"`, `"m2"`, `"ha"`), or `"site"` for counts
#'     of a whole site. Write what was sampled: the individuals in a trap of
#'     30 x 10 cm are per `"300 cm2"`, those of two such traps counted
#'     together per `"600 cm2"`, those in 10 L of filtered water per
#'     `"10 L"`, a density per litre per `"L"`. The biomass is given per this
#'     basis. To add up groups sampled differently, or to get the totals of a
#'     site from the share of it that was sampled, put them on one basis
#'     first with [convert_abundance()].
#'   Other columns (treatment, habitat...) are kept. The life stage, the sex
#'   and the counting unit belong to the taxon unit, in `taxa`.
#' @param taxa A data frame with one row per taxon unit (a taxon, or a taxon
#'   in one life stage and sex): `taxon_id` (unique), `broad_group`
#'   (`"phytoplankton"`, `"zooplankton"`, `"macroinvertebrate"`, `"fish"` or
#'   `"amphibian"`) and the classification in the columns `kingdom`,
#'   `phylum`, `class`, `order`, `family`, `genus` and `species` (the full
#'   binomial), `NA` below the rank identified, as returned by
#'   [standardize_taxa()]. Optional columns:
#'   * `col_id`: the identifier of the taxon in the Catalogue of Life (COL
#'     XR), used first to find it in the reference.
#'   * `life_stage`: `"adult"`, `"juvenile"`, `"copepodite"`, `"nauplius"`,
#'     `"larva"`, `"nymph"` or `"pupa"`; `NA` when not stated.
#'   * `sex`: `"female"` or `"male"` (adults); `NA` when not stated.
#'   * `count_unit`: `"cell"` or `"natural_unit"` (colony or filament as
#'     counted) for phytoplankton, where it is required; `"individual"` or
#'     `NA` for animals.
#'   * `cells_per_unit` (>= 1): for phytoplankton counted in natural units,
#'     the mean number of cells of a unit, where it was counted; the unit is
#'     then estimated as that many cells (see Details). `NA` otherwise.
#'   Other columns, such as `scientific_name`, are kept, and
#'   `scientific_name` is added to the result. Mean lengths in `taxa`
#'   (`length_mean`, `length_unit`, `length_type`, `length_sd`, `length_cv`)
#'   are deprecated: they are still used, after any row of `lengths`, with a
#'   warning; give them in `lengths`.
#' @param lengths `NULL` (default) or a data frame with the lengths of the
#'   taxon units, one row per mean length or per measured length:
#'   * `taxon_id`: key of the taxon unit in `taxa`.
#'   * `site` and `date` (optional): where the lengths apply. A row without
#'     `site` applies to all the samples of the unit, one with `site` and
#'     without `date` to all the dates of that site, and one with both to that
#'     sample. For each sample the most specific lengths are used.
#'   * Either `length_mean` (> 0), the mean length of the unit, the usual case,
#'     with its spread in `length_sd` or `length_cv` when known; or `length`
#'     (> 0), the length of measured individuals, with `n` (default 1) the
#'     number of individuals of that length: one row per individual, or per
#'     length class with its count. A unit and sample have one or the other.
#'   * `length_unit` (`"um"`, `"mm"` or `"cm"`) and `length_type`
#'     (`"body_length"`, `"shell_length"`, `"shell_height"`, `"fork_length"`,
#'     `"total_length"`, `"standard_length"` or `"snout_vent_length"`), one
#'     type per unit and sample.
#' @param currency Currency of all the rows of the result: `"wet_mass"`
#'   (default; the only one that all the groups reach), `"dry_mass"`,
#'   `"afdm"`, `"carbon"` or `"energy"`. Rows that cannot be converted to it
#'   (phytoplankton to dry mass or AFDM, amphibians to AFDM or carbon,
#'   phytoplankton counted in natural units without `cells_per_unit` to carbon
#'   or energy) are not estimated, with the reason. See [convert_currency()].
#' @param unit `NULL` (default: mg for masses, J for energy) or the unit of
#'   all the rows: `"pg"`, `"ng"`, `"ug"`, `"mg"` or `"g"` for masses;
#'   `"mJ"`, `"J"`, `"kJ"` or `"MJ"` for energy.
#' @param habitat `NULL` (default: unknown) or the habitat of the organisms,
#'   `"freshwater"`, `"marine"` or `"transitional"`, for the conversion models
#'   that have a habitat effect (see [convert_currency()]).
#' @param bias_correction Back-transformation of the length-mass relations,
#'   required when lengths are used: `"none"` (the relation gives the median
#'   mass of the individuals of a given length) or `"lognormal"` (their mean:
#'   multiplied by \eqn{\exp(\sigma^2/2)}{exp(sigma^2/2)}, Sprugel 1983, with
#'   the typical residual standard deviation of the relations fitted from
#'   individual data; `sigma_ind` in [get_reference()]). No default:
#'   there is no general criterion to choose; for biomass, which needs the
#'   mean mass, `"lognormal"` is the consistent choice.
#' @param level Level of the prediction intervals (default 0.95).
#' @param extrapolate If `FALSE` (default), a length outside the lengths of
#'   the reference relations of the taxon and of all its relatives is not
#'   estimated; if `TRUE`, it is estimated from the closest relative and
#'   marked `length_status = "outside"`.
#' @param reference Reference tables, from [get_reference()] (default).
#'
#' @return An object of class `bmassa_estimate`: a list with `estimates`, a
#'   data frame with the rows of `abundance` (same order) and these columns
#'   added (use `as.data.frame()` to get it):
#'   * `broad_group`, `scientific_name` (if in `taxa`), `life_stage`, `sex`
#'     and `count_unit` of the taxon unit, and `domain` (`"phytoplankton"`,
#'     `"invertebrate"`, `"fish"`, `"amphibian"`).
#'   * `method` (`"length"` or `"typical"`) and `model` (see
#'     [get_reference()]).
#'   * `stage_used`: stage group of the model of masses without lengths:
#'     for insects `"larva"`, `"pupa"`, `"adult"`, or `"larva (assumed)"`
#'     when the stage is not given; for copepods `"copepod_adult"`,
#'     `"copepod_copepodite"`, `"copepod_nauplius"`, or
#'     `"copepod_adult (assumed)"` when the stage is not given (a copepod
#'     given as `"juvenile"`, `"larva"`, `"nymph"` or `"pupa"` is not
#'     estimated); `"non_insect"` for the other invertebrates; `"adult"` for
#'     amphibians; `NA` otherwise (and for estimates from lengths, whose
#'     relations have no stage effect; the lengths of copepod nauplii use
#'     relations of their own).
#'   * `length_mm`, `length_type`, `length_basis` and `length_n`: the length
#'     used (mm), its type, whether it is a mean length (`"mean"`) or comes
#'     from measured individuals (`"individuals"`: then `length_mm` is the
#'     length of their mean mass and `length_n` their number), and
#'     `length_status` (`"within"` the lengths of the reference relations
#'     used, `"outside"`, or `NA`).
#'   * `match_level`, `match_name`: the node of the reference the estimate
#'     comes from (the lowest taxonomic level shared with the reference), and
#'     `n_reference` (records or equations below it) and `sources`.
#'   * `currency` and `unit`, the same in all the rows, and, for the rows
#'     whose models are in another currency, `converted_from`, `conversion`
#'     (the factors used) and `conversion_sources`, as in [convert_currency()].
#'   * `ind_mean`, `ind_median`, `ind_lower`, `ind_upper`: mean, median and
#'     prediction interval of the mean individual mass of a population of the
#'     taxon, in `unit`.
#'   * `biomass_mean`, `biomass_median`, `biomass_lower`, `biomass_upper`:
#'     `abundance` times the above, in `biomass_unit` (`unit` per `per`).
#'   * `level`, `interval` (`"prediction"`), `uncertainty` (components
#'     included), `corrections` (`"jensen"`, `"lognormal"`).
#'   * `status` (`"estimated"` or `"not_estimated"`), `reason` (why not) and
#'     `note` (assumptions and substitutions of an estimate).
#'
#'   The list also holds `settings` (the arguments used and the version of the
#'   reference) and `diagnostics` (counts by method, level and reason), and
#'   what [summarize_biomass()] needs to simulate totals. A zero-row
#'   `abundance` gives zero rows.
#'
#' @details
#' **Model.** For a taxon matched to node \eqn{j} of the reference, the log of
#' the mean individual mass of a new population is predicted as
#' \deqn{\ln m \sim N(\eta_j,\ \tau^2), \quad \tau^2 = V_j + \sum_{\ell\ \mathrm{missing}} \sigma^2_\ell + \sigma^2_{pop} + \sigma^2_{source}}{ln m ~ N(eta_j, tau^2), tau^2 = V_j + sum(sigma2_l, missing levels) + sigma2_pop + sigma2_source}
#' where \eqn{\eta_j} is the estimate of the path from the root of the model to
#' the node (fixed effects plus the random effects of the node and its
#' ancestors), \eqn{V_j} its conditional variance, the sum runs over the
#' taxonomic levels below the node (the taxon's own genus or species when they
#' are not in the reference), \eqn{\sigma^2_{pop}} is the variation between
#' populations (records or equations of the same species; in the model of
#' fish without lengths, the variation between populations plus that between
#' samplings of one population, separated by the populations sampled many
#' times) and \eqn{\sigma^2_{source}} that between sources or studies (a
#' user's data are a new one). The models of invertebrates and fish without
#' lengths have no effect of source, so there the variation between sources is
#' part of that between populations. The length-mass relations of
#' invertebrates are fitted without an effect of study, so that each species
#' keeps the level of its own equations, and \eqn{\sigma^2_{source}} is the
#' variance between studies of a second fit with that effect, taken out of the
#' variance between equations: \eqn{\tau^2} is that of the fit without study,
#' and in [summarize_biomass()] the part of a new study is shared by all the
#' rows of the model. With a length \eqn{L} (mm), \eqn{\eta_j = \alpha_j + \beta_j (\ln L -
#' c)}{eta_j = alpha_j + beta_j (ln L - c)} and the variances are evaluated at
#' \eqn{\ln L}. The median is \eqn{\exp(\eta)}{exp(eta)}, the mean
#' \eqn{\exp(\eta + \tau^2/2)}{exp(eta + tau^2/2)} and the interval
#' \eqn{\exp(\eta \pm z\tau)}{exp(eta +/- z tau)}. Totals over taxa or sites
#' need the correlations between rows: use [summarize_biomass()].
#'
#' **Stages.** Without lengths, the model of invertebrates has random effects
#' of stage group for insects (larva, pupa, adult; without a stage, taken as
#' larva) and one group for the other invertebrates. Copepods are in that
#' group, and their copepodites and nauplii add a fixed contrast against the
#' adults (a copepod without a stage is taken as an adult): its estimate, its
#' variance and twice its covariance with the path enter \eqn{\eta}{eta} and
#' \eqn{V_j}{V_j}. Most records of nauplii are of an order, so their contrast
#' is wide. A count of a species that includes its copepodites is therefore
#' overestimated when given without a stage (1.4 to 3.6 times for the
#' copepods of Tuesday Lake, whose means seem to include them): give the
#' copepodites as rows of their own. With lengths there is no stage effect, but the body length of
#' copepod nauplii has relations of its own (length type
#' `"nauplius_body_length"`, used for rows of copepods with `life_stage =
#' "nauplius"` and `length_type = "body_length"`).
#'
#' **Natural units of phytoplankton.** Without `cells_per_unit`, a natural
#' unit is predicted from the natural units of the reference, which are
#' those of its sources: in colonial genera often single cells or small
#' groups, not whole colonies (a note says so). With `cells_per_unit` = k, the
#' biovolume of the unit is that of a cell of the taxon times k:
#' \eqn{\eta}{eta} of the model of cells plus \eqn{\ln k}{ln k}, with the
#' variances of the cell (k is taken as known).
#'
#' **Fish without lengths.** The model of fish without lengths is fitted to
#' the mean wet mass of whole populations (all the sizes caught): species and
#' site of an electrofishing survey of Iberian streams (Arranz et al. 2022),
#' species in the catch of a river fishery (Koning et al. 2024), and species
#' in each sampling (season and year) and gear of the streams and lakes of the
#' National Ecological Observatory Network of the United States (NEON,
#' electrofishing and fyke and gill nets, 2016-2024), whose samplings of a
#' species at a site are one population with a random effect. The mean of a
#' population depends on its size structure, which the variation between
#' populations and between samplings takes up, so the intervals are wide; a
#' mean length gives a closer estimate. A fish of an order absent from these
#' data is predicted from the fish of its class, and one of a class absent
#' too from all their fish (the root of the model), with the variances of the
#' levels below; a fish with `life_stage` or `sex` is not estimated without a
#' length.
#'
#' **Matching.** A taxon with a `col_id` is first looked for by it (the
#' identifier follows a taxon whose name or genus changed). Otherwise, from
#' the lowest rank of the taxon upwards, the first node of the model with the
#' same name at the same rank (homonyms: the one whose higher ranks agree most
#' with the taxon). A taxon that matches no node is not estimated. With a length, the node used is the lowest one whose
#' reference relations cover that length.
#'
#' **Lengths.** If \eqn{L} is a mean length and the spread of lengths is
#' given, the mean mass is multiplied by \eqn{(1 + CV^2)^{\beta(\beta-1)/2}}{(1 + CV^2)^(beta (beta - 1)/2)}
#' (lognormal lengths; `corrections = "jensen"`); without the spread the mean
#' mass is underestimated (by 29 % with \eqn{\beta = 3}{beta = 3} and CV = 0.3),
#' which is noted. With the lengths \eqn{L_i} of \eqn{n_i} measured
#' individuals, the mean mass is the mean of their masses, which under the
#' relation is that at the length
#' \eqn{L_* = (\sum n_i L_i^\beta / \sum n_i)^{1/\beta}}{L* = (sum n_i L_i^beta / sum n_i)^(1/beta)}
#' (`length_mm`), with the uncertainty of the relation at that length; the
#' node used is the lowest one whose relations cover all their lengths.
#'
#' **Currency.** The models give phytoplankton biovolume, dry mass of
#' invertebrates and wet mass of fish and amphibians; every row is then
#' converted to `currency` and `unit` as [convert_currency()] does, so that
#' all the rows can be added up and compared. Biovolume goes to wet mass with
#' the density of water (1 mm\out{<sup>3</sup>} = 1 mg).
#'
#' @section Missing values and zeros:
#' An abundance of 0 gives a biomass of 0 when the mass is estimated; a
#' missing abundance gives a missing biomass (`status = "not_estimated"`,
#' `reason = "abundance is missing"`) but the individual mass is kept. Rows
#' without an estimate keep `NA`, with the reason, and one warning reports how
#' many there are. Amphibians given as larvae or juveniles (with or without
#' lengths), amphibians without lengths that are not declared adults, and
#' fish without lengths given a life stage or a sex, are not estimated: the
#' reference has no masses for them.
#'
#' @section Assumptions:
#' The interval is for the mean individual mass of a population, not for one
#' individual, and includes the uncertainty of the model estimates, of the
#' taxonomic substitution and the variation between populations and between
#' sources. It does not include the error of the abundances, the uncertainty
#' of the structure of the models, nor that of the conversion factors applied
#' to some reference data before fitting (e.g. wet to dry mass; only their
#' central values were used). Masses without
#' lengths vary much more between populations than masses at a given length:
#' in cross-validation, the median 95 % interval of a new population of a
#' species of the reference spans a factor of 7 to 38 either way without
#' lengths, against 1.2 to 3.3 with them (`dev/validation` in the source
#' repository): measure lengths where possible. Insects without lengths and
#' without a stage are taken as larvae, and copepods as adults. The
#' length-mass relations of insects come mostly from larvae. The masses of
#' rotifers come from direct weighings of dry mass (Dumont et al. 1975;
#' Bottrell et al. 1976, Table VIII) and from volumes (wet mass, converted to
#' dry mass with the convention of 10 %, 3.9 % in *Asplanchna*). The weighings
#' are a fixed contrast in the model (about four times the masses from
#' volumes), and rotifers are predicted as from volumes, the usual method of
#' field studies. The model of invertebrates without lengths has no effect of
#' source: its sources measure different groups of taxa by different methods,
#' and each taxon keeps the level of the sources that measured it, except for
#' one difference of method. The benthic invertebrates of Moreau et al.
#' (2021) are masses of individual prey compiled from studies of diet and
#' calorimetry, often of large or adult animals, and weigh about 6 times
#' (95 % interval 3.6-11) the means of all the individuals of field samples of
#' Hrycik et al. (2024) for the same taxa; they are a fixed contrast in the
#' model, and benthic invertebrates are predicted as population means.
#'
#' @references
#' Arranz I, Brucet S, Bartrons M, et al. (2022). Individual body mass and
#' length dataset for over 12,000 fish from Iberian streams. *Data in Brief*,
#' 108248. \doi{10.1016/j.dib.2022.108248}
#'
#' Koning A, Pin K, Chhuoy S, et al. (2024). Weight and length measurements
#' for 37 Mekong River fish species of the 'Dai' fishery. Dryad.
#' \doi{10.5061/dryad.dr7sqvb5h}
#'
#' NEON (National Ecological Observatory Network) (2026). Fish
#' electrofishing, gill netting, and fyke netting counts (DP1.20107.001),
#' RELEASE-2026. \doi{10.48443/s1e6-df79}
#'
#' Sprugel DG (1983). Correcting for bias in log-transformed allometric
#' equations. *Ecology*, 64, 209-210.
#'
#' The sources of the reference data are listed by [cite_sources()].
#'
#' @seealso [check_coverage()] to see beforehand at which level each taxon
#'   will be estimated; [summarize_biomass()] for totals with intervals;
#'   [cite_sources()] for the references to cite.
#' @export
#' @examples
#' est <- estimate_biomass(example_survey$abundance, example_survey$taxa,
#'                         example_survey$lengths, bias_correction = "lognormal")
#' est
#' head(as.data.frame(est)[, c("site", "taxon_id", "method", "match_level",
#'                             "ind_median", "ind_lower", "ind_upper", "unit")])
estimate_biomass <- function(abundance, taxa, lengths = NULL, currency = "wet_mass", unit = NULL,
                             habitat = NULL, bias_correction, level = 0.95, extrapolate = FALSE,
                             reference = get_reference()) {
  check_data_frame(abundance)
  check_data_frame(taxa)
  if (!is.null(lengths)) check_data_frame(lengths)
  check_conversion_args(currency, habitat, unit, to_arg = "currency")
  check_level(level)
  check_flag(extrapolate)
  check_reference(reference)
  tx <- check_taxa_table(taxa)
  ab <- check_abundance_table(abundance, tx)
  ln <- if (is.null(lengths)) NULL else check_lengths_table(lengths, tx, ab)
  warn_taxa_lengths(tx)
  lg <- bind_length_groups(group_lengths(taxa_length_rows(tx), "taxa"),
                           group_lengths(ln, "lengths"))
  rows <- resolve_rows(ab, tx, lg)
  if (any(!is.na(rows$length_mm))) {
    if (missing(bias_correction)) {
      abort_input(paste(
        "`bias_correction` is required when lengths are used: \"none\" (median mass of the",
        "individuals of a given length) or \"lognormal\" (their mean, needed for biomass)."
      ))
    }
    check_choice(bias_correction, bias_correction_levels)
  } else if (missing(bias_correction)) {
    bias_correction <- NA_character_
  } else {
    check_choice(bias_correction, bias_correction_levels)
  }
  pred <- predict_rows(rows, tx, reference, bias_correction, extrapolate, ind = lg$ind)
  out <- assemble_estimates(ab, tx, rows, pred, level, reference)
  est <- new_bmassa_estimate(out, rows, pred, reference,
                             lineage = tx[rows$ti, c(main_ranks, "col_id")],
                             settings = list(level = level, bias_correction = bias_correction,
                                             extrapolate = extrapolate, reference = reference$info),
                             lengths_ind = lg$ind)
  # All the rows in one currency and unit (version 0.5.0).
  est <- convert_estimate(est, currency, habitat, unit)
  attr(est, "n_failed") <- NULL
  n_not <- sum(est$estimates$status == "not_estimated")
  if (n_not > 0L) {
    warn_not_estimable(sprintf(
      "%d of %d rows were not estimated; see the `reason` column of the result.",
      n_not, nrow(est$estimates)
    ))
  }
  est
}

# Input tables ----------------------------------------------------------------

check_level <- function(x, arg = arg_name(substitute(x)), call = sys.call(-1)) {
  force(arg)
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    abort_input(sprintf("`%s` must be a single number between 0 and 1, not %s.", arg,
                        describe_value(x)), call = call)
  }
  invisible(x)
}

# Rows of `data` that repeat a key.
check_unique_rows <- function(data, cols, arg, call = sys.call(-1)) {
  key <- make_key(data, cols)
  dup <- which(duplicated(key) | duplicated(key, fromLast = TRUE))
  if (length(dup) > 0L) {
    abort_input(sprintf(
      "`%s` has more than one row for the same %s (%s). Combine them or tell them apart.",
      arg, format_items(quote_items(cols), max = Inf), format_rows(dup)
    ), rows = dup, call = call)
  }
  invisible(data)
}

check_text_column <- function(data, column, arg, allow_na = FALSE, call = sys.call(-1)) {
  x <- data[[column]]
  if (!is.character(x) && !is.factor(x) && !all(is.na(x))) {
    abort_input(sprintf("Column '%s' of `%s` must be text, not %s.", column, arg, class(x)[1L]),
                call = call)
  }
  if (!allow_na && anyNA(x)) {
    abort_input(sprintf("Column '%s' of `%s` has missing values in %s.", column, arg,
                        format_rows(which(is.na(x)))), rows = which(is.na(x)), call = call)
  }
  invisible(data)
}

check_taxa_table <- function(taxa, require_count_unit = TRUE, call = sys.call(-1)) {
  check_columns(taxa, c("taxon_id", "broad_group", main_ranks), arg = "taxa", call = call)
  tx <- as_plain_data_frame(taxa)
  check_text_column(tx, "taxon_id", "taxa", call = call)
  check_text_column(tx, "broad_group", "taxa", call = call)
  for (r in main_ranks) {
    check_text_column(tx, r, "taxa", allow_na = TRUE, call = call)
    v <- trimws(gsub("\\s+", " ", as_text(tx[[r]])))
    v[!is.na(v) & !nzchar(v)] <- NA_character_
    tx[[r]] <- v
  }
  tx$taxon_id <- as_text(tx$taxon_id)
  tx$broad_group <- as_text(tx$broad_group)
  if ("col_id" %in% names(tx)) {
    check_text_column(tx, "col_id", "taxa", allow_na = TRUE, call = call)
    v <- trimws(as_text(tx$col_id))
    v[!is.na(v) & !nzchar(v)] <- NA_character_
    tx$col_id <- v
  } else {
    tx$col_id <- rep(NA_character_, nrow(tx))
  }
  check_unique_rows(tx, "taxon_id", "taxa", call = call)
  check_values_in(tx, "broad_group", broad_group_levels, allow_na = FALSE, arg = "taxa", call = call)
  empty <- which(rowSums(!is.na(as.matrix(tx[main_ranks]))) == 0L)
  if (length(empty) > 0L) {
    abort_input(sprintf("`taxa` has no classification (all of %s are NA) in %s.",
                        format_items(main_ranks, max = Inf), format_rows(empty)),
                rows = empty, call = call)
  }
  # The taxon unit: life stage, sex and counting unit.
  if ("life_stage" %in% names(tx)) {
    old <- which(as_text(tx$life_stage) %in% c("female", "female_ovigerous", "male"))
    if (length(old) > 0L) {
      abort_input(sprintf(paste(
        "Column 'life_stage' of `taxa` has a sex in %s: give it in `sex` (\"female\" or",
        "\"male\") and the stage, if any, in `life_stage` (\"adult\", \"juvenile\"...)."),
        format_rows(old)), rows = old, call = call)
    }
    check_values_in(tx, "life_stage", unit_stage_levels, arg = "taxa", call = call)
    tx$life_stage <- as_text(tx$life_stage)
  } else {
    tx$life_stage <- rep(NA_character_, nrow(tx))
  }
  if ("sex" %in% names(tx)) {
    check_values_in(tx, "sex", sex_levels, arg = "taxa", call = call)
    tx$sex <- as_text(tx$sex)
  } else {
    tx$sex <- rep(NA_character_, nrow(tx))
  }
  if ("count_unit" %in% names(tx)) {
    check_values_in(tx, "count_unit", count_unit_levels, arg = "taxa", call = call)
    tx$count_unit <- as_text(tx$count_unit)
  } else {
    tx$count_unit <- rep(NA_character_, nrow(tx))
  }
  phyto <- tx$broad_group == "phytoplankton"
  bad <- which(phyto & !tx$count_unit %in% c("cell", "natural_unit") &
                 (require_count_unit | !is.na(tx$count_unit)))
  if (length(bad) > 0L) {
    abort_input(sprintf(paste(
      "Phytoplankton taxa need `count_unit` \"cell\" or \"natural_unit\" in `taxa` (%s):",
      "the reference has biovolumes per cell and per natural unit."), format_rows(bad)),
      rows = bad, call = call)
  }
  bad <- which(!phyto & !is.na(tx$count_unit) & tx$count_unit != "individual")
  if (length(bad) > 0L) {
    abort_input(sprintf("Animal taxa must have `count_unit` \"individual\" or NA in `taxa` (%s).",
                        format_rows(bad)), rows = bad, call = call)
  }
  tx$count_unit[!phyto] <- "individual"
  # Cells per natural unit (colonies, filaments), when they were counted.
  if ("cells_per_unit" %in% names(tx)) {
    check_numeric_column(tx, "cells_per_unit", lower = 1, allow_na = TRUE, arg = "taxa",
                         call = call)
    bad <- which(!is.na(tx$cells_per_unit) & !(phyto & tx$count_unit %in% "natural_unit"))
    if (length(bad) > 0L) {
      abort_input(sprintf(paste(
        "Column 'cells_per_unit' of `taxa` is for phytoplankton counted in natural units",
        "(`count_unit` \"natural_unit\"); give NA in %s."), format_rows(bad)),
        rows = bad, call = call)
    }
    tx$cells_per_unit <- as.numeric(tx$cells_per_unit)
  } else {
    tx$cells_per_unit <- rep(NA_real_, nrow(tx))
  }
  if (any(length_input_columns %in% names(tx))) {
    check_length_columns(tx, "taxa", required = FALSE, call = call)
  }
  tx
}

# The mean lengths of a table (`taxa`, optional; `lengths`, required): a row
# with a length needs its unit and type; the spread is a SD or a CV, not both.
check_length_columns <- function(d, arg, required, call = sys.call(-1)) {
  force(call)
  if (!required && !"length_mean" %in% names(d)) {
    stray <- intersect(length_input_columns, names(d))
    abort_input(sprintf("`%s` has %s but no `length_mean`.", arg,
                        format_items(quote_items(stray), max = Inf)), call = call)
  }
  check_columns(d, "length_mean", arg = arg, call = call)
  if (!required && all(is.na(d$length_mean))) {
    return(invisible(d))
  }
  check_numeric_column(d, "length_mean", lower = 0, include_lower = FALSE,
                       allow_na = !required, arg = arg, call = call)
  check_columns(d, c("length_unit", "length_type"), arg = arg, call = call)
  has <- !is.na(d$length_mean)
  for (col in c("length_unit", "length_type")) {
    rows <- which(has & is.na(d[[col]]))
    if (length(rows) > 0L) {
      abort_input(sprintf("Column '%s' of `%s` is missing in %s, which have a length.", col, arg,
                          format_rows(rows)), rows = rows, call = call)
    }
  }
  check_values_in(d, "length_unit", names(length_unit_mm), arg = arg, call = call)
  check_values_in(d, "length_type", length_type_levels, arg = arg, call = call)
  for (col in c("length_sd", "length_cv")) {
    if (col %in% names(d)) {
      check_numeric_column(d, col, lower = 0, allow_na = TRUE, arg = arg, call = call)
    }
  }
  if (all(c("length_sd", "length_cv") %in% names(d))) {
    both <- which(!is.na(d$length_sd) & !is.na(d$length_cv))
    if (length(both) > 0L) {
      abort_input(sprintf("`%s` gives both `length_sd` and `length_cv` in %s; keep one.", arg,
                          format_rows(both)), rows = both, call = call)
    }
  }
  invisible(d)
}

check_abundance_table <- function(abundance, tx, call = sys.call(-1)) {
  check_columns(abundance, c("site", "taxon_id", "abundance", "per"), arg = "abundance",
                call = call)
  moved <- intersect(c("life_stage", "sex", "count_unit"), names(abundance))
  if (length(moved) > 0L) {
    abort_input(sprintf(paste(
      "`abundance` has %s: these belong to the taxon unit and go in `taxa`, with one",
      "`taxon_id` per taxon, life stage and sex (and counting unit in phytoplankton)."),
      format_items(quote_items(moved), max = Inf)), call = call)
  }
  check_no_clash(abundance, c(estimate_biomass_columns, conversion_columns), arg = "abundance",
                 call = call)
  ab <- as_plain_data_frame(abundance)
  check_text_column(ab, "site", "abundance", call = call)
  check_text_column(ab, "taxon_id", "abundance", call = call)
  check_text_column(ab, "per", "abundance", call = call)
  check_numeric_column(ab, "abundance", lower = 0, allow_na = TRUE, arg = "abundance",
                       call = call)
  if ("date" %in% names(ab) && (is.list(ab$date) || !is.atomic(ab$date))) {
    abort_input("Column 'date' of `abundance` must be a vector of dates, days or text.",
                call = call)
  }
  missing_ids <- which(!as_text(ab$taxon_id) %in% tx$taxon_id)
  if (length(missing_ids) > 0L) {
    abort_input(sprintf(
      "Some `taxon_id` of `abundance` are not in `taxa` (%s in %s).",
      format_items(quote_items(unique(as_text(ab$taxon_id)[missing_ids]))), format_rows(missing_ids)
    ), rows = missing_ids, call = call)
  }
  check_unique_rows(ab, intersect(c("site", "date", "taxon_id"), names(ab)), "abundance",
                    call = call)
  ab
}

# The table of lengths (version 0.5.0): each row is either the mean length of
# a taxon unit (`length_mean`, with its spread if known) or the length of
# measured individuals (`length`, with `n` individuals of that length, 1 by
# default: one row per individual, or per length class with its count).
# `site` and `date` say where the lengths apply: a row without `site` applies
# to all the samples of the unit, one with `site` and without `date` to all
# the dates of the site, and one with both to that sample.
check_lengths_table <- function(lengths, tx, ab = NULL, call = sys.call(-1)) {
  force(call)
  check_columns(lengths, c("taxon_id", "length_unit", "length_type"), arg = "lengths", call = call)
  if (!any(c("length_mean", "length") %in% names(lengths))) {
    abort_input(paste(
      "`lengths` needs `length_mean` (the mean length of a taxon unit) or `length` (the",
      "lengths of measured individuals), or both."), call = call)
  }
  if ("life_stage" %in% names(lengths)) {
    abort_input(paste(
      "`lengths` has `life_stage`: the stage belongs to the taxon unit, in `taxa`;",
      "give each stage its own `taxon_id` and its lengths under it."), call = call)
  }
  ln <- as_plain_data_frame(lengths)
  nl <- nrow(ln)
  # Columns left empty (all NA, read as logical) are numeric.
  for (col in intersect(c("length_mean", "length", "n", "length_sd", "length_cv"), names(ln))) {
    if (is.logical(ln[[col]]) && all(is.na(ln[[col]]))) ln[[col]] <- as.numeric(ln[[col]])
  }
  check_text_column(ln, "taxon_id", "lengths", call = call)
  ln$taxon_id <- as_text(ln$taxon_id)
  if ("site" %in% names(ln)) {
    check_text_column(ln, "site", "lengths", allow_na = TRUE, call = call)
    ln$site <- as_text(ln$site)
  } else {
    ln$site <- rep(NA_character_, nl)
  }
  if ("date" %in% names(ln)) {
    if (is.list(ln$date) || !is.atomic(ln$date)) {
      abort_input("Column 'date' of `lengths` must be a vector of dates, days or text.", call = call)
    }
    if (!is.null(ab) && !"date" %in% names(ab)) {
      abort_input("`lengths` has `date` but `abundance` does not: the samples have no dates.",
                  call = call)
    }
    nosite <- which(!is.na(ln$date) & is.na(ln$site))
    if (length(nosite) > 0L) {
      abort_input(sprintf(paste(
        "`lengths` has a `date` without `site` in %s: a date belongs to the sample of a site",
        "(leave both empty for lengths that apply to all the samples)."), format_rows(nosite)),
        rows = nosite, call = call)
    }
  }
  for (col in c("length_mean", "length")) {
    if (col %in% names(ln)) {
      check_numeric_column(ln, col, lower = 0, include_lower = FALSE, allow_na = TRUE,
                           arg = "lengths", call = call)
      ln[[col]] <- as.numeric(ln[[col]])
    } else {
      ln[[col]] <- rep(NA_real_, nl)
    }
  }
  both <- which(!is.na(ln$length_mean) & !is.na(ln$length))
  if (length(both) > 0L) {
    abort_input(sprintf(paste(
      "`lengths` gives both `length_mean` and `length` in %s: a row is either the mean length",
      "of a taxon unit or the length of measured individuals."), format_rows(both)),
      rows = both, call = call)
  }
  none <- which(is.na(ln$length_mean) & is.na(ln$length))
  if (length(none) > 0L) {
    abort_input(sprintf("`lengths` has no length (`length_mean` or `length`) in %s; drop them.",
                        format_rows(none)), rows = none, call = call)
  }
  for (col in c("length_unit", "length_type")) {
    rows <- which(is.na(ln[[col]]))
    if (length(rows) > 0L) {
      abort_input(sprintf("Column '%s' of `lengths` is missing in %s.", col, format_rows(rows)),
                  rows = rows, call = call)
    }
  }
  check_values_in(ln, "length_unit", names(length_unit_mm), arg = "lengths", call = call)
  check_values_in(ln, "length_type", length_type_levels, arg = "lengths", call = call)
  ind <- !is.na(ln$length)
  for (col in c("length_sd", "length_cv")) {
    if (col %in% names(ln)) {
      check_numeric_column(ln, col, lower = 0, allow_na = TRUE, arg = "lengths", call = call)
      bad <- which(ind & !is.na(ln[[col]]))
      if (length(bad) > 0L) {
        abort_input(sprintf(paste(
          "`lengths` gives `%s` for lengths of individuals in %s: the spread is for a mean",
          "length (`length_mean`); the individuals give their own."), col, format_rows(bad)),
          rows = bad, call = call)
      }
    }
  }
  if (all(c("length_sd", "length_cv") %in% names(ln))) {
    both <- which(!is.na(ln$length_sd) & !is.na(ln$length_cv))
    if (length(both) > 0L) {
      abort_input(sprintf("`lengths` gives both `length_sd` and `length_cv` in %s; keep one.",
                          format_rows(both)), rows = both, call = call)
    }
  }
  if ("n" %in% names(ln)) {
    check_numeric_column(ln, "n", lower = 0, include_lower = FALSE, allow_na = TRUE,
                         arg = "lengths", call = call)
    bad <- which(!ind & !is.na(ln$n))
    if (length(bad) > 0L) {
      abort_input(sprintf(paste(
        "`lengths` gives `n` for a mean length in %s: `n` is the number of individuals of a",
        "row of `length` (a length class)."), format_rows(bad)), rows = bad, call = call)
    }
    ln$n <- as.numeric(ln$n)
    ln$n[ind & is.na(ln$n)] <- 1
  } else {
    ln$n <- ifelse(ind, 1, NA_real_)
  }
  missing_ids <- which(!ln$taxon_id %in% tx$taxon_id)
  if (length(missing_ids) > 0L) {
    abort_input(sprintf(
      "Some `taxon_id` of `lengths` are not in `taxa` (%s in %s).",
      format_items(quote_items(unique(ln$taxon_id[missing_ids]))), format_rows(missing_ids)
    ), rows = missing_ids, call = call)
  }
  # One kind of length and one length type per unit and sample, and one mean.
  date <- if ("date" %in% names(ln)) as_text(ln$date) else rep(NA_character_, nl)
  key <- paste(ln$taxon_id, ln$site, date, sep = "\u001f")
  mixed <- which(key %in% key[ind] & key %in% key[!ind])
  if (length(mixed) > 0L) {
    abort_input(sprintf(paste(
      "`lengths` has both a mean length and lengths of individuals for the same taxon unit",
      "and sample in %s: give one or the other."), format_rows(mixed)), rows = mixed, call = call)
  }
  n_type <- tapply(as_text(ln$length_type), key, function(x) length(unique(x)))
  multi <- which(n_type[key] > 1L)
  if (length(multi) > 0L) {
    abort_input(sprintf(paste(
      "`lengths` has more than one length type for the same taxon unit and sample in %s:",
      "use one type of length per unit and sample."), format_rows(multi)), rows = multi,
      call = call)
  }
  dup <- which(!ind & (duplicated(key) | duplicated(key, fromLast = TRUE)))
  if (length(dup) > 0L) {
    abort_input(sprintf(paste(
      "`lengths` has more than one mean length for the same taxon unit and sample in %s",
      "(`taxon_id`, `site` and `date`). Combine them, or give the lengths of the",
      "individuals in `length`."), format_rows(dup)), rows = dup, call = call)
  }
  ln
}

# Mean lengths of `taxa` (deprecated in version 0.5.0) as rows of `lengths`
# without site: they apply to all the samples of the unit, after any row of
# `lengths`.
taxa_length_rows <- function(tx) {
  has <- if ("length_mean" %in% names(tx)) which(!is.na(tx$length_mean)) else integer(0)
  if (length(has) == 0L) {
    return(NULL)
  }
  out <- data.frame(taxon_id = tx$taxon_id[has], site = NA_character_,
                    length_mean = as.numeric(tx$length_mean[has]), length = NA_real_, n = NA_real_,
                    length_unit = as_text(tx$length_unit)[has],
                    length_type = as_text(tx$length_type)[has], stringsAsFactors = FALSE)
  for (col in c("length_sd", "length_cv")) {
    if (col %in% names(tx)) out[[col]] <- as.numeric(tx[[col]][has])
  }
  out
}

warn_taxa_lengths <- function(tx, call = sys.call(-1)) {
  if ("length_mean" %in% names(tx) && any(!is.na(tx$length_mean))) {
    warn_deprecated(paste(
      "Lengths in `taxa` are deprecated and will not be read in a future version: give them",
      "in `lengths`, where a row without `site` applies to all the samples of the taxon unit",
      "(rows of `lengths` take precedence)."), call = call)
  }
}

# The lengths of each taxon unit and sample: one group per taxon unit, site
# and date of `lengths` (or of `taxa`), with its kind (`"mean"` or
# `"individuals"`), the length type, the mean length (mm) and its coefficient
# of variation, or the number of individuals and the range of their lengths;
# and the individual lengths (mm) of each group.
group_lengths <- function(ln, origin) {
  groups <- data.frame(taxon_id = character(0), site = character(0), date = character(0),
                       source = character(0), basis = character(0), length_type = character(0),
                       length_mm = numeric(0), length_cv = numeric(0), length_n = numeric(0),
                       length_lo = numeric(0), length_hi = numeric(0), stringsAsFactors = FALSE)
  ind <- data.frame(group = integer(0), length_mm = numeric(0), n = numeric(0))
  if (is.null(ln) || nrow(ln) == 0L) {
    return(list(groups = groups, ind = ind))
  }
  date <- if ("date" %in% names(ln)) as_text(ln$date) else rep(NA_character_, nrow(ln))
  key <- paste(ln$taxon_id, ln$site, date, sep = "\u001f")
  gid <- match(key, unique(key))
  first <- !duplicated(gid)
  to_mm <- unname(length_unit_mm[as_text(ln$length_unit)])
  is_ind <- !is.na(ln$length)
  cv <- rep(NA_real_, nrow(ln))
  if ("length_cv" %in% names(ln)) cv <- ln$length_cv
  if ("length_sd" %in% names(ln)) cv <- ifelse(is.na(cv), ln$length_sd / ln$length_mean, cv)
  groups <- data.frame(taxon_id = ln$taxon_id[first], site = ln$site[first], date = date[first],
                       source = origin, basis = ifelse(is_ind[first], "individuals", "mean"),
                       length_type = as_text(ln$length_type)[first],
                       length_mm = ln$length_mean[first] * to_mm[first], length_cv = cv[first],
                       length_n = NA_real_, length_lo = NA_real_, length_hi = NA_real_,
                       stringsAsFactors = FALSE)
  groups$length_lo <- groups$length_hi <- groups$length_mm
  ind <- data.frame(group = gid[is_ind], length_mm = ln$length[is_ind] * to_mm[is_ind],
                    n = ln$n[is_ind])
  if (nrow(ind) > 0L) {
    g <- sort(unique(ind$group))
    w <- tapply(ind$n, ind$group, sum)[as.character(g)]
    groups$length_n[g] <- w
    groups$length_lo[g] <- tapply(ind$length_mm, ind$group, min)[as.character(g)]
    groups$length_hi[g] <- tapply(ind$length_mm, ind$group, max)[as.character(g)]
    groups$length_mm[g] <- tapply(ind$length_mm * ind$n, ind$group, sum)[as.character(g)] / w
  }
  list(groups = groups, ind = ind)
}

# Both sets of groups in one: those of `taxa` first.
bind_length_groups <- function(a, b) {
  off <- nrow(a$groups)
  bi <- b$ind
  bi$group <- bi$group + off
  list(groups = rbind(a$groups, b$groups), ind = rbind(a$ind, bi))
}

# Code of the stage of a unit in the reference data: the stage, or the sex of
# an adult or unstaged unit ("female", "male").
unit_stage <- function(life_stage, sex) {
  out <- life_stage
  by_sex <- (is.na(life_stage) | life_stage == "adult") & !is.na(sex)
  out[by_sex] <- sex[by_sex]
  out
}


# Empty length columns of the rows given to predict_rows().
add_length_columns <- function(rows) {
  n <- nrow(rows)
  rows$length_mm <- rep(NA_real_, n)
  rows$length_type <- rep(NA_character_, n)
  rows$length_cv <- rep(NA_real_, n)
  rows$length_basis <- rep(NA_character_, n)
  rows$length_n <- rep(NA_real_, n)
  rows$length_lo <- rep(NA_real_, n)
  rows$length_hi <- rep(NA_real_, n)
  rows$length_group <- rep(NA_integer_, n)
  rows
}

# Fill the length columns of rows `i` from groups `g` of group_lengths().
set_row_lengths <- function(rows, i, g, groups) {
  for (col in c("length_mm", "length_type", "length_cv", "length_basis", "length_n",
                "length_lo", "length_hi")) {
    src <- if (col == "length_basis") "basis" else col
    rows[[col]][i] <- groups[[src]][g]
  }
  rows$length_group[i] <- g
  rows
}

# One row per abundance row with what the prediction needs: taxon, domain,
# counting unit, stage, and the lengths that apply, if any: the most specific
# group of `lengths` (the sample, the site in all its dates, or all the
# samples), otherwise the mean length of `taxa` (deprecated).
resolve_rows <- function(ab, tx, lg) {
  n <- nrow(ab)
  ti <- match(as_text(ab$taxon_id), tx$taxon_id)
  domain <- unname(broad_group_domain[tx$broad_group[ti]])
  rows <- data.frame(row = seq_len(n), taxon_id = as_text(ab$taxon_id), ti = ti,
                     site = as_text(ab$site), domain = domain,
                     count_unit = tx$count_unit[ti], cells_per_unit = tx$cells_per_unit[ti],
                     stage = unit_stage(tx$life_stage, tx$sex)[ti], stringsAsFactors = FALSE)
  rows <- add_length_columns(rows)
  g <- lg$groups
  if (n == 0L || nrow(g) == 0L) {
    return(rows)
  }
  a_date <- if ("date" %in% names(ab)) as_text(ab$date) else rep(NA_character_, n)
  sep <- "\u001f"
  g_key <- paste(g$taxon_id, g$site, g$date, g$source, sep = sep)
  choose <- rep(NA_integer_, n)
  # From the least to the most specific: `taxa`, all the samples, the site in
  # all its dates, the sample.
  keys <- list(paste(rows$taxon_id, NA, NA, "taxa", sep = sep),
               paste(rows$taxon_id, NA, NA, "lengths", sep = sep),
               paste(rows$taxon_id, rows$site, NA, "lengths", sep = sep),
               paste(rows$taxon_id, rows$site, a_date, "lengths", sep = sep))
  for (k in keys) {
    m <- match(k, g_key)
    choose[!is.na(m)] <- m[!is.na(m)]
  }
  has <- which(!is.na(choose))
  if (length(has) > 0L) rows <- set_row_lengths(rows, has, choose[has], g)
  rows
}

# Length of the mean mass of a set of individual lengths (mm, with weights n)
# for a relation of slope b: the mean of L^b, back to a length.
mean_mass_length <- function(mm, n, b) {
  if (abs(b) < 1e-8) {
    return(exp(sum(n * log(mm)) / sum(n)))
  }
  (sum(n * mm^b) / sum(n))^(1 / b)
}

# Prediction ----------------------------------------------------------------------

# Prediction for each row: model, node, linear predictor and variance parts.
# `ind` holds the individual lengths of the length groups (group_lengths()).
# Rows of the same taxon unit with the same lengths (or none) have the same
# prediction, whatever their site or date: it is computed once for each and
# copied (a survey repeats each taxon in many samples).
predict_rows <- function(rows, tx, reference, bias_correction, extrapolate, ind = NULL) {
  if (is.null(rows$length_basis)) rows <- add_length_columns(rows)
  cols <- intersect(c("ti", "domain", "count_unit", "cells_per_unit", "stage", "length_mm",
                      "length_type", "length_cv", "length_basis", "length_n", "length_lo",
                      "length_hi", "length_group"), names(rows))
  key <- do.call(paste, c(unname(as.list(rows[cols])), sep = "\u001f"))
  first <- which(!duplicated(key))
  if (length(first) == nrow(rows)) {
    return(predict_unique_rows(rows, tx, reference, bias_correction, extrapolate, ind))
  }
  res <- predict_unique_rows(rows[first, , drop = FALSE], tx, reference, bias_correction,
                             extrapolate, ind)
  idx <- match(key, key[first])
  pred <- res$pred[idx, , drop = FALSE]
  rownames(pred) <- NULL
  list(pred = pred, missing_keys = res$missing_keys[idx])
}

predict_unique_rows <- function(rows, tx, reference, bias_correction, extrapolate, ind = NULL) {
  n <- nrow(rows)
  nodes <- reference$nodes
  comps <- reference$components
  models <- reference$models
  chr <- rep(NA_character_, n)
  num <- rep(NA_real_, n)
  int <- rep(NA_integer_, n)
  pr <- data.frame(model = chr, method = chr, stage_used = chr, length_status = chr, node = int,
                   stage_node = int, contrast = chr, match_level = chr, match_name = chr, n_reference = int,
                   sources = chr, length_used = num, xc = rep(0, n), offset = rep(0, n), eta = num,
                   v_param = num, v_missing = num,
                   v_population = num, v_source = num, corrections = chr, uncertainty = chr,
                   status = rep("not_estimated", n), reason = chr, note = chr,
                   stringsAsFactors = FALSE)
  if (is.null(rows$length_basis)) rows <- add_length_columns(rows)
  basis <- rows$length_basis
  basis[is.na(basis) & !is.na(rows$length_mm)] <- "mean"
  lo <- ifelse(is.na(rows$length_lo), rows$length_mm, rows$length_lo)
  hi <- ifelse(is.na(rows$length_hi), rows$length_mm, rows$length_hi)
  missing_keys <- vector("list", n)
  if (n == 0L) {
    return(list(pred = pr, missing_keys = missing_keys))
  }
  comp_value <- function(m, l, p) {
    v <- comps$value[comps$model == m & comps$level == l & comps$parameter == p]
    if (length(v)) v else 0
  }
  insect <- tx$class[rows$ti] %in% "Insecta"
  copepod <- tx$class[rows$ti] %in% "Copepoda" | tx$order[rows$ti] %in% copepod_orders
  cells <- if ("cells_per_unit" %in% names(rows)) rows$cells_per_unit else rep(NA_real_, n)
  colonial <- tx$genus[rows$ti] %in% reference$colonial_genera
  # Length type of the model: the body length of copepod nauplii has its own
  # relations (their shape differs from that of copepodites and adults).
  lt_m <- rows$length_type
  lt_m[copepod & rows$stage %in% "nauplius" & lt_m %in% "body_length"] <- "nauplius_body_length"
  # Model and stage of each row.
  for (i in seq_len(n)) {
    d <- rows$domain[i]
    st <- rows$stage[i]
    if (!is.na(rows$length_mm[i])) {
      pr$method[i] <- "length"
      pr$model[i] <- paste0(d, "_lwr")
      lt <- reference$lwr_types
      if (!any(lt$model == pr$model[i] & lt$length_type == lt_m[i])) {
        pr$reason[i] <- sprintf("no reference length-mass relations of %s for %s",
                                gsub("_", " ", lt_m[i]), d)
        pr$model[i] <- NA_character_
      } else if (d == "amphibian" && st %in% c("larva", "nymph", "juvenile")) {
        pr$reason[i] <- "the length-mass relations of amphibians are for animals after metamorphosis"
        pr$model[i] <- NA_character_
      }
      next
    }
    pr$method[i] <- "typical"
    if (d == "phytoplankton") {
      if (rows$count_unit[i] == "natural_unit" && !is.na(cells[i])) {
        # A natural unit of known cells: the biovolume of a cell times the cells.
        pr$model[i] <- "phytoplankton_cell"
        pr$offset[i] <- log(cells[i])
        pr$note[i] <- append_note(pr$note[i], sprintf(
          "natural unit of %s cells: the biovolume of one cell times the number of cells",
          format(signif(cells[i], 4))))
      } else {
        pr$model[i] <- paste0("phytoplankton_", rows$count_unit[i])
        if (rows$count_unit[i] == "natural_unit" && colonial[i]) {
          pr$note[i] <- append_note(pr$note[i], paste(
            "colonial genus counted in natural units: the reference has the natural units of",
            "its sources (often single cells or small groups), not whole colonies; give",
            "cells_per_unit for colonies counted whole"))
        }
      }
    } else if (d == "invertebrate") {
      pr$model[i] <- "invertebrate_individual"
      if (insect[i]) {
        pr$stage_used[i] <- if (is.na(st)) "larva (assumed)" else insect_stage_group[[st]]
      } else if (copepod[i]) {
        if (is.na(st)) {
          pr$stage_used[i] <- "copepod_adult (assumed)"
        } else if (st %in% names(copepod_stage_group)) {
          pr$stage_used[i] <- copepod_stage_group[[st]]
        } else {
          pr$model[i] <- NA_character_
          pr$reason[i] <- sprintf(paste(
            "life stage \"%s\" of a copepod: give \"adult\", \"copepodite\" or \"nauplius\"",
            "(their masses differ about tenfold)"), st)
        }
      } else {
        pr$stage_used[i] <- "non_insect"
        if (!is.na(st)) pr$note[i] <- "life stage not modelled for animals other than insects"
      }
    } else if (d == "fish") {
      if (is.na(st)) {
        pr$model[i] <- "fish_individual"
        pr$note[i] <- append_note(pr$note[i], paste(
          "fish without length: mean mass of a population of all the sizes caught in the",
          "reference (electrofishing of streams, nets in lakes, a river fishery), which may",
          "differ from the sizes of this population; a mean length gives a closer estimate"))
      } else {
        pr$reason[i] <- paste(
          "fish without length are estimated as whole populations (all the sizes caught in",
          "the reference), not by life stage or sex: give a length, or leave life_stage and",
          "sex empty")
      }
    } else if (d == "amphibian") {
      if (st %in% c("adult", "female", "female_ovigerous", "male")) {
        pr$model[i] <- "amphibian_adult"
        pr$stage_used[i] <- "adult"
      } else {
        pr$reason[i] <- if (is.na(st)) {
          paste("amphibians without length and without a stage are taken as larvae, which have",
                "no reference masses; give life_stage = \"adult\" for adults, or a length")
        } else {
          "amphibian larvae and juveniles have no reference masses without lengths"
        }
      }
    }
  }
  # Matching, once per taxon, model and length type.
  key <- paste(rows$ti, pr$model, lt_m, sep = "\u001f")
  for (u in unique(key[!is.na(pr$model)])) {
    idx <- which(key == u & !is.na(pr$model))
    i <- idx[1]
    m <- pr$model[i]
    lv <- strsplit(models$levels[models$model == m], ",")[[1]]
    cand <- which(nodes$model == m & nodes$level %in% lv &
                    (is.na(nodes$length_type) | nodes$length_type %in% lt_m[i]))
    lin <- unlist(tx[rows$ti[i], main_ranks])
    mt <- match_lineage(lin, nodes[cand, ], lv, col_id = tx$col_id[rows$ti[i]])
    root <- if (is.null(mt$row) && m == "fish_individual") {
      which(nodes$model == m & nodes$level == "root")
    } else {
      integer(0)
    }
    if (is.null(mt$row) && length(root) != 1L) {
      pr$reason[idx] <- sprintf("no taxon of the reference model %s matches the classification%s",
                                m, if (is.na(mt$note)) "" else paste0(" (", mt$note, ")"))
      next
    }
    if (is.null(mt$row)) {
      # Fish without length of a class (or order, the highest level of the
      # model) absent from the reference: the mean of all its fish, with the
      # variances of all the levels.
      j <- root
      mt$note <- append_note(mt$note, paste(lv[1], "not in the reference: estimated from all",
                                            "the fish of the reference"))
    } else {
      j <- cand[mt$row]
    }
    for (i in idx) {
      jj <- j
      if (pr$method[i] == "length") {
        cov <- find_covering_node(j, lo[i], hi[i], lt_m[i], nodes, comps)
        if (is.na(cov) && !extrapolate) {
          rt <- nodes[nodes$model == m & nodes$level == "root" & nodes$node == lt_m[i], ]
          what <- if (basis[i] == "individuals") {
            sprintf("lengths of individuals %s-%s mm", format(signif(lo[i], 3)), format(signif(hi[i], 3)))
          } else {
            sprintf("length %s mm", format(signif(rows$length_mm[i], 3)))
          }
          pr$reason[i] <- sprintf(paste(
            "%s outside the lengths of the reference relations of %s and its",
            "relatives (%s-%s mm for %s); use extrapolate = TRUE to estimate it anyway"),
            what, nodes$node[j], format(signif(rt$length_min, 3)),
            format(signif(rt$length_max, 3)), gsub("_", " ", lt_m[i]))
          next
        }
        if (is.na(cov)) {
          pr$length_status[i] <- "outside"
        } else {
          pr$length_status[i] <- "within"
          if (cov != j) {
            pr$note[i] <- append_note(pr$note[i], sprintf(
              "length outside the relations of %s: estimated from %s", nodes$node[j], nodes$node[cov]))
          }
          jj <- cov
        }
      }
      pr$node[i] <- jj
      pr$note[i] <- append_note(pr$note[i], mt$note)
      pr$note[i] <- append_note(pr$note[i], lineage_note(lin, nodes$path[jj], nodes$level[jj]))
    }
  }
  # Linear predictor and variances.
  for (i in which(!is.na(pr$node))) {
    j <- pr$node[i]
    m <- pr$model[i]
    lv <- strsplit(models$levels[models$model == m], ",")[[1]]
    lwr <- pr$method[i] == "length"
    lt_row <- if (lwr) reference$lwr_types[reference$lwr_types$model == m &
                                              reference$lwr_types$length_type == lt_m[i], ]
    if (lwr) {
      # A set of individual lengths enters as the length of their mean mass
      # under the relation of the node.
      L <- rows$length_mm[i]
      if (basis[i] == "individuals") {
        v <- ind[ind$group == rows$length_group[i], , drop = FALSE]
        L <- mean_mass_length(v$length_mm, v$n, nodes$path_slope[j])
      }
      pr$length_used[i] <- L
    }
    xc <- if (lwr) log(pr$length_used[i]) - lt_row$centre else 0
    eta <- nodes$path_intercept[j] + pr$offset[i]
    vp <- nodes$path_var[j]
    if (lwr) {
      eta <- eta + nodes$path_slope[j] * xc
      vp <- vp + 2 * xc * nodes$path_cov[j] + xc^2 * nodes$path_var_slope[j]
    }
    miss <- if (nodes$level[j] == "root") lv else lv[seq_along(lv) > match(nodes$level[j], lv)]
    vm <- 0
    for (l in miss) {
      vm <- vm + comp_value(m, l, "intercept")
      if (lwr) vm <- vm + 2 * xc * comp_value(m, l, "intercept_slope") + xc^2 * comp_value(m, l, "slope")
    }
    lin <- unlist(tx[rows$ti[i], main_ranks])
    missing_keys[[i]] <- vapply(miss, function(l) {
      nm <- lin[[l]]
      paste(m, if (lwr) lt_m[i] else "", l,
            if (is.na(nm)) paste0("<", rows$taxon_id[i], ">") else paste(lin[main_ranks[seq_len(match(l, main_ranks))]], collapse = "|"),
            sep = "\u001f")
    }, character(1))
    if (!is.na(pr$stage_used[i])) {
      sg <- sub(" (assumed)", "", pr$stage_used[i], fixed = TRUE)
      # Copepods are in the group of the other invertebrates; their
      # copepodites and nauplii add a fixed contrast, with its variance and
      # twice its covariance with the path of the node.
      ck <- if (startsWith(sg, "copepod_")) unname(copepod_contrast[sg]) else NA_character_
      if (startsWith(sg, "copepod_")) sg <- "non_insect"
      ks <- which(nodes$model == m & nodes$level == "stage" & nodes$node == sg)
      if (length(ks) == 1L) {
        pr$stage_node[i] <- ks
        eta <- eta + nodes$effect[ks]
        vp <- vp + nodes$effect_var[ks]
      } else if (any(nodes$model == m & nodes$level == "stage")) {
        abort_internal(sprintf("stage group %s is not in the model %s", sg, m))
      }
      if (!is.na(ck)) {
        fx <- reference$fixed[reference$fixed$model == m & reference$fixed$term == paste0("c_", ck), ]
        cc <- paste0("contrast_cov_", ck)
        if (nrow(fx) != 1L || !cc %in% names(nodes)) {
          abort_internal(sprintf("contrast %s is not in the model %s", ck, m))
        }
        pr$contrast[i] <- ck
        eta <- eta + fx$estimate
        vp <- vp + fx$se^2 + 2 * nodes[[cc]][j]
      }
    }
    corr <- character(0)
    if (lwr) {
      b <- nodes$path_slope[j]
      if (basis[i] == "individuals") {
        pr$note[i] <- append_note(pr$note[i], sprintf(
          "mean mass of %s individuals measured (%s-%s mm)", format(signif(rows$length_n[i], 4)),
          format(signif(lo[i], 3)), format(signif(hi[i], 3))))
      } else if (!is.na(rows$length_cv[i])) {
        eta <- eta + b * (b - 1) / 2 * log1p(rows$length_cv[i]^2)
        corr <- c(corr, "jensen")
      } else {
        pr$note[i] <- append_note(pr$note[i], paste(
          "mean length without its spread (length_sd or length_cv): the mean mass is underestimated"))
      }
      if (identical(bias_correction, "lognormal")) {
        eta <- eta + lt_row$sigma_ind^2 / 2
        corr <- c(corr, "lognormal")
      }
    }
    pr$xc[i] <- xc
    pr$eta[i] <- eta
    pr$v_param[i] <- vp
    pr$v_missing[i] <- vm
    pr$v_population[i] <- comp_value(m, "residual", "intercept") + comp_value(m, "equation", "intercept") +
      comp_value(m, "population", "intercept")
    pr$v_source[i] <- comp_value(m, "source", "intercept") + comp_value(m, "study", "intercept")
    pr$corrections[i] <- if (length(corr)) paste(corr, collapse = "; ") else NA_character_
    has_src <- any(comps$model == m & comps$level %in% c("source", "study"))
    pr$uncertainty[i] <- paste(c("parameter", if (length(miss)) "taxonomic_imputation",
                                 "between_population",
                                 if (has_src) if (lwr) "between_study" else "between_source"),
                               collapse = "; ")
    pr$match_level[i] <- nodes$level[j]
    pr$match_name[i] <- nodes$node[j]
    pr$n_reference[i] <- nodes$n_units[j]
    pr$sources[i] <- nodes$sources[j]
    pr$status[i] <- "estimated"
    if (identical(pr$stage_used[i], "larva (assumed)")) {
      pr$note[i] <- append_note(pr$note[i], "life stage not given: taken as larva")
    } else if (identical(pr$stage_used[i], "copepod_adult (assumed)")) {
      pr$note[i] <- append_note(pr$note[i], paste(
        "life stage not given: taken as adult (too high if the count includes copepodites)"))
    }
  }
  list(pred = pr, missing_keys = missing_keys)
}

# Deepest node of a model whose name and rank match the lineage of a taxon.
# A taxon with the identifier of the COL (`col_id`) is looked for by it first:
# the identifier follows the taxon when its name or genus changes. Homonyms
# are resolved by the number of higher ranks that agree.
match_lineage <- function(lin, cand, levels, col_id = NA_character_) {
  if (!is.na(col_id) && "col_id" %in% names(cand)) {
    k <- which(cand$col_id %in% col_id & cand$level %in% levels)
    if (length(k) == 1L) {
      own <- lin[[cand$level[k]]]
      note <- if (identical(unname(own), cand$node[k])) NA_character_ else
        sprintf("matched by col_id to %s %s", cand$level[k], cand$node[k])
      return(list(row = k, note = note))
    }
  }
  notes <- character(0)
  for (l in rev(levels)) {
    name <- lin[[l]]
    if (is.na(name)) next
    k <- which(cand$level == l & cand$node == name)
    if (length(k) > 1L) {
      upper <- lin[main_ranks[seq_len(match(l, main_ranks) - 1L)]]
      score <- vapply(k, function(j) {
        p <- strsplit(cand$path[j], "|", fixed = TRUE)[[1]]
        sum(p[seq_along(upper)] == upper, na.rm = TRUE)
      }, numeric(1))
      best <- k[score == max(score)]
      if (length(best) == 1L) {
        return(list(row = best, note = sprintf("homonym %s told apart by the higher ranks", name)))
      }
      notes <- c(notes, sprintf("%s is ambiguous in the reference", name))
      next
    }
    if (length(k) == 1L) {
      return(list(row = k, note = if (length(notes)) paste(notes, collapse = "; ") else NA_character_))
    }
  }
  list(row = NULL, note = if (length(notes)) paste(notes, collapse = "; ") else NA_character_)
}

# The node itself, or its closest ancestor, whose reference relations cover a
# length; NA if none does, up to the root of the length type. At a level
# without random slopes (amphibian species, mostly one mean length each) the
# shape of the relation comes from the closest ancestor with slopes, so its
# range of lengths is the one that counts.
find_covering_node <- function(j, length_lo, length_hi, length_type, nodes, comps) {
  m <- nodes$model[j]
  slope_lv <- unique(comps$level[comps$model == m & comps$parameter == "slope"])
  parent_of <- function(k) {
    parent_path <- nodes$parent[k]
    nxt <- if (is.na(parent_path)) {
      which(nodes$model == m & nodes$level == "root" & nodes$node == length_type)
    } else {
      which(nodes$model == m & nodes$length_type %in% length_type & nodes$path == parent_path &
              nodes$level %in% main_ranks)
    }
    if (length(nxt) != 1L) {
      abort_internal(sprintf("parent of node %s of model %s not found", nodes$node[k], m))
    }
    nxt
  }
  cur <- j
  repeat {
    rng <- cur
    while (nodes$level[rng] != "root" && !nodes$level[rng] %in% slope_lv) rng <- parent_of(rng)
    if (!is.na(nodes$length_min[rng]) && length_lo >= nodes$length_min[rng] &&
        length_hi <= nodes$length_max[rng]) {
      return(cur)
    }
    if (nodes$level[cur] == "root") {
      return(NA_integer_)
    }
    cur <- parent_of(cur)
  }
}

# Higher ranks of the taxon that differ from those of the reference node.
lineage_note <- function(lin, path, level) {
  if (level == "root") return(NA_character_)
  p <- strsplit(path, "|", fixed = TRUE)[[1]]
  ranks <- main_ranks[seq_along(p)]
  diff <- which(!is.na(lin[ranks]) & lin[ranks] != p & ranks %in% c("phylum", "class", "order", "family"))
  if (length(diff) == 0L) return(NA_character_)
  paste0("the reference places it in ", paste(sprintf("%s %s (%s in taxa)", ranks[diff], p[diff],
                                                      lin[ranks[diff]]), collapse = ", "))
}

# Output ------------------------------------------------------------------------

assemble_estimates <- function(ab, tx, rows, pred, level, reference) {
  pr <- pred$pred
  n <- nrow(ab)
  models <- reference$models
  km <- match(pr$model, models$model)
  tau2 <- pr$v_param + pr$v_missing + pr$v_population + pr$v_source
  z <- stats::qnorm(1 - (1 - level) / 2)
  ok <- pr$status == "estimated"
  ind_median <- ifelse(ok, exp(pr$eta), NA_real_)
  ind_mean <- ifelse(ok, exp(pr$eta + tau2 / 2), NA_real_)
  ind_lower <- ifelse(ok, exp(pr$eta - z * sqrt(tau2)), NA_real_)
  ind_upper <- ifelse(ok, exp(pr$eta + z * sqrt(tau2)), NA_real_)
  a <- ab$abundance
  status <- pr$status
  reason <- pr$reason
  na_ab <- ok & is.na(a)
  status[na_ab] <- "not_estimated"
  reason[na_ab] <- "abundance is missing"
  out <- ab
  out$broad_group <- tx$broad_group[rows$ti]
  if ("scientific_name" %in% names(tx) && !"scientific_name" %in% names(ab)) {
    out$scientific_name <- as_text(tx$scientific_name)[rows$ti]
  }
  out$life_stage <- tx$life_stage[rows$ti]
  out$sex <- tx$sex[rows$ti]
  out$count_unit <- rows$count_unit
  out$domain <- rows$domain
  out$method <- pr$method
  out$model <- pr$model
  out$stage_used <- pr$stage_used
  out$length_mm <- ifelse(is.na(pr$length_used), rows$length_mm, pr$length_used)
  out$length_type <- rows$length_type
  out$length_basis <- rows$length_basis
  out$length_n <- rows$length_n
  out$length_status <- pr$length_status
  out$match_level <- pr$match_level
  out$match_name <- pr$match_name
  out$n_reference <- pr$n_reference
  out$sources <- pr$sources
  out$currency <- ifelse(ok, models$currency[km], NA_character_)
  out$unit <- ifelse(ok, models$unit[km], NA_character_)
  out$ind_mean <- ind_mean
  out$ind_median <- ind_median
  out$ind_lower <- ind_lower
  out$ind_upper <- ind_upper
  out$biomass_mean <- a * ind_mean
  out$biomass_median <- a * ind_median
  out$biomass_lower <- a * ind_lower
  out$biomass_upper <- a * ind_upper
  out$biomass_unit <- ifelse(ok, paste0(out$unit, "/", as_text(ab$per)), NA_character_)
  out$level <- rep(level, n)
  out$interval <- ifelse(ok, "prediction", NA_character_)
  out$uncertainty <- ifelse(ok, pr$uncertainty, "not_quantified")
  out$corrections <- pr$corrections
  out$status <- status
  out$reason <- reason
  out$note <- pr$note
  rownames(out) <- NULL
  out
}

new_bmassa_estimate <- function(estimates, rows, pred, reference, lineage, settings,
                                lengths_ind = NULL) {
  pr <- pred$pred
  sim <- data.frame(row = rows$row, taxon_id = rows$taxon_id, site = rows$site,
                    stage_used = pr$stage_used, model = pr$model, method = pr$method,
                    length_type = rows$length_type, length_used = pr$length_used,
                    length_cv = rows$length_cv, length_basis = rows$length_basis,
                    length_group = rows$length_group, node = pr$node, stage_node = pr$stage_node,
                    contrast = pr$contrast, cells_per_unit = rows$cells_per_unit,
                    xc = pr$xc, offset = pr$offset, eta = pr$eta, v_param = pr$v_param, v_missing = pr$v_missing,
                    v_population = pr$v_population, v_source = pr$v_source,
                    stringsAsFactors = FALSE)
  st <- estimates$status
  diagnostics <- list(
    n_rows = nrow(estimates),
    n_estimated = sum(st == "estimated"),
    by_method_level = as.data.frame(table(method = estimates$method[st == "estimated"],
                                          match_level = factor(estimates$match_level[st == "estimated"],
                                                               levels = c(rev(main_ranks), "root"))),
                                    stringsAsFactors = FALSE),
    reasons = as.data.frame(table(reason = estimates$reason[st != "estimated"]),
                            stringsAsFactors = FALSE)
  )
  diagnostics$by_method_level <- diagnostics$by_method_level[diagnostics$by_method_level$Freq > 0, ]
  if (is.null(lengths_ind)) lengths_ind <- data.frame(group = integer(0), length_mm = numeric(0),
                                                       n = numeric(0))
  structure(list(estimates = estimates, settings = settings, diagnostics = diagnostics,
                 simulation = list(rows = sim, missing_keys = pred$missing_keys,
                                   lineage = `rownames<-`(lineage, NULL), lengths = lengths_ind),
                 reference = reference),
            class = "bmassa_estimate")
}

#' @export
as.data.frame.bmassa_estimate <- function(x, ...) {
  x$estimates
}

#' @export
print.bmassa_estimate <- function(x, ...) {
  e <- x$estimates
  ok <- e$status == "estimated"
  cat(sprintf("<bmassa_estimate> %d rows: %d estimated (%d from lengths, %d without lengths), %d not estimated\n",
              nrow(e), sum(ok), sum(ok & e$method %in% "length"), sum(ok & e$method %in% "typical"),
              sum(!ok)))
  if (any(ok)) {
    tab <- table(factor(e$match_level[ok], levels = c(rev(main_ranks), "root")))
    tab <- tab[tab > 0]
    cat("Lowest level shared with the reference: ",
        paste(sprintf("%s %d", names(tab), as.integer(tab)), collapse = ", "), "\n", sep = "")
    cat(sprintf("Interval: %s %g %%, for the mean individual mass of a population%s\n",
                "prediction", 100 * x$settings$level,
                if (is.na(x$settings$bias_correction)) "" else
                  paste0("; bias_correction = \"", x$settings$bias_correction, "\"")))
  }
  cv <- x$settings$conversion
  if (!is.null(cv)) {
    unit <- if (is.null(cv$unit)) default_unit(cv$to) else cv$unit
    cat(sprintf("Currency: %s (%s); %d rows converted from the currency of their models; habitat %s\n",
                gsub("_", " ", cv$to), unit, sum(!is.na(e$converted_from)),
                if (is.null(cv$habitat)) "unknown" else cv$habitat))
  }
  if (any(!ok)) {
    r <- table(e$reason[!ok])
    cat("Not estimated:\n")
    for (i in seq_along(r)) cat(sprintf("  %3d  %s\n", as.integer(r[i]), names(r)[i]))
  }
  cat("Use as.data.frame() for the table and summarize_biomass() for totals.\n")
  invisible(x)
}
