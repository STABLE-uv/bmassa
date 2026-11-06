# Currency conversion of estimates with the hierarchical conversion models -------------

conversion_targets <- c("wet_mass", "dry_mass", "afdm", "carbon", "energy")
mass_unit_g <- c(pg = 1e-12, ng = 1e-9, ug = 1e-6, mg = 1e-3, g = 1)
energy_unit_J <- c(mJ = 1e-3, J = 1, kJ = 1e3, MJ = 1e6)
# Wet mass of phytoplankton from biovolume: the density of water, 1 g/cm3 =
# 1 pg/um3 (EN 16695:2015, section 6.6), source "cen2015" of the reference.
phyto_density_pg_um3 <- 1
# Dry mass of rotifers as a fraction of their wet mass: the convention of US
# EPA (2016, LG403, section 8.2.2), source "usepa2016", and 0.039 for
# Asplanchna (Dumont et al. 1975, p. 94; also in US EPA 2016).
rotifer_dry_wet <- 0.1
asplanchna_dry_wet <- 0.039
# Amphibians (after metamorphosis): the values of frogs and toads of US EPA
# (1993, Table 4-1, p. 4-13), source "usepa1993", each the mean (SD) of 3
# studies: 85 % water (4.7), so a dry mass of 15 % of the wet mass, and a
# gross energy of 4.6 kcal/g of dry mass (0.45), in J/ug (4.184 J/cal).
amphibian_dry_wet <- 0.15
amphibian_dry_wet_sd <- 0.047
amphibian_energy_dry <- 4.6 * 4.184e-3
amphibian_energy_dry_sd <- 0.45 * 4.184e-3
amphibian_n_studies <- 3
# Energy of phytoplankton from its carbon: the general factor of Platt and
# Irwin (1973, p. 309), 11.40 cal/mg C, the mean of their Table 3 (10 samples
# of net phytoplankton of a marine spring bloom, gross energy by bomb
# calorimetry; SD 1.108 cal/mg C, computed from the table), source
# "platt1973", in J per pg of carbon (4.184 J/cal).
phyto_energy_carbon <- 11.40 * 4.184e-9
phyto_energy_carbon_sd <- 1.108 * 4.184e-9
phyto_energy_n <- 10
habitat_levels <- c("freshwater", "marine", "transitional")
# Classes of diatoms (carbon-volume equation of Menden-Deuer and Lessard 2000
# for diatoms; the protist one for the rest).
diatom_classes <- c(
  "Bacillariophyceae", "Coscinodiscophyceae", "Fragilariophyceae", "Mediophyceae"
)
prokaryote_names <- c("Bacteria", "Bacillati", "Monera", "Cyanobacteriota", "Cyanobacteria",
                      "Cyanophyta", "Cyanophycota", "Cyanoprokaryota")
conversion_columns <- c("converted_from", "conversion", "conversion_sources")

#' Convert estimates to another currency, with the uncertainty of the factors
#'
#' `convert_currency()` expresses the mean individual mass and the biomass of
#' each row of an estimate in another currency: wet mass, dry mass, ash-free
#' dry mass (AFDM), carbon or energy. Invertebrates and fish are converted with
#' hierarchical models of the conversion factors fitted to the data of Brey et
#' al. (2010), Hrycik et al. (2024), Coughlan et al. (2021), Cummins and
#' Wuycheck (1971), Evans-White et al. (2005) and the CADSES data sets: a
#' taxon gets the factor of its closest relatives, and the interval widens
#' with the uncertainty of the factor. Rotifers go between dry and wet mass
#' with the conventional ratio of US EPA (2016): dry mass 10 % of the wet
#' mass, 3.9 % in *Asplanchna*. Amphibians go between wet and dry mass and to
#' energy with the values of frogs and toads compiled by US EPA (1993): dry
#' mass 15 % (SD 4.7) of the wet mass and 4.6 kcal/g (SD 0.45) of dry mass.
#' Phytoplankton biovolumes are converted to wet mass with the conventional
#' density of 1 g/cm3 (EN 16695) and, for single cells (or natural units with
#' a known number of cells), to carbon with the carbon-volume equations of
#' Menden-Deuer and Lessard (2000), and from carbon to energy with the factor
#' of Platt and Irwin (1973), 11.40 cal/mg C (SD 1.108).
#'
#' @param data Result of [estimate_biomass()] (or of `convert_currency()` on
#'   it).
#' @param to Target currency: `"wet_mass"`, `"dry_mass"`, `"afdm"`,
#'   `"carbon"` or `"energy"`.
#' @param habitat `NULL` (default: unknown) or the habitat of the organisms,
#'   `"freshwater"`, `"marine"` or `"transitional"`. Several conversion models
#'   have a habitat effect: if it is unknown its variance is added to the
#'   uncertainty; if it is given its estimate is used.
#' @param unit `NULL` (default: mg for invertebrates, g for fish, pg per
#'   counting unit for the wet mass and carbon of phytoplankton and J for
#'   energy) or a unit for all the
#'   rows in `to`: for masses `"pg"`, `"ng"`, `"ug"`, `"mg"` or `"g"`; for
#'   energy `"mJ"`, `"J"`, `"kJ"` or `"MJ"`. A common unit is needed to
#'   compare or add up groups (with the same `per`).
#' @param ... Not used; an error if anything is given.
#'
#' @return An object of the same class
#'   with the same rows in the same order, in which the rows converted to `to`
#'   have new `currency`, `unit`, `ind_*`, `biomass_*` and `biomass_unit`,
#'   and three columns are added: `converted_from` (the currency of the
#'   estimate, `NA` if not converted), `conversion` (steps, models and the
#'   node of each, e.g. `"dry_mass->carbon: invertebrate dry_mass_to_carbon,
#'   family Chironomidae"`) and `conversion_sources`. `uncertainty` gains
#'   `"conversion"` when the uncertainty of the factors is included. Rows
#'   that cannot be converted get `NA`, `status = "not_estimated"` and the
#'   `reason`, and one warning reports how many there are. The conversion
#'   always starts from the currency of the estimate, so converting a
#'   converted result again does not add up steps. [summarize_biomass()]
#'   simulates the totals with the errors of the factors.
#'
#' @details
#' **Factors.** For a row whose log mass in its own currency \eqn{A} is
#' normal with mean \eqn{\eta_A} and variance \eqn{\tau_A^2}, the log mass in
#' \eqn{B} is
#' \deqn{\ln m_B = \ln m_A + s \ln f, \qquad \ln f \sim N(\eta_f, \tau_f^2)}{ln m_B = ln m_A + s ln f, ln f ~ N(eta_f, tau_f^2)}
#' where \eqn{\ln f} is predicted for the taxon by the conversion model of its
#' domain from \eqn{A} to \eqn{B} (\eqn{s = 1}) or, if there is none, from
#' \eqn{B} to \eqn{A} (\eqn{s = -1}); without either, the conversion goes
#' through dry mass (for fish, carbon comes from wet to dry mass and from dry
#' mass to carbon). As in [estimate_biomass()], \eqn{\eta_f} is the estimate of
#' the path to the deepest node of the model shared with the taxon (matched by
#' `col_id` or by name and rank), and \eqn{\tau_f^2} adds its conditional
#' variance, the variances of the missing levels, the variation between
#' records of a species (a new population) and between habitats (when
#' unknown). The models have no effect of source: most of their factors come
#' from the data bank of Brey et al. (2010), which compiles many studies, so
#' the differences between studies are part of the variation between records,
#' and two to four sources are too few to estimate a variance between them.
#' A source with its own methods may still deviate as a whole: in
#' cross-validation, with each source left out in turn, 66 % of the wet to
#' dry mass factors and 61 % of the dry mass to carbon factors of invertebrates
#' fell inside their 95 % intervals. The errors of the factor and of the mass are
#' independent: \eqn{\tau_B^2 = \tau_A^2 + \tau_f^2}{tau_B^2 = tau_A^2 + tau_f^2}.
#'
#' **Energy.** The factor is an energy density in J per ug of the mass of the
#' estimate, from the models of energy per dry mass (invertebrates) and per wet
#' mass (fish), or the convention of amphibians (below); a mass in mg or g is
#' first expressed in ug. The energy
#' densities of the data bank of Brey et al. (2010), of Cummins and Wuycheck
#' (1971) (measured by combustion on whole aquatic animals, from original
#' sources not in the data bank) and of a few other articles were screened
#' before the fit: values per AFDM outside 15-45 J/mg (the
#' combustion energy of carbohydrate, protein and lipid sets 17.5-39.5 J/mg;
#' Gnaiger and Bitterlich 1984), values inconsistent with another energy
#' density of the same material, values more than 4 robust standard
#' deviations from those of their class, and values with a standardized
#' residual beyond 4 in the fitted model were left out. There is no model of
#' energy per carbon (the carbon of some sources is total rather than
#' organic), nor of the energy of masses with shell or skeleton.
#'
#' **Rotifer wet mass.** The conversion models have no wet to dry mass
#' factors of rotifers, so their dry mass is taken as 10 % of the wet mass,
#' and as 3.9 % in *Asplanchna*, the convention of US EPA (2016, section
#' 8.2.2). EPA attributes the 10 % to Doohan (1973), whose article gives dry
#' weights of *Brachionus plicatilis* but no dry to wet mass ratio; the 3.9 %
#' is the dry weight to volume ratio of *Asplanchna brightwelli* at a density
#' of 1 (Dumont et al. 1975, p. 94). The uncertainty of the ratio is not
#' included (the rows say so in `note`), although it varies among rotifers.
#' Their AFDM, carbon and energy come from the conversion models.
#'
#' **Amphibian dry mass and energy.** The reference data have too few
#' conversion factors of amphibians for a model, so they take the values of
#' frogs and toads compiled by US EPA (1993, Table 4-1, p. 4-13), each the
#' mean (SD) of three studies: 85 % water (4.7), that is a dry mass of 15 % of
#' the wet mass, and a gross energy of 4.6 kcal/g (0.45) of dry mass, 19.2
#' J/mg (4.184 J/cal). A factor \eqn{f} with mean \eqn{\mu} and standard
#' deviation \eqn{s} over \eqn{n = 3} studies is taken as lognormal,
#' \deqn{\ln f \sim N(\ln \mu - v/2, v), \qquad v = \sigma^2 + \sigma^2/n, \qquad \sigma^2 = \ln(1 + s^2/\mu^2)}{ln f ~ N(ln mu - v/2, v), v = sigma^2 + sigma^2/n, sigma^2 = ln(1 + s^2/mu^2)}
#' where \eqn{\sigma^2} is the variation between populations and
#' \eqn{\sigma^2/n} the uncertainty of the mean of the studies, so that the
#' mean of a new population is \eqn{\mu}. The mean dry mass is then 15 % of
#' the mean wet mass (median 14.1 %, 95 % range 7.0-28 %), and the energy
#' 19.2 J per mg of dry mass on average (median 19.1 J/mg, 95 % range
#' 15.3-23.9 J/mg). [summarize_biomass()] shares the variation between
#' populations among the rows of a taxon as `population_effect` says, and the
#' uncertainty of the mean among all the amphibian rows. The energy per wet
#' mass of the table (1.2 kcal/g) comes from a single study and is not used.
#' Salamanders, newts and caecilians get the same values, marked as an
#' extrapolation in `conversion` and `note`. Amphibians are not converted to
#' AFDM or carbon.
#'
#' **Phytoplankton wet mass.** The wet mass is the biovolume times the
#' density of water, 1 g/cm3 = 1 pg/um3, the convention of EN 16695 (CEN 2015,
#' section 6.6: 1 mm3/L of biovolume is 1 mg/L of wet mass). The conversion is
#' linear, so cells and natural units (colonies, filaments) are converted
#' alike, and the interval keeps its width. The uncertainty of the density is
#' not included (the rows say so in `note`): measured densities of freshwater
#' algae are up to 10 % higher, 1.04-1.10 g/cm3 (Oliver et al. 1981, p. 289).
#'
#' **Phytoplankton carbon.** For rows counted in cells,
#' \deqn{\log_{10} C = \log_{10} a + b \log_{10} V}{log10(C) = log10(a) + b * log10(V)}
#' with \eqn{C} in pg C per cell and \eqn{V} in um3 per cell; for diatoms
#' \eqn{\log_{10} a = -0.541}{log10(a) = -0.541} and \eqn{b = 0.811}, for the
#' other protists \eqn{-0.665} and \eqn{0.939} (Menden-Deuer and Lessard 2000,
#' Table 4). The log carbon is then normal with mean
#' \eqn{\ln 10 \log_{10} a + b \eta_V}{ln(10) log10(a) + b eta_V} and variance
#' \eqn{b^2 \tau_V^2}{b^2 tau_V^2}. The uncertainty of the equation itself
#' is not included: Table 4 gives confidence intervals for \eqn{\log_{10} a}
#' and \eqn{b} but not their covariance nor the residual variance; the rows
#' say so in `note`. Cyanobacteria use the equation of protists, marked as an
#' extrapolation (the authors note that two cyanobacterial values fall close
#' to it). A natural unit of k cells (`cells_per_unit` in the `taxa` of
#' [estimate_biomass()]) has the carbon of k cells of its biovolume divided
#' by k: \eqn{\ln C = \ln 10 \log_{10} a + (1 - b) \ln k + b \ln V}{ln C = ln(10) log10(a) + (1 - b) ln k + b ln V}.
#'
#' **Phytoplankton energy.** From the carbon, with the general factor of
#' Platt and Irwin (1973, p. 309), 11.40 cal (47.7 J) per mg of carbon, the
#' mean of the 10 samples of their Table 3 (SD 1.108 cal/mg C, computed from
#' the table; coefficient of variation 10 %, half that of the energy per dry
#' mass). It is lognormal with that mean, as the conventions of amphibians:
#' variance between populations \eqn{\ln(1 + CV^2)}{ln(1 + CV^2)} plus a tenth
#' of it for the mean of the 10 samples. The samples were net phytoplankton
#' (mesh 76 um) of the spring bloom of 1969 in St. Margaret's Bay (Nova
#' Scotia), a marine coastal bay, and the energy is gross energy by bomb
#' calorimetry: its use for freshwater phytoplankton is an extrapolation (the
#' rows say so in `note`).
#'
#' @section Missing values and zeros:
#' Rows without an estimate stay `NA`. Rows that cannot be converted get
#' `NA`: amphibians to AFDM or carbon (the reference data have too few
#' conversion factors for them), phytoplankton to dry mass or AFDM (no
#' factors from biovolume), phytoplankton counted in natural units without
#' `cells_per_unit` to carbon or energy (the equation is for single cells,
#' and with \eqn{b < 1} it would underestimate the carbon of a colony),
#' and taxa that match no node of the conversion model (a phylum without
#' factors in it). An abundance of 0 keeps a biomass of 0.
#'
#' @section Assumptions:
#' The composition of the organisms is that of their relatives in the
#' reference data; mollusc masses are without shell, as in
#' [estimate_biomass()]. Factors are independent of the mass estimated and,
#' in a chain, of each other. The biovolume of a counting unit of
#' phytoplankton is that of its cells, with the cell wall and without
#' mucilage, lorica or spines, as defined by EN 16695 (section 3.3). The
#' carbon of the cell of mean volume stands for
#' the mean carbon per cell (with \eqn{b} between 0.8 and 0.94 and a
#' coefficient of variation of cell volumes of 0.3, the difference is below 1
#' %). Energy densities are gross energies of combustion (bomb calorimetry) of
#' whole animals, gut contents included as in most sources; they are not the
#' energy assimilated by a consumer. Carbon is as measured by each source,
#' total carbon in some of them (in crustaceans it then includes the
#' carbonate of the cuticle). Factors measured on animals kept in ethanol or
#' formalin were left out where the source says so, because preservation
#' changes the wet and dry mass by amounts that depend on the preservative,
#' the time and the species (Dermott and Paterson 1974); the data bank of Brey
#' et al. (2010) does not record the preservation.
#'
#' @references
#' Brey T, Müller-Wiegmann C, Zittier ZMC, Hagen W (2010). Body composition in
#' aquatic organisms - a global data bank of relationships between mass,
#' elemental composition and energy content. *Journal of Sea Research*, 64,
#' 334-340. \doi{10.1016/j.seares.2010.05.002}
#'
#' CEN (2015). *EN 16695:2015. Water quality - Guidance on the estimation of
#' phytoplankton biovolume*. European Committee for Standardization, Brussels.
#'
#' Cummins KW, Wuycheck JC (1971). Caloric equivalents for investigations in
#' ecological energetics. *Mitteilungen der Internationalen Vereinigung für
#' Theoretische und Angewandte Limnologie*, 18, 1-158.
#' \doi{10.1080/05384680.1971.11903918}
#'
#' Dermott RM, Paterson CG (1974). Determining dry weight and percentage dry
#' matter of chironomid larvae. *Canadian Journal of Zoology*, 52, 1243-1250.
#' \doi{10.1139/z74-165}
#'
#' Doohan M (1973). An energy budget for adult *Brachionus plicatilis* Muller
#' (Rotatoria). *Oecologia*, 13, 351-362.
#'
#' Dumont HJ, Van de Velde I, Dumont S (1975). The dry weight estimate of
#' biomass in a selection of Cladocera, Copepoda and Rotifera from the
#' plankton, periphyton and benthos of continental waters. *Oecologia*, 19,
#' 75-97. \doi{10.1007/BF00377592}
#'
#' Evans-White MA, Stelzer RS, Lamberti GA (2005). Taxonomic and regional
#' patterns in benthic macroinvertebrate elemental composition in streams.
#' *Freshwater Biology*, 50, 1786-1799.
#' \doi{10.1111/j.1365-2427.2005.01455.x}
#'
#' Gnaiger E, Bitterlich G (1984). Proximate biochemical composition and
#' caloric content calculated from elemental CHN analysis: a stoichiometric
#' concept. *Oecologia*, 62, 289-298. \doi{10.1007/BF00384259}
#'
#' Menden-Deuer S, Lessard EJ (2000). Carbon to volume relationships for
#' dinoflagellates, diatoms, and other protist plankton. *Limnology and
#' Oceanography*, 45, 569-579. \doi{10.4319/lo.2000.45.3.0569}
#'
#' Oliver RL, Kinnear AJ, Ganf GG (1981). Measurements of cell density of three
#' freshwater phytoplankters by density gradient centrifugation. *Limnology
#' and Oceanography*, 26, 285-294. \doi{10.4319/lo.1981.26.2.0285}
#'
#' Platt T, Irwin B (1973). Caloric content of phytoplankton. *Limnology and
#' Oceanography*, 18, 306-310. \doi{10.4319/lo.1973.18.2.0306}
#'
#' US EPA (1993). *Wildlife Exposure Factors Handbook*. EPA/600/R-93/187.
#' U.S. Environmental Protection Agency, Washington, DC.
#'
#' US EPA (2016). *Standard Operating Procedure for Zooplankton Analysis*.
#' LG403, Revision 07. U.S. Environmental Protection Agency, Great Lakes
#' National Program Office, Chicago.
#'
#' The other sources of the conversion factors are listed by [cite_sources()].
#'
#' @seealso [get_reference()] for the conversion models (`conversion_models`,
#'   `conversion_nodes`) and the carbon equations.
#' @export
#' @examples
#' est <- estimate_biomass(example_survey$abundance, example_survey$taxa,
#'                         example_survey$lengths, bias_correction = "lognormal")
#' carbon <- convert_currency(est, to = "carbon", unit = "ug")
#' head(as.data.frame(carbon)[, c("site", "taxon_id", "converted_from", "ind_median",
#'                                "ind_lower", "ind_upper", "unit", "status")])
#' summarize_biomass(carbon, by = c("site", "broad_group"), na_rm = TRUE, seed = 1)
#'
#' # Energy, in J per individual and per unit of area or volume
#' energy <- convert_currency(est, to = "energy", unit = "J")
#' head(as.data.frame(energy)[, c("site", "taxon_id", "converted_from", "ind_median",
#'                                "ind_lower", "ind_upper", "unit", "status")])
#'
#' # Wet mass in mg; phytoplankton at the density of water (1 mm3/L = 1 mg/L)
#' wet <- convert_currency(est, to = "wet_mass", unit = "mg")
#' head(as.data.frame(wet)[, c("site", "taxon_id", "converted_from", "ind_median",
#'                             "unit", "conversion")])
convert_currency <- function(data, ...) {
  UseMethod("convert_currency")
}

#' @rdname convert_currency
#' @export
convert_currency.default <- function(data, ...) {
  abort_input(sprintf("`data` must be the result of estimate_biomass(), not %s.",
                      describe_value(data)))
}

#' @rdname convert_currency
#' @export
convert_currency.bmassa_estimate <- function(data, to, habitat = NULL, unit = NULL, ...) {
  check_dots_empty(...)
  check_choice(to, conversion_targets)
  if (!is.null(habitat)) check_choice(habitat, habitat_levels)
  if (!is.null(unit)) {
    check_choice(unit, names(if (identical(to, "energy")) energy_unit_J else mass_unit_g))
  }
  native <- if (is.null(data$native)) data$estimates else data$native
  plan <- plan_conversions(native, data$simulation, data$reference, to, habitat, unit)
  out <- apply_conversions(native, plan, to, data$settings$level)
  failed <- !is.na(native$ind_median) & is.na(out$ind_median)
  data$native <- native
  data$estimates <- out
  data$simulation$conversion <- plan[c("steps", "missing_keys", "shift", "slope")]
  data$settings$conversion <- list(to = to, habitat = habitat, unit = unit)
  data$diagnostics$conversion <- as.data.frame(table(
    result = ifelse(failed, "not converted", ifelse(is.na(out$converted_from),
                                                    ifelse(is.na(out$ind_median), "no estimate",
                                                           "already in the currency"),
                                                    "converted"))), stringsAsFactors = FALSE)
  if (any(failed)) {
    warn_not_estimable(sprintf(
      "%d of %d rows could not be converted to %s; see the `reason` column of the result.",
      sum(failed), nrow(out), to))
  }
  data
}

# Plan of the conversion of each row: the linear map of its log mass (shift
# and slope), the extra variance, the steps through the conversion models and
# what to report.
plan_conversions <- function(native, simulation, reference, to, habitat, unit) {
  n <- nrow(native)
  sim <- simulation$rows
  lin <- simulation$lineage
  cm <- reference$conversion_models
  cn <- reference$conversion_nodes
  cc <- reference$conversion_components
  ce <- reference$carbon_equations
  comp_value <- function(m, l) {
    v <- cc$value[cc$model == m & cc$level == l & cc$parameter == "intercept"]
    if (length(v)) v else 0
  }
  shift <- rep(0, n)
  slope <- rep(1, n)
  v_add <- rep(0, n)
  new_unit <- native$unit
  converted <- rep(FALSE, n)
  quantified <- rep(TRUE, n)
  label <- rep(NA_character_, n)
  csrc <- rep(NA_character_, n)
  reason <- rep(NA_character_, n)
  note <- rep(NA_character_, n)
  steps <- list()
  missing_keys <- list()
  route <- function(dom, from, target) {
    direct <- paste(dom, paste0(from, "_to_", target))
    inverse <- paste(dom, paste0(target, "_to_", from))
    if (direct %in% cm$model) return(list(list(model = direct, sign = 1)))
    if (inverse %in% cm$model) return(list(list(model = inverse, sign = -1)))
    if (from != "dry_mass" && target != "dry_mass") {
      a <- route(dom, from, "dry_mass")
      b <- route(dom, "dry_mass", target)
      if (!is.null(a) && !is.null(b)) return(c(a, b))
    }
    NULL
  }
  # Planned once per taxon, model and currency; the other rows of the same
  # taxon (other sites) copy the plan.
  key <- paste(sim$taxon_id, sim$model, native$currency, sep = "\u001f")
  has_eta <- !is.na(sim$eta)
  reps <- which(has_eta)[!duplicated(key[has_eta])]
  for (i in reps) {
    from <- native$currency[i]
    dom <- native$domain[i]
    if (from == to) next
    if (dom == "amphibian") {
      if (!to %in% c("wet_mass", "dry_mass", "energy")) {
        reason[i] <- sprintf(paste("the reference data have too few conversion factors for",
                                   "amphibians to %s (by convention, they go only between wet",
                                   "and dry mass and to energy)"), gsub("_", " ", to))
        next
      }
      # Conventional factors of frogs and toads, lognormal (convention_lognormal()):
      # wet to dry mass, or its inverse, and dry mass to energy.
      ext <- if (lin$order[i] %in% "Anura") "" else "; extrapolated: not a frog or toad"
      conv <- list()
      if (from == "wet_mass" || to == "wet_mass") {
        conv[[length(conv) + 1L]] <- c(
          convention_lognormal(amphibian_dry_wet, amphibian_dry_wet_sd, amphibian_n_studies),
          list(sign = if (from == "wet_mass") 1 else -1,
               model = "amphibian wet_mass_to_dry_mass convention",
               label = sprintf("%s: amphibians, dry mass 15 %% (SD 4.7) of wet mass",
                               if (from == "wet_mass") "wet_mass->dry_mass" else
                                 "dry_mass->wet_mass")))
      }
      if (to == "energy") {
        conv[[length(conv) + 1L]] <- c(
          convention_lognormal(amphibian_energy_dry, amphibian_energy_dry_sd,
                               amphibian_n_studies),
          list(sign = 1, model = "amphibian dry_mass_to_energy convention",
               label = "dry_mass->energy: amphibians, 4.6 kcal/g (SD 0.45) of dry mass"))
      }
      for (cv in conv) {
        shift[i] <- shift[i] + cv$sign * cv$eta
        v_add[i] <- v_add[i] + cv$v_population + cv$v_mean
        steps[[length(steps) + 1L]] <- data.frame(
          row = i, type = "convention", model = cv$model, node = NA_integer_, sign = cv$sign,
          eta = cv$eta, v_population = cv$v_population, v_source = cv$v_mean, v_habitat = 0,
          habitat_node = NA_integer_, stringsAsFactors = FALSE)
        missing_keys[[length(missing_keys) + 1L]] <- character(0)
      }
      if (to == "energy") {
        shift[i] <- shift[i] + log(mass_unit_g[[native$unit[i]]] / mass_unit_g[["ug"]])
        new_unit[i] <- "J"
      }
      converted[i] <- TRUE
      label[i] <- paste(sprintf("%s (frogs and toads, US EPA 1993%s)",
                                vapply(conv, `[[`, character(1), "label"), ext), collapse = "; ")
      csrc[i] <- "usepa1993"
      if (nzchar(ext)) {
        note[i] <- "values of frogs and toads, extrapolated to other amphibians"
      }
      next
    }
    if (dom == "phytoplankton") {
      if (to == "wet_mass") {
        # Linear in the biovolume, so natural units are converted too (their
        # volume is that of their cells, without mucilage, as in EN 16695).
        shift[i] <- log(phyto_density_pg_um3)
        new_unit[i] <- "pg"
        converted[i] <- TRUE
        quantified[i] <- FALSE
        label[i] <- "biovolume->wet_mass: density 1 g/cm3 (EN 16695)"
        csrc[i] <- "cen2015"
        note[i] <- paste("wet mass with the conventional density of 1 g/cm3 (uncertainty not",
                         "quantified; measured densities of algal cells are 1.04-1.10 g/cm3)")
        steps[[length(steps) + 1L]] <- data.frame(
          row = i, type = "density", model = NA_character_, node = NA_integer_, sign = 1,
          eta = NA_real_, v_population = 0, v_source = 0, v_habitat = 0, habitat_node = NA_integer_,
          stringsAsFactors = FALSE)
        missing_keys[[length(missing_keys) + 1L]] <- character(0)
        next
      }
      if (!to %in% c("carbon", "energy")) {
        reason[i] <- sprintf("no conversion from biovolume to %s in the reference", gsub("_", " ", to))
        next
      }
      if (!identical(sim$model[i], "phytoplankton_cell")) {
        reason[i] <- paste("the carbon-volume equation is for single cells: rows counted in",
                           "natural units are converted only with `cells_per_unit` in `taxa`")
        next
      }
      li <- lin[i, ]
      diatom <- li$phylum %in% "Bacillariophyta" | li$class %in% diatom_classes
      if (!diatom && is.na(li$phylum) && is.na(li$class)) {
        reason[i] <- "taxonomic group unknown (diatom or not) for the carbon-volume equation"
        next
      }
      eq <- ce[ce$equation_id == if (diatom) "mdl2000_diatoms" else "mdl2000_protist", ]
      # Per cell; a natural unit of k cells has k times the carbon of a cell of
      # its biovolume divided by k: ln C = ln10 log10(a) + (1 - b) ln k + b ln V.
      k <- if (!is.null(sim$cells_per_unit) && !is.na(sim$cells_per_unit[i])) sim$cells_per_unit[i] else 1
      shift[i] <- log(10) * eq$log10_a + (1 - eq$b) * log(k)
      slope[i] <- eq$b
      new_unit[i] <- "pg"
      converted[i] <- TRUE
      quantified[i] <- FALSE
      prok <- li$kingdom %in% prokaryote_names | li$phylum %in% prokaryote_names
      label[i] <- sprintf("biovolume->carbon: Menden-Deuer and Lessard 2000, %s%s",
                          if (diatom) "diatoms" else "protists other than diatoms",
                          if (prok) " (extrapolated: prokaryote)" else "")
      csrc[i] <- eq$source_id
      note[i] <- "uncertainty of the carbon-volume equation not quantified"
      steps[[length(steps) + 1L]] <- data.frame(
        row = i, type = "carbon", model = NA_character_, node = NA_integer_, sign = 1,
        eta = NA_real_, v_population = 0, v_source = 0, v_habitat = 0, habitat_node = NA_integer_,
        stringsAsFactors = FALSE)
      missing_keys[[length(missing_keys) + 1L]] <- character(0)
      if (to == "energy") {
        # Energy from the carbon, lognormal (convention_lognormal()).
        cv <- convention_lognormal(phyto_energy_carbon, phyto_energy_carbon_sd, phyto_energy_n)
        shift[i] <- shift[i] + cv$eta
        v_add[i] <- v_add[i] + cv$v_population + cv$v_mean
        new_unit[i] <- "J"
        quantified[i] <- TRUE
        label[i] <- paste0(label[i], "; carbon->energy: 11.40 cal/mg C (SD 1.11; Platt and Irwin 1973)")
        csrc[i] <- paste(eq$source_id, "platt1973", sep = ";")
        note[i] <- paste("energy from carbon with the factor of marine net phytoplankton of Platt",
                         "and Irwin (1973), extrapolated; gross energy; uncertainty of the",
                         "carbon-volume equation not quantified")
        steps[[length(steps) + 1L]] <- data.frame(
          row = i, type = "convention", model = "phytoplankton carbon_to_energy convention",
          node = NA_integer_, sign = 1, eta = cv$eta, v_population = cv$v_population,
          v_source = cv$v_mean, v_habitat = 0, habitat_node = NA_integer_,
          stringsAsFactors = FALSE)
        missing_keys[[length(missing_keys) + 1L]] <- character(0)
      }
      next
    }
    if (dom == "invertebrate" && lin$phylum[i] %in% "Rotifera" &&
        setequal(c(from, to), c("dry_mass", "wet_mass"))) {
      # Rotifers: the conversion models have no wet to dry mass factors of
      # rotifers; a conventional ratio instead.
      asp <- lin$genus[i] %in% "Asplanchna"
      r <- if (asp) asplanchna_dry_wet else rotifer_dry_wet
      shift[i] <- if (to == "wet_mass") -log(r) else log(r)
      converted[i] <- TRUE
      quantified[i] <- FALSE
      label[i] <- sprintf("%s->%s: rotifers, dry mass %s %% of wet mass (%s)", from, to,
                          if (asp) "3.9" else "10",
                          if (asp) "Asplanchna, Dumont et al. 1975" else "US EPA 2016")
      csrc[i] <- if (asp) "dumont1975;usepa2016" else "usepa2016"
      note[i] <- "conventional dry to wet mass ratio of rotifers (uncertainty not quantified)"
      steps[[length(steps) + 1L]] <- data.frame(
        row = i, type = "convention", model = NA_character_, node = NA_integer_, sign = 1,
        eta = NA_real_, v_population = 0, v_source = 0, v_habitat = 0, habitat_node = NA_integer_,
        stringsAsFactors = FALSE)
      missing_keys[[length(missing_keys) + 1L]] <- character(0)
      next
    }
    rt <- route(dom, from, to)
    if (is.null(rt)) {
      reason[i] <- sprintf("no conversion model from %s to %s for %s", gsub("_", " ", from),
                           gsub("_", " ", to), dom)
      next
    }
    li <- unlist(lin[i, main_ranks])
    ok <- TRUE
    labs <- character(0)
    srcs <- character(0)
    row_steps <- list()
    row_keys <- list()
    eta_add <- 0
    var_add <- 0
    for (st in rt) {
      m <- st$model
      lv <- strsplit(cm$levels[cm$model == m], ",")[[1]]
      cand <- which(cn$model == m & cn$level %in% lv)
      mt <- match_lineage(li, cn[cand, ], lv, col_id = lin$col_id[i])
      if (is.null(mt$row)) {
        reason[i] <- sprintf("no taxon of the conversion model %s matches the classification", m)
        ok <- FALSE
        break
      }
      j <- cand[mt$row]
      eta <- cn$path_intercept[j]
      vp <- cn$path_var[j]
      miss <- lv[seq_along(lv) > match(cn$level[j], lv)]
      vm <- sum(vapply(miss, function(l) comp_value(m, l), numeric(1)))
      keys <- vapply(miss, function(l) {
        nm <- li[[l]]
        paste(m, l, if (is.na(nm)) paste0("<", sim$taxon_id[i], ">") else
          paste(li[main_ranks[seq_len(match(l, main_ranks))]], collapse = "|"), sep = "\u001f")
      }, character(1))
      v_pop <- comp_value(m, "residual")
      v_src <- comp_value(m, "source")
      v_hab <- 0
      h_node <- NA_integer_
      if (any(cc$model == m & cc$level == "habitat")) {
        if (is.null(habitat)) {
          v_hab <- comp_value(m, "habitat")
        } else {
          k <- which(cn$model == m & cn$level == "habitat" & cn$node == habitat)
          if (length(k) == 1L) {
            eta <- eta + cn$effect[k]
            v_hab <- cn$effect_var[k]
            h_node <- k
          } else {
            v_hab <- comp_value(m, "habitat")
            note[i] <- append_note(note[i], sprintf("no %s factors in %s: habitat taken as unknown",
                                                    habitat, m))
          }
        }
      }
      eta_add <- eta_add + st$sign * eta
      var_add <- var_add + vp + vm + v_pop + v_src + v_hab
      ends <- unlist(cm[cm$model == m, c("from", "to")])
      if (st$sign < 0) ends <- rev(ends)
      labs <- c(labs, sprintf("%s->%s: %s%s, %s %s", ends[1], ends[2],
                              if (st$sign < 0) "inverse of " else "", m, cn$level[j], cn$node[j]))
      srcs <- c(srcs, strsplit(cn$sources[j], ";", fixed = TRUE)[[1]])
      if (!is.na(mt$note)) note[i] <- append_note(note[i], mt$note)
      row_steps[[length(row_steps) + 1L]] <- data.frame(
        row = i, type = "factor", model = m, node = j, sign = st$sign, eta = eta,
        v_population = v_pop, v_source = v_src, v_habitat = v_hab, habitat_node = h_node,
        stringsAsFactors = FALSE)
      row_keys[[length(row_keys) + 1L]] <- unname(keys)
    }
    if (!ok) next
    shift[i] <- eta_add
    if (to == "energy") {
      # Energy densities are in J/ug: the mass of the estimate is first
      # expressed in micrograms.
      shift[i] <- shift[i] + log(mass_unit_g[[native$unit[i]]] / mass_unit_g[["ug"]])
      new_unit[i] <- "J"
    }
    v_add[i] <- var_add
    converted[i] <- TRUE
    label[i] <- paste(labs, collapse = "; ")
    csrc[i] <- paste(sort(unique(srcs)), collapse = ";")
    steps <- c(steps, row_steps)
    missing_keys <- c(missing_keys, row_keys)
  }
  rep_of <- reps[match(key, key[reps])]
  copy <- which(has_eta & rep_of != seq_len(n))
  if (length(copy)) {
    src <- rep_of[copy]
    shift[copy] <- shift[src]
    slope[copy] <- slope[src]
    v_add[copy] <- v_add[src]
    new_unit[copy] <- new_unit[src]
    converted[copy] <- converted[src]
    quantified[copy] <- quantified[src]
    label[copy] <- label[src]
    csrc[copy] <- csrc[src]
    reason[copy] <- reason[src]
    note[copy] <- note[src]
  }
  # Unit: the factors are ratios of masses, so a converted mass keeps the unit
  # of the estimate (mg, g); carbon of cells is in pg and energy in J. `unit`
  # rescales every row in `to`.
  in_to <- !is.na(sim$eta) & (converted | native$currency %in% to)
  if (!is.null(unit) && to == "energy") {
    k <- in_to & new_unit %in% names(energy_unit_J)
    shift[k] <- shift[k] + log(energy_unit_J[new_unit[k]] / energy_unit_J[[unit]])
    new_unit[k] <- unit
  } else if (!is.null(unit)) {
    k <- in_to & new_unit %in% names(mass_unit_g)
    shift[k] <- shift[k] + log(mass_unit_g[new_unit[k]] / mass_unit_g[[unit]])
    new_unit[k] <- unit
  }
  # Steps of the representative rows, then one copy per row of the taxon.
  tmpl <- if (length(steps)) do.call(rbind, steps) else
    data.frame(row = integer(0), type = character(0), model = character(0), node = integer(0),
               sign = numeric(0), eta = numeric(0), v_population = numeric(0),
               v_source = numeric(0), v_habitat = numeric(0), habitat_node = integer(0),
               stringsAsFactors = FALSE)
  rr <- which(has_eta)
  idx_list <- lapply(rr, function(r) which(tmpl$row == rep_of[r]))
  idx <- unlist(idx_list)
  steps <- tmpl[idx, , drop = FALSE]
  steps$row <- rep(rr, lengths(idx_list))
  rownames(steps) <- NULL
  missing_keys <- missing_keys[idx]
  list(shift = shift, slope = slope, v_add = v_add, unit = new_unit, converted = converted,
       in_to = in_to, quantified = quantified, label = label, sources = csrc, reason = reason,
       note = note, steps = steps, missing_keys = missing_keys, sim = sim)
}

# Log-scale parameters of a conventional factor given as the mean and SD of n
# studies: lognormal with that mean and SD between populations
# (`v_population`), plus the uncertainty of the mean of the studies
# (`v_mean`, shared by all the rows); `eta` keeps the mean of a new
# population at `mean`.
convention_lognormal <- function(mean, sd, n) {
  v_population <- log1p((sd / mean)^2)
  v_mean <- v_population / n
  list(eta = log(mean) - (v_population + v_mean) / 2, v_population = v_population,
       v_mean = v_mean)
}

apply_conversions <- function(native, plan, to, level) {
  out <- native
  sim <- plan$sim
  tau2 <- sim$v_param + sim$v_missing + sim$v_population + sim$v_source
  eta <- plan$shift + plan$slope * sim$eta
  tau2 <- plan$slope^2 * tau2 + plan$v_add
  z <- stats::qnorm(1 - (1 - level) / 2)
  ok <- plan$in_to
  fail <- !is.na(plan$reason)
  out$currency[ok] <- to
  out$unit[ok] <- plan$unit[ok]
  a <- native$abundance
  set_ind <- function(col, value) {
    v <- out[[col]]
    v[ok] <- value[ok]
    v[fail] <- NA_real_
    v
  }
  out$ind_median <- set_ind("ind_median", exp(eta))
  out$ind_mean <- set_ind("ind_mean", exp(eta + tau2 / 2))
  out$ind_lower <- set_ind("ind_lower", exp(eta - z * sqrt(tau2)))
  out$ind_upper <- set_ind("ind_upper", exp(eta + z * sqrt(tau2)))
  for (col in c("mean", "median", "lower", "upper")) {
    out[[paste0("biomass_", col)]] <- a * out[[paste0("ind_", col)]]
  }
  out$biomass_unit[ok] <- paste0(out$unit[ok], "/", as_text(native$per)[ok])
  out$currency[fail] <- NA_character_
  out$unit[fail] <- NA_character_
  out$biomass_unit[fail] <- NA_character_
  out$interval[fail] <- NA_character_
  out$uncertainty[fail] <- "not_quantified"
  q <- plan$converted & plan$quantified & !fail
  out$uncertainty[q] <- paste0(out$uncertainty[q], "; conversion")
  was <- out$status == "estimated" & fail
  out$status[fail] <- "not_estimated"
  out$reason[fail] <- ifelse(was[fail], plan$reason[fail],
                             append_note(out$reason[fail], plan$reason[fail]))
  out$note <- append_note(out$note, plan$note)
  out$converted_from <- ifelse(plan$converted & !fail, native$currency, NA_character_)
  out$conversion <- ifelse(plan$converted & !fail, plan$label, NA_character_)
  out$conversion_sources <- ifelse(plan$converted & !fail, plan$sources, NA_character_)
  rownames(out) <- NULL
  out
}
