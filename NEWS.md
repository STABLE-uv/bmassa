# bmassa 0.2.0

## Changes that alter results

* The package is renamed from biomassa to bmassa (version 0.1.0 was
  biomassa): `library(bmassa)`, and the classes of the results and
  conditions start with `bmassa_` (`bmassa_estimate`, `bmassa_error_input`,
  ...).
* Masses are in mg (invertebrates, dry mass; molluscs without shell) and g
  (fish and amphibians, wet mass), not in µg as in 0.1.0; biovolumes are
  still in µm³. The `unit` column says which.
* The estimates come from hierarchical models fitted to all the reference
  data, not from the reference values of the taxon or its relatives: a taxon
  with data of its own is shrunk towards its relatives, and one without data
  is predicted from them, with a wider interval.
* `convert_currency()` uses hierarchical models of the
  conversion factors of Brey et al. (2010), Hrycik et al. (2024), Coughlan et
  al. (2021), Evans-White et al. (2005; carbon of 46 taxa of stream
  invertebrates) and the CADSES data sets, with their uncertainty, instead of
  factors given by the user. Factors measured on animals kept in ethanol or
  formalin are left out where the source says so (Dermott and Paterson
  1974). The energy per dry mass of the 19 records of Salonen et al. (1976)
  in the data bank of Brey et al. (2010), computed there with the carbon
  content instead of the ash-free fraction, is not used.
* `summarize_biomass()` draws all the rows with a mass, also those with a
  missing abundance, so that its draws are those of `simulate()`. With a
  seed, its medians and intervals can change within the Monte Carlo error
  when such rows bring nodes or populations of their own (not in the example
  data). Its result has a new column, `population_by`.
* The zooplankton of the reference gains six sources (phase 4k): stage
  weights and length-dry mass regressions of copepods and cladocerans of
  Culver et al. (1985); length-dry mass regressions compiled by Bottrell et
  al. (1976) and McCauley (1984), each equation once (from its original
  compilation, not again from a later one); individual weights of rotifers
  of Bottrell et al. (1976); equations of the wet mass of 222 rotifers and
  standard weights of 250 species of Ejsmont-Karabin (1998), converted to dry
  mass with the convention of 10 % (3.9 % in *Asplanchna*); and the mean dry
  mass of 139 freshwater cladocerans and copepods from their typical lengths
  (Hébert et al. 2016). Masses without lengths change: copepod nauplii and
  copepodites no longer take the mass of adults (nauplii of Copepoda, 9.8 to
  0.22 ug of dry mass), rotifers are mostly 1.1 to 3 times lighter (direct
  weighings of Dumont et al. 1975 are about seven times the masses from
  volumes for the same species; the model takes the difference as variation
  between sources), and more rotifers and cladocerans are estimated from
  their own species.
* Copepods without lengths have stages: copepodites and nauplii are fixed
  contrasts against the adults in the model of invertebrates (`stage_used`
  `"copepod_adult"`, `"copepod_copepodite"`, `"copepod_nauplius"`); a
  copepod without a stage is taken as an adult (`"copepod_adult
  (assumed)"`), and one given as `"juvenile"` is not estimated. With
  lengths, copepod nauplii use relations of their own (length type
  `"nauplius_body_length"`, chosen from `life_stage`).
* Orders and classes that the Catalogue of Life does not give are taken from
  its nearest rank above them (a superorder for an order) before the
  classification of the sources, so that sister families share their node
  in the models (the families of Hygrophila were split among Basommatophora,
  Pulmonata and Hygrophila). Copepods take the orders Calanoida, Cyclopoida
  and Harpacticoida, and the orders of free-living flatworms the class
  Turbellaria; a taxon given as Calanoida now matches its order.

## New features

* `estimate_biomass()`: mean individual mass of a population of each taxon
  (mean, median and prediction interval) and biomass, for phytoplankton,
  zooplankton, macroinvertebrates, fish and amphibians, from two tables: the
  taxon units (a taxon, or a taxon in one life stage and sex, with its
  classification, its counting unit and, where measured, its mean length)
  and their abundances by site (and date, if any); a third, optional table
  holds lengths measured in particular samples.
* `summarize_biomass()` works on the results of `estimate_biomass()` and
  gives totals with prediction intervals by simulation, with shared draws for
  the rows that share nodes of the reference, missing taxa, populations
  (`population_effect`) and the source (the mean of a total is the exact sum
  of the means of its rows; the median and the interval come from the
  draws).
* `simulate()` for results of `estimate_biomass()` and `convert_currency()`:
  joint draws of the mean individual mass and the biomass of every row, with
  the errors shared as in the totals of `summarize_biomass()`, for models that
  use the biomass of every taxon, such as the energy fluxes of a food web;
  `as.data.frame()` gives one row per row and draw.
* `convert_abundance()`: the abundances of the abundance table per a common
  volume, area or whole site, from the size of each sample in `per` ("L",
  "5 L", "0.25 m2"...) and, when needed, the depth, volume or area of each
  site, so that the biomasses of groups sampled per litre and per square
  metre can be added up; the factor and the assumption of each row are kept.
* `population_by` in `summarize_biomass()` and `simulate()`: the populations
  of a taxon are shared only within the groups of some columns of the
  abundance table, such as habitats stocked from different water bodies.
* `convert_currency()` works on the results of `estimate_biomass()` and
  converts to wet mass, dry mass, ash-free dry mass, carbon (phytoplankton
  cells with Menden-Deuer and Lessard 2000) or energy, with the uncertainty of
  the factors in the intervals and in the totals. Phytoplankton biovolumes
  (cells and natural units) go to wet mass with the conventional density of
  1 g/cm3 of EN 16695 (1 mm3/L = 1 mg/L), without uncertainty; measured
  densities of algal cells are up to 10 % higher (Oliver et al. 1981).
  Rotifers go between dry and wet mass with the convention of US EPA (2016):
  dry mass 10 % of the wet mass, 3.9 % in *Asplanchna* (Dumont et al. 1975).
  Amphibians go between wet and dry mass and to energy with the values of
  frogs and toads compiled by US EPA (1993, Table 4-1; means of three
  studies): dry mass 15 % of the wet mass (85 % water, SD 4.7) and 4.6
  kcal/g of dry mass (19.2 J/mg, SD 0.45), with their uncertainty in the
  intervals and the totals; they are not converted to ash-free dry mass or
  carbon.
* `convert_currency(to = "energy")` (units `"mJ"`, `"J"`, `"kJ"`, `"MJ"`)
  uses hierarchical models of energy per dry mass, ash-free dry mass and wet
  mass of invertebrates and fish, fitted to the energy densities of Brey et al.
  (2010), Cummins and Wuycheck (1971) and other articles after a screening:
  values per ash-free dry mass outside 15-45 J/mg (the combustion energy of
  carbohydrate, protein and lipid; Gnaiger and Bitterlich 1984), values
  inconsistent with another value of the same material, and values far from
  those of their class or from the fitted model are left out. There is no
  model of energy per carbon; amphibians use the convention above, and
  phytoplankton is not converted to energy.
* `standardize_taxa()`: classification and identifier (`col_id`) of a list of
  taxa in the Catalogue of Life (COL XR), through the GBIF API. Needs the
  suggested packages curl and jsonlite, and internet; when the service cannot
  be reached or gives an answer that cannot be read, it stops with an error
  of class `bmassa_error_network`.
* `estimate_biomass()` and `check_coverage()` look a taxon up by its `col_id`
  first, when `taxa` has one.
* `check_coverage()`: the taxonomic level from which each taxon will be
  estimated, and the width of its interval.
* `cite_sources()`: the sources behind an estimate and its conversions, with
  citation, DOI and licence; `get_reference()`: the reference tables. The
  package holds the parameters of the models fitted to the data of the
  sources (and the equations published in them), not the data.
* Example data set `example_survey` in the input format of
  `estimate_biomass()`.

## Removed

* The functions of 0.1.0: `estimate_mass_typical()`, `estimate_mass_length()`
  and `calc_biomass()` (replaced by `estimate_biomass()`), `match_taxa()`
  (`check_coverage()`), `calc_phyto_carbon()` (`convert_currency(to =
  "carbon")`) and `read_brey_databank()` (the factors of Brey et al. 2010 are
  in the conversion models), with the methods of `summarize_biomass()` and
  `convert_currency()` for their results: both now need the result of
  `estimate_biomass()` and give an error otherwise.
* The data sets of 0.1.0: `typical_mass`, `length_mass`, `carbon_volume` and
  `taxonomy` (their data are in the models of the reference), and
  `example_counts` and `example_lengths` (replaced by `example_survey`).
* readxl is no longer a suggested package.

# bmassa 0.1.0

First version, published as biomassa.

## Changes that alter results

* None (first version).

## New features

* `estimate_mass_typical()`: size of one counting unit (individual, cell or
  natural unit) from typical sizes: phytoplankton biovolumes of Kremer et al.
  (2014) and dry masses of Dumont et al. (1975) and Moreau et al. (2021)
  (`typical_mass`).
* `estimate_mass_length()`: individual dry mass from body length with the
  length-dry mass regressions of van Leeuwen et al. (2025) and Dumont et al.
  (1975) (`length_mass`), within their calibration range unless
  `extrapolate = TRUE`.
* `match_taxa()`: the reference records that represent each taxon (taxon
  itself, taxa within it, and relatives only up to `fallback_rank`), with a
  taxonomy table (`taxonomy`) that can be replaced or extended.
* `calc_biomass()`, `calc_phyto_carbon()` (carbon-volume regressions of
  Menden-Deuer and Lessard 2000, `carbon_volume`), `convert_currency()` and
  `read_brey_databank()` (factors from a local copy of the data bank of Brey
  et al. 2010), and `summarize_biomass()`.
* Example data sets `example_counts` and `example_lengths`, and a vignette.
