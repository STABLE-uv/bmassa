# bmassa 0.6.0

## Changes that alter results

* The carbon of invertebrates and fish comes from more data, those of the
  nitrogen and phosphorus below: the models of carbon per dry mass have 1336
  factors of invertebrates (563 of freshwater animals) instead of 911 (144),
  and 132 of fish (62) instead of 97 (27). Of the 61 taxon units of animals
  of the vignette `azoxystrobin`, the carbon per dry mass changes by less
  than 5 % in 54 and by less than 10 % in all but the ostracods (from 38 % to
  33 %); planarians (*Dugesia*), without factors before, now go to carbon.
* The model of invertebrates without lengths gains the typical mass of the
  harpacticoid copepods of Hébert et al. (2016), left out before for lack of
  a classification (Harpacticoida is not in the COL XR; now placed in
  Podoplea, as for the nitrogen and phosphorus below). Taxa estimated from
  their species or genus do not change; some estimated from their order or
  class do: in the vignette `azoxystrobin`, the copepod nauplii (class
  Copepoda) 0.73 times their previous mass, the ostracods 0.83, the
  cyclopoid and calanoid copepods (orders) 0.85 and 0.90, *Dugastella
  valentina* 1.04; the other 57 units of invertebrates change by less than
  2 %.
* Amphibians given as larvae, or without a stage and without a snout-vent
  length, are tadpoles if they are frogs or toads, and are estimated (below);
  before, they were not. `bias_correction` is therefore required when a
  survey has tadpoles without lengths. Larvae of salamanders and newts, and
  juveniles, are still not estimated.
* The conversion models of invertebrates gain the soft tissue of freshwater
  molluscs from articles and theses (Ravera et al. 2007, Benke et al. 1999,
  Doherty 1986 and nine others, listed by `cite_sources()`): 34 factors of
  17 species of dry to wet mass, AFDM to dry mass, carbon and energy (and
  three shares of the shell in the dry mass). The edible flesh of snails,
  animals exposed to diuron or to high calcium in experiments and clams
  starved in the laboratory are left out. The zebra and quagga mussels of
  the Szczecin Lagoon of Coughlan et al. (2021) are now transitional instead
  of freshwater. The fit of the carbon per dry mass of invertebrates now
  reaches the better of its two optima (the previous fit had stopped at one
  without variance between phyla, with a worse REML criterion for the same
  data), which changes most the taxa matched only by their phylum. In the
  vignette `azoxystrobin`, the AFDM of *Melanopsis*, *Theodoxus* and
  *Physella acuta* (estimated from the class Gastropoda) rises by 2 %, and
  the carbon of *Physella acuta* by 3 %, of the ostracods by 5 % and of the
  water mites by 3 %; the other estimates change by less than 1 %.
* Snails from their shell height are about twice as heavy as before. The
  snails of the length-mass data (the CADSES lagoons) were weighed with their
  shell, taken out before with the ratio of the dry mass without shell to the
  ash-free dry mass with shell of freshwater bivalves, which left them 0.05 to
  0.19 of their dry mass. Now it is taken out with the share of the soft
  tissue in the dry mass of snails, from articles and a thesis on their shells
  (Kołodziejczyk 1984, Sokolowski et al. 2003, Kamler 1969, James et al.
  1998, listed by `cite_sources()`), by their family, order or class: 0.17 to
  0.20. As that share does not follow the shell of each species, it is capped
  by the ash of the taxon: the soft tissue cannot hold more organic matter
  than the whole snail, which lowers it to 0.10-0.16 in four snails with
  heavy shells (*Tritia neritea*, *Bittium reticulatum*, *Theodoxus
  fluviatilis* and *Pirenella*). A snail estimated from the class Gastropoda
  (such as *Physella acuta* or *Melanopsis*) is 1.82 times as heavy at the
  same shell height, and *Theodoxus* 2.26 times; those with thin shells are
  probably still underestimated (see the known limitations in the README).
  Bivalves change by less than 5 %.
* The ash-free dry mass of the soft tissue of snails is 0.75 of their dry
  mass instead of 0.68 (class Gastropoda; 3 to 8 % more in the orders and
  species with factors): four snails of the CADSES leaf packs without an
  order in the data set had been taken as molluscs without shell, so their
  ash with the shell entered the model of the soft tissue (and a sea slug had
  been taken as shelled). The factors of Coughlan et al. (2021) pooled over
  the sites of a species are no longer used along with those of each site,
  which they repeat (as with the length-mass equations, they would only stand
  in for species without sites of three or more individuals), and the dry
  mass of the soft tissue over the wet mass with shell gains 60 snails weighed
  one by one (Crichton 2003): with both, the molluscs of Hrycik et al.
  (2024), weighed wet with their shell, have 3 % more dry mass.
* All the models are fitted from four starting points (`fit_model()` in
  `data-raw/fit_models.R`), keeping the fit with the best REML criterion. The
  length-mass models of fish and of adult amphibians reach better optima than
  their previous fits (REML criterion -2441.2 instead of -2378.9, and -26.2
  instead of -15.0): fish known by their family get narrower intervals (95 %
  factor 1.4 instead of 2.0) and those known only by their order wider ones
  (2.8 instead of 2.0); in the example survey, *Squalius pyrenaicus* (by its
  genus) is 3 % heavier and *Pelophylax perezi* 0.8 % lighter. The other
  models have no better optimum.
* In the vignette `azoxystrobin`, without lengths: *Melanopsis* and
  *Theodoxus* 1.5 % heavier and *Physella acuta* 0.7 % (more dry mass in the
  molluscs of Hrycik et al. 2024), and their AFDM 10 to 11 % more; the AFDM
  of the other animals changes by -1.9 to +1.1 %, and the rest by less than
  0.1 %.
* The conversion models of invertebrates gain the ash and water of freshwater
  zooplankton from articles and a report (Riccardi and Mangoni 1999, Birge
  and Juday 1922, Watanabe et al. 1978 and Hasan et al. 2023, listed by
  `cite_sources()`): 29 factors of AFDM to dry mass of 14 taxa of cladocerans
  and copepods, and 4 of dry to wet mass of *Daphnia* and *Moina*. The energy
  densities of Riccardi and Mangoni (1999), calculated from the elemental
  composition, and the moisture of zooplankton air-dried on filter paper and
  stored before it was weighed (Ovie and Ovie 2006) are not used. The dry
  mass of *Daphnia* is 0.132 of its wet mass instead of 0.145: in the
  vignette `azoxystrobin`, *Daphnia* has 9.7 % more wet mass (as *Daphnia
  magna* in the example survey), *Ceriodaphnia*, *Simocephalus* and
  *Scapholeberis* 4.6 % and the other cladocerans 2.1 %; the AFDM of the
  animals changes by -0.2 to +0.5 %, and the rest by less than 0.1 %.
* The models of carbon, nitrogen and phosphorus per dry mass gain the aquatic
  animals of the Limno-STOICH database (Corman et al. 2026; release of
  2025-12-11), averaged by population: 9091 factors of 65 of its data
  sources, some unpublished. The models of invertebrates have 4999, 5005 and
  1707 factors instead of 1337, 1346 and 581, those of fish 343, 553 and 191
  instead of 132, 342 and 93, and those of tadpoles gain the wood frog
  (*Boreorana sylvatica*, *Lithobates sylvaticus* in the database) and
  *Pseudacris maculata* (carbon and nitrogen, Gosner stages 25 to 42). Left out: data sets already used, parts of animals (also when
  entered as whole organisms), eggs and parasites, terrestrial animals and
  the emerged adults of aquatic insects, fish analysed as muscle (Kiffney et
  al. 2002; Volk and Kiffney 2012), molluscs with the carbon and nitrogen of a
  shell, and amphibians other than those tadpoles. The intervals of carbon and
  nitrogen narrow (invertebrates known by their family: 95 % factor 1.33
  instead of 1.48 for carbon) and those of phosphorus widen (2.27 instead of
  1.89). In the vignette `azoxystrobin`, the carbon of the animals changes by
  -13 to +6 % (54 of 62 within 5 %), their nitrogen by -21 to +12 % and their
  phosphorus by -33 to +25 % (median -10 %); in the example survey, the carbon
  of *Limnodrilus hoffmeisteri* is 15 % lower and the phosphorus of
  *Gammarus* 21 % higher. The masses do not change.
* The conversion models gain the ash and energy of four triclads (Calow and
  Jennings 1974) and the carbon of 14 species of fish of Polish lakes
  (Penczak et al. 1985): 22 factors. Farmed animals fed formulated diets
  (bass, bullfrog tadpoles), triclads starved in the laboratory and trout
  farmed in cages are left out, and so are a leech already in the data bank of
  Brey et al. (2010) and the water mites and triclads of Morante et al. (2012).
  The estimates of the vignette `azoxystrobin` and of the example survey
  change by less than 1 %: the AFDM and energy of *Dugesia* by +0.6 %, the
  carbon of the fish of the example survey by up to +0.9 %.

## New features

* Nitrogen and phosphorus: `currency` in `estimate_biomass()` and `to` in
  `convert_currency()` take `"nitrogen"` and `"phosphorus"`, with models of
  their content per dry mass of invertebrates and fish (issue #3). They are
  fitted to the data bank of Brey et al. (2010), the compilation of Hébert
  et al. (2016) for crustacean zooplankton, the appendix of Evans-White et
  al. (2005) for stream invertebrates and data sets and tables of the
  elemental composition of freshwater animals (Beck et al. 2022; González et
  al. 2018, 2025; Moody et al. 2018; May et al. 2026; and others, listed by
  `cite_sources()`) and the Limno-STOICH database (Corman et al. 2026, below),
  averaged by population. Phytoplankton and adult
  amphibians are not converted to them, nor rotifers to phosphorus.
* Tadpoles: estimated in dry mass (without the gut) from their total length,
  with a length-mass model fitted to tadpoles of 11 species of four wetlands
  of Georgia, USA, measured one by one before the metamorphic climax (Knapp
  et al. 2021). Without lengths, from the typical length of their taxon, a
  model of the mean total length of the tadpoles of each sample of the same
  data (`method = "typical_length"`, `length_basis = "typical"`), with its
  uncertainty. They go to carbon, nitrogen and phosphorus with models of
  their own (13 species, with Vanni et al. 2002 through González et al.
  2025), but not to wet mass, ash-free dry mass or energy (no factors of
  tadpoles in the reference).

# bmassa 0.5.0

## Changes that alter results

* `estimate_biomass()` gives all its rows in one currency and unit: wet mass
  in mg by default (the only currency that all the groups reach), or the one
  chosen with the new argument `currency` (`"dry_mass"`, `"afdm"`,
  `"carbon"`, `"energy"`), with `unit` and `habitat`. The rows are converted
  from the currency of their models as `convert_currency()` does, with the
  uncertainty of the factors; rows that cannot reach the currency are not
  estimated, with the reason. Before, the result had phytoplankton in
  biovolume, invertebrates in dry mass and vertebrates in wet mass (issue #1).
* `convert_currency()` gives every row in mg (or J for energy) unless `unit`
  says otherwise; before, each group kept the unit of its estimate.

## New features

* The table `lengths` holds all the lengths (issue #5): the mean length of a
  taxon unit (`length_mean`, with `length_sd` or `length_cv`), the usual case,
  or the lengths of measured individuals (`length`, with `n` individuals of
  each length: one row per individual or per length class). `site` and `date`
  are optional: a row without `site` applies to all the samples of the unit,
  one with `site` to all the dates of the site, and one with both to that
  sample; each sample takes the most specific lengths. With individual
  lengths the mean mass is the mean of their masses (the mass at the length
  of their mean mass, reported in `length_mm`), without assuming a
  distribution; the new columns `length_basis` and `length_n` say which
  lengths were used.
* `summarize_size()` describes the size structure of each sample or group
  (issue #4): the mean and the median size of the individuals, weighted by
  their abundances, and, with `breaks`, the shares of the individuals and of
  the biomass in size classes, with prediction intervals from the draws of
  the masses. The size is the individual mass, in the currency of the
  estimate, or the length. It uses the lengths of measured individuals and
  the spread of mean lengths where `lengths` gives them; individuals known
  only by the mean size of their taxon make the median and the classes
  approximations (`approximation`, `share_taxon_mean`), and the share of
  individuals without a size is reported (`share_excluded`).
* Two vignettes with field data (issue #2): `vignette("tuesday-lake")`
  compares the predicted masses and biomasses with those measured in Tuesday
  Lake (91 % of the measured masses inside the 95 % intervals), and
  `vignette("azoxystrobin")` goes from the raw counts of a mesocosm experiment
  (Amador et al. 2024) to the effect of the fungicide on biomass and size
  structure with a generalized linear mixed model of the R package *glmmTMB*
  (Brooks et al. 2017, <https://doi.org/10.32614/RJ-2017-066>; now in
  Suggests), with the uncertainty of the masses carried by refitting to the
  draws of `simulate()`. Their data are in `inst/extdata` (`?azoxystrobin`,
  CC BY 4.0; `?tuesday_lake`, from the R package *cheddar*, Hudson et al.
  2013, <https://doi.org/10.1111/2041-210X.12005>, BSD 2-Clause; see the file
  `COPYRIGHTS`).
* `convert_abundance()` says in `per_assumption` what share of the site a
  sample is when it gives the totals of a site (a trap of 300 cm² in a tank
  of 1.54 m²: 1.95 %), and the help of `estimate_biomass()` and
  `convert_abundance()` explains how to give in `per` the size of each sample,
  the sampling effort (`"300 cm2"` for a trap of 30 × 10 cm, `"600 cm2"` for
  two such traps counted together, `"10 L"` for 10 L of filtered water).

## Bug fixes

* `standardize_taxa()`: a name whose ending implies its rank (*Orthocladiinae*,
  a subfamily) takes the usage of that rank before one of another rank (the
  COL XR also has a placeholder genus *Orthocladiinae*, which was taken).

## Performance

* `estimate_biomass()` and `check_coverage()` predict each taxon unit (with
  its lengths) once, however many samples it has: about 20 times faster with
  12 000 rows of 104 taxa. `summarize_biomass()`, `simulate()` and the new
  `summarize_size()` draw the deviations of missing taxonomic levels faster.
  The results are the same.

## Deprecated

* Lengths in `taxa` (`length_mean`, `length_unit`, `length_type`,
  `length_sd`, `length_cv`) still work in this version, after any row of
  `lengths`, with a warning: give them in `lengths`, in rows without `site`.

## Documentation

* New version of the logo, with the hexagon filled.
* R packages are named in italics in the vignettes, the README and the help
  pages, and the external ones are cited, with their references.

## Data

* `example_survey$taxa` has no lengths; `example_survey$lengths` holds the
  mean lengths of the units (rows without site), those of the chub in each
  site and the lengths of the *Gammarus* measured in L1, by length class.

# bmassa 0.4.0

## Changes that alter results

* The reference of fish includes the fish of the National Ecological
  Observatory Network of the United States (NEON; Fish electrofishing, gill
  netting, and fyke netting counts, DP1.20107.001, RELEASE-2026), measured in
  2016-2024 in 27 streams and lakes: the mean wet mass of 94 species in each
  sampling and gear (1,180 records), for fish without lengths, and 120
  length-mass equations by total length of 75 species. The model of fish
  without lengths goes from 57 to 145 species and from 21 to 31 families, and
  that of their length-mass relations from 84 to 157 species, lampreys
  included.
  The samplings of a species at a site are one population with a random
  effect (level `"population"` of `get_reference()$components`), so that the
  variation between populations and that between samplings of a population
  both enter the interval of a new population.
* The model of fish without lengths has no effect of source, like that of
  invertebrates since 0.3.0: its three sources come from different regions
  and share few species, and the effect of source put the species of NEON
  1.5 times above their own data, with wider intervals and the same coverage
  in cross-validation. The 95 % interval of a new population of a species of
  the reference narrows from a factor of about 8.5 to 7.3, and the masses of
  the fish species already in 0.3.0 change little (median ×0.94). In
  cross-validation, the 95 % intervals of fish without lengths cover 96 % of
  new populations and of new species (0.3.0: 94 % and 91 %). In the test
  with Tuesday Lake, the fish go from 8.3 to 1.2 times their measured masses:
  the dace and the mudminnow are now estimated from their own genus and
  species instead of from Cyprinidae and all fish.
* A fish without length of an order absent from the reference is estimated
  from the fish of its class (Teleostei) instead of from all fish.
* The taxonomy of the reference is that of the COL26.10 XR (it was the
  COL26.6 XR): the classification of 1,085 of the 8,059 names of the reference
  data changed, mostly the classes and orders of algae (Zygnematophyceae for
  Conjugatophyceae, Volvocales for Chlamydomonadales), with small changes in
  the models of phytoplankton and in the conversion factors.
* The length-mass relations of invertebrates are fitted without an effect of
  study, so that each species keeps the level of its own equations. With the
  effect, which predicted for an average study, the species measured by a
  single study moved away from them: the rotifers of Ejsmont-Karabin (1998)
  1.26 times above, the species of CADSES 1.33 to 1.37 times above, the
  bivalves of Coughlan et al. (2021) 0.80 times below. Masses from lengths of
  rotifers and of the benthos of transitional waters go down by about a
  fifth, and those of insects and microcrustaceans change little. In
  cross-validation, the 95 % interval of a new population of a species of
  the reference narrows from a factor of 3.7 to 3.3, with the same coverage
  (95 %; new species 96 %, new studies 94 %). The variance between studies,
  from a second fit with that effect (SD 0.34 in ln), is still the deviation
  of a new study that `summarize_biomass()` shares among all the
  invertebrates with lengths, because the equations of a study left out
  deviate together. It is taken out of the variance between equations
  (levels `"study"` and `"equation"` of `get_reference()$components`), so the
  interval of each row is that of the fit without study.
* The conversion models have no effect of source. They have two to four
  sources, most of their factors come from the data bank of Brey et al.
  (2010), which compiles many studies, and the effect was estimated as zero in
  8 of the 13 models with more than one source; where it was not, it put the
  species of the data bank below their own factors (0.94 times in dry mass to
  AFDM and 0.95 in wet to dry mass of invertebrates) with the same coverage
  in cross-validation. The conversions of invertebrates between wet mass, dry
  mass and AFDM change by a few per cent, and so do the masses of Hrycik et
  al. (2024) that the reference converts to dry mass (median ×0.98).

## Bug fixes

* *Thryssocypris tonlesapensis* (Koning et al. 2024) was missing from the
  model of fish without lengths.

## Documentation

* Copepods without a stage are still taken as adults, and the documentation
  and the note of those rows warn that a count of a species that includes
  its copepodites is then overestimated (1.4 to 3.6 times for the copepods
  of Tuesday Lake, whose means seem to include them): give the copepodites
  as rows of their own.
* The README shows the logo of the package.

# bmassa 0.3.0

## Changes that alter results

* The model of the masses of invertebrates without lengths has no effect of
  source. Its sources measure different groups of taxa by different methods,
  and the effect of source, which predicted for an average source, moved
  whole groups away from their own data: rotifers 2 to 4 times above their
  masses from volumes, benthic species measured only by Moreau et al. (2021)
  3 to 5 times below them. Comparing the sources only within rotifers,
  microcrustaceans and the other invertebrates did not remove the shift. Each
  taxon now keeps the level of the sources that measured it, the weighings of
  rotifers are corrected by their contrast (`c_rotifer_weighing`, now about
  four times the masses from volumes), and in cross-validation the 95 %
  intervals keep their coverage and are narrower. Masses without lengths
  change by group (median for species measured by a single source): rotifers
  about ×0.5, cladocerans and copepods ×1.3 to 1.7, and benthic species
  measured only by Hrycik et al. (2024) about ×0.6.
* The benthic records of Moreau et al. (2021) are a fixed contrast in the
  model of invertebrates without lengths (`c_moreau_benthos`): they are
  masses of individual prey compiled from studies of diet and calorimetry,
  often of large or adult animals, and weigh about 6 times (95 % interval
  3.6 to 11) the means of all the individuals of field samples of Hrycik et
  al. (2024) for the same taxa. Benthic invertebrates are predicted as
  population means, so benthic species measured only by Moreau et al. (2021)
  are about 0.4 times their mass in 0.2.0. In the test with Tuesday Lake,
  *Chaoborus* goes from 4.2 to 1.9 times the measured mass.

## Bug fixes

* `standardize_taxa()` prefers, among homonyms, the taxa expected in the
  broad group: in the zooplankton, rotifers, branchiopods, copepods and
  ostracods, then other arthropods; in the macroinvertebrates, the phyla of
  freshwater invertebrates. The zooplankton genus *Latona* was taken as the
  homonymous bivalve. A genus without family in the COL XR is now accepted
  when it has an order (the cladoceran *Sida* was not found).

# bmassa 0.2.0

## Changes that alter results

* The package is renamed from *biomassa* to *bmassa* (version 0.1.0 was
  *biomassa*): `library(bmassa)`, and the classes of the results and
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
* The zooplankton of the reference gains six sources: stage
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
  0.22 ug of dry mass), and more rotifers and cladocerans are estimated from
  their own species.
* Rotifers are predicted as from volumes, the usual method of field studies:
  direct weighings of their dry mass (Dumont et al. 1975; Bottrell et al.
  1976, Table VIII) are a fixed contrast in the model (`c_rotifer_weighing`,
  about twice the masses from volumes), not variation between sources. Their
  masses without lengths are about half of what the weighings would give, and
  the intervals of all the invertebrates without lengths are narrower (the
  variation between sources is smaller). A test with the fresh masses of
  Tuesday Lake (Jonsson et al. 2005) prompted the change.
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
  model of energy per carbon of animals; amphibians use the convention
  above, and phytoplankton the factor of Platt and Irwin (1973, below).
* Fish without lengths are estimated as whole populations (model
  `fish_individual`), from the mean wet mass of the populations of an
  electrofishing survey of Iberian streams (Arranz et al. 2022) and of the
  catch of a river fishery (Koning et al. 2024), all sizes together; the
  intervals are wide (the mean mass of a population depends on its sizes). A
  fish of an order absent from those data takes the mean of all of them; a
  fish given a life stage or a sex still needs a length.
* `cells_per_unit` in the `taxa` of `estimate_biomass()` and
  `check_coverage()`: for phytoplankton counted in natural units, the number
  of cells of a unit, which is then estimated as that many cells (and goes to
  carbon and energy). Without it, natural units of colonial genera get a
  note: in the reference they are those of its sources, often single cells or
  small groups, not whole colonies.
* `convert_currency(to = "energy")` converts phytoplankton cells (and natural
  units with `cells_per_unit`) from their carbon with the factor of Platt and
  Irwin (1973), 11.40 cal per mg of carbon (SD 1.1 among their 10 samples of
  marine net phytoplankton), an extrapolation for freshwater that the rows
  note.
* `standardize_taxa()`: classification and identifier (`col_id`) of a list of
  taxa in the Catalogue of Life (COL XR), through the GBIF API. Needs the
  suggested packages *curl* (Ooms 2023,
  <https://CRAN.R-project.org/package=curl>) and *jsonlite* (Ooms 2014,
  <https://arxiv.org/abs/1403.2805>), and internet; when the service cannot be
  reached or gives an answer that cannot be read, it stops with an error of
  class `bmassa_error_network`.
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
* *readxl* (Wickham and Bryan 2023,
  <https://CRAN.R-project.org/package=readxl>) is no longer a suggested
  package.

# bmassa 0.1.0

First version, published as *biomassa*.

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
