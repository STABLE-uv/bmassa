
<!-- README.md is generated from README.Rmd. Edit README.Rmd and render it with rmarkdown::render("README.Rmd", output_format = rmarkdown::github_document(html_preview = FALSE)) after installing the package. -->

# bmassa

<!-- badges: start -->

[![R-CMD-check](https://github.com/STABLE-uv/bmassa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/STABLE-uv/bmassa/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

bmassa estimates the mean individual mass and the biomass of aquatic
taxa (phytoplankton, zooplankton, macroinvertebrates, fish and
amphibians) from abundance data, with prediction intervals. Hierarchical
models fitted to published compilations and data sets predict each taxon
from its own reference data or, when it has none, from its closest
relatives, and the interval widens with the substitution. Masses are in
the native currency of each group: biovolume (µm³) for phytoplankton,
dry mass (mg; without shell in molluscs) for invertebrates and wet mass
(g) for fish and amphibians.

## Installation

From GitHub, with its dependencies:

``` r
# install.packages("remotes")
remotes::install_github("STABLE-uv/bmassa", dependencies = TRUE)
```

bmassa needs R (\>= 3.6.0) and the packages stats and withr. curl and
jsonlite are needed only by `standardize_taxa()`, which queries the
Catalogue of Life through the GBIF API.

## Example

Two tables: the taxon units (a taxon, or a taxon in one life stage and
sex, with its classification, its counting unit and its mean length
where it was measured), and their abundances by site (and date, if any).
A third, optional table holds lengths measured in particular samples.
Constructed example data shipped with the package:

``` r
library(bmassa)
tx <- example_survey$taxa
ab <- example_survey$abundance
ln <- example_survey$lengths
tx[c(1, 8, 19, 20), c("taxon_id", "broad_group", "genus", "species", "life_stage",
                      "count_unit", "length_mean", "length_unit", "length_type")]
#>    taxon_id   broad_group      genus                 species life_stage count_unit length_mean
#> 1       t01 phytoplankton Cyclotella Cyclotella meneghiniana       <NA>       cell          NA
#> 8       t08   zooplankton    Daphnia           Daphnia magna       <NA> individual         2.2
#> 19      t19     amphibian Pelophylax       Pelophylax perezi      adult individual        58.0
#> 20      t20     amphibian Pelophylax       Pelophylax perezi      larva individual          NA
#>    length_unit       length_type
#> 1         <NA>              <NA>
#> 8           mm       body_length
#> 19          mm snout_vent_length
#> 20        <NA>              <NA>
head(ab, 3)
#>   site taxon_id abundance per
#> 1   L1      t01      1800  mL
#> 2   L1      t02     25000  mL
#> 3   L1      t03       300  mL
ln
#>   taxon_id site length_mean length_unit length_type length_sd
#> 1      t17   L1         110          mm fork_length        20
#> 2      t17   L2          85          mm fork_length        15
```

For your own list of names, `standardize_taxa()` returns this
classification and the identifier of each taxon (`col_id`) from the
Catalogue of Life (COL XR); it needs internet.

Before estimating, the taxonomic level from which each taxon will come,
and the factor by which the median is divided and multiplied to get its
95 % interval:

``` r
check_coverage(tx, ln)
#> <bmassa_coverage> 20 taxa, 22 rows (14 without lengths, 8 with lengths)
#>   phytoplankton      typical species 4, genus 1; median factor 12
#>   zooplankton        typical species 3; median factor 12
#>   zooplankton        length  species 1; median factor 3
#>   macroinvertebrate  typical species 1, genus 1, family 1; median factor 14
#>   macroinvertebrate  length  species 2, family 1; median factor 3.5
#>   fish               typical species 1, genus 1; median factor 11
#>   fish               length  species 1, genus 2; median factor 1.4
#>   amphibian          typical none, not estimable 1; median factor -
#>   amphibian          length  species 1; median factor 1.7
#> Factor: the 95 % prediction interval is the median divided and multiplied by it.
```

The mean individual mass of a population of each taxon and the biomass
of each row. With lengths, the back-transformation of the length-mass
relations has to be chosen: `"lognormal"` gives the mean mass of the
individuals of each length, which biomass needs:

``` r
est <- estimate_biomass(ab, tx, ln, bias_correction = "lognormal")
#> Warning in estimate_biomass(ab, tx, ln, bias_correction = "lognormal"): 3 of 40 rows were not
#> estimated; see the `reason` column of the result.
est
#> <bmassa_estimate> 40 rows: 37 estimated (13 from lengths, 24 without lengths), 3 not estimated
#> Lowest level shared with the reference: species 27, genus 6, family 4
#> Interval: prediction 95 %, for the mean individual mass of a population; bias_correction = "lognormal"
#> Not estimated:
#>     1  abundance is missing
#>     2  amphibian larvae and juveniles have no reference masses without lengths
#> Use as.data.frame() for the table and summarize_biomass() for totals.
e <- as.data.frame(est)
print(e[e$site == "L1", c("taxon_id", "method", "match_level", "ind_median", "ind_lower",
                          "ind_upper", "unit")], digits = 3)
#>    taxon_id  method match_level ind_median ind_lower ind_upper unit
#> 1       t01 typical     species   2.79e+03  2.48e+02  3.14e+04  um3
#> 2       t02 typical     species   1.10e+02  1.01e+01  1.20e+03  um3
#> 3       t03 typical       genus   1.00e+03  4.76e+01  2.12e+04  um3
#> 4       t04 typical     species   5.48e+01  4.49e+00  6.68e+02  um3
#> 5       t05 typical     species   1.01e+02  8.43e+00  1.21e+03  um3
#> 6       t06 typical     species   2.18e-05  1.87e-06  2.54e-04   mg
#> 7       t07 typical     species   8.53e-05  7.25e-06  1.00e-03   mg
#> 8       t08  length     species   5.71e-02  1.87e-02  1.74e-01   mg
#> 9       t09 typical     species   6.35e-03  5.35e-04  7.54e-02   mg
#> 10      t10  length     species   5.19e-01  1.68e-01  1.61e+00   mg
#> 11      t11  length     species   6.07e-01  1.74e-01  2.12e+00   mg
#> 12      t12  length      family   1.39e+01  2.97e+00  6.54e+01   mg
#> 13      t13 typical       genus   4.34e-01  3.21e-02  5.86e+00   mg
#> 14      t14 typical      family   4.00e-01  1.22e-02  1.31e+01   mg
#> 15      t15 typical     species   1.90e-01  1.55e-02  2.33e+00   mg
#> 16      t16  length     species   5.55e+01  4.56e+01  6.76e+01    g
#> 17      t17  length       genus   2.14e+01  1.48e+01  3.10e+01    g
#> 18      t18 typical     species   7.14e-01  5.55e-02  9.18e+00    g
#> 19      t19  length     species   2.31e+01  1.33e+01  3.99e+01    g
#> 20      t20 typical        <NA>         NA        NA        NA <NA>
```

In another currency, here carbon in µg for all the groups, with the
uncertainty of the conversion factors (amphibians have too few factors
in the reference data: they go to dry mass and energy with the values of
frogs and toads of US EPA 1993, but not to carbon):

``` r
est_c <- convert_currency(est, to = "carbon", unit = "ug")
#> Warning in convert_currency.bmassa_estimate(est, to = "carbon", unit = "ug"): 2 of 40 rows could
#> not be converted to carbon; see the `reason` column of the result.
ec <- as.data.frame(est_c)
print(ec[ec$site == "L1", c("taxon_id", "converted_from", "ind_median", "ind_lower",
                            "ind_upper", "unit")], digits = 3)
#>    taxon_id converted_from ind_median ind_lower ind_upper unit
#> 1       t01      biovolume   1.79e-04  2.52e-05  1.28e-03   ug
#> 2       t02      biovolume   1.79e-05  1.90e-06  1.68e-04   ug
#> 3       t03      biovolume   1.42e-04  8.13e-06  2.49e-03   ug
#> 4       t04      biovolume   9.28e-06  8.86e-07  9.71e-05   ug
#> 5       t05      biovolume   1.65e-05  1.60e-06  1.69e-04   ug
#> 6       t06       dry_mass   8.28e-03  6.64e-04  1.03e-01   ug
#> 7       t07       dry_mass   3.05e-02  2.50e-03  3.72e-01   ug
#> 8       t08       dry_mass   2.60e+01  7.83e+00  8.62e+01   ug
#> 9       t09       dry_mass   2.96e+00  2.36e-01  3.72e+01   ug
#> 10      t10       dry_mass   2.30e+02  6.89e+01  7.69e+02   ug
#> 11      t11       dry_mass   2.81e+02  7.20e+01  1.10e+03   ug
#> 12      t12       dry_mass   6.38e+03  1.17e+03  3.47e+04   ug
#> 13      t13       dry_mass   1.58e+02  1.12e+01  2.23e+03   ug
#> 14      t14       dry_mass   1.67e+02  4.98e+00  5.62e+03   ug
#> 15      t15       dry_mass   7.84e+01  5.02e+00  1.23e+03   ug
#> 16      t16       wet_mass   5.86e+06  3.45e+06  9.94e+06   ug
#> 17      t17       wet_mass   2.21e+06  1.16e+06  4.20e+06   ug
#> 18      t18       wet_mass   8.26e+04  6.00e+03  1.14e+06   ug
#> 19      t19           <NA>         NA        NA        NA <NA>
#> 20      t20           <NA>         NA        NA        NA <NA>
```

Totals by site and group, by simulation with shared draws: two taxa
estimated from the same family, or the same taxon in two sites, have
correlated errors. Currencies and units are never added together:

``` r
tot <- summarize_biomass(est_c, by = c("site", "broad_group"), na_rm = TRUE, seed = 1)
print(tot[, c("site", "broad_group", "biomass_unit", "biomass_mean", "biomass_lower",
              "biomass_upper", "n_missing")], digits = 3)
#>    site       broad_group biomass_unit biomass_mean biomass_lower biomass_upper n_missing
#> 1    L1         amphibian         <NA>           NA            NA            NA         2
#> 2    L1              fish    ug/100 m2     5.43e+07      3.22e+07      9.54e+07         0
#> 3    L1 macroinvertebrate        ug/m2     3.63e+05      1.17e+05      9.88e+05         0
#> 4    L1     phytoplankton        ug/mL     1.71e+00      2.78e-01      6.31e+00         0
#> 5    L1       zooplankton         ug/L     3.11e+02      6.52e+01      1.11e+03         0
#> 6    L2         amphibian         <NA>           NA            NA            NA         2
#> 7    L2              fish    ug/100 m2     5.18e+07      2.81e+07      1.13e+08         0
#> 8    L2 macroinvertebrate        ug/m2     3.40e+05      8.86e+04      1.16e+06         0
#> 9    L2     phytoplankton        ug/mL     1.98e+00      2.61e-01      8.25e+00         0
#> 10   L2       zooplankton         ug/L     1.08e+02      8.83e+00      4.86e+02         1
```

The draws themselves, for models that use the biomass of every taxon
(such as the energy fluxes of a food web): the mean individual mass and
the biomass of every row in every draw, with the same shared errors as
the totals. Here all the groups per square metre of each site
(`convert_abundance()`, from the depth of each site) and in wet mass, so
that they can be added up, with the populations of each habitat apart
(`population_by`):

``` r
sites <- data.frame(site = c("L1", "L2"), depth = c(1.5, 0.8), depth_unit = "m")
ab_m2 <- convert_abundance(ab, to = "m2", sites = sites)
ab_m2$habitat <- ifelse(ab_m2$site == "L1", "lake", "spring")
est_w <- convert_currency(estimate_biomass(ab_m2, tx, ln, bias_correction = "lognormal"),
                          to = "wet_mass", unit = "mg")
sim <- simulate(est_w, nsim = 1000, seed = 1, population_by = "habitat")
d <- as.data.frame(sim, draws = 1:2)
print(head(d[, c("draw", "site", "taxon_id", "habitat", "ind_mass", "biomass",
                 "biomass_unit")]), digits = 3)
#>   draw site taxon_id habitat ind_mass biomass biomass_unit
#> 1    1   L1      t01    lake 1.39e-06  3754.3        mg/m2
#> 2    1   L1      t02    lake 2.20e-08   823.3        mg/m2
#> 3    1   L1      t03    lake 3.62e-07   162.9        mg/m2
#> 4    1   L1      t04    lake 1.37e-08   185.4        mg/m2
#> 5    1   L1      t05    lake 2.59e-08    23.3        mg/m2
#> 6    1   L1      t06    lake 1.44e-03   757.2        mg/m2
```

## Methods and decisions

- The models were fitted with lme4 to the reference data (one per domain
  and kind of data): masses without lengths with nested random effects
  by taxonomic level, length-mass relations with random intercepts and
  slopes by level and length type, and conversion factors between
  currencies with nested random effects. The package predicts from their
  tables with base R; `get_reference()` returns them.
- A taxon is matched by its `col_id`, if given, and otherwise by name
  and rank, from its lowest rank upwards, to the taxonomy of the
  Catalogue of Life (COL XR) used for the reference.
- The interval is for the mean individual mass of a population of the
  taxon: it includes the uncertainty of the model, of the taxonomic
  substitution and the variation between populations and between
  sources, but not the error of the abundances.
- The totals and `simulate()` share the errors of rows estimated from
  the same node of the reference, of the same population and of the same
  conversion factors; `simulate()` returns the draws of the totals, row
  by row.
- Insects without lengths and without a stage are taken as larvae, and
  copepods as adults (their copepodites and nauplii are fixed
  contrasts); fish without lengths are estimated as whole populations
  (all sizes), with wide intervals; amphibians only after metamorphosis
  (larvae are not estimated), and without lengths only when declared
  adults.
- Phytoplankton counted in natural units takes the natural units of the
  reference, often single cells or small groups in colonial genera; give
  `cells_per_unit` for colonies counted whole.
- Lengths outside those of the reference relations of a taxon use its
  closest relative that covers them, and are not extrapolated unless
  `extrapolate = TRUE`.
- The vignette (`vignette("bmassa")`) goes through the workflow and the
  decisions; the cross-validation of the intervals is in the source
  repository (`dev/validacion/validacion_modelos.md`).

## Data sources

``` r
s <- cite_sources()
knitr::kable(s[, c("source_id", "citation", "licence")])
```

| source_id          | citation                                                                                                                                                                                                                                                                                                                                | licence                                                                                                                                 |
|:-------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------|
| arranz2022         | Arranz I, Brucet S, Bartrons M, García-Comas C, Alcaraz C, Bardina M, Navarro Barquero P, Casals F, Caiola N, Duran MC, García-Berthou E, Maceda-Veiga A, Munné A, Rodríguez-Pérez MJ, Solà C, de Sostoa A, Benejam L (2022). Individual body mass and length dataset for over 12,000 fish from Iberian streams. Data in Brief, 108248. | Zenodo record: Other (Open), with no further description; article in Data in Brief                                                      |
| bottrell1976       | Bottrell HH, Duncan A, Gliwicz ZM, Grygierek E, Herzig A, Hillbricht-Ilkowska A, Kurasawa H, Larsson P, Weglenska T (1976). A review of some problems in zooplankton production studies. Norwegian Journal of Zoology 24, 419-456.                                                                                                      | Values transcribed from the paper                                                                                                       |
| brey2010           | Brey T, Müller-Wiegmann C, Zittier ZMC, Hagen W (2010). Body composition in aquatic organisms - a global data bank of relationships between mass, elemental composition and energy content. Journal of Sea Research 64, 334-340.                                                                                                        | No licence stated                                                                                                                       |
| cadses_boxcorer    | Rosati I, Basset A (2020). Benthic macroinvertebrates communities data collected in 15 coastal lagoons of Mediterranean and Black sea. LifeWatch Italy.                                                                                                                                                                                 | CC BY 4.0 (LifeWatch Italy; EML metadata of the download)                                                                               |
| cadses_leafpack    | Rosati I, Basset A (2021). Benthic macroinvertebrates communities data collected in 15 coastal lagoons of Mediterranean and Black sea by using the leaf pack technique. LifeWatch Italy.                                                                                                                                                | CC BY-NC 4.0 (LifeWatch Italy)                                                                                                          |
| cen2015            | CEN (2015). EN 16695:2015 Water quality - Guidance on the estimation of phytoplankton biovolume. European Committee for Standardization, Brussels.                                                                                                                                                                                      | Convention of the standard; no data used                                                                                                |
| ciancio2006        | Ciancio JE, Pascual MA (2006). Energy density of freshwater Patagonian organisms. Ecología Austral 16, 91-94.                                                                                                                                                                                                                           | Values transcribed from the paper                                                                                                       |
| coughlan2021       | Coughlan NE, Cunningham EM, Cuthbert RN, et al. (2021). Biometric conversion factors as a unifying platform for comparative assessment of invasive freshwater bivalves. Journal of Applied Ecology 58, 1945-1956.                                                                                                                       | CC0 1.0 (Dryad)                                                                                                                         |
| culver1985         | Culver DA, Boucherle MM, Bean DJ, Fletcher JW (1985). Biomass of freshwater crustacean zooplankton from length-weight regressions. Canadian Journal of Fisheries and Aquatic Sciences 42, 1380-1390.                                                                                                                                    | Values transcribed from the paper                                                                                                       |
| cummins1971        | Cummins KW, Wuycheck JC (1971). Caloric equivalents for investigations in ecological energetics. Mitteilungen der Internationalen Vereinigung für Theoretische und Angewandte Limnologie 18, 1-158.                                                                                                                                     | Values transcribed from the paper                                                                                                       |
| dumont1975         | Dumont HJ, Van de Velde I, Dumont S (1975). The dry weight estimate of biomass in a selection of Cladocera, Copepoda and Rotifera from the plankton, periphyton and benthos of continental waters. Oecologia 19, 75-97.                                                                                                                 | Values transcribed from the paper                                                                                                       |
| ejsmontkarabin1998 | Ejsmont-Karabin J (1998). Empirical equations for biomass calculation of planktonic rotifers. Polskie Archiwum Hydrobiologii 45, 513-522.                                                                                                                                                                                               | Values transcribed from the paper                                                                                                       |
| evanswhite2005     | Evans-White MA, Stelzer RS, Lamberti GA (2005). Taxonomic and regional patterns in benthic macroinvertebrate elemental composition in streams. Freshwater Biology 50, 1786-1799.                                                                                                                                                        | Values transcribed from the paper                                                                                                       |
| hebert2016         | Hébert M-P, Beisner BE, Maranger R (2016). A compilation of quantitative functional traits for marine and freshwater crustacean zooplankton. Ecology 97, 1081.                                                                                                                                                                          | Copyright and proprietary restrictions: none (metadata of the data paper); cite the paper                                               |
| hrycik2024         | Hrycik AR, Burlakova LE, Karatayev AY, Daniel SE, Dermott R, Tarbell M, Hinchey EK (2024). A dataset of individual wet weights of benthic macroinvertebrates. Limnology and Oceanography Letters 9, 696-715.                                                                                                                            | CC0 1.0 (Dryad data); CC BY 4.0 (article, Table 3)                                                                                      |
| jensen2006         | Jensen TC, Anderson TR, Daufresne M, Hessen DO (2006). Does excess carbon affect respiration of the rotifer Brachionus calyciflorus Pallas? Freshwater Biology 51, 2320-2333.                                                                                                                                                           | Values transcribed from the paper                                                                                                       |
| koning2024         | Koning A, Pin K, Chhuoy S, Uy Sophorn S, Everest E, Ngor PB, Chandra S, Hogan Z (2024). Weight and length measurements for 37 Mekong River fish species of the ‘Dai’ fishery. Dryad.                                                                                                                                                    | CC0 1.0 (Dryad)                                                                                                                         |
| kremer2014         | Kremer CT, Gillette JP, Rudstam LG, Brettum P, Ptacnik R (2014). A compendium of cell and natural unit biovolumes for \>1200 freshwater phytoplankton species. Ecology 95, 2984.                                                                                                                                                        | No copyright restrictions (Ecological Archives metadata); cite the paper                                                                |
| mccauley1984       | McCauley E (1984). The estimation of the abundance and biomass of zooplankton in samples. In: Downing JA, Rigler FH (eds), A Manual on Methods for the Assessment of Secondary Productivity in Fresh Waters, 2nd edn, IBP Handbook 17. Blackwell Scientific Publications, Oxford, pp. 228-265. ISBN 0-632-00616-1.                      | Values transcribed from the book chapter                                                                                                |
| mdl2000            | Menden-Deuer S, Lessard EJ (2000). Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. Limnology and Oceanography 45, 569-579.                                                                                                                                                                     | Values transcribed from the paper                                                                                                       |
| moreau2021         | Moreau A, Dupuy C, Bocher P, Farau S (2021). Morphological, calorific and nutritive characteristics of 656 freshwater invertebrates taxa. Biodiversity Data Journal 9, e70214.                                                                                                                                                          | CC0 1.0 (Dryad)                                                                                                                         |
| platt1973          | Platt T, Irwin B (1973). Caloric content of phytoplankton. Limnology and Oceanography 18, 306-310.                                                                                                                                                                                                                                      | Values transcribed from the paper                                                                                                       |
| rimet2018          | Rimet F, Druart J-C (2018). A trait database for Phytoplankton of temperate lakes. Annales de Limnologie - International Journal of Limnology 54, 18.                                                                                                                                                                                   | CC BY 4.0 (Zenodo)                                                                                                                      |
| santini2018        | Santini L, Benítez-López A, Ficetola GF, Huijbregts MAJ (2018). Length–mass allometries in amphibians. Integrative Zoology 13, 36-45.                                                                                                                                                                                                   | Publisher copyright (article: (c) 2017 ISZS, IOZ CAS and John Wiley & Sons Australia); no licence stated for the supplementary material |
| theilacker1971     | Theilacker GH, McMaster MF (1971). Mass culture of the rotifer Brachionus plicatilis and its evaluation as a food for larval anchovies. Marine Biology 10, 183-188.                                                                                                                                                                     | Values transcribed from the paper                                                                                                       |
| theilacker1984     | Theilacker GH, Kimball AS (1984). Comparative quality of rotifers and copepods as foods for larval fishes. California Cooperative Oceanic Fisheries Investigations Reports 25, 80-86.                                                                                                                                                   | Values transcribed from the paper                                                                                                       |
| trochet2014        | Trochet A, Moulherat S, Calvez O, Stevens V, Clobert J, Schmeller D (2014). A database of life-history traits of European amphibians. Biodiversity Data Journal 2, e4123.                                                                                                                                                               | CC0 (data); CC BY 4.0 (article)                                                                                                         |
| usepa1993          | US EPA (1993). Wildlife Exposure Factors Handbook. EPA/600/R-93/187. U.S. Environmental Protection Agency, Washington, DC.                                                                                                                                                                                                              | US government work; values used, no data                                                                                                |
| usepa2016          | US EPA (2016). Standard Operating Procedure for Zooplankton Analysis. LG403, Revision 07, July 2016. U.S. Environmental Protection Agency, Great Lakes National Program Office, Chicago.                                                                                                                                                | Convention of the standard operating procedure; no data used                                                                            |
| vanleeuwen2025     | van Leeuwen CHA, Csabai Z, Szloboda A, Móra A, Declerck SAJ (2025). Length-dry mass relationships of aquatic insects: geographic and taxonomic variation in a digital database. Freshwater Biology 70, e70056.                                                                                                                          | CC BY 4.0 (Figshare)                                                                                                                    |

## Status

Version 0.3.0, not on CRAN. In cross-validation, the 95 % intervals
cover 93-97 % of new populations of species in the reference and 91-97 %
of species left out. Masses without lengths are much more uncertain (the
median is divided and multiplied by a factor of about 6 to 40 for the 95
% interval of a new population) than masses from lengths (a factor of
1.2 to 4): measure lengths where possible.

## Known limitations

- Benthic invertebrates without lengths are predicted as population
  means: the masses of individual prey of Moreau et al. (2021) are
  corrected by a factor of about 6 estimated from the few families they
  share with the population means of Hrycik et al. (2024), so the
  estimates of benthic taxa measured only by Moreau depend on that
  correction (95 % interval 3.6 to 11).
- Phytoplankton counted in natural units goes to carbon and energy only
  with `cells_per_unit`; without it, it stays in biovolume and wet mass.
- Amphibian larvae are not estimated, with or without lengths.
- Fish without lengths: the mean masses of populations come from Iberian
  streams (Arranz et al. 2022) and the Mekong (Koning et al. 2024);
  orders absent from them take the mean of all fish.
- There is no dry mass or AFDM of phytoplankton, no AFDM or carbon of
  amphibians, and no energy of masses with shell.
- Detritus, periphyton and macrophytes are outside the scope of the
  package.
