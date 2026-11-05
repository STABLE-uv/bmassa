# data-raw

Scripts and hand transcriptions that build the data sets of bmassa. This
folder is excluded from the package tarball (`.Rbuildignore`).

The package holds only the parameters of the fitted models, by taxon, with
the sources to cite (`R/sysdata.rda`, `get_reference()`, `cite_sources()`).
The public repository holds the scripts of this folder and this file; the
transcriptions, the source registry and the tables in `derived/` are kept in
the development repository of the package.

## Scripts

| Script | Builds | Run from the package root |
|---|---|---|
| `references.R` | The reference tables of version 0.1.0, no longer in the package but read by `taxonomy_colxr.R` and `model_data.R`: `derived/taxonomy_0.1.0.csv`, `derived/typical_mass_0.1.0.csv` and `derived/length_mass_0.1.0.csv` (numbers with 17 significant digits, read back exactly), and the build log `references_log.md` | `Rscript data-raw/references.R` (package readxl; helpers of `names.R`) |
| `examples.R` | `data/example_survey.rda` (constructed data, not field data; `col_id` of the taxa from the COL26.6 XR) | `Rscript data-raw/examples.R` (package usethis) |
| `fit_lwr.R` | `derived/lwr_fitted.csv` (length–mass equations fitted from individual measurements), `derived/conversions_cadses.csv` (dry mass to ash-free dry mass), `derived/fit_lwr_log.md` | `Rscript data-raw/fit_lwr.R` (package readxl; originals under `BMASSA_DATABASES`) |
| `typical_records.R` | `derived/typical_records.csv` (mean masses or biovolumes of individuals, cells or colonies, with the mean length where given), `derived/typical_records_log.md` | `Rscript data-raw/typical_records.R` (package readxl) |
| `published_lwr.R` | `derived/lwr_published.csv` (length–mass equations transcribed from articles), `derived/lwr_published_log.md` (checks against the typical records) | `Rscript data-raw/published_lwr.R`, after `typical_records.R` (package readxl) |
| `conversions.R` | `derived/conversions.csv` (factors between wet, dry and ash-free dry mass, carbon and energy, with or without shell; `sample` links the factors of the same material), `derived/conversions_log.md` | `Rscript data-raw/conversions.R`, after `fit_lwr.R` (package readxl; reads also `energy_literature.csv` and `cummins1971_*.csv`, and the data bank of Brey et al. 2010 with `brey.R`) |
| `taxonomy_colxr.R` | `derived/taxa_colxr.csv` (each name of the reference tables matched to the Catalogue of Life eXtended Release, with its classification and flags), `derived/colxr_matches.csv` (the answers of the GBIF API, also used as cache), `derived/colxr_version.csv` (version of the COL XR of those answers), `derived/taxa_colxr_log.md` | `Rscript data-raw/taxonomy_colxr.R`, after the other scripts (packages curl, jsonlite and readxl; internet only for names not yet in the cache, or for all with `BMASSA_COLXR_REFRESH = "true"`) |
| `hierarchy.R` | `derived/ref_taxa.csv` (one row per source and name of the reference tables: COL taxon, lowest rank, domain and the seven main ranks, with where each rank comes from), `derived/ref_taxa_log.md` | `Rscript data-raw/hierarchy.R`, after `taxonomy_colxr.R` |
| `fit_conversions.R` | `derived/conversion_nodes.csv`, `derived/conversion_components.csv`, `derived/conversion_fixed.csv` (hierarchical models of the conversion factors, one per pair of currencies and domain, including energy per dry mass, ash-free dry mass and wet mass), `derived/conversion_used.csv` (the factors of each model), `derived/conversion_screening.csv` (the energy densities left out by the screening, and why), `derived/conversion_fit_log.md` | `Rscript data-raw/fit_conversions.R`, after `hierarchy.R` (packages lme4 and Matrix) |
| `model_data.R` | `derived/model_typical.csv` (typical masses and biovolumes by model, in the native currency, with their lineage), `derived/model_lwr.csv` (length–mass equations and amphibian population means by model), `derived/model_data_log.md` | `Rscript data-raw/model_data.R`, after `fit_conversions.R` |
| `fit_models.R` | `derived/model_nodes.csv`, `derived/model_components.csv`, `derived/model_fixed.csv` (the fitted hierarchical models: estimates and covariances of each node of the reference, variance components, fixed effects), `derived/model_fit_log.md` | `Rscript data-raw/fit_models.R`, after `model_data.R` (packages lme4 and Matrix; about 1 minute) |
| `model_tables.R` | `R/sysdata.rda`: the internal object `ref_tables` with the fitted models (nodes, components, fixed effects, ranges of lengths and sources by node, typical residual standard deviation of the length-mass equations fitted from individual data), the conversion models between mass currencies and to energy of invertebrates and fish, the carbon-volume equations of Menden-Deuer and Lessard (2000), the identifier of each node in the COL XR (`col_id`, from the cache of `taxonomy_colxr.R`) and the source registry, returned by `get_reference()` | `Rscript data-raw/model_tables.R`, after `fit_models.R` and `fit_conversions.R` (package usethis) |
| `check_permissions.R` | Report of the permissions of the sources and of the MD5 checks; with `--release`, an error if an included source was refused (see below) | `Rscript data-raw/check_permissions.R --release` |
| `sources.R`, `brey.R`, `names.R` | Helpers sourced by the scripts: the source registry and the original files; the reader of the data bank of Brey et al. (2010); the normalisation of taxon names of version 0.1.0 | — |

The tables in `derived/` are intermediate: they keep one row per published
or fitted value, with its source, original reference, permission and the
reason when a row is not to be used (`use_in_model`, `flag`), and are the
input of the taxonomic model of version 0.2.0. Run the scripts in a UTF-8
locale (for instance `LANG=C.UTF-8`).

## Source registry (version 0.2.0)

`sources.csv` lists every source, used or not: citation, DOI, data URL,
licence, `status` (`included`, `planned`, `excluded`), `permission`
(`not_needed`, `pending`, `unverified`, `granted`, `refused`,
`not_requested`; see `sources.R`), the original files and their MD5
checksums. The scripts for version 0.2.0 read the originals from the folder
given by `BMASSA_DATABASES` (default `data-raw/databases`), with one
subfolder per source, for instance:

```r
Sys.setenv(BMASSA_DATABASES = "E:/biomassa/databases")
```

The records of the sources are used to fit the models and the sources are
cited (`cite_sources()`), whatever their permission, unless it was refused
(decision of 2026-11-05): the package and its public repository hold only
the parameters of the models. Sources whose permission is `pending` or
`unverified` are left out of the tables with
`BMASSA_INCLUDE_PENDING = "false"`. Before a release,
`Rscript data-raw/check_permissions.R --release` must pass.

`references.R` reads the original files from the folder given by the
environment variable `BMASSA_SOURCES` (default `data-raw/sources`), checks
their MD5 checksums and stops if any file differs. The original files are
not kept in the repository. Run it in a UTF-8 locale (for instance
`LANG=C.UTF-8`).

## Original files

| File | Source | Where to get it | Licence | MD5 |
|---|---|---|---|---|
| `taxa_table_030614.csv` | Kremer et al. (2014) | Ecological Archives E095-257, <http://esapubs.org/archive/ecol/E095/257/> | "Copyright restrictions: None"; the paper must be cited (metadata page) | `147cbba2ba0fd96d33dd0eb9ad32a0e4` |
| `bvd_species_ag_030614.csv` | Kremer et al. (2014) | as above | as above | `b3638fd46f5d1518b8206bf96ea6bbaf` |
| `bvd_genus_ag_030614.csv` | Kremer et al. (2014) | as above | as above | `3a6f1dc5c9be57c1d4db5731205cba50` |
| `bvd_species_raw_030614.csv` | Kremer et al. (2014) | as above | as above | `1e3e815263e9d5a65dc5f91506ca30a3` |
| `bvd_genus_raw_030614.csv` | Kremer et al. (2014) | as above | as above | `9d8889b1d03b470bc60c812f7177f758` |
| `Database_invertebrates.txt` | Moreau et al. (2021) | Dryad, <https://doi.org/10.5061/dryad.j3tx95xfg> | CC0 (Dryad) | `b14b453a45a39b91552f03485c743f03` |
| `Database_references.txt` | Moreau et al. (2021) | as above | CC0 (Dryad) | `3c52052294fbeeebfb5a2c1b88bb2680` |
| `vanLeeuwen2025_TableS2.xlsx` | van Leeuwen et al. (2025), Table S2 | Figshare, <https://doi.org/10.6084/m9.figshare.29085596> | CC BY 4.0 | `eeb61729f6fb47065b292871c9c8b442` |

The Kremer et al. biovolumes are in log10(um3); `references.R` checks that
the published taxon means are the means of the log10 raw values before
using the raw values.

## Checks and rules applied by `references.R`

- Kremer et al.: species rows (all records of each accepted name, including
  those identified to variety or form) and genus rows (all records of the
  genus, flagged `pooled`) reproduce the numbers of measurements and the
  log10 means of the published species and genus products; the script stops
  otherwise. Per-cell values: in genera where at least 5% of the records
  with an explicit cell count report several cells per natural unit, only
  those records are used.
- van Leeuwen et al.: equations recomputed from the original intercept and
  slope when both are given; equations whose prediction differs by more
  than a factor of 20 from the median of the other equations of their family
  are left out (four, listed in `references_log.md`).
- Dumont et al.: continent from the origin codes of the paper.
- Every decision and every row left out is written to `references_log.md`.

## Hand transcriptions

| File | Content | Source |
|---|---|---|
| `dumont1975_length_mass.csv` | 39 length-dry mass equations of cladocerans, copepods and rotifers, with the printed coefficients and units, the life stage, the length range and the reason when an equation is left out | Dumont et al. (1975), Tables 2 and 6 and Fig. 7 (page images) |
| `dumont1975_typical.csv` | 112 mean dry masses of animals weighed in groups | Dumont et al. (1975), Tables 3, 7 and 8 |
| `mdl2000_carbon_volume.csv` | 11 carbon-volume regressions of protist plankton | Menden-Deuer and Lessard (2000), Table 4 |
| `lineage_manual.csv` | Lineages of the Dumont et al. taxa not covered by the other sources, and documented corrections of the classifications | WoRMS (<https://www.marinespecies.org>), consulted on 2026-11-01 |
| `hrycik2024_table3.csv` | 145 rows of mean DW/WW ratios of freshwater benthic invertebrates (shell-free, SFDW/SFWW and with shell), with the cited authors | Hrycik et al. (2024), Table 3, pp. 704-708 (page images) |
| `santini2018_table3.csv` | 13 snout-vent length–mass equations of frogs and salamanders (log10 coefficients, standard errors, adjusted R², number of species, RMSE, MAPE) | Santini et al. (2018), Table 3, p. 41 (page image) |
| `energy_literature.csv` | Energy densities and mass ratios (% water, % dry mass, % ash, % AFDM, % carbon) of 24 taxa from articles: value and unit as printed, error, n, stage and parts of the animals; rows kept but not used say why in `flag` | Yúfera et al. (1997), p. 85 (kept, not used: in the data bank of Brey et al. 2010); Mitchell (2001), Table 3; Ciancio and Pascual (2006), Table 1; Nelson and Kruuk (1997), pp. 76-77 (kept, not used: animals prepared as prey); Theilacker and McMaster (1971), p. 186; Theilacker and Kimball (1984), Tables 2-3 (the value per ash-free dry mass, computed with the ash of the 1971 study, kept but not used); Jensen et al. (2006), p. 2328; Méthot et al. (2012), Table 2, p. 754, and Cianciara (1980), Tables III and VI, pp. 384 and 390 (kept, not used: animals preserved in ethanol and in formalin); all checked against the page images |
| `cummins1971_table3.csv` | Energy per gram of dry mass, ash-free dry mass and wet mass (with the error term and its code), % ash, % water, number of samples and codes (season, stage, method, source, general notes) of the 279 printed rows of aquatic macroconsumers, as printed; added columns: `taxon` (abbreviations, repeat marks and misprints resolved), `habitat` (only where it differs from that of the source), `curation_flag` (reason to leave the row out), `transcription_note` (doubtful digits, second error terms) and `edit_note` (what was edited and why) | Cummins and Wuycheck (1971), Table 3, pp. 91-119 (PDF pp. 92-120), transcribed from the page images at 200-600 dpi in five page ranges; the values were checked against each other (energy per dry and ash-free dry mass and % ash), against the OCR layer of the PDF and, for five pages, against an independent second transcription |
| `cummins1971_sources.csv` | The 36 original sources cited by Cummins and Wuycheck (1971) for those rows: number, reference, habitat, whether their energy densities (`brey_energy`) or their water or ash contents (`brey_mass`) are in the data bank of Brey et al. (2010), and a note | Cummins and Wuycheck (1971), notations V and VI; reference sheet and records of Brey et al. (2010) |
| `evanswhite2005_appendix.csv` | Carbon as % of dry mass of 46 taxa of stream invertebrates (30 insects, 5 crustaceans, 11 molluscs without shell): mean, SD and range among streams, number of streams and region, as printed; added columns: `taxon` (two misprints corrected, `edit_note`), `group` and the phylum, class and family used as hints for the Catalogue of Life | Evans-White et al. (2005), Appendix, pp. 1798-1799: read from the text layer of the publisher's PDF and checked against the page images; the rows of two streams have their mean at the midpoint of the range and their SD equal to the range over the square root of 2 (within rounding), and the means weighted by the number of streams reproduce those of Table 1 for insects, molluscs and crustaceans (48.6, 41.8 and 34.8 % against 48.3, 42.2 and 34.8 %) |
| `taxa_manual.csv` | Corrections to the matching of reference names with the Catalogue of Life (COL XR): another name or rank to query, an answer accepted although flagged, or a name left out, each with its reason | Checked by hand against the answers of the GBIF API (2026-11-02) |
| `culver1985_table2.csv` | Lengths and dry weights of the stages of the copepods of Lake Erie (NI-NVI, CI-CV, adult females and males; eggs), as printed, with the name looked up (`taxon`, `taxon_note`) and the phylum and class used as hints | Culver et al. (1985), Table 2, p. 1383 (page images at 400 dpi) |
| `culver1985_table3.csv` | 30 regressions W = a L^b (ug dry mass, mm) of copepods (to the base of the caudal rami) and cladocerans (Daphnia by standard and total length), with n, the confidence limit of b, var L, RMS and r² | Culver et al. (1985), Table 3, p. 1384 |
| `ejsmontkarabin1998_table2.csv` | 251 species of rotifers: state (contracted or living), formula of the body volume (mm3 = mg wet weight) from the body length or width, appendices (% of the weight), mucus sheaths, range of lengths, standard length (with what is measured) and standard wet weight | Ejsmont-Karabin (1998), Table II, pp. 516-521 (scans at 300 dpi, read at 400 dpi) |
| `ejsmontkarabin1998_table1.csv` | The geometric formulas of the body volume of 20 genera (Ruttner-Kolisko 1977) | Ejsmont-Karabin (1998), Table I, p. 515 |
| `mccauley1984_table7_2.csv` | 88 regressions ln W = ln a + b ln L (ug dry mass, mm) of cladocerans and copepods, with the length measured (code and definition), preservation, n, F or r, range and the original study | McCauley (1984), Table 7.2, book pp. 234-237 (PDF pp. 262-265) |
| `bottrell1976_length_mass.csv` | 65 regressions ln W = ln a + b ln L (ug dry mass, mm) of rotifers (Table VII), cladocerans (Tables X and XI) and copepods (Table XII), with the stages, locality, contributor, n and fit statistics | Bottrell et al. (1976), pp. 432-443 (scan read at 400 dpi) |
| `bottrell1976_typical.csv` | 75 individual weights of Protozoa and rotifers: dry weights compiled from several studies (Table VIII) and wet weights from calculated volumes, as the range of the means of IBP sites (Table IX; printed in mg, read as ug) | Bottrell et al. (1976), Tables VIII and IX, pp. 433-434 |

The zooplankton transcriptions of phase 4k were read twice from separate
crops of the page images (Culver et al., Ejsmont-Karabin, Bottrell et al.,
McCauley) and checked as follows. Culver et al.: each regression gives the
weights of Table 2 at their lengths (roughly) and its range covers them; a
and b of three regressions are also quoted in the text. Ejsmont-Karabin: the
standard weight of 248 of the 251 species is its formula at its standard
length (with the appendices) within 5 %; of the other three, Asplanchnopus
multiceps (printed weight a tenth of its formula) is left out of the
typical masses, and two Trichocerca keep their printed weights (those of a
slightly shorter length); Lindia truncata is printed with the body width
(b^3) but computed with its length. McCauley: the 48 regressions taken from
Bottrell et al. (1976) and Burns (1969) coincide in ln a and b with the
independent transcription of Bottrell et al.; one range differs
(Bosmina longispina, 0.38-0.50 in McCauley, 0.30-0.50 in the page of
Bottrell et al., which is used). Bottrell et al.: the regressions of Table
VII are for lengths in um, not mm as captioned (Dumont et al. 1975; left
out), and the wet weights of Table IX are in ug, not mg. The original file
of Hébert et al. (2016) is read from `BMASSA_DATABASES`
(`Hebert2016/zooplankton_traits.csv`, `references.csv`).

The transcription of Dumont et al. was checked against the representative
values of their Tables 4 and 6 (`tests/testthat/fixtures/dumont1975_representative.csv`,
used by the tests): 30 of the 33 equations with representative values
reproduce them within 5.3%; the other three are left out (see
`exclude_reason`).

The transcription of Hrycik et al. was compared row by row with the page
images; 17 of its 23 values taken only from Brey et al. (2010) are
reproduced by the mean of the matching records of the data bank
(`derived/conversions_log.md`). The equations of Santini et al. predict the
adult masses of their own species means with a median ratio of 1.01 and those
of Trochet et al. (2014) with 1.01 (`derived/lwr_published_log.md`).

How the columns of Hrycik et al. Table 3 are read for molluscs (their wet
masses include the shell), and which rows are not used (summaries of Brey et
al. records, values repeated for many taxa, implausible values), is
explained at the start of `conversions.R` and listed in its log.

## Taxonomy (version 0.2.0)

`taxonomy_colxr.R` matches every name of the reference tables (the package
data of 0.1.0 and the tables in `derived/`) with the Catalogue of Life
eXtended Release through the GBIF name-matching API v2
(`checklistKey = 7ddf754f-d193-4cc9-b351-99906754a03b`). The classification
given by each source is sent as hints to tell homonyms apart, the names are
cleaned first (diacritics, cf., sp. 1, descriptions of morphotypes), and the
answers are checked: Animalia for animals and not for algae, the phylum of
the source, the same name or a close spelling variant, and the rank the name
should have. Names without an acceptable answer are sent again with their
rank (taxonRank) and with fewer hints, and species that are not in the COL
go to their genus. `taxa_manual.csv` holds the corrections made by hand. The
answers are kept in `derived/colxr_matches.csv`, so the script runs again
without internet; `BMASSA_COLXR_REFRESH = "true"` sends all the queries
again (the COL XR changes every month: the log records the version).

## Hierarchical models (version 0.2.0)

`hierarchy.R` gives every reference name the seven main ranks of the COL
classification; ranks the COL leaves empty are filled from related COL taxa,
then (orders and classes) from the nearest COL rank above them (a superorder
for a missing order, so that sister families share their node), then from the
classification of the source and, failing that, with placeholders, so that
every name has one parent per rank. Copepods, for which the COL gives only
the superorders Gymnoplea and Podoplea, take the orders Calanoida, Cyclopoida
and Harpacticoida (WoRMS), and the orders of free-living flatworms the class
Turbellaria. `model_data.R` joins those lineages to the typical masses and
the length–mass equations, in the native currency of each domain:
biovolume (um3) per cell or per natural unit for phytoplankton, dry mass (mg;
shell-free for molluscs) for invertebrates, wet mass (g) for fish and
amphibians. Records identified above species get nodes of their own below
their lowest rank. Wet masses of Hrycik et al. (2024) and the dry masses with
shell of the CADSES molluscs are converted with the models of
`fit_conversions.R`. Equations without a length range get the range of their
relatives (genus, family or order).

`fit_models.R` fits one model per domain and kind of data with lme4 (REML):
typical masses with nested random intercepts by taxonomic level, and
length–mass relations from three points per equation within its length range,
with random intercepts and slopes by taxonomic level within length type and
random effects of equation (populations) and study. For each node of the
reference it writes the estimate of the path from the root and its
conditional covariance, computed from the mixed model equations so that the
effects of a node, of its ancestors and of the fixed effects are not treated as
independent; the package predicts from these tables with base R. The
cross-validation of the intervals is in `dev/validacion/validacion_modelos.md`.

The order of the scripts for version 0.2.0 is `fit_lwr.R`,
`typical_records.R`, `published_lwr.R`, `conversions.R`, `taxonomy_colxr.R`,
`hierarchy.R`, `fit_conversions.R`, `model_data.R`, `fit_models.R` and
`model_tables.R`; `examples.R` builds the example data sets.

## References

Brey T, Müller-Wiegmann C, Zittier ZMC, Hagen W (2010). Body composition in
aquatic organisms - a global data bank of relationships between mass,
elemental composition and energy content. *Journal of Sea Research* 64,
334-340. <https://doi.org/10.1016/j.seares.2010.05.002>

Coughlan NE, Cunningham EM, Cuthbert RN, et al. (2021). Biometric conversion
factors as a unifying platform for comparative assessment of invasive
freshwater bivalves. *Journal of Applied Ecology* 58, 1945-1956.
<https://doi.org/10.1111/1365-2664.13941>

Dumont HJ, Van de Velde I, Dumont S (1975). The dry weight estimate of
biomass in a selection of Cladocera, Copepoda and Rotifera from the
plankton, periphyton and benthos of continental waters. *Oecologia* 19,
75-97. <https://doi.org/10.1007/BF00377592>

Hrycik AR, Burlakova LE, Karatayev AY, Daniel SE, Dermott R, Tarbell M,
Hinchey EK (2024). A dataset of individual wet weights of benthic
macroinvertebrates. *Limnology and Oceanography Letters* 9, 696-715.
<https://doi.org/10.1002/lol2.10428>

Kremer CT, Gillette JP, Rudstam LG, Brettum P, Ptacnik R (2014). A
compendium of cell and natural unit biovolumes for >1200 freshwater
phytoplankton species. *Ecology* 95, 2984. <https://doi.org/10.1890/14-0603.1>

Menden-Deuer S, Lessard EJ (2000). Carbon to volume relationships for
dinoflagellates, diatoms, and other protist plankton. *Limnology and
Oceanography* 45, 569-579. <https://doi.org/10.4319/lo.2000.45.3.0569>

Moreau A, Dupuy C, Bocher P, Farau S (2021). Morphological, calorific and
nutritive characteristics of 656 freshwater invertebrates taxa.
*Biodiversity Data Journal* 9, e70214. <https://doi.org/10.3897/BDJ.9.e70214>

Rimet F, Druart J-C (2018). A trait database for phytoplankton of temperate
lakes. *Annales de Limnologie - International Journal of Limnology* 54, 18.
<https://doi.org/10.1051/limn/2018009>

Santini L, Benítez-López A, Ficetola GF, Huijbregts MAJ (2018). Length–mass
allometries in amphibians. *Integrative Zoology* 13, 36-45.
<https://doi.org/10.1111/1749-4877.12268>

Trochet A, Moulherat S, Calvez O, Stevens V, Clobert J, Schmeller D (2014). A
database of life-history traits of European amphibians. *Biodiversity Data
Journal* 2, e4123. <https://doi.org/10.3897/BDJ.2.e4123>

van Leeuwen CHA, Csabai Z, Szloboda A, Móra A, Declerck SAJ (2025).
Length-dry mass relationships of aquatic insects: geographic and taxonomic
variation in a digital database. *Freshwater Biology* 70, e70056.
<https://doi.org/10.1111/fwb.70056>
