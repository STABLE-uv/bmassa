# The whole workflow of a user of bmassa (runs with Rscript --vanilla in an
# empty directory). It is the example of the README, of the vignette and of
# tests/testthat/test-workflow.R.
library(bmassa)

# 1. Two tables: the taxon units (a taxon, or a taxon in one life stage and
#    sex, with its broad group, its classification and its counting unit) and
#    their abundances by site (and date, if any). A third, optional table holds
#    the lengths: usually a mean length per taxon unit (a row without site
#    applies to all the samples), or the mean of a site or a sample, or the
#    lengths of measured individuals (one row per individual or length class,
#    with `n`).
tx <- example_survey$taxa
ab <- example_survey$abundance
ln <- example_survey$lengths
head(tx[, c("taxon_id", "broad_group", "scientific_name", "life_stage", "count_unit")])
head(ab)
ln

# 0. (With internet) The classification and the col_id of your own list of
#    names come from standardize_taxa(), which queries the Catalogue of Life
#    (COL XR):
#    my_taxa <- data.frame(taxon_id = c("a", "b"), broad_group = "macroinvertebrate",
#                          scientific_name = c("Baetis rhodani", "Ephemerella ignita"))
#    standardize_taxa(my_taxa)

# 2. Before estimating: the taxonomic level from which each taxon will come and
#    how wide its interval will be (an intermediate result that conditions the
#    estimate).
cov <- check_coverage(tx, ln)
cov
as.data.frame(cov)[, c("taxon_id", "method", "match_level", "match_name", "factor", "status")]

# 3. Mean individual mass of a new population of each taxon (mean, median and
#    95 % prediction interval) and biomass of each row, all in one currency
#    and unit (wet mass in mg by default; `currency` chooses another). With
#    lengths, the back-transformation has to be chosen: "lognormal" gives the
#    mean of the individuals of each length; "none", their median.
est <- estimate_biomass(ab, tx, ln, bias_correction = "lognormal")
est
e <- as.data.frame(est)
e[, c("site", "taxon_id", "method", "length_basis", "match_level", "ind_median", "ind_lower",
      "ind_upper", "unit", "status")]
e[e$status != "estimated", c("site", "taxon_id", "life_stage", "reason")]

# 4. Another currency afterwards, with the uncertainty of the factors: carbon
#    in µg for all the groups (phytoplankton per cell with Menden-Deuer and
#    Lessard 2000; amphibians without conversion). The same comes from
#    estimate_biomass(ab, tx, ln, currency = "carbon", unit = "ug", ...).
est_c <- convert_currency(est, to = "carbon", unit = "ug")
ec <- as.data.frame(est_c)
ec[ec$site == "L1", c("taxon_id", "converted_from", "ind_median", "ind_lower", "ind_upper",
                      "unit", "status")]

# 5. Totals by site and group with shared draws (two species estimated from
#    the same family have correlated errors); different units of abundance
#    (per L, per m2) are never added together.
tot <- summarize_biomass(est_c, by = c("site", "broad_group"), na_rm = TRUE, seed = 1)
tot

# 6. Joint draws for other packages (for instance the energy fluxes of a food
#    web in bwebs): mean individual mass and biomass of every row in every
#    draw, with the same correlations as the totals. With population_by, each
#    habitat has its own populations. To add up groups sampled per litre and
#    per m2, all the abundances on the same basis before estimating:
#    convert_abundance() with the size of each sample in `per` and the depth,
#    volume or area of each site (per m2 here; per tank with to = "site").
sites <- data.frame(site = c("L1", "L2"), depth = c(1.5, 0.8), depth_unit = "m")
ab_h <- convert_abundance(ab, to = "m2", sites = sites)
unique(ab_h[, c("per_original", "per", "per_factor", "per_assumption")])
ab_h$habitat <- ifelse(ab_h$site == "L1", "lake", "spring")
est_w <- estimate_biomass(ab_h, tx, ln, bias_correction = "lognormal")
sim <- simulate(est_w, nsim = 1000, seed = 1, population_by = "habitat")
sim
d <- as.data.frame(sim, draws = 1:2)
d[d$draw == 1, c("site", "taxon_id", "habitat", "ind_mass", "biomass", "biomass_unit")]

# 7. Size structure of each sample and group: mean and median wet mass of the
#    individuals and shares of individuals and biomass by mass class; the
#    lengths of the macroinvertebrates (the Gammarus of L1 were measured).
sz <- summarize_size(est_w, by = c("site", "broad_group"), breaks = c(0.01, 1, 100), seed = 1)
sz
sz$classes[sz$classes$site == "L1", c("broad_group", "class", "share_individuals",
                                      "share_biomass", "approximation")]
sl <- summarize_size(est_w, by = c("site", "broad_group"), size = "length",
                     breaks = seq(4, 20, by = 4), seed = 1)
as.data.frame(sl)[, c("site", "broad_group", "size_type", "mean_size", "median_size",
                      "share_excluded", "share_measured")]

# 8. The sources to cite.
cite_sources(est_c)[, c("source_id", "citation", "licence")]

# 9. Export (never to the working directory by default).
out <- tempfile(fileext = ".csv")
utils::write.csv(ec, out, row.names = FALSE)
