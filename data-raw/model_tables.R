# Reference tables of version 0.2.0 inside the package (R/sysdata.rda)
#
# Reads the fitted models (data-raw/derived/model_*.csv, from fit_models.R),
# their data (model_typical.csv, model_lwr.csv, from model_data.R), the
# equations fitted from individual data (lwr_fitted.csv) and the source
# registry (sources.csv), and writes the internal object `ref_tables` to
# R/sysdata.rda:
#   models      one row per model: domain, kind, counting unit, currency, unit
#               and taxonomic levels
#   lwr_types   one row per length-mass model and length type: centre of the
#               log lengths, range of lengths, and the typical residual
#               standard deviation (ln) of the equations fitted from
#               individual data (for the lognormal back-transformation)
#   nodes       one row per node of each model (records identified above
#               species keep their own nodes in the fit, but those cannot be
#               matched by name and are left out here): path estimates and
#               covariances (see fit_models.R), number of records or
#               equations, sources and, in length-mass models, the range of
#               lengths of the equations below the node; plus the stage
#               effects of the invertebrate typical-mass model
#   components  variance components by model and level
#   fixed       fixed effects
#   conversion_models, conversion_nodes, conversion_components,
#               conversion_fixed
#               the models of the conversion factors between mass currencies,
#               and from them to energy, of invertebrates and fish
#               (fit_conversions.R), with the same
#               columns as the tables above (pairs with shells, and the
#               domain "other" of algae and plants, are left out: no model
#               of the package gives those currencies)
#   carbon_equations
#               the carbon-volume equations of Menden-Deuer and Lessard
#               (2000, Table 4) for phytoplankton cells
#   sources     the source registry: citation, DOI and licence (the package
#               holds only parameters fitted to the data of the sources or
#               published in them)
#   info        versions: COL release, date of the fits, lme4
# Every node also gets the identifier of its name and rank in the COL XR
# (`col_id`), from the classifications cached by taxonomy_colxr.R.
# Run after fit_models.R and fit_conversions.R:
#   Rscript data-raw/model_tables.R  (package usethis)

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
Sys.unsetenv("BMASSA_FIT_ONLY_FUNCTIONS")

rd <- function(f) read.csv(file.path(out_dir, f), stringsAsFactors = FALSE, encoding = "UTF-8")
nodes <- rd("model_nodes.csv")
components <- rd("model_components.csv")
fixed <- rd("model_fixed.csv")
dat <- read_model_data()
lf <- rd("lwr_fitted.csv")
colv <- rd("colxr_version.csv")
src <- read.csv("data-raw/sources.csv", stringsAsFactors = FALSE, encoding = "UTF-8")

# Models -------------------------------------------------------------------------------
specs <- model_specs(dat$typ, dat$lwr)
domain <- sub("_(cell|natural_unit|individual|adult|lwr)$", "", specs$model)
models <- data.frame(
  model = specs$model, domain = domain, kind = specs$type,
  count_unit = ifelse(grepl("_cell$", specs$model), "cell",
                      ifelse(grepl("_natural_unit$", specs$model), "natural_unit", "individual")),
  life_stage = ifelse(specs$model == "amphibian_adult", "adult", NA_character_),
  currency = c(phytoplankton = "biovolume", invertebrate = "dry_mass", fish = "wet_mass",
               amphibian = "wet_mass")[domain],
  unit = c(phytoplankton = "um3", invertebrate = "mg", fish = "g", amphibian = "g")[domain],
  levels = vapply(specs$model, function(m) {
    paste(intersect(levels_all, components$level[components$model == m]), collapse = ",")
  }, character(1)),
  stringsAsFactors = FALSE)
rownames(models) <- NULL

# Sources and length ranges below each node ---------------------------------------------
per_node <- list()
for (i in seq_len(nrow(specs))) {
  m <- specs$model[i]
  type <- specs$type[i]
  d <- model_rows(dat, m, type)
  lv <- strsplit(models$levels[models$model == m], ",")[[1]]
  pre <- if (type == "lwr") paste0(d$length_type, "|") else ""
  unit <- if (type == "lwr") d$equation_id else d$record_id
  lo <- if (type == "lwr") ifelse(d$kind == "mean", d$length_mean, d$length_min) else NA_real_
  hi <- if (type == "lwr") ifelse(d$kind == "mean", d$length_mean, d$length_max) else NA_real_
  keys <- c(list(root = if (type == "lwr") d$length_type else rep("all", nrow(d))),
            stats::setNames(lapply(lv, function(l) paste0(pre, node_path(d, l))), lv))
  for (l in names(keys)) {
    k <- keys[[l]]
    per_node[[paste(m, l)]] <- data.frame(
      model = m, level = l, key = sort(unique(k)),
      sources = vapply(split(d$source_id, k), function(z) paste(sort(unique(z)), collapse = ";"),
                       character(1))[sort(unique(k))],
      length_min = if (type == "lwr") vapply(split(lo, k), min, numeric(1))[sort(unique(k))] else NA_real_,
      length_max = if (type == "lwr") vapply(split(hi, k), max, numeric(1))[sort(unique(k))] else NA_real_,
      stringsAsFactors = FALSE)
  }
}
per_node <- do.call(rbind, per_node)
rownames(per_node) <- NULL
key <- ifelse(nodes$level == "root", nodes$node,
              ifelse(is.na(nodes$length_type), nodes$path, paste(nodes$length_type, nodes$path, sep = "|")))
key[nodes$level == "root" & !grepl("_lwr$", nodes$model)] <- "all"
k <- match(paste(nodes$model, nodes$level, key), paste(per_node$model, per_node$level, per_node$key))
is_stage <- nodes$level == "stage"
stopifnot(!anyNA(k[!is_stage]))
nodes$sources <- per_node$sources[k]
nodes$length_min <- per_node$length_min[k]
nodes$length_max <- per_node$length_max[k]
nodes$sources[is_stage] <- NA_character_
# Nodes of records identified above species: their names carry the source and
# the original name, and no user taxon can match them.
pseudo <- grepl(" not identified; ", nodes$node, fixed = TRUE)
nodes <- nodes[!pseudo, ]
rownames(nodes) <- NULL

# Identifiers of the COL XR for the nodes ----------------------------------------------------
# From the classifications of the answers cached by taxonomy_colxr.R: the key of
# each name at each main rank, within its kingdom. Names with two keys at the
# same rank and kingdom (homonyms) get none.
cache <- rd("colxr_matches.csv")
cache <- cache[!is.na(cache$classification) & !is.na(cache$classification_keys), ]
pairs <- do.call(rbind, lapply(seq_len(nrow(cache)), function(i) {
  nm <- strsplit(cache$classification[i], "|", fixed = TRUE)[[1]]
  ky <- strsplit(cache$classification_keys[i], "|", fixed = TRUE)[[1]]
  stopifnot(length(nm) == length(ky))
  rank <- sub(":.*$", "", nm)
  king <- sub("^KINGDOM:", "", nm[rank == "KINGDOM"])
  if (length(king) != 1L) return(NULL)
  keep <- rank %in% toupper(levels_all)
  data.frame(kingdom = king, level = tolower(rank[keep]), node = sub("^[A-Z_]+:", "", nm[keep]),
             col_id = sub("^[A-Z_]+:", "", ky[keep]), stringsAsFactors = FALSE)
}))
pairs <- unique(pairs)
dup <- duplicated(pairs[, c("kingdom", "level", "node")]) |
  duplicated(pairs[, c("kingdom", "level", "node")], fromLast = TRUE)
pairs <- pairs[!dup, ]
node_col_id <- function(nd) {
  king <- sub("\\|.*$", "", nd$path)
  k <- match(paste(king, nd$level, nd$node), paste(pairs$kingdom, pairs$level, pairs$node))
  out <- pairs$col_id[k]
  out[!nd$level %in% levels_all] <- NA_character_
  out
}
nodes$col_id <- node_col_id(nodes)

# Conversion models -------------------------------------------------------------------
# Only the pairs between the currencies of the package (wet mass, dry mass,
# ash-free dry mass and carbon, without shells) of invertebrates and fish, and
# from those masses to energy (fit_conversions.R fits no model from carbon to
# energy).
cn <- rd("conversion_nodes.csv")
cc <- rd("conversion_components.csv")
cf <- rd("conversion_fixed.csv")
mass_cur <- c("wet_mass", "dry_mass", "afdm", "carbon")
pair <- sub("^[a-z]+ ", "", unique(cn$model))
conv_keep <- unique(cn$model)[sub(" .*$", "", unique(cn$model)) %in% c("invertebrate", "fish") &
                                sub("_to_.*$", "", pair) %in% mass_cur &
                                sub("^.*_to_", "", pair) %in% c(mass_cur, "energy")]
cn <- cn[cn$model %in% conv_keep, ]
cc <- cc[cc$model %in% conv_keep, ]
cf <- cf[cf$model %in% conv_keep, ]
conversion_models <- data.frame(
  model = sort(conv_keep), domain = sub(" .*$", "", sort(conv_keep)),
  from = sub("_to_.*$", "", sub("^[a-z]+ ", "", sort(conv_keep))),
  to = sub("^.*_to_", "", sub("^[a-z]+ ", "", sort(conv_keep))),
  levels = vapply(sort(conv_keep), function(m) {
    paste(intersect(levels_all, cc$level[cc$model == m]), collapse = ",")
  }, character(1)),
  context = vapply(sort(conv_keep), function(m) {
    paste(intersect(c("habitat", "source"), cc$level[cc$model == m]), collapse = ",")
  }, character(1)),
  stringsAsFactors = FALSE)
rownames(conversion_models) <- NULL
# Sources below each node, from the factors used in the fits (listed by
# fit_conversions.R, after the screening of the energy densities).
rt <- rd("ref_taxa.csv")
cvd <- rd("conversions.csv")
used <- rd("conversion_used.csv")
shell_cur <- function(x, shell) paste0(x, ifelse(shell, "+shell", ""))
cvd$pair <- paste0(shell_cur(cvd$from, cvd$from_shell), "_to_", shell_cur(cvd$to, cvd$to_shell))
cu <- cvd[match(used$conversion_id, cvd$conversion_id), ]
stopifnot(!anyNA(cu$conversion_id))
cd <- data.frame(record_id = cu$conversion_id, source_id = cu$source_id, taxon = cu$taxon,
                 pair = cu$pair, stringsAsFactors = FALSE)
cd <- add_lineage(cd, rt)
cd$model <- paste(cd$domain, cd$pair)
stopifnot(identical(cd$model, used$model))
cd <- cd[cd$model %in% conv_keep, ]
cper <- list()
for (m in conv_keep) {
  d <- cd[cd$model == m, ]
  lv <- strsplit(conversion_models$levels[conversion_models$model == m], ",")[[1]]
  keys <- c(list(root = rep("all", nrow(d))), stats::setNames(lapply(lv, function(l) node_path(d, l)), lv))
  for (l in names(keys)) {
    k <- keys[[l]]
    cper[[paste(m, l)]] <- data.frame(
      model = m, level = l, key = sort(unique(k)),
      sources = vapply(split(d$source_id, k), function(z) paste(sort(unique(z)), collapse = ";"),
                       character(1))[sort(unique(k))], stringsAsFactors = FALSE)
  }
}
cper <- do.call(rbind, cper)
ckey <- ifelse(cn$level == "root", "all", cn$path)
k <- match(paste(cn$model, cn$level, ckey), paste(cper$model, cper$level, cper$key))
in_tree <- cn$level %in% c("root", levels_all)
stopifnot(!anyNA(k[in_tree]))
cn$sources <- ifelse(in_tree, cper$sources[k], NA_character_)
cn <- cn[!grepl(" not identified; ", cn$node, fixed = TRUE), ]
# Same columns as the nodes of the mass models, so that the same code draws
# their errors.
for (col in setdiff(names(nodes), names(cn))) cn[[col]] <- NA
cn$col_id <- node_col_id(cn)
conversion_nodes <- cn[, names(nodes)]
rownames(conversion_nodes) <- NULL

# Carbon of phytoplankton cells (Menden-Deuer and Lessard 2000, Table 4): the
# equations recommended by the authors for mixed protist plankton.
mdl <- read.csv("data-raw/mdl2000_carbon_volume.csv", stringsAsFactors = FALSE, encoding = "UTF-8")
carbon_equations <- mdl[match(c("mdl2000_diatoms", "mdl2000_protist"), mdl$equation_id),
                        c("equation_id", "group", "log10_a", "log10_a_ci95", "b", "b_ci95", "r2", "n",
                          "source_detail")]
carbon_equations$source_id <- "mdl2000"
rownames(carbon_equations) <- NULL
stopifnot(nrow(carbon_equations) == 2L, !anyNA(carbon_equations$log10_a))

# Length types: centre, range and residual standard deviation ----------------------------
lt <- fixed[grepl("^i_", fixed$term), c("model", "term", "centre")]
lt$length_type <- sub("^i_", "", lt$term)
used <- dat$lwr[dat$lwr$kind == "equation" & is.finite(dat$lwr$sigma_ln), ]
lwr_types <- data.frame(model = lt$model, length_type = lt$length_type, centre = lt$centre,
                        stringsAsFactors = FALSE)
root <- nodes[nodes$level == "root" & grepl("_lwr$", nodes$model), ]
kr <- match(paste(lwr_types$model, lwr_types$length_type), paste(root$model, root$node))
lwr_types$length_min <- root$length_min[kr]
lwr_types$length_max <- root$length_max[kr]
lwr_types$sigma_ind <- NA_real_
lwr_types$n_sigma <- 0L
for (i in seq_len(nrow(lwr_types))) {
  s <- used$sigma_ln[paste0(used$domain, "_lwr") == lwr_types$model[i] &
                       used$length_type == lwr_types$length_type[i]]
  lwr_types$n_sigma[i] <- length(s)
  if (length(s)) lwr_types$sigma_ind[i] <- stats::median(s)
}
# A length type without equations fitted from individual data takes the
# median of its domain.
for (i in which(lwr_types$n_sigma == 0L)) {
  s <- used$sigma_ln[paste0(used$domain, "_lwr") == lwr_types$model[i]]
  lwr_types$sigma_ind[i] <- stats::median(s)
}
stopifnot(all(is.finite(lwr_types$sigma_ind)))

# Sources ---------------------------------------------------------------------------------
# Conventions of R/conversion.R: cen2015 (density of 1 g/cm3 for the wet mass of
# phytoplankton), usepa2016 and dumont1975 (dry to wet mass of rotifers),
# usepa1993 (dry to wet mass and energy of amphibians).
used_src <- sort(unique(c(unlist(strsplit(nodes$sources[!is.na(nodes$sources)], ";")),
                          unlist(strsplit(conversion_nodes$sources[!is.na(conversion_nodes$sources)], ";")),
                          carbon_equations$source_id, "cen2015", "usepa2016", "dumont1975",
                          "usepa1993")))
sources <- src[src$source_id %in% used_src, c("source_id", "citation", "doi", "licence")]
rownames(sources) <- NULL
stopifnot(setequal(sources$source_id, used_src))

fit_log <- readLines(file.path(out_dir, "model_fit_log.md"), encoding = "UTF-8")
fitted_on <- regmatches(fit_log, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", fit_log))[1]
conv_log <- readLines(file.path(out_dir, "conversion_fit_log.md"), encoding = "UTF-8")
conversions_fitted_on <- regmatches(conv_log, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", conv_log))[1]
info <- list(
  data_version = "0.2.0",
  col_release = colv$dataset_alias[1], col_dataset_key = as.character(colv$clb_dataset_key[1]),
  col_index_built = colv$index_built[1], col_checked_on = colv$checked_on[1],
  fitted_on = fitted_on, conversions_fitted_on = conversions_fitted_on,
  lme4 = as.character(utils::packageVersion("lme4")))
stopifnot(!anyNA(unlist(info)))

ref_tables <- list(models = models, lwr_types = lwr_types, nodes = nodes,
                   components = components, fixed = fixed,
                   conversion_models = conversion_models, conversion_nodes = conversion_nodes,
                   conversion_components = cc, conversion_fixed = cf,
                   carbon_equations = carbon_equations, sources = sources, info = info)
usethis::use_data(ref_tables, internal = TRUE, overwrite = TRUE, compress = "xz")
cat(sprintf("nodes with col_id: %d of %d (models), %d of %d (conversions)\n",
            sum(!is.na(nodes$col_id)), sum(nodes$level %in% levels_all),
            sum(!is.na(conversion_nodes$col_id)), sum(conversion_nodes$level %in% levels_all)))
cat(sprintf("models %d, nodes %d (record nodes left out: %d), components %d, conversion models %d (nodes %d), sources %d; %s\n",
            nrow(models), nrow(nodes), sum(pseudo), nrow(components), nrow(conversion_models),
            nrow(conversion_nodes), nrow(sources),
            format(structure(file.size("R/sysdata.rda"), class = "object_size"), units = "Kb")))
