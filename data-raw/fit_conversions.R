# Hierarchical models of the conversion factors (version 0.2.0)
#
# Reads data-raw/derived/conversions.csv (conversions.R) and ref_taxa.csv
# (hierarchy.R) and writes, in data-raw/derived/:
#   conversion_nodes.csv       one row per model and node, as model_nodes.csv
#   conversion_components.csv  variance components by model and level
#   conversion_fixed.csv       fixed effects
#   conversion_used.csv        the factors used by each fitted model
#   conversion_screening.csv   the energy densities left out by the screening, and why
#   conversion_fit_log.md      what was fitted, and diagnostics
# A model per pair of currencies and domain (and, from version 0.6.0, for the
# tadpoles apart: domain "amphibian_larva"): y = ln of the factor (a mass
# ratio, or J/ug for energy); a fixed intercept; random intercepts of each nested
# taxonomic level that splits the data and of habitat (freshwater, marine,
# transitional) when two or more habitats have 10 or more factors. No effect
# of source (decision of 2026-11-06; see fit_one()). The functions are those
# of fit_models.R. Pairs with fewer than 30 factors or fewer than 5 species nodes
# are not fitted. model_data.R uses these tables to convert the wet masses of
# Hrycik et al. (2024) and the dry masses with shell of the CADSES molluscs; the
# package uses them in convert_currency().
#
# Energy densities (decision of 2026-11-03): only
# per dry mass, ash-free dry mass (AFDM) and wet mass without shell or
# skeleton (no model from carbon: carbon may be total rather than organic,
# and 80 of 232 values per carbon of invertebrates and fish are outside 35-55
# J/mg). Screening, in this order:
#   1. per AFDM: outside 15-45 J/mg. The combustion energy of carbohydrate,
#      protein and lipid is 17.5, 23.9 and 39.5 kJ/g (Gnaiger and Bitterlich
#      1984, Table 1), so organic matter lies within 17.5-39.5 J/mg AFDM; the
#      margin allows for measurement error.
#   2. same material (`sample`): the energy per dry mass cannot exceed that per
#      AFDM, nor that per wet mass the one per dry mass (5 % tolerance); both
#      values of an inconsistent pair are left out.
#   3. robust z, by model, against the values of the class of the taxon (or of
#      its phylum, or of the whole model, when the class or the phylum has
#      fewer than 10 values): |ln value - median| / max(1.4826 MAD, 0.1) > 4.
#   4. standardized residuals of the fitted model (conditional residual over
#      the residual SD) beyond 4: left out and the model fitted again, once.
# The mass ratio models are not screened (their factors were checked in
# conversions.R).
#
# Run from the package root, in a UTF-8 locale:  Rscript data-raw/fit_conversions.R

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
Sys.unsetenv("BMASSA_FIT_ONLY_FUNCTIONS")

rt <- read.csv(file.path(out_dir, "ref_taxa.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
cv <- read.csv(file.path(out_dir, "conversions.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
log_lines <- c("# Fits of the conversion factors", "",
               paste0("Built by `data-raw/fit_conversions.R` on ", format(Sys.Date()), " (", R.version.string,
                      ", lme4 ", as.character(utils::packageVersion("lme4")), ", REML)."), "")
say <- function(...) log_lines <<- c(log_lines, sprintf(...))

currency <- function(x, shell) paste0(x, ifelse(shell, "+shell", ""))
cv$pair <- paste0(currency(cv$from, cv$from_shell), "_to_", currency(cv$to, cv$to_shell))
energy_pairs <- c("dry_mass_to_energy", "afdm_to_energy", "wet_mass_to_energy")
is_energy <- cv$pair %in% energy_pairs
u <- cv[cv$use_in_model & (cv$factor_unit == "ratio" | is_energy) & is.finite(cv$factor) &
          cv$factor > 0, ]
say("- %d of %d factors used: mass ratios with `use_in_model` and energy densities per dry mass, AFDM and wet mass without shell (screened below; the models of algae and other organisms are fitted too, but the package uses only those of invertebrates and fish); %d energy densities per carbon or of masses with shell left out.",
    nrow(u), nrow(cv), sum(cv$to == "energy" & !is_energy))

# Screening of the energy densities -------------------------------------------------
screen <- list()
drop_ids <- function(ids, step, reason) {
  if (length(ids)) {
    screen[[length(screen) + 1L]] <<- data.frame(conversion_id = ids, step = step, reason = reason,
                                                 stringsAsFactors = FALSE)
  }
  u <<- u[!u$conversion_id %in% ids, ]
}
J_mg <- function(x) x$factor * 1000
e <- u[u$pair %in% energy_pairs, ]
drop_ids(e$conversion_id[e$pair == "afdm_to_energy" & (J_mg(e) < 15 | J_mg(e) > 45)], 1,
         "outside 15-45 J/mg AFDM")
e <- u[u$pair %in% energy_pairs & !is.na(u$sample), ]
by_sample <- split(e, e$sample)
bad <- unlist(lapply(by_sample, function(s) {
  v <- function(p) s$factor[s$pair == p][1]
  ids <- character(0)
  dw <- v("dry_mass_to_energy")
  af <- v("afdm_to_energy")
  ww <- v("wet_mass_to_energy")
  if (!is.na(dw) && !is.na(af) && dw > 1.05 * af) {
    ids <- c(ids, s$conversion_id[s$pair %in% c("dry_mass_to_energy", "afdm_to_energy")])
  }
  if (!is.na(ww) && !is.na(dw) && ww > 1.05 * dw) {
    ids <- c(ids, s$conversion_id[s$pair %in% c("wet_mass_to_energy", "dry_mass_to_energy")])
  }
  if (!is.na(ww) && is.na(dw) && !is.na(af) && ww > 1.05 * af) {
    ids <- c(ids, s$conversion_id[s$pair %in% c("wet_mass_to_energy", "afdm_to_energy")])
  }
  unique(ids)
}), use.names = FALSE)
drop_ids(unique(bad), 2, "inconsistent with another energy density of the same material")

d <- data.frame(record_id = u$conversion_id, source_id = u$source_id, taxon = u$taxon, pair = u$pair,
                habitat = u$habitat, stage = "unspecified", y = log(u$factor), stringsAsFactors = FALSE)
d <- add_lineage(d, rt)
stopifnot(attr(d, "lost") == 0)
# Tadpoles (decision of 2026-11-09): their factors form models of their own,
# domain "amphibian_larva", apart from those of adult amphibians (whose few
# factors give no model).
tad <- d$domain == "amphibian" & u$life_stage[match(d$record_id, u$conversion_id)] %in% "larva"
d$domain[tad] <- "amphibian_larva"
d$model <- paste(d$domain, d$pair)
robust_out <- character(0)
for (m in unique(d$model[d$pair %in% energy_pairs])) {
  x <- d[d$model == m, ]
  n_class <- table(x$class)
  n_phylum <- table(x$phylum)
  by_class <- !is.na(x$class) & as.vector(n_class[x$class]) >= 10
  by_phylum <- !by_class & !is.na(x$phylum) & as.vector(n_phylum[x$phylum]) >= 10
  ref <- ifelse(by_class, paste("class", x$class), ifelse(by_phylum, paste("phylum", x$phylum), "all"))
  for (g in unique(ref)) {
    k <- ref == g
    pool <- if (startsWith(g, "class ")) x$y[x$class %in% sub("^class ", "", g)] else
      if (startsWith(g, "phylum ")) x$y[x$phylum %in% sub("^phylum ", "", g)] else x$y
    spread <- max(stats::mad(pool), 0.1)
    z <- abs(x$y[k] - stats::median(pool)) / spread
    robust_out <- c(robust_out, x$record_id[k][z > 4])
  }
}
drop_ids(robust_out, 3, "robust z beyond 4 within its taxonomic group")
d <- d[!d$record_id %in% robust_out, ]

n_sp <- tapply(node_path(d, "species"), d$model, function(z) length(unique(z)))
n_row <- table(d$model)
fit_these <- names(n_row)[n_row >= 30 & n_sp[names(n_row)] >= 5]
say("- %d models (pair and domain) with 30 or more factors and 5 or more species nodes; not fitted: %s.",
    length(fit_these), paste(sprintf("%s (%d)", setdiff(names(n_row), fit_these),
                                     as.integer(n_row[setdiff(names(n_row), fit_these)])), collapse = ", "))
say("")
fit_one <- function(x) {
  levels <- split_levels(x)
  opt <- model_options(x, "typical")
  # No random effect of source (decision of 2026-11-06). The models have two
  # to four sources, one of which (Brey et al. 2010, a compilation of many
  # studies whose differences are part of the residual) holds most of the
  # factors: too few to estimate a variance between sources. It was estimated
  # as 0 in 8 of the 13 models with two or more sources; in the others it put
  # the species of the larger source below their own factors (0.94 times in
  # dry mass to AFDM, 0.95 in wet to dry mass) and, in wet mass with shell to
  # dry mass, those of Hrycik et al. (2024) at 0.78 times theirs (0.90
  # without it), with the same coverage in cross-validation and intervals as
  # wide or wider (dev/validation/source_effect_conversions.md).
  opt$source <- FALSE
  x <- model_frame(x, "typical", levels)
  warn <- character(0)
  fit <- withCallingHandlers(fit_model(x, "typical", levels, opt),
                             warning = function(w) {
                               warn <<- c(warn, conditionMessage(w))
                               invokeRestart("muffleWarning")
                             },
                             message = function(w) {
                               warn <<- c(warn, trimws(conditionMessage(w)))
                               invokeRestart("muffleMessage")
                             })
  list(fit = fit, x = x, levels = levels, opt = opt, warn = warn)
}
out <- list(nodes = list(), components = list(), fixed = list())
used <- list()
for (m in sort(fit_these)) {
  x <- d[d$model == m, ]
  f <- fit_one(x)
  n_resid <- 0L
  if (sub("^[a-z_]+ ", "", m) %in% energy_pairs) {
    r <- stats::residuals(f$fit) / stats::sigma(f$fit)
    out_ids <- f$x$record_id[abs(r) > 4]
    if (length(out_ids)) {
      n_resid <- length(out_ids)
      drop_ids(out_ids, 4, "standardized residual beyond 4")
      x <- x[!x$record_id %in% out_ids, ]
      f <- fit_one(x)
    }
  }
  fit <- f$fit
  x <- f$x
  levels <- f$levels
  opt <- f$opt
  warn <- f$warn
  used[[m]] <- data.frame(conversion_id = x$record_id, model = m, stringsAsFactors = FALSE)
  tabs <- model_tables(fit, x, m, "typical", levels)
  for (k in names(out)) out[[k]][[m]] <- tabs[[k]]
  cp <- tabs$components
  v <- cp[cp$quantity == "variance", ]
  vv <- stats::setNames(v$value, v$level)
  pop <- sum(vv[intersect(names(vv), c("residual", "source"))])
  below <- rev(cumsum(rev(vv[levels])))
  fac <- exp(1.96 * sqrt(pop + c(below[-1], 0)))
  ctx <- tabs$nodes[tabs$nodes$level == "habitat", ]
  energy_model <- sub("^[a-z_]+ ", "", m) %in% energy_pairs
  say("- **%s**: %d factors of %d species nodes (%s); levels %s%s; median %s %.3g; SD (ln): %s; 95 %% factor (x/÷) for a new population by the lowest level known: %s%s%s%s.",
      m, nrow(x), length(unique(node_path(x, "species"))),
      paste(sprintf("%s %d", names(table(x$source_id)), as.integer(table(x$source_id))), collapse = ", "),
      paste(levels, collapse = ", "),
      if (length(opt$context)) paste0("; plus ", paste(opt$context, collapse = ", ")) else "",
      if (energy_model) "(J/mg)" else "factor",
      exp(stats::median(x$y)) * if (energy_model) 1000 else 1,
      paste(sprintf("%s %.2f", v$level, v$sd_cor), collapse = ", "),
      paste(sprintf("%s %.2f", rev(levels), rev(fac)), collapse = ", "),
      if (nrow(ctx)) paste0("; habitat effects (ln) ", paste(sprintf("%s %.2f", ctx$node, ctx$effect), collapse = ", ")) else "",
      if (n_resid) sprintf("; %d values with standardized residual beyond 4 left out and the model refitted", n_resid) else "",
      if (length(warn)) paste0(" (", paste(unique(warn), collapse = "; "), ")") else "")
}

# Screening summary --------------------------------------------------------------------
scr <- if (length(screen)) do.call(rbind, screen) else
  data.frame(conversion_id = character(0), step = integer(0), reason = character(0))
k <- match(scr$conversion_id, cv$conversion_id)
scr$source_id <- cv$source_id[k]
scr$taxon <- cv$taxon[k]
scr$pair <- cv$pair[k]
scr$habitat <- cv$habitat[k]
scr$J_mg <- cv$factor[k] * 1000
scr$original_reference <- cv$original_reference[k]
scr <- scr[order(scr$step, scr$pair, scr$taxon, method = "radix"), ]
say("")
say("## Screening of the energy densities")
say("")
e_all <- cv[cv$use_in_model & is_energy & is.finite(cv$factor) & cv$factor > 0, ]
say("- %d energy densities per dry mass, AFDM or wet mass without shell entered the screening (%s).",
    nrow(e_all), paste(sprintf("%s %d", names(table(e_all$source_id)), as.integer(table(e_all$source_id))),
                       collapse = ", "))
for (s in 1:4) {
  z <- scr[scr$step == s, ]
  say("- Step %d (%s): %d left out (%s).", s, c("bounds per AFDM", "same material", "robust z",
                                               "model residuals")[s], nrow(z),
      if (nrow(z)) paste(sprintf("%s %d", names(table(z$pair)), as.integer(table(z$pair))), collapse = ", ") else "none")
}
say("- Values left out are listed in `conversion_screening.csv`.")

nodes <- do.call(rbind, out$nodes)
components <- do.call(rbind, out$components)
fixed <- do.call(rbind, out$fixed)
used <- do.call(rbind, used)
rownames(nodes) <- rownames(components) <- rownames(fixed) <- rownames(used) <- NULL
write.csv(nodes, file.path(out_dir, "conversion_nodes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(components, file.path(out_dir, "conversion_components.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(fixed, file.path(out_dir, "conversion_fixed.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(used, file.path(out_dir, "conversion_used.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(scr, file.path(out_dir, "conversion_screening.csv"), row.names = FALSE, fileEncoding = "UTF-8")
writeLines(log_lines, file.path(out_dir, "conversion_fit_log.md"), useBytes = TRUE)
cat(log_lines, sep = "\n")
