# Study effect in the length-mass model of invertebrates (outside the check)
#
# The length-mass models (data-raw/fit_models.R) have a random effect of study
# when they have five or more studies: only that of invertebrates (46 studies;
# fish have 3 and amphibians 2). A prediction is then for a new study (the
# effects of the studies left out, their variance added). Many of the studies
# cover a single group (the rotifers of Ejsmont-Karabin 1998, the zooplankton
# of Dumont et al. 1975, the benthos of streams of Benke et al. 1999), so the
# effect of a study may take part of the effect of its taxa. Compares the fit
# with the study effect (A, as in 0.3.0) and without it (N):
#  - the components (SD);
#  - the median, by study, of the ratio between the prediction for a species
#    node at the centre of the length range of an equation and the value of
#    the equation there, in the species with equations of a single study
#    (1 = the prediction keeps the level of its data);
#  - a cross-validation with the folds of dev/validation/cross_validation.R
#    (same seed): 10 groups of equations (`unit`), 10 groups of species
#    (`species`) and each study left out in turn (`study`), each equation
#    evaluated at the centre of its length range; coverage of the 90 and 95 %
#    intervals, standard deviation of z, the median ×/÷ factor of the 95 %
#    interval, bias and median absolute error (ln), overall and by study.
# Partial results in dev/validation/cv/study_effect/ (not in git; reused if
# they exist). Writes dev/validation/study_effect_lwr.md.
# From the package root: Rscript dev/validation/study_effect_lwr.R  (about 45 min with 2 cores)

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
dat <- read_model_data()
m <- "invertebrate_lwr"
type <- "lwr"
d <- model_rows(dat, m, type)
stopifnot(length(unique(d$study)) >= 5)
cores <- min(2L, parallel::detectCores())
k_folds <- 10L
seed <- 20261102L

fit_variant <- function(train, v) {
  levels <- split_levels(train)
  opt <- model_options(train, type)
  opt$study <- v == "A" && length(unique(train$study)) >= 5
  train <- model_frame(train, type, levels)
  fit <- suppressMessages(suppressWarnings(fit_model(train, type, levels, opt)))
  list(tabs = model_tables(fit, train, m, type, levels), levels = levels, fit = fit)
}
predict_centre <- function(a, nd) {
  nd <- nd[nd$point == 2L, ]
  pr <- predict_tables(a$tabs, nd, type, a$levels)
  data.frame(equation_id = nd$equation_id, study = nd$study, species = node_path(nd, "species"),
             y = nd$y, pr, stringsAsFactors = FALSE)
}

# Folds as in cross_validation.R.
unit <- d$equation_id
sp <- node_path(d, "species")
set.seed(seed)
u <- unique(unit)
fold_unit <- stats::setNames(sample(rep(seq_len(k_folds), length.out = length(u))), u)
s <- unique(sp)
fold_sp <- stats::setNames(sample(rep(seq_len(k_folds), length.out = length(s))), s)
folds <- c(lapply(seq_len(k_folds), function(f) list(scheme = "unit", fold = as.character(f), test = fold_unit[unit] == f)),
           lapply(seq_len(k_folds), function(f) list(scheme = "species", fold = as.character(f), test = fold_sp[sp] == f)),
           lapply(sort(unique(d$study)), function(z) list(scheme = "study", fold = z, test = d$study == z)))
jobs <- c(lapply(c("A", "N"), function(v) list(variant = v, scheme = "full", fold = "full", test = NULL)),
          unlist(lapply(c("A", "N"), function(v) lapply(folds, function(f) c(list(variant = v), f))),
                 recursive = FALSE))

part_dir <- "dev/validation/cv/study_effect"
dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
part_file <- function(j) file.path(part_dir, sprintf("%s_%s_%s.rds", j$variant, j$scheme,
                                                     gsub("[^A-Za-z0-9]+", "_", j$fold)))
t0 <- Sys.time()
invisible(parallel::mclapply(jobs, function(j) {
  f <- part_file(j)
  if (file.exists(f)) return(NULL)
  out <- tryCatch({
    if (is.null(j$test)) {
      a <- fit_variant(d, j$variant)
      list(pred = predict_centre(a, d), components = a$tabs$components)
    } else {
      a <- fit_variant(d[!j$test, ], j$variant)
      list(pred = predict_centre(a, d[j$test, ]))
    }
  }, error = function(e) list(error = conditionMessage(e)))
  saveRDS(c(j[c("variant", "scheme", "fold")], out), f)
  NULL
}, mc.cores = cores, mc.preschedule = FALSE))
message(sprintf("%d fits in %.1f min", length(jobs), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
parts <- lapply(jobs, function(j) readRDS(part_file(j)))
err <- Filter(function(p) !is.null(p$error), parts)
if (length(err)) message("Fits that failed: ", paste(vapply(err, function(p)
  paste(p$variant, p$scheme, p$fold, p$error), character(1)), collapse = "; "))

# Components ----------------------------------------------------------------------
full <- Filter(function(p) p$scheme == "full", parts)
comp <- do.call(rbind, lapply(full, function(p) {
  x <- p$components[p$components$quantity == "variance" & p$components$parameter == "intercept", ]
  data.frame(variant = p$variant, level = x$level, sd = round(x$sd_cor, 3), stringsAsFactors = FALSE)
}))
lv_order <- c(levels_all, "equation", "study", "residual")
comp <- comp[order(match(comp$level, lv_order), comp$variant), ]
comp <- stats::reshape(comp, idvar = "level", timevar = "variant", direction = "wide")
names(comp) <- sub("^sd\\.", "sd_", names(comp))
slopes <- do.call(rbind, lapply(full, function(p) {
  x <- p$components[p$components$quantity == "variance" & p$components$parameter == "slope", ]
  if (!nrow(x)) return(NULL)
  data.frame(variant = p$variant, level = x$level, sd = round(x$sd_cor, 3), stringsAsFactors = FALSE)
}))

# Species of a single study: prediction / value of the equation at its centre.
pf <- do.call(rbind, lapply(full, function(p) cbind(variant = p$variant, p$pred)))
pf <- pf[!grepl("not identified", pf$species, fixed = TRUE), ]
n_st <- tapply(pf$study, pf$species, function(z) length(unique(z)))
single <- pf[n_st[pf$species] == 1, ]
single$ratio <- exp(single$eta - single$y)
sm <- stats::aggregate(ratio ~ study + variant, single, stats::median)
sm <- stats::reshape(sm, idvar = "study", timevar = "variant", direction = "wide")
ns <- stats::aggregate(list(equations = single$equation_id[single$variant == "A"]),
                       by = list(study = single$study[single$variant == "A"]), FUN = length)
nsp <- stats::aggregate(list(species = single$species[single$variant == "A"]),
                        by = list(study = single$study[single$variant == "A"]),
                        FUN = function(z) length(unique(z)))
sm <- merge(merge(ns, nsp, by = "study"), sm, by = "study")
sm$ratio.A <- round(sm$ratio.A, 2)
sm$ratio.N <- round(sm$ratio.N, 2)
sm <- sm[order(-sm$equations, sm$study), ]
names(sm) <- c("study", "equations", "species", "ratio_A", "ratio_N")
all_single <- stats::aggregate(ratio ~ variant, single, function(r) stats::quantile(r, c(0.1, 0.5, 0.9)))

# Cross-validation -----------------------------------------------------------------
cvr <- do.call(rbind, lapply(Filter(function(p) p$scheme != "full" && is.null(p$error), parts),
                             function(p) cbind(variant = p$variant, scheme = p$scheme, fold = p$fold, p$pred)))
cvr$z <- (cvr$y - cvr$eta) / sqrt(cvr$v_total)
n_na <- sum(!is.finite(cvr$z))
cvr <- cvr[is.finite(cvr$z), ]
summ <- function(x) {
  data.frame(n = nrow(x), cov90 = round(100 * mean(abs(x$z) < stats::qnorm(0.95))),
             cov95 = round(100 * mean(abs(x$z) < stats::qnorm(0.975))), z_sd = round(stats::sd(x$z), 2),
             factor95 = round(stats::median(exp(stats::qnorm(0.975) * sqrt(x$v_total))), 2),
             bias = round(mean(x$y - x$eta), 3), mdae = round(stats::median(abs(x$y - x$eta)), 3))
}
overall <- do.call(rbind, lapply(split(cvr, list(cvr$scheme, cvr$variant), drop = TRUE), function(x)
  cbind(x[1, c("scheme", "variant")], summ(x))))
overall <- overall[order(match(overall$scheme, c("unit", "species", "study")), overall$variant), ]
big <- names(which(table(d$study[d$point == 2L]) >= 20))
by_study <- do.call(rbind, lapply(split(cvr[cvr$study %in% big, ], list(cvr$scheme[cvr$study %in% big],
                                                                       cvr$variant[cvr$study %in% big],
                                                                       cvr$study[cvr$study %in% big]), drop = TRUE),
                                  function(x) cbind(x[1, c("study", "scheme", "variant")], summ(x))))
by_study <- by_study[order(by_study$study, match(by_study$scheme, c("unit", "species", "study")), by_study$variant), ]

# Deviation shared by the equations of a study left out: SD of the mean of
# y - eta by study, in the equations of species measured by other studies
# (predicted at the species level), against the SD expected if their errors
# were independent (square root of the mean of v_total / n).
sh <- cvr[cvr$scheme == "study" & cvr$matched_level == "species", ]
shared <- do.call(rbind, lapply(c("A", "N"), function(v) {
  x <- sh[sh$variant == v, ]
  mr <- tapply(x$y - x$eta, x$study, mean)
  n <- tapply(x$y, x$study, length)
  vt <- tapply(x$v_total, x$study, mean)
  do.call(rbind, lapply(c(5, 10), function(k) {
    j <- n >= k
    s_obs <- stats::sd(mr[j])
    s_ind <- sqrt(mean(vt[j] / n[j]))
    data.frame(variant = v, min_equations = k, studies = sum(j), equations = sum(n[j]),
               sd_of_means = round(s_obs, 2), sd_if_independent = round(s_ind, 2),
               shared_sd = round(sqrt(max(s_obs^2 - s_ind^2, 0)), 2))
  }))
}))

# Agreement of variant A with the cache of cross_validation.R (same folds).
cache <- "dev/validation/cv/invertebrate_lwr.csv"
if (file.exists(cache)) {
  cc <- read.csv(cache, stringsAsFactors = FALSE, encoding = "UTF-8")
  cc$scheme[cc$scheme == "source"] <- "study"
  a <- cvr[cvr$variant == "A", ]
  k <- match(paste(a$scheme, a$fold, a$equation_id), paste(cc$scheme, cc$fold, cc$unit))
  message(sprintf("Variant A against cross_validation.R: %d of %d units matched, max |eta difference| %.2g",
                  sum(!is.na(k)), nrow(a), max(abs(a$eta - cc$eta[k]), na.rm = TRUE)))
}

md_table <- function(x) {
  x[] <- lapply(x, function(v) if (is.numeric(v)) format(v, trim = TRUE) else v)
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep("---|", ncol(x)), collapse = "")),
    apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
q <- function(v) {
  r <- all_single$ratio[all_single$variant == v, ]
  sprintf("%.2f (10 %%: %.2f; 90 %%: %.2f)", r[2], r[1], r[3])
}
lines <- c("# Study effect in the length-mass model of invertebrates", "",
           paste0("Generated by `dev/validation/study_effect_lwr.R` on ", format(Sys.Date()), " (",
                  R.version.string, ", lme4 ", as.character(utils::packageVersion("lme4")), "; seed ", seed,
                  ", the folds of `cross_validation.R`)."), "",
           sprintf("Model `%s`: %d equations of %d studies (%d species nodes). Variants: `A`, random effect of study (as fitted in 0.3.0); `N`, without it. Both keep the random effect of equation (variation between populations).",
                   m, length(unique(d$equation_id)), length(unique(d$study)), length(unique(sp))), "",
           "## Components (SD of the intercepts, ln)", "", md_table(comp), "",
           if (!is.null(slopes)) c(paste0("SD of the slopes: ", paste(sprintf("%s %s %.3f", slopes$variant, slopes$level,
                                                                              slopes$sd), collapse = "; "), "."), "") else NULL,
           "## Species of a single study: prediction / equation at the centre of its range (median)", "",
           sprintf("All equations of species of a single study: A %s; N %s.", q("A"), q("N")), "",
           md_table(sm), "",
           "## Cross-validation", "",
           "`unit`: 10 groups of equations (a new population of a species of the reference); `species`: 10 groups of species nodes (a new species); `study`: each study left out in turn (a new study). Each equation is evaluated at the centre of its length range. `cov90`, `cov95`: % of values inside the 90 and 95 % prediction intervals; `z_sd`: SD of z (1 if the variance is well calibrated); `factor95`: median ×/÷ factor of the 95 % interval; `bias`, `mdae`: mean of y − η and median of |y − η| (ln).", "",
           md_table(overall), "",
           sprintf("Units without prediction (length type without data in the fit): %d.", n_na), "",
           "### Deviation shared by the equations of a study left out", "",
           "Scheme `study`, equations of species measured by other studies (predicted at the species level), in the studies with at least `min_equations` of them: SD of the mean of y − η by study (`sd_of_means`), the SD expected if the errors of the equations were independent (`sd_if_independent`, square root of the mean of v_total / n) and the SD of a deviation shared by the equations of a study implied by the difference (`shared_sd`). The intervals of each row do not depend on it, but a total over many taxa of a new study does.", "",
           md_table(shared), "",
           "### Studies with 20 or more equations", "",
           md_table(by_study))
writeLines(lines, "dev/validation/study_effect_lwr.md", useBytes = TRUE)
cat(lines, sep = "\n")
