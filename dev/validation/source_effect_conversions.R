# Source effect in the conversion models (outside the check)
#
# The conversion models (data-raw/fit_conversions.R) have a random effect of
# source when they have two or more sources; a prediction is then for a new
# source (the effects of the sources left out, their variance added). For
# each such model, compares the fit with it (A, as in 0.3.0) and without it
# (N):
#  - the components (SD of source, of the residual and of the levels);
#  - the median, by source, of the ratio between the prediction for a species
#    node and the mean of its factors, in the species measured by a single
#    source (1 = the prediction keeps the level of its data);
#  - a cross-validation: 10 groups of factors (`unit`), 10 groups of species
#    (`species`) and each source left out in turn (`source`, when what is left
#    has 30 or more factors of 5 or more species), with the coverage of the
#    95 % interval, the median ×/÷ factor of that interval and the bias (mean
#    of y − η, ln) by source.
# The data are those of the fitted models (data-raw/derived/conversion_used.csv,
# after the screening of fit_conversions.R).
# Writes dev/validation/source_effect_conversions.md.
# From the package root: Rscript dev/validation/source_effect_conversions.R  (about 2 h)

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
rt <- read.csv(file.path(out_dir, "ref_taxa.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
cv0 <- read.csv(file.path(out_dir, "conversions.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
used <- read.csv(file.path(out_dir, "conversion_used.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
currency <- function(x, shell) paste0(x, ifelse(shell, "+shell", ""))
cv0$pair <- paste0(currency(cv0$from, cv0$from_shell), "_to_", currency(cv0$to, cv0$to_shell))
u <- cv0[match(used$conversion_id, cv0$conversion_id), ]
d_all <- data.frame(record_id = u$conversion_id, source_id = u$source_id, taxon = u$taxon,
                    pair = u$pair, habitat = u$habitat, stage = "unspecified", y = log(u$factor),
                    stringsAsFactors = FALSE)
d_all <- add_lineage(d_all, rt)
stopifnot(attr(d_all, "lost") == 0)
# Tadpoles: models of their own, as in fit_conversions.R.
d_all$domain[d_all$domain == "amphibian" & u$life_stage %in% "larva"] <- "amphibian_larva"
d_all$model <- used$model
stopifnot(identical(d_all$model, paste(d_all$domain, d_all$pair)))

# Models with two or more sources; which are used by the package
# (convert_currency()) or by data-raw/model_data.R.
package_models <- c(paste("fish", c("afdm_to_carbon", "afdm_to_energy", "dry_mass_to_afdm",
                                     "dry_mass_to_carbon", "dry_mass_to_energy", "dry_mass_to_nitrogen",
                                     "dry_mass_to_phosphorus", "wet_mass_to_afdm",
                                     "wet_mass_to_dry_mass", "wet_mass_to_energy")),
                    paste("invertebrate", c("afdm_to_carbon", "afdm_to_energy", "dry_mass_to_afdm",
                                             "dry_mass_to_carbon", "dry_mass_to_energy",
                                             "dry_mass_to_nitrogen", "dry_mass_to_phosphorus",
                                             "wet_mass_to_afdm", "wet_mass_to_dry_mass",
                                             "wet_mass_to_energy")),
                    paste("amphibian_larva", c("dry_mass_to_carbon", "dry_mass_to_nitrogen",
                                               "dry_mass_to_phosphorus")))
model_data_models <- paste("invertebrate", c("wet_mass_to_dry_mass", "wet_mass+shell_to_dry_mass",
                                             "dry_mass+shell_to_afdm+shell", "dry_mass+shell_to_dry_mass",
                                             "dry_mass_to_afdm"))
used_by <- function(m) {
  u <- c(if (m %in% package_models) "package", if (m %in% model_data_models) "model_data.R")
  if (length(u)) paste(u, collapse = ", ") else "—"
}
n_src <- tapply(d_all$source_id, d_all$model, function(s) length(unique(s)))
models <- sort(names(n_src)[n_src >= 2])

fit_variant <- function(d, m, v) {
  levels <- split_levels(d)
  opt <- model_options(d, "typical")
  opt$source <- v == "A" && length(unique(d$source_id)) >= 2
  f <- model_frame(d, "typical", levels)
  fit <- suppressMessages(suppressWarnings(fit_model(f, "typical", levels, opt)))
  list(tabs = model_tables(fit, f, m, "typical", levels), levels = levels, fit = fit)
}
# A habitat without factors in the fit is taken as unknown, as in
# convert_currency(): no habitat effect, and the habitat variance added.
predict_rows <- function(a, nd) {
  tabs <- a$tabs
  hab <- tabs$nodes$level == "habitat"
  if (any(hab) && any(!nd$habitat %in% tabs$nodes$node[hab])) {
    extra <- tabs$nodes[which(hab)[1], ]
    extra$node <- extra$path <- "<unknown>"
    extra$effect <- 0
    extra$effect_var <- tabs$components$value[tabs$components$level == "habitat" &
                                                tabs$components$quantity == "variance"]
    tabs$nodes <- rbind(tabs$nodes, extra)
    nd$habitat[!nd$habitat %in% tabs$nodes$node[hab]] <- "<unknown>"
  }
  predict_tables(tabs, nd, "typical", a$levels)
}

set.seed(20261106)
one_model <- function(m) {
  d0 <- d_all[d_all$model == m, ]
  rownames(d0) <- NULL
  u_id <- unique(d0$record_id)
  fold_u <- stats::setNames(sample(rep(1:10, length.out = length(u_id))), u_id)
  sp <- node_path(d0, "species")
  u_sp <- unique(sp)
  fold_s <- stats::setNames(sample(rep(1:10, length.out = length(u_sp))), u_sp)
  src <- sort(unique(d0$source_id))
  jobs <- c(lapply(1:10, function(f) list(scheme = "unit", test = fold_u[d0$record_id] == f)),
            lapply(1:10, function(f) list(scheme = "species", test = fold_s[sp] == f)),
            lapply(src, function(s) list(scheme = "source", test = d0$source_id == s)))
  jobs <- Filter(function(j) {
    rest <- d0[!j$test, ]
    nrow(rest) >= 30 && length(unique(node_path(rest, "species"))) >= 5
  }, jobs)
  res <- list()
  cv <- list()
  param <- list()
  for (v in c("A", "N")) {
    a <- fit_variant(d0, m, v)
    pr <- predict_rows(a, d0)
    res[[v]] <- data.frame(model = m, variant = v, species = sp, source = d0$source_id, y = d0$y,
                           eta = pr$eta, stringsAsFactors = FALSE)
    vc <- as.data.frame(lme4::VarCorr(a$fit))
    sdc <- function(g) { s <- vc$sdcor[vc$grp == g]; if (length(s)) s else NA_real_ }
    param[[v]] <- data.frame(model = m, variant = v, sd_source = sdc("source_id"),
                             sd_residual = sdc("Residual"), stringsAsFactors = FALSE)
    cv[[v]] <- do.call(rbind, lapply(jobs, function(j) {
      b <- tryCatch(fit_variant(d0[!j$test, ], m, v), error = function(e) NULL)
      if (is.null(b)) return(NULL)
      p2 <- tryCatch(predict_rows(b, d0[j$test, ]), error = function(e) NULL)
      if (is.null(p2)) return(NULL)
      data.frame(model = m, variant = v, scheme = j$scheme, source = d0$source_id[j$test],
                 y = d0$y[j$test], eta = p2$eta, vt = p2$v_total, stringsAsFactors = FALSE)
    }))
  }
  list(res = do.call(rbind, res), cv = do.call(rbind, cv), param = do.call(rbind, param),
       n = data.frame(model = m, factors = nrow(d0),
                      sources = paste(sprintf("%s %d", names(table(d0$source_id)),
                                              as.integer(table(d0$source_id))), collapse = ", "),
                      used_by = used_by(m),
                      stringsAsFactors = FALSE))
}
# Each model is saved on its own (dev/validation/cv/conversions/, not
# tracked; BMASSA_CV_REDO=true computes them again), so an interrupted run
# goes on from the models it had finished.
cache_dir <- file.path("dev", "validation", "cv", "conversions")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
if (identical(Sys.getenv("BMASSA_CV_REDO"), "true")) unlink(file.path(cache_dir, "*.rds"))
cache_file <- function(m) file.path(cache_dir, paste0(gsub("[^A-Za-z0-9_.+-]", "_", m), ".rds"))
todo <- models[!file.exists(cache_file(models))]
t0 <- Sys.time()
invisible(parallel::mclapply(todo, function(m) {
  saveRDS(one_model(m), cache_file(m))
  NULL
}, mc.cores = min(2L, parallel::detectCores()), mc.preschedule = FALSE))
out <- lapply(models, function(m) readRDS(cache_file(m)))
message(sprintf("%d models (%d run now) in %.1f min", length(models), length(todo),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
res <- do.call(rbind, lapply(out, `[[`, "res"))
cv <- do.call(rbind, lapply(out, `[[`, "cv"))
param <- do.call(rbind, lapply(out, `[[`, "param"))
param$sd_source <- round(param$sd_source, 3)
param$sd_residual <- round(param$sd_residual, 3)
info <- do.call(rbind, lapply(out, `[[`, "n"))
info$used_by <- vapply(info$model, used_by, "")

# Species measured by a single source: prediction of their node / mean of their factors.
res <- res[!grepl("not identified", res$species, fixed = TRUE), ]
key <- paste(res$model, res$species)
n_s <- tapply(res$source, key, function(s) length(unique(s)))
single <- res[n_s[key] == 1, ]
sm <- stats::aggregate(cbind(y, eta) ~ model + variant + source + species, single, mean)
sm$ratio <- exp(sm$eta - sm$y)
ratio <- stats::aggregate(ratio ~ model + source + variant, sm, stats::median)
n_sp <- stats::aggregate(species ~ model + source, sm[sm$variant == "A", ], length)
ratio <- stats::reshape(ratio, idvar = c("model", "source"), timevar = "variant", direction = "wide")
ratio <- merge(n_sp, ratio, by = c("model", "source"))
names(ratio)[3] <- "species"
ratio <- ratio[order(ratio$model, -ratio$species), ]

z <- (cv$y - cv$eta) / sqrt(cv$vt)
cv$inside <- abs(z) <= stats::qnorm(0.975)
cv$factor <- exp(stats::qnorm(0.975) * sqrt(cv$vt))
summ <- function(x) data.frame(n = nrow(x), cov95 = round(100 * mean(x$inside)),
                               factor95 = round(stats::median(x$factor), 2),
                               bias = round(mean(x$y - x$eta), 3))
overall <- do.call(rbind, lapply(split(cv, list(cv$model, cv$scheme, cv$variant), drop = TRUE), function(x)
  cbind(x[1, c("model", "scheme", "variant")], summ(x))))
overall <- overall[order(overall$model, overall$scheme, overall$variant), ]
by_source <- do.call(rbind, lapply(split(cv, list(cv$model, cv$scheme, cv$variant, cv$source), drop = TRUE),
                                   function(x) cbind(x[1, c("model", "scheme", "variant", "source")], summ(x))))
by_source <- by_source[by_source$scheme != "source" & by_source$source != "brey2010", ]
by_source <- by_source[order(by_source$model, by_source$source, by_source$scheme, by_source$variant), ]

md_table <- function(x) {
  x[] <- lapply(x, function(v) if (is.numeric(v)) format(v, digits = 3, trim = TRUE) else v)
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep("---|", ncol(x)), collapse = "")),
    apply(x, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
lines <- c("# Source effect in the conversion models", "",
           paste0("Generated by `dev/validation/source_effect_conversions.R` on ", format(Sys.Date()), " (",
                  R.version.string, "; *lme4* ", as.character(utils::packageVersion("lme4")),
                  ", Bates et al. 2015, <https://doi.org/10.18637/jss.v067.i01>",
                  "; seed 20261106)."), "",
           "Variants: `A`, random effect of source (as fitted in `data-raw/fit_conversions.R`); `N`, without it. Only the models with two or more sources. `used_by`: `package`, used by `convert_currency()`; `model_data.R`, used to convert reference data before the fits of the mass models.", "",
           "## Models", "", md_table(info), "",
           "## Components (SD, ln)", "", md_table(param), "",
           "## Species of a single source: prediction / mean of their factors (median)", "",
           md_table(ratio), "",
           "## Cross-validation", "",
           "`unit`: 10 groups of factors; `species`: 10 groups of species; `source`: each source left out in turn (when the rest has 30 or more factors of 5 or more species). `cov95`: % of values inside the 95 % interval; `factor95`: median ×/÷ factor of the interval; `bias`: mean of y − η (ln).", "",
           md_table(overall), "",
           "### Sources other than Brey et al. (2010), schemes `unit` and `species`", "",
           md_table(by_source))
writeLines(lines, "dev/validation/source_effect_conversions.md", useBytes = TRUE)
cat(lines, sep = "\n")
