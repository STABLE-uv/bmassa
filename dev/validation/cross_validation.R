# Validation of the hierarchical models (outside the check)
#
# Cross-validation of the prediction intervals:
#  - "unit": 10 groups of records (typical masses) or of equations (length-mass):
#    a new population of a species of the reference; the samplings of a
#    population (fish of NEON, column `population`) go together;
#  - "species": 10 groups of species nodes: a species that is not in the
#    reference, predicted from the lowest level left;
#  - "source": each source (typical masses) or study (length-mass) left out in
#    turn.
# For each unit left out (record; centre of the length range of each
# equation; population mean): z = (y - eta) / sqrt(v_total); coverage of the
# 50, 90 and 95 % intervals, bias and median absolute error (ln), by model,
# scheme and level of the prediction.
# Writes dev/validation/cross_validation.csv (one row per unit left out) and
# dev/validation/cross_validation.md.
# From the package root: Rscript dev/validation/cross_validation.R  (about 30 min)

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
dat <- read_model_data()
specs <- model_specs(dat$typ, dat$lwr)
cores <- min(2L, parallel::detectCores())
k_folds <- 10L
seed <- 20261102L

one_fold <- function(d_all, type, test, model, scheme, fold) {
  train <- d_all[!test, ]
  levels <- split_levels(train)
  opt <- model_options(train, type)
  train <- model_frame(train, type, levels)
  fit <- suppressMessages(suppressWarnings(fit_model(train, type, levels, opt)))
  tabs <- model_tables(fit, train, model, type, levels)
  nd <- d_all[test, ]
  if (type == "lwr") nd <- nd[nd$point %in% c(0L, 2L), ]
  pr <- predict_tables(tabs, nd, type, levels)
  unit <- if (type == "lwr") nd$equation_id else nd$record_id
  data.frame(model = model, scheme = scheme, fold = fold, unit = unit, source_id = nd$source_id,
             taxon = nd$taxon, lowest_rank = nd$lowest_rank,
             length_type = if (type == "lwr") nd$length_type else NA_character_,
             y = nd$y, pr, stringsAsFactors = FALSE)
}

# Partial results by model in dev/validation/cv/ (reused if they exist;
# BMASSA_CV_REDO=true computes them again).
cv_dir <- "dev/validation/cv"
dir.create(cv_dir, showWarnings = FALSE)
redo <- identical(Sys.getenv("BMASSA_CV_REDO"), "true")
res <- list()
for (i in seq_len(nrow(specs))) {
  m <- specs$model[i]
  type <- specs$type[i]
  part <- file.path(cv_dir, paste0(m, ".csv"))
  if (file.exists(part) && !redo) {
    res[[m]] <- read.csv(part, stringsAsFactors = FALSE, encoding = "UTF-8")
    next
  }
  d <- model_rows(dat, m, type)
  unit <- if (type == "lwr") d$equation_id else d$record_id
  # Groups by population: the records of a population sampled several times
  # are left out together.
  if (type != "lwr" && "population" %in% names(d)) {
    unit <- ifelse(is.na(d$population), unit, d$population)
  }
  sp <- node_path(d, "species")
  src <- if (type == "lwr") d$study else d$source_id
  set.seed(seed)
  u <- unique(unit)
  fold_unit <- stats::setNames(sample(rep(seq_len(k_folds), length.out = length(u))), u)
  s <- unique(sp)
  fold_sp <- stats::setNames(sample(rep(seq_len(k_folds), length.out = length(s))), s)
  jobs <- c(lapply(seq_len(k_folds), function(f) list(scheme = "unit", fold = f, test = fold_unit[unit] == f)),
            lapply(seq_len(k_folds), function(f) list(scheme = "species", fold = f, test = fold_sp[sp] == f)),
            if (length(unique(src)) > 1)
              lapply(sort(unique(src)), function(z) list(scheme = "source", fold = z, test = src == z)))
  t0 <- Sys.time()
  out <- parallel::mclapply(jobs, function(j) {
    tryCatch(one_fold(d, type, j$test, m, j$scheme, as.character(j$fold)),
             error = function(e) {
               data.frame(model = m, scheme = j$scheme, fold = as.character(j$fold), unit = NA_character_,
                          source_id = NA_character_, taxon = NA_character_, lowest_rank = NA_character_,
                          length_type = NA_character_, y = NA_real_, eta = NA_real_, v_param = NA_real_,
                          v_missing = NA_real_, v_population = NA_real_, v_total = NA_real_,
                          matched_level = paste("error:", conditionMessage(e)), stringsAsFactors = FALSE)
             })
  }, mc.cores = cores)
  res[[m]] <- do.call(rbind, out)
  write.csv(res[[m]], part, row.names = FALSE, fileEncoding = "UTF-8")
  message(sprintf("%s: %d fits in %.1f min", m, length(jobs),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, res)
rownames(res) <- NULL
res$z <- (res$y - res$eta) / sqrt(res$v_total)
write.csv(res, "dev/validation/cross_validation.csv", row.names = FALSE, fileEncoding = "UTF-8")

# Summary --------------------------------------------------------------------------
ok <- is.finite(res$z)
summ <- function(x) {
  data.frame(n = nrow(x), cov50 = mean(abs(x$z) < stats::qnorm(0.75)), cov90 = mean(abs(x$z) < stats::qnorm(0.95)),
             cov95 = mean(abs(x$z) < stats::qnorm(0.975)), z_mean = mean(x$z), z_sd = stats::sd(x$z),
             bias = mean(x$y - x$eta), mdae = stats::median(abs(x$y - x$eta)),
             factor95 = stats::median(exp(1.96 * sqrt(x$v_total))))
}
r <- res[ok, ]
lv <- c("root", levels_all)
by_level <- do.call(rbind, lapply(split(r, list(r$model, r$scheme, r$matched_level), drop = TRUE), function(x)
  cbind(x[1, c("model", "scheme", "matched_level")], summ(x))))
by_level <- by_level[order(by_level$model, by_level$scheme, match(by_level$matched_level, lv)), ]
overall <- do.call(rbind, lapply(split(r, list(r$model, r$scheme), drop = TRUE), function(x)
  cbind(x[1, c("model", "scheme")], summ(x))))
overall <- overall[order(overall$model, overall$scheme), ]
fmt <- function(x) {
  x$cov50 <- sprintf("%.0f %%", 100 * x$cov50)
  x$cov90 <- sprintf("%.0f %%", 100 * x$cov90)
  x$cov95 <- sprintf("%.0f %%", 100 * x$cov95)
  for (j in c("z_mean", "z_sd", "bias", "mdae")) x[[j]] <- sprintf("%.2f", x[[j]])
  x$factor95 <- sprintf("%.1f", x$factor95)
  x
}
md_table <- function(x) {
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep("---|", ncol(x)), collapse = "")),
    apply(x, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |")))
}
not_pred <- res[!ok, ]
not_pred$reason <- ifelse(grepl("^error:", not_pred$matched_level),
                          paste("fit impossible:", sub("^error: ", "", not_pred$matched_level)),
                          ifelse(!is.na(not_pred$length_type) & not_pred$matched_level == "root",
                                 "length type without data in the fit",
                                 "stage without data in the fit"))
np <- if (nrow(not_pred)) {
  a <- stats::aggregate(list(n = rep(1L, nrow(not_pred))),
                        by = not_pred[, c("model", "scheme", "fold", "reason")], FUN = sum)
  a[order(a$model, a$scheme, a$fold), ]
} else NULL
lines <- c(paste0("# Cross-validation of the hierarchical models (version ",
                  read.dcf("DESCRIPTION", fields = "Version")[1, 1], ")"), "",
           paste0("Generated by `dev/validation/cross_validation.R` on ", format(Sys.Date()), " (",
                  R.version.string, ", lme4 ", as.character(utils::packageVersion("lme4")), "; seed ", seed, ")."),
           "",
           "Schemes: `unit`, a new population of a species of the reference (10 groups of records or equations; the samplings of a population of fish of NEON go together); `species`, a new species (10 groups of species); `source`, a new source or study (each one left out in turn). Equations are evaluated at the centre of their length range. `cov50`, `cov90`, `cov95`: proportion of observed values inside the 50, 90 and 95 % prediction intervals; `z_sd`: standard deviation of z (1 if the variance is well calibrated); `bias` and `mdae`: mean of y − η and median of |y − η| (ln); `factor95`: median of the ×/÷ factor of the 95 % interval.",
           "", "## Overall", "", md_table(fmt(overall)), "", "## By level of the prediction", "",
           md_table(fmt(by_level)), "",
           "## Reading", "",
           "- `unit` and `species` measure what the model promises: a new population of a known species and a new species of a known higher taxon. The 95 % coverage should lie between 90 and 97 %.",
           "- `source` is harder than real use will be: with two sources (phytoplankton, amphibians), leaving one out leaves no variance between sources to estimate, so the intervals do not include it; in invertebrates the sources cover different groups (zooplankton of Dumont, benthos of Hrycik and Moreau) and in fish different regions (Iberian streams, the Mekong, the United States), and leaving one out means predicting from orders, classes or phyla with another mix of sizes. In the full model the variance between sources does enter the prediction, except in invertebrates and fish without lengths, which have no source effect. In the length-mass models `source` leaves out each study; in that of invertebrates, the variance of a new study is part of that between equations (fit without study; `study_effect_lwr.md`).",
           "",
           sprintf("## Units without prediction (%d)", nrow(not_pred)), "",
           if (is.null(np)) "None." else md_table(np))
writeLines(lines, "dev/validation/cross_validation.md", useBytes = TRUE)
cat(lines, sep = "\n")
