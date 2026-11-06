# Source effect in the model of fish without lengths (outside the check)
#
# With the fish of NEON, the model fish_individual has three sources (Arranz
# et al. 2022, Iberian streams; Koning et al. 2024, a Mekong fishery; NEON,
# streams and lakes of the United States) that share few species. Compares two
# ways of handling them, both with the population effect (the samplings of a
# species at a NEON site):
#  - A: random effect of source (rule of 0.3.0 for the models with two or more
#    sources except that of invertebrates);
#  - N: no source effect, as in the model of invertebrates since 0.3.0.
# For each one: the median, by source, of the ratio between the prediction and
# the mean of the records of the species measured by a single source, and a
# cross-validation (10 groups of populations, 10 groups of species and each
# source left out in turn) with the coverage of the 95 % interval, the median
# ×/÷ factor of that interval and the bias (mean of y − η, ln) by source.
# Writes dev/validation/source_effect_fish.md.
# From the package root: Rscript dev/validation/source_effect_fish.R  (about 1 min)

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
dat <- read_model_data()
m <- "fish_individual"
d0 <- model_rows(dat, m, "typical")
variants <- c(A = "with source effect", N = "without source effect")

fit_variant <- function(d, v) {
  levels <- split_levels(d)
  opt <- model_options(d, "typical")
  opt$source <- v == "A"
  f <- model_frame(d, "typical", levels)
  fit <- suppressMessages(suppressWarnings(fit_model(f, "typical", levels, opt)))
  list(tabs = model_tables(fit, f, m, "typical", levels), levels = levels, fit = fit)
}

pop <- ifelse(is.na(d0$population), d0$record_id, d0$population)
set.seed(20261106)
up <- unique(pop)
fold_p <- stats::setNames(sample(rep(1:10, length.out = length(up))), up)
us <- unique(d0$species)
fold_s <- stats::setNames(sample(rep(1:10, length.out = length(us))), us)
jobs <- c(lapply(1:10, function(f) list(scheme = "population", test = fold_p[pop] == f)),
          lapply(1:10, function(f) list(scheme = "species", test = fold_s[d0$species] == f)),
          lapply(sort(unique(d0$source_id)), function(s) list(scheme = "source", test = d0$source_id == s)))
res <- list()
cv <- list()
param <- list()
for (v in names(variants)) {
  a <- fit_variant(d0, v)
  pr <- predict_tables(a$tabs, d0, "typical", a$levels)
  res[[v]] <- data.frame(variant = v, species = d0$species, source = d0$source_id, y = d0$y,
                         eta = pr$eta, stringsAsFactors = FALSE)
  vc <- as.data.frame(lme4::VarCorr(a$fit))
  sdc <- function(g) { s <- vc$sdcor[vc$grp == g]; if (length(s)) round(s, 2) else NA }
  param[[v]] <- data.frame(variant = v, sd_source = sdc("source_id"), sd_population = sdc("population"),
                           sd_residual = sdc("Residual"), sd_genus = sdc("g_genus"),
                           sd_family = sdc("g_family"), sd_order = sdc("g_order"))
  cv[[v]] <- do.call(rbind, lapply(jobs, function(j) {
    b <- tryCatch(fit_variant(d0[!j$test, ], v), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    p2 <- tryCatch(predict_tables(b$tabs, d0[j$test, ], "typical", b$levels), error = function(e) NULL)
    if (is.null(p2)) return(NULL)
    data.frame(variant = v, scheme = j$scheme, source = d0$source_id[j$test], y = d0$y[j$test],
               eta = p2$eta, vt = p2$v_total, stringsAsFactors = FALSE)
  }))
}
res <- do.call(rbind, res)
cv <- do.call(rbind, cv)
param <- do.call(rbind, param)

# Species measured by a single source: prediction of their node / mean of their records.
n_src <- tapply(d0$source_id, d0$species, function(s) length(unique(s)))
single <- res[res$species %in% names(n_src)[n_src == 1], ]
sp <- aggregate(cbind(y, eta) ~ variant + source + species, single, mean)
sp$ratio <- exp(sp$eta - sp$y)
ratio <- aggregate(ratio ~ variant + source, sp, stats::median)

z <- (cv$y - cv$eta) / sqrt(cv$vt)
cv$inside <- abs(z) <= stats::qnorm(0.975)
cv$factor <- exp(stats::qnorm(0.975) * sqrt(cv$vt))
cv$res <- cv$y - cv$eta
summ <- function(x) data.frame(n = nrow(x), cov95 = round(100 * mean(x$inside)),
                               factor95 = round(stats::median(x$factor), 1),
                               bias = round(mean(x$res), 2))
overall <- do.call(rbind, lapply(split(cv, list(cv$variant, cv$scheme), drop = TRUE), function(x)
  cbind(x[1, c("variant", "scheme")], summ(x))))
by_source <- do.call(rbind, lapply(split(cv, list(cv$variant, cv$scheme, cv$source), drop = TRUE), function(x)
  cbind(x[1, c("variant", "scheme", "source")], summ(x))))
overall <- overall[order(overall$scheme, overall$variant), ]
by_source <- by_source[order(by_source$scheme, by_source$source, by_source$variant), ]

md_table <- function(d) {
  c(paste0("| ", paste(names(d), collapse = " | "), " |"),
    paste0("|", paste(rep("---|", ncol(d)), collapse = "")),
    apply(d, 1, function(r) paste0("| ", paste(trimws(format(r)), collapse = " | "), " |")))
}
out <- c("# Source effect in the model of fish without lengths", "",
         paste0("Generated by `dev/validation/source_effect_fish.R` on ", format(Sys.Date()),
                " with the data of `data-raw/derived/model_typical.csv` (", nrow(d0), " records of ",
                length(unique(d0$species)), " species; sources: ",
                paste(sprintf("%s %d", names(table(d0$source_id)), as.integer(table(d0$source_id))),
                      collapse = ", "), ")."), "",
         "Variants: A, random effect of source; N, without it. Both with the population effect (the samplings of a species at a NEON site).", "",
         "## Components (standard deviations, ln)", "", md_table(param), "",
         "## Species of a single source: prediction / mean of their records (median)", "",
         md_table(stats::reshape(ratio, idvar = "source", timevar = "variant", direction = "wide")), "",
         "## Cross-validation", "",
         "`population`: 10 groups of populations (all the samplings of a species at a NEON site go together); `species`: 10 groups of species; `source`: each source left out in turn. `cov95`: % of values inside the 95 % interval; `factor95`: median ×/÷ factor of the interval; `bias`: mean of y − η (ln; negative if the model predicts above).", "",
         md_table(overall), "", "### By source", "", md_table(by_source))
writeLines(out, "dev/validation/source_effect_fish.md", useBytes = TRUE)
cat(out, sep = "\n")
