# Source effect in the model of invertebrates without lengths (outside the check)
#
# Compares four ways of handling the sources in invertebrate_individual:
#  - A: random effect of source crossed with the whole hierarchy (version 0.2.0);
#  - C: source within group of taxa (rotifers, microcrustaceans and the other
#    invertebrates);
#  - N: no source effect;
#  - M: no source effect, with the fixed contrast of the benthic records of
#    Moreau et al. (2021) (`moreau_benthos`, model_data.R): version 0.3.0,
#    decisions of 2026-11-05.
# A, C and N are fitted without that contrast.
# For each one: the median, by source, of the ratio between the prediction and
# the mean of the records of the species measured by a single source, and a
# quick cross-validation (10 groups of records and each source left out in
# turn) with the coverage of the 95 % interval and the mean factor
# data/prediction by group of taxa.
# Writes dev/validation/source_effect_invertebrates.md.
# From the package root: Rscript dev/validation/source_effect_invertebrates.R  (about 2 min)

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
dat <- read_model_data()
m <- "invertebrate_individual"
d0 <- model_rows(dat, m, "typical")
group <- ifelse(d0$phylum %in% "Rotifera", "rotifers",
                ifelse(d0$class %in% c("Branchiopoda", "Copepoda", "Ostracoda"),
                       "microcrustaceans", "other invertebrates"))
variants <- c(A = "crossed", C = "within group", N = "no source",
              M = "no source, with the contrast of Moreau in the benthos")
without_moreau <- function(d, v) {
  if (v != "M") d$contrast[d$contrast %in% "moreau_benthos"] <- NA_character_
  d
}

fit_variant <- function(d, g, v) {
  d <- without_moreau(d, v)
  levels <- split_levels(d)
  opt <- model_options(d, "typical")
  opt$source <- v %in% c("A", "C")
  f <- model_frame(d, "typical", levels)
  if (v == "C") f$source_id <- paste(f$source_id, g, sep = "|")
  fit <- suppressMessages(suppressWarnings(fit_model(f, "typical", levels, opt)))
  list(tabs = model_tables(fit, f, m, "typical", levels), levels = levels, fit = fit)
}

set.seed(20261105)
u <- unique(d0$record_id)
fold <- stats::setNames(sample(rep(1:10, length.out = length(u))), u)
res <- list()
cv <- list()
param <- list()
for (v in names(variants)) {
  a <- fit_variant(d0, group, v)
  pr <- predict_tables(a$tabs, without_moreau(d0, v), "typical", a$levels)
  res[[v]] <- data.frame(variant = v, species = d0$species, source = d0$source_id, y = d0$y,
                         eta = pr$eta, stringsAsFactors = FALSE)
  vc <- as.data.frame(lme4::VarCorr(a$fit))
  s <- vc$sdcor[vc$grp == "source_id"]
  fx <- lme4::fixef(a$fit)
  se <- sqrt(diag(as.matrix(stats::vcov(a$fit))))
  param[[v]] <- data.frame(variant = v, sd_source = if (length(s)) round(s, 2) else NA,
                           sd_residual = round(vc$sdcor[vc$grp == "Residual"], 2),
                           weighings = round(exp(fx[["c_rotifer_weighing"]]), 1),
                           moreau_benthos = if ("c_moreau_benthos" %in% names(fx))
                             sprintf("%.1f (%.1f-%.1f)", exp(fx[["c_moreau_benthos"]]),
                                     exp(fx[["c_moreau_benthos"]] - 1.96 * se[["c_moreau_benthos"]]),
                                     exp(fx[["c_moreau_benthos"]] + 1.96 * se[["c_moreau_benthos"]]))
                           else "—")
  jobs <- c(lapply(1:10, function(f) list(scheme = "unit", test = fold[d0$record_id] == f)),
            lapply(sort(unique(d0$source_id)), function(s) list(scheme = "source", test = d0$source_id == s)))
  cv[[v]] <- do.call(rbind, lapply(jobs, function(j) {
    b <- tryCatch(fit_variant(d0[!j$test, ], group[!j$test], v), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    p2 <- tryCatch(predict_tables(b$tabs, without_moreau(d0[j$test, ], v), "typical", b$levels),
                   error = function(e) NULL)
    if (is.null(p2)) return(NULL)
    data.frame(variant = v, scheme = j$scheme, group = group[j$test], y = d0$y[j$test],
               eta = p2$eta, vt = p2$v_total, stringsAsFactors = FALSE)
  }))
}
res <- do.call(rbind, res)
cv <- do.call(rbind, cv)
param <- do.call(rbind, param)

# Species of a single source (without those not identified to species).
res <- res[!grepl("not identified", res$species, fixed = TRUE), ]
nsrc <- tapply(res$source, res$species, function(s) length(unique(s)))
r1 <- res[nsrc[res$species] == 1, ]
sm <- stats::aggregate(cbind(y, eta) ~ variant + species + source, r1, mean)
sm$ratio <- exp(sm$eta - sm$y)
t1 <- stats::aggregate(ratio ~ source + variant, sm, stats::median)
t1 <- stats::reshape(t1, idvar = "source", timevar = "variant", direction = "wide")
n_sp <- table(sm$source[sm$variant == "A"])
t1 <- data.frame(source = t1$source, species = as.integer(n_sp[t1$source]),
                 A = sprintf("%.2f", t1$ratio.A), C = sprintf("%.2f", t1$ratio.C),
                 N = sprintf("%.2f", t1$ratio.N), M = sprintf("%.2f", t1$ratio.M))

cv$z <- (cv$y - cv$eta) / sqrt(cv$vt)
cv <- cv[is.finite(cv$z), ]
t2 <- do.call(rbind, lapply(split(cv, list(cv$scheme, cv$group, cv$variant), drop = TRUE), function(x) {
  data.frame(scheme = x$scheme[1], group = x$group[1], variant = x$variant[1], n = nrow(x),
             cov95 = sprintf("%.0f %%", 100 * mean(abs(x$z) < 1.96)),
             data_prediction = sprintf("%.2f", exp(mean(x$y - x$eta))))
}))
t2 <- t2[order(t2$scheme, t2$group, t2$variant), ]
t3 <- do.call(rbind, lapply(split(cv, list(cv$scheme, cv$variant), drop = TRUE), function(x) {
  data.frame(scheme = x$scheme[1], variant = x$variant[1], n = nrow(x),
             cov95 = sprintf("%.0f %%", 100 * mean(abs(x$z) < 1.96)),
             factor95 = sprintf("%.1f", stats::median(exp(1.96 * sqrt(x$vt)))))
}))
t3 <- t3[order(t3$scheme, t3$variant), ]

md_table <- function(x) {
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep("---|", ncol(x)), collapse = "")),
    apply(x, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |")))
}
lines <- c("# Source effect in the model of invertebrates without lengths", "",
           paste0("Generated by `dev/validation/source_effect_invertebrates.R` on ", format(Sys.Date()), " (",
                  R.version.string, ", lme4 ", as.character(utils::packageVersion("lme4")),
                  "; seed 20261105)."), "",
           "Variants: `A`, crossed source effect (version 0.2.0); `C`, source within group of taxa (rotifers, microcrustaceans and the rest); `N`, no source effect; `M`, no source effect and a fixed contrast for the benthic records of Moreau et al. (2021) (version 0.3.0). A, C and N are fitted without that contrast. The predictions of the tables carry the contrasts of each record (weighings of rotifers, stages of copepods and, in M, that of Moreau).", "",
           "## Parameters", "",
           "`weighings`: factor of the contrast of the weighings of rotifers against the masses from volumes. `moreau_benthos`: factor (95 % CI) of the contrast of the benthic records of Moreau against the population means of Hrycik; in M, a benthic species measured only by Moreau is predicted at its mass divided by that factor.", "",
           md_table(param), "",
           "## Prediction against the data of each species", "",
           "Median, by source, of prediction/mean of the records in the species measured by a single source (1 = the prediction matches their data).", "",
           md_table(t1), "",
           "## Cross-validation by group", "",
           "`unit`: 10 groups of records; `source`: each source left out in turn. `data_prediction`: geometric mean of data/prediction.", "",
           md_table(t2), "",
           "## Overall cross-validation", "",
           md_table(t3), "",
           "## Reading", "",
           "- With a source effect (A and C), the prediction is for an average source, and the species of a single source shift towards it: the rotifers of Ejsmont-Karabin and Bottrell end up 2-3 times above their data, and the species of Moreau 3-5 times below. C does not fix it: the differences lie within each group (weighings against volumes in rotifers; Moreau against Hrycik in the benthos), and in C the microcrustaceans end up above theirs.",
           "- Without a source effect (N), each species stays at its data, the contrast of the weighings takes up the difference of method of the rotifers, the coverage is kept and the intervals are narrower. The difference between Moreau and Hrycik in the benthos is left uncorrected (lower coverage of the `source` scheme in that group).",
           "- Moreau et al. (2021) compile masses of individual prey from studies of the diet of birds and fish and of calorimetry (mostly Driver et al. 1974 and Driver 1981), often of large or adult animals; Hrycik et al. (2024) give means of all the individuals of field samples, which is what a biomass needs. M estimates that difference with the families they share (few, hence its interval) and predicts for population means. With Hrycik left out (`source` scheme), M cannot calibrate the contrast and its benthic predictions are of no use: it depends on that source.")
writeLines(lines, "dev/validation/source_effect_invertebrates.md", useBytes = TRUE)
cat(lines, sep = "\n")
