# Coverage of a list of taxa by the reference ----------------------------------

#' At which taxonomic level each taxon will be estimated
#'
#' `check_coverage()` goes through the taxon units of a survey (and their
#' lengths) before any abundance is used and says, for each one, from which
#' node of the reference its mass would be predicted, and how wide the
#' prediction interval of the mean individual mass of a population would be.
#' Use it to see which taxa rest on their own data and which on relatives,
#' and where measuring lengths would pay off.
#'
#' @inheritParams estimate_biomass
#' @param level Level of the intervals whose width is reported (default 0.95).
#'
#' @return A data frame of class `bmassa_coverage` with one row per taxon unit
#'   of `taxa` (from its length, if it has one, or without it; phytoplankton
#'   without `count_unit`, twice: per cell and per natural unit) and one per
#'   row of `lengths`. Columns: `taxon_id`, `broad_group`, `scientific_name` (if
#'   in `taxa`), `method` (`"typical"` or `"length"`), `count_unit`, `site`
#'   (and `date`, if `lengths` has it) of the rows of `lengths`, `life_stage`,
#'   `sex`, `length_mm`, `length_type`, `model`, `stage_used`, `match_level`,
#'   `match_name`, `n_reference`, `factor` (the interval is the median divided
#'   and multiplied by this factor), `level`, `status` (`"estimable"` or
#'   `"not_estimable"`), `reason` and `note`. Its `print()` method summarises
#'   the levels by broad group.
#'
#' @details
#' Each unit is checked as [estimate_biomass()] would estimate it: with its
#' life stage, sex and counting unit, insects without a stage as larvae, and
#' amphibians without a length estimable only as adults. The width does not
#' depend on the back-transformation chosen in [estimate_biomass()], which
#' only shifts the estimate.
#'
#' @export
#' @examples
#' cov <- check_coverage(example_survey$taxa, example_survey$lengths)
#' cov
#' as.data.frame(cov)[, c("taxon_id", "method", "match_level", "factor")]
check_coverage <- function(taxa, lengths = NULL, level = 0.95, reference = get_reference()) {
  check_data_frame(taxa)
  if (!is.null(lengths)) check_data_frame(lengths)
  check_level(level)
  check_reference(reference)
  tx <- check_taxa_table(taxa, require_count_unit = FALSE)
  ln <- if (is.null(lengths)) NULL else check_lengths_table(lengths, tx)
  dom <- unname(broad_group_domain[tx$broad_group])
  nt <- nrow(tx)
  base <- data.frame(ti = seq_len(nt), count_unit = tx$count_unit, site = rep(NA_character_, nt),
                     date = rep(NA_character_, nt), length_mm = rep(NA_real_, nt),
                     length_type = rep(NA_character_, nt), length_cv = rep(NA_real_, nt),
                     stringsAsFactors = FALSE)
  has <- if ("length_mean" %in% names(tx)) which(!is.na(tx$length_mean)) else integer(0)
  if (length(has) > 0L) {
    v <- length_values(tx, has)
    base$length_mm[has] <- v$mm
    base$length_type[has] <- v$type
    base$length_cv[has] <- v$cv
  }
  # Phytoplankton without a counting unit: both, per cell and per natural unit.
  open_cu <- which(is.na(base$count_unit))
  if (length(open_cu) > 0L) {
    nu <- base[open_cu, , drop = FALSE]
    base$count_unit[open_cu] <- "cell"
    nu$count_unit <- "natural_unit"
    base <- rbind(base, nu)
  }
  if (!is.null(ln) && nrow(ln) > 0L) {
    v <- length_values(ln, seq_len(nrow(ln)))
    li <- match(as_text(ln$taxon_id), tx$taxon_id)
    lrows <- data.frame(ti = li, count_unit = tx$count_unit[li], site = as_text(ln$site),
                        date = if ("date" %in% names(ln)) as_text(ln$date) else NA_character_,
                        length_mm = v$mm, length_type = v$type, length_cv = v$cv,
                        stringsAsFactors = FALSE)
    base <- rbind(base, lrows)
  }
  base$count_unit[is.na(base$count_unit)] <- "cell"
  base <- base[order(base$ti, !is.na(base$site), base$count_unit, method = "radix"), ]
  stage <- unit_stage(tx$life_stage, tx$sex)
  rows <- data.frame(row = seq_len(nrow(base)), taxon_id = tx$taxon_id[base$ti], ti = base$ti,
                     site = base$site, domain = dom[base$ti], count_unit = base$count_unit,
                     cells_per_unit = ifelse(base$count_unit %in% "natural_unit",
                                             tx$cells_per_unit[base$ti], NA_real_),
                     stage = stage[base$ti], length_mm = base$length_mm,
                     length_type = base$length_type, length_cv = base$length_cv,
                     stringsAsFactors = FALSE)
  pr <- predict_rows(rows, tx, reference, "none", extrapolate = FALSE)$pred
  z <- stats::qnorm(1 - (1 - level) / 2)
  tau2 <- pr$v_param + pr$v_missing + pr$v_population + pr$v_source
  ok <- pr$status == "estimated"
  out <- data.frame(taxon_id = rows$taxon_id, broad_group = tx$broad_group[rows$ti],
                    stringsAsFactors = FALSE)
  if ("scientific_name" %in% names(tx)) out$scientific_name <- as_text(tx$scientific_name)[rows$ti]
  out$method <- pr$method
  out$count_unit <- rows$count_unit
  out$site <- rows$site
  if (!is.null(ln) && "date" %in% names(ln)) out$date <- base$date
  out$life_stage <- tx$life_stage[rows$ti]
  out$sex <- tx$sex[rows$ti]
  out$length_mm <- rows$length_mm
  out$length_type <- rows$length_type
  out$model <- pr$model
  out$stage_used <- pr$stage_used
  out$match_level <- pr$match_level
  out$match_name <- pr$match_name
  out$n_reference <- pr$n_reference
  out$factor <- ifelse(ok, exp(z * sqrt(tau2)), NA_real_)
  out$level <- rep(level, nrow(out))
  out$status <- ifelse(ok, "estimable", "not_estimable")
  out$reason <- pr$reason
  out$note <- pr$note
  rownames(out) <- NULL
  class(out) <- c("bmassa_coverage", "data.frame")
  out
}

#' @export
print.bmassa_coverage <- function(x, ...) {
  d <- as.data.frame(x)
  cat(sprintf("<bmassa_coverage> %d taxa, %d rows (%d without lengths, %d with lengths)\n",
              length(unique(d$taxon_id)), nrow(d), sum(d$method %in% "typical"),
              sum(d$method %in% "length")))
  lv <- c(rev(main_ranks), "root")
  for (g in broad_group_levels[broad_group_levels %in% d$broad_group]) {
    s <- d[d$broad_group == g, ]
    for (m in c("typical", "length")) {
      sm <- s[s$method %in% m, ]
      if (nrow(sm) == 0L) next
      ok <- sm$status == "estimable"
      tab <- table(factor(sm$match_level[ok], levels = lv))
      tab <- tab[tab > 0]
      cat(sprintf("  %-18s %-7s %s%s; median factor %s\n", g, m,
                  if (length(tab)) paste(sprintf("%s %d", names(tab), as.integer(tab)), collapse = ", ") else "none",
                  if (any(!ok)) sprintf(", not estimable %d", sum(!ok)) else "",
                  if (any(ok)) format(signif(stats::median(sm$factor[ok]), 2)) else "-"))
    }
  }
  cat(sprintf("Factor: the %g %% prediction interval is the median divided and multiplied by it.\n",
              100 * d$level[1]))
  invisible(x)
}

#' @export
as.data.frame.bmassa_coverage <- function(x, ...) {
  class(x) <- "data.frame"
  x
}
