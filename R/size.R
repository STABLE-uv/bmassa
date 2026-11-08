# Size structure of the community by group ----------------------------------------

size_levels <- c("mass", "length")

size_summary_columns <- c(
  "per", "size", "size_type", "size_unit", "abundance", "share_excluded", "share_measured",
  "share_lognormal", "share_taxon_mean", "mean_size", "mean_size_lower", "mean_size_upper",
  "median_size", "median_size_lower", "median_size_upper", "approximation", "level", "interval",
  "n_rows", "n_sized", "n_missing", "status", "n_draws"
)

size_class_columns <- c(
  "class", "lower", "upper", "share_individuals", "share_individuals_lower",
  "share_individuals_upper", "share_biomass", "share_biomass_lower", "share_biomass_upper"
)

#' Size structure of the community by sample or any other grouping
#'
#' `summarize_size()` describes the sizes of the individuals of each group of
#' rows of an estimate (each site, by default): the mean and the median size
#' of the individuals, weighted by their abundances, and, with `breaks`, the
#' share of the individuals and of the biomass in each size class. The size is
#' the mass of an individual, in the currency and unit of the estimate, or its
#' length. It uses the lengths of measured individuals (or length classes)
#' where `lengths` gives them; where only the mean size of a taxon is known,
#' the indicators that depend on the spread of sizes are flagged as
#' approximations.
#'
#' @param data Result of [estimate_biomass()] or [convert_currency()].
#' @param by Character vector with the names of the grouping columns, e.g.
#'   `"site"` (default) or `c("site", "date")`. Abundances of different units
#'   (`per`) are never pooled: each group has one row per unit (put all the
#'   abundances on one basis with [convert_abundance()] to pool them).
#' @param size `"mass"` (default), the mass of each individual in the currency
#'   and unit of `data`, or `"length"`, its length in mm, only for the rows
#'   estimated from lengths. Lengths of different types (body length, fork
#'   length...) are not pooled: a group with more than one is an error.
#' @param breaks `NULL` (default: no classes) or the increasing limits of the
#'   size classes, 0 or more, in the unit of the sizes (that of `data` for
#'   masses, mm for lengths). Classes are closed on the left and open on the
#'   right, \eqn{[a, b)}; the sizes below the first limit and from the last one
#'   upwards are classes of their own, so that the shares of each group add up
#'   to 1.
#' @param na_rm If `FALSE` (default), a group with rows whose abundance is
#'   missing gets `NA` indicators and `status = "incomplete"`: the composition
#'   of the group is not known. If `TRUE`, those rows are left out and
#'   `n_missing` says how many.
#' @inheritParams summarize_biomass
#'
#' @return An object of class `bmassa_size`, a list with:
#'   * `summary`: one row per combination of the `by` columns and the unit of
#'     the abundances, sorted by them (C locale order), with columns `by...`,
#'     `per` (the unit of the abundances), `size` (`"mass"` or `"length"`),
#'     `size_type` (the currency of the masses, or the length type),
#'     `size_unit`, `abundance` (individuals of the group, per `per`),
#'     `share_excluded` (share of them without a size: rows not estimated
#'     and, for lengths, rows without lengths), `share_measured`,
#'     `share_lognormal` and `share_taxon_mean` (shares of the individuals with
#'     a size whose sizes come from the lengths of measured individuals, from a
#'     mean length with its spread, or from the mean size of their taxon unit:
#'     they add up to 1), `mean_size`, `mean_size_lower`, `mean_size_upper`
#'     (mean size of the individuals and prediction interval), `median_size`,
#'     `median_size_lower`, `median_size_upper`, `approximation` (`TRUE` when
#'     some individuals have the mean size of their taxon: the median and the
#'     shares of the classes are then approximations), `level`, `interval`
#'     (`"prediction"`), `n_rows`, `n_sized` (rows with a size), `n_missing`
#'     (rows without an abundance), `status` (`"complete"` or
#'     `"incomplete"`) and `n_draws` (`NA` when nothing was drawn).
#'   * `classes`: `NULL` without `breaks`; otherwise one row per group and
#'     class, with the columns of the group, `per`, `size`, `size_type`,
#'     `size_unit`, `class` (label), `lower`, `upper` (limits),
#'     `share_individuals`, `share_individuals_lower`,
#'     `share_individuals_upper`, `share_biomass`, `share_biomass_lower`,
#'     `share_biomass_upper` (shares of the individuals with a size and of their
#'     biomass: mean of the draws and prediction interval), `approximation` and
#'     `level`.
#'   * `settings`: the arguments used and the class limits.
#'
#'   `as.data.frame()` returns `summary`. A zero-row input gives zero-row
#'   tables.
#'
#' @details
#' **Sizes of the individuals.** Each row with an estimated mass and a known
#' abundance contributes its individuals with the sizes that its data give:
#' * the lengths of measured individuals (`length_basis = "individuals"`):
#'   each length with its share of the individuals of the row (`n` in
#'   `lengths`), and the mass of the length-mass relation of the row at that
#'   length;
#' * a mean length with its spread (`length_sd` or `length_cv`): lengths
#'   lognormal with that mean and coefficient of variation, so masses
#'   lognormal with a standard deviation of the log mass \eqn{|\beta|} times
#'   that of the log length (\eqn{\beta}, the slope of the relation of the
#'   row);
#' * a mean length without spread, or a mass estimated without lengths: all
#'   the individuals of the row at that length, or at the mean mass of the
#'   population (`share_taxon_mean`).
#'
#' The mass of an individual is that of the relation at its length, without
#' the scatter of masses between individuals of the same length. With
#' `bias_correction = "lognormal"` in [estimate_biomass()] these are mean
#' masses at each length, and the masses of the individuals of a row add up
#' to its biomass.
#'
#' **Indicators.** `mean_size` is the mean size of the individuals with a
#' size, weighted by their abundances; for masses, the biomass of the group
#' divided by its individuals with a mass. `median_size` is the smallest size
#' at which the cumulative share of the individuals reaches one half. The
#' share of a class is that of the individuals with a size (or of their
#' biomass) between its limits. When individuals have the mean size of their
#' taxon unit (`share_taxon_mean` above 0), the median and the shares of the
#' classes are approximations based on those means (`approximation = TRUE`):
#' the spread of sizes within those taxa is unknown, and their individuals all
#' fall in the class of their mean. The mean is not affected.
#'
#' **Intervals.** The mean individual masses of the rows are drawn as in
#' [summarize_biomass()] (the draws of [simulate()][simulate.bmassa_estimate()]
#' with the same `seed` and `nsim = n_draws`, populations shared); in each draw
#' the masses of the individuals of a row move with the mass of the row, and
#' the indicators are computed. The estimate of the mean and the median size
#' is the median of their draws, and that of a share the mean of its draws,
#' so that the shares of a group add up to 1 (a share can be bimodal: a taxon
#' at its mean size falls in one class or another depending on the draw); the
#' limits are quantiles of the draws. The lengths and the abundances are taken
#' as known: with `size = "length"`, the mean and median lengths and the
#' shares of individuals by class have no uncertainty (the limits equal the
#' estimate), and only the shares of biomass are drawn.
#'
#' @section Missing values and zeros:
#' Rows whose mass was not estimated (and, for lengths, rows without lengths)
#' have no size: their individuals are counted in `share_excluded` and left
#' out of the indicators. Rows with a missing abundance make the group
#' incomplete (see `na_rm`). Rows with a zero abundance add nothing. A group
#' without individuals with a size has `NA` indicators.
#'
#' @section Reproducibility:
#' With the same `seed`, the draws are the same on the same platform.
#'
#' @seealso [summarize_biomass()] for the biomass by group.
#' @export
#' @examples
#' est <- estimate_biomass(example_survey$abundance, example_survey$taxa,
#'                         example_survey$lengths, bias_correction = "lognormal")
#' # Mean and median wet mass (mg) of the individuals of each group, by site.
#' sz <- summarize_size(est, by = c("site", "broad_group"), seed = 1)
#' sz
#' as.data.frame(sz)[, c("site", "broad_group", "mean_size", "median_size",
#'                       "share_excluded", "share_taxon_mean")]
#'
#' # Lengths of the macroinvertebrates of L1, with classes of 4 mm.
#' l1 <- summarize_size(est, by = c("site", "broad_group"), size = "length",
#'                      breaks = seq(0, 20, by = 4), seed = 1)
#' cl <- l1$classes
#' cl[cl$site == "L1" & cl$broad_group == "macroinvertebrate",
#'    c("class", "share_individuals", "share_biomass", "approximation")]
summarize_size <- function(data, ...) {
  UseMethod("summarize_size")
}

#' @rdname summarize_size
#' @export
summarize_size.default <- function(data, ...) {
  abort_input(sprintf(
    "`data` must be the result of estimate_biomass() or convert_currency(), not %s.",
    describe_value(data)
  ))
}

#' @rdname summarize_size
#' @export
summarize_size.bmassa_estimate <- function(data, by = "site", size = "mass", breaks = NULL,
                                           level = 0.95, n_draws = 1000, na_rm = FALSE,
                                           seed = NULL, ...) {
  check_dots_empty(...)
  check_by(by)
  check_choice(size, size_levels)
  breaks <- check_breaks(breaks)
  check_level(level)
  check_n_draws(n_draws)
  check_flag(na_rm)
  check_seed(seed)
  e <- data$estimates
  check_columns(e, by, arg = "data$estimates")
  clash <- intersect(by, c(size_summary_columns, size_class_columns))
  if (length(clash) > 0L) {
    abort_input(sprintf("`by` cannot contain %s, which are columns of the result.",
                        format_items(quote_items(clash), max = Inf)))
  }
  bounds <- if (is.null(breaks)) NULL else
    unique(c(if (breaks[1] > 0) 0, breaks, if (is.finite(breaks[length(breaks)])) Inf))
  n <- nrow(e)
  known <- !is.na(e$abundance)
  has_mass <- !is.na(e$ind_median)
  sized <- has_mass & known & (size == "mass" | e$method %in% "length")
  # One currency and unit for the masses (always, from version 0.5.0).
  cu <- unique(e[has_mass, c("currency", "unit"), drop = FALSE])
  if (nrow(cu) > 1L) {
    abort_input(sprintf(paste(
      "`data` has masses in more than one currency or unit (%s): convert them to one with",
      "convert_currency()."), format_items(quote_items(paste(cu$currency, cu$unit)))))
  }
  # Groups (C locale order) and, within them, the unit of the abundances.
  per <- as_text(e$per)
  key <- if (n) paste(make_key(e, by), per, sep = "\r") else character(0)
  ukey <- unique(key)
  first <- match(ukey, key)
  ord <- do.call(order, c(unname(lapply(e[first, by, drop = FALSE], as_text)),
                          list(per[first], method = "radix")))
  ukey <- ukey[ord]
  first <- first[ord]
  gid <- match(key, ukey)
  ng <- length(ukey)
  size_unit <- if (size == "length") "mm" else if (nrow(cu)) cu$unit else NA_character_
  size_type <- rep(if (size == "mass" && nrow(cu)) cu$currency else NA_character_, ng)
  if (size == "length") {
    for (g in seq_len(ng)) {
      lt <- unique(as_text(e$length_type[sized & gid == g]))
      if (length(lt) > 1L) {
        abort_input(sprintf(paste(
          "The group %s has lengths of more than one type (%s), which are not comparable:",
          "add a column that separates them to `by` (such as \"broad_group\" or",
          "\"length_type\"), or subset the estimate."),
          group_label(e[first[g], by, drop = FALSE], per[first[g]]),
          format_items(quote_items(lt))))
      }
      if (length(lt) == 1L) size_type[g] <- lt
    }
  }
  # Draws of the masses of the rows, as in simulate(): needed for masses and
  # for the shares of biomass of length classes.
  use_draws <- size == "mass" || !is.null(bounds)
  drawn <- which(!is.na(e$ind_mean))
  ld <- NULL
  if (use_draws && length(drawn) > 0L) {
    ld <- if (is.null(seed)) {
      draw_log_mass(data, drawn, n_draws, "shared", NULL)
    } else {
      withr::with_seed(seed, draw_log_mass(data, drawn, n_draws, "shared", NULL))
    }
  }
  p <- c((1 - level) / 2, 1 - (1 - level) / 2)
  contributing <- which(sized & e$abundance > 0)
  comp <- size_components(data, contributing, size)
  nb <- if (is.null(bounds)) 0L else length(bounds) - 1L
  num <- rep(NA_real_, ng)
  res <- data.frame(abundance = num, share_excluded = num, share_measured = num,
                    share_lognormal = num, share_taxon_mean = num, mean_size = num,
                    mean_size_lower = num, mean_size_upper = num, median_size = num,
                    median_size_lower = num, median_size_upper = num,
                    approximation = rep(NA, ng), n_rows = tabulate(gid, nbins = ng),
                    n_sized = tabulate(gid[sized], nbins = ng),
                    n_missing = tabulate(gid[!known], nbins = ng))
  cls <- array(NA_real_, dim = c(ng, max(nb, 1L), 6L))
  levels_g <- seq_len(ng)
  known_of <- split(which(known), factor(gid[known], levels = levels_g))
  rows_of <- split(contributing, factor(gid[contributing], levels = levels_g))
  comp_of <- split(seq_len(nrow(comp)), factor(gid[comp$row], levels = levels_g))
  for (g in seq_len(ng)) {
    if (res$n_missing[g] > 0L && !na_rm) next
    res$abundance[g] <- sum(e$abundance[known_of[[g]]])
    rg <- rows_of[[g]]
    cg <- comp[comp_of[[g]], , drop = FALSE]
    w_rows <- e$abundance[rg]
    total <- sum(w_rows)
    if (res$abundance[g] > 0) res$share_excluded[g] <- 1 - total / res$abundance[g]
    if (length(rg) == 0L || total <= 0) next
    a_c <- e$abundance[cg$row] * cg$share
    for (k in c("measured", "lognormal", "taxon_mean")) {
      res[[paste0("share_", k)]][g] <- sum(a_c[cg$kind == k]) / total
    }
    res$approximation[g] <- any(cg$kind == "taxon_mean" & a_c > 0)
    # Log shifts of the masses of the rows in each draw (columns), or none.
    shift <- if (is.null(ld)) NULL else
      t(ld[, match(cg$row, drawn), drop = FALSE]) - log(e$ind_median[cg$row])
    out <- size_indicators(cg, a_c, shift, size, bounds, log(e$ind_median[cg$row]), p)
    res[g, c("mean_size", "mean_size_lower", "mean_size_upper")] <- out$mean
    res[g, c("median_size", "median_size_lower", "median_size_upper")] <- out$median
    if (nb > 0L) cls[g, , ] <- out$classes
  }
  status <- ifelse(res$n_missing > 0L, "incomplete", "complete")
  grp <- e[first, by, drop = FALSE]
  rownames(grp) <- NULL
  has_est <- !is.na(res$mean_size)
  summary <- data.frame(grp, per = per[first], size = rep(size, ng), size_type = size_type,
                        size_unit = rep(size_unit, ng), res[c("abundance", "share_excluded",
                                                              "share_measured", "share_lognormal",
                                                              "share_taxon_mean", "mean_size",
                                                              "mean_size_lower", "mean_size_upper",
                                                              "median_size", "median_size_lower",
                                                              "median_size_upper",
                                                              "approximation")],
                        level = rep(level, ng),
                        interval = ifelse(has_est, "prediction", NA_character_),
                        res[c("n_rows", "n_sized", "n_missing")], status = status,
                        n_draws = rep(if (use_draws) as.integer(n_draws) else NA_integer_, ng),
                        stringsAsFactors = FALSE)
  rownames(summary) <- NULL
  classes <- NULL
  if (nb > 0L) {
    lab <- format_class_labels(bounds)
    gi <- rep(seq_len(ng), each = nb)
    ci <- rep(seq_len(nb), times = ng)
    v <- function(k) as.vector(t(matrix(cls[, , k], nrow = ng, ncol = nb)))
    classes <- data.frame(grp[gi, , drop = FALSE], per = per[first][gi], size = rep(size, ng * nb),
                          size_type = size_type[gi], size_unit = rep(size_unit, ng * nb),
                          class = lab[ci], lower = bounds[ci], upper = bounds[ci + 1L],
                          share_individuals = v(1), share_individuals_lower = v(2),
                          share_individuals_upper = v(3), share_biomass = v(4),
                          share_biomass_lower = v(5), share_biomass_upper = v(6),
                          approximation = summary$approximation[gi], level = rep(level, ng * nb),
                          stringsAsFactors = FALSE)
    rownames(classes) <- NULL
  }
  structure(list(summary = summary, classes = classes,
                 settings = list(by = by, size = size, size_unit = size_unit, breaks = breaks,
                                 bounds = bounds, level = level,
                                 n_draws = if (use_draws) as.integer(n_draws) else NA_integer_,
                                 na_rm = na_rm, seed = seed,
                                 bias_correction = data$settings$bias_correction)),
            class = "bmassa_size")
}

#' @export
as.data.frame.bmassa_size <- function(x, ...) {
  x$summary
}

#' @export
print.bmassa_size <- function(x, ...) {
  s <- x$summary
  st <- x$settings
  what <- if (st$size == "mass") {
    sprintf("mass of the individuals (%s, %s)",
            if (nrow(s) && any(!is.na(s$size_type))) gsub("_", " ", s$size_type[!is.na(s$size_type)][1])
            else "no masses", st$size_unit)
  } else {
    "length of the individuals (mm)"
  }
  cat(sprintf("<bmassa_size> %s by %s: %d groups%s\n", what, paste(st$by, collapse = ", "),
              nrow(s), if (any(s$status == "incomplete")) sprintf(", %d incomplete",
                                                                  sum(s$status == "incomplete")) else ""))
  if (nrow(s) == 0L) {
    return(invisible(x))
  }
  ok <- !is.na(s$mean_size)
  if (any(ok)) {
    ex <- s$share_excluded[ok]
    cat(sprintf("Individuals without a size: %s of each group (median %s)\n",
                if (diff(range(ex)) < 0.005) format_percent(ex[1]) else
                  paste(format_percent(min(ex), FALSE), format_percent(max(ex)), sep = "-"),
                format_percent(stats::median(ex))))
    cat(sprintf("Approximate median and classes (individuals at the mean size of their taxon): %d of %d groups\n",
                sum(s$approximation[ok]), sum(ok)))
    cat(sprintf("Interval: prediction %g %%%s\n", 100 * st$level,
                if (is.na(st$n_draws)) " (lengths and abundances taken as known)" else
                  sprintf(", from %d draws", st$n_draws)))
  }
  if (!is.null(st$bounds)) {
    cat(sprintf("Classes: %s\n", paste(format_class_labels(st$bounds), collapse = " ")))
  }
  cols <- c(st$by, "per", "mean_size", "median_size", "share_excluded", "approximation")
  show <- s[seq_len(min(nrow(s), 10L)), cols, drop = FALSE]
  print(show, digits = 3, row.names = FALSE)
  if (nrow(s) > 10L) cat(sprintf("... %d more groups\n", nrow(s) - 10L))
  cat("Use as.data.frame() for all the indicators and x$classes for the size classes.\n")
  invisible(x)
}

# Components of the sizes of the individuals of the rows `r` (rows with a mass
# and a positive abundance): one per measured length of a row, or one per row
# (a lognormal from a mean length with its spread, or a single value). `mu`
# and `sigma` are on the log scale of the size (sigma 0 for a single value),
# `share` is the share of the individuals of the row, `tilt` the exponent of
# the size in the mass of an individual (1 for masses; the slope of the
# length-mass relation of the row for lengths), `logq` the log of the mass of
# its individuals relative to the mean mass of the row (measured lengths) and
# `kind` the origin of the sizes.
size_components <- function(data, r, size) {
  empty <- data.frame(row = integer(0), mu = numeric(0), sigma = numeric(0), share = numeric(0),
                      tilt = numeric(0), logq = numeric(0), kind = character(0),
                      stringsAsFactors = FALSE)
  if (length(r) == 0L) {
    return(empty)
  }
  e <- data$estimates
  sim <- data$simulation$rows
  ind <- data$simulation$lengths
  nodes <- data$reference$nodes
  cv <- data$simulation$conversion
  conv_slope <- if (is.null(cv)) rep(1, nrow(e)) else cv$slope
  lm <- log(e$ind_median[r])
  from_length <- sim$method[r] %in% "length"
  b <- rep(1, length(r))
  b[from_length] <- conv_slope[r[from_length]] * nodes$path_slope[sim$node[r[from_length]]]
  tilt <- if (size == "mass") rep(1, length(r)) else b
  measured <- from_length & sim$length_basis[r] %in% "individuals"
  lognormal <- from_length & !measured & !is.na(sim$length_cv[r])
  single <- !measured & !lognormal
  len <- log(sim$length_used[r])
  parts <- list(empty)
  if (any(single)) {
    parts[[2]] <- data.frame(row = r[single], mu = if (size == "mass") lm[single] else len[single],
                             sigma = 0, share = 1, tilt = tilt[single], logq = 0,
                             kind = "taxon_mean", stringsAsFactors = FALSE)
  }
  if (any(lognormal)) {
    sl <- sqrt(log1p(sim$length_cv[r[lognormal]]^2))
    sm <- if (size == "mass") abs(b[lognormal]) * sl else sl
    parts[[3]] <- data.frame(row = r[lognormal],
                             mu = (if (size == "mass") lm[lognormal] else len[lognormal]) - sm^2 / 2,
                             sigma = sm, share = 1, tilt = tilt[lognormal], logq = 0,
                             kind = "lognormal", stringsAsFactors = FALSE)
  }
  for (k in which(measured)) {
    i <- r[k]
    v <- ind[ind$group == sim$length_group[i], , drop = FALSE]
    logq <- b[k] * (log(v$length_mm) - len[k])
    parts[[length(parts) + 1L]] <- data.frame(
      row = i, mu = if (size == "mass") lm[k] + logq else log(v$length_mm), sigma = 0,
      share = v$n / sum(v$n), tilt = tilt[k], logq = logq, kind = "measured",
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, parts)
  out <- out[order(match(out$row, r), method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Indicators of one group from its components `cg`, with `a_c` individuals
# each: the mean and the median size and the shares of the classes between
# `bounds`, in every draw. `shift` holds the log shifts of the masses of the
# rows of the components in each draw (one column per draw; NULL: none), and
# `lm` the log mean mass of their rows. Returns the estimate (the median of
# the draws for the sizes, their mean for the shares) and the limits of each
# indicator.
size_indicators <- function(cg, a_c, shift, size, bounds, lm, p) {
  sig <- cg$sigma
  w <- a_c / sum(a_c)
  # Log locations of the sizes and log masses of the components in each draw.
  mu0 <- matrix(cg$mu, nrow = nrow(cg), ncol = 1L)
  mass_log <- if (size == "mass") {
    if (is.null(shift)) mu0 + sig^2 / 2 else cg$mu + sig^2 / 2 + shift
  } else {
    if (is.null(shift)) matrix(lm + cg$logq, ncol = 1L) else lm + cg$logq + shift
  }
  mu <- if (size == "mass" && !is.null(shift)) cg$mu + shift else mu0
  summarise <- function(x, centre = stats::median) {
    c(centre(x), stats::quantile(x, probs = p, names = FALSE))
  }
  mean_d <- colSums(w * exp(mu + sig^2 / 2))
  median_d <- exp(mixture_quantile(0.5, mu, sig, w))
  out <- list(mean = summarise(mean_d), median = summarise(median_d), classes = NULL)
  if (!is.null(bounds)) {
    nb <- length(bounds) - 1L
    bw <- a_c * exp(mass_log - max(mass_log))
    bw <- sweep(bw, 2L, colSums(bw), "/")
    f_ind <- lapply(bounds, function(x) colSums(w * below(x, mu, sig, 0)))
    # With lengths, the shares below a limit are the same in every draw.
    f_bio <- lapply(bounds, function(x) colSums(bw * as.vector(below(x, mu, sig, cg$tilt))))
    # Shares below 1e-12 (far tails of the lognormals, rounding) are 0.
    share <- function(hi, lo) {
      d <- hi - lo
      d[d < 1e-12] <- 0
      d
    }
    cl <- matrix(NA_real_, nrow = nb, ncol = 6L)
    for (j in seq_len(nb)) {
      cl[j, 1:3] <- summarise(share(f_ind[[j + 1L]], f_ind[[j]]), mean)
      cl[j, 4:6] <- summarise(share(f_bio[[j + 1L]], f_bio[[j]]), mean)
    }
    out$classes <- cl
  }
  out
}

# Share of each component below `x` (sizes < x), one column per draw: for
# single values, whether they are below it (a value within 1e-9 of a limit, on
# the log scale, is taken as at the limit); for lognormals, their distribution
# function, tilted by `tilt` (the share of the biomass when the mass is
# proportional to the size to the power `tilt`).
below <- function(x, mu, sigma, tilt) {
  lx <- log(x)
  z <- lx - mu
  out <- (z > 1e-9) * 1
  cont <- sigma > 0
  if (any(cont)) {
    tilt <- rep_len(tilt, length(sigma))
    out[cont, ] <- stats::pnorm(z[cont, , drop = FALSE] / sigma[cont] - tilt[cont] * sigma[cont])
  }
  out
}

# Smallest value with a share of at least `prob` of the individuals at or
# below it, for mixtures of single values (sigma 0) and lognormals with
# weights `w`, on the log scale; `mu` has one column per draw. Single values
# only: from the cumulative weights of the sorted values; otherwise bisection.
mixture_quantile <- function(prob, mu, sigma, w, iter = 50L) {
  target <- prob - 1e-10
  cont <- sigma > 0
  if (!any(cont)) {
    nr <- nrow(mu)
    if (nr == 1L) {
      return(mu[1L, ])
    }
    nd <- ncol(mu)
    # All the draws sorted at once (positions in `mu`), and the cumulative
    # weights down each column.
    o <- matrix(order(rep(seq_len(nd), each = nr), mu), nrow = nr)
    ws <- matrix(w[(o - 1L) %% nr + 1L], nrow = nr)
    cw <- if (nr <= 200L) lower.tri(diag(nr), diag = TRUE) %*% ws else apply(ws, 2L, cumsum)
    k <- pmin(colSums(cw < target) + 1L, nr)
    return(mu[o[cbind(k, seq_len(nd))]])
  }
  spread <- ifelse(cont, 10 * sigma, 0)
  lo <- apply(mu - spread, 2L, min) - 1
  hi <- apply(mu + spread, 2L, max) + 1
  cdf <- function(y) {
    z <- matrix(y, nrow = nrow(mu), ncol = ncol(mu), byrow = TRUE) - mu
    pr <- (z >= -1e-9) * 1
    if (any(cont)) pr[cont, ] <- stats::pnorm(z[cont, , drop = FALSE] / sigma[cont])
    colSums(w * pr)
  }
  for (k in seq_len(iter)) {
    mid <- (lo + hi) / 2
    up <- cdf(mid) >= target
    hi[up] <- mid[up]
    lo[!up] <- mid[!up]
  }
  hi
}

format_class_labels <- function(bounds) {
  f <- vapply(bounds, function(x) format(x, digits = 6, trim = TRUE), character(1))
  sprintf("[%s, %s)", f[-length(f)], f[-1L])
}

format_percent <- function(x, sign = TRUE) {
  paste0(format(signif(100 * x, 2), trim = TRUE), if (sign) " %")
}

# A group described by its columns and the unit of its abundances.
group_label <- function(g, per) {
  paste(c(sprintf("%s = %s", names(g), vapply(g, function(v) as_text(v)[1], character(1))),
          sprintf("per = %s", per)), collapse = ", ")
}

check_by <- function(by, call = sys.call(-1)) {
  if (!is.character(by) || length(by) == 0L || anyNA(by)) {
    abort_input("`by` must be a character vector of column names.", call = call)
  }
  invisible(by)
}

check_n_draws <- function(n_draws, call = sys.call(-1)) {
  if (!is.numeric(n_draws) || length(n_draws) != 1L || is.na(n_draws) || n_draws < 100 ||
      n_draws != round(n_draws)) {
    abort_input(sprintf("`n_draws` must be a whole number of at least 100, not %s.",
                        describe_value(n_draws)), call = call)
  }
  invisible(n_draws)
}

check_breaks <- function(breaks, call = sys.call(-1)) {
  if (is.null(breaks)) {
    return(NULL)
  }
  if (!is.numeric(breaks) || length(breaks) == 0L || anyNA(breaks) || any(breaks < 0) ||
      is.unsorted(breaks, strictly = TRUE)) {
    abort_input(sprintf(paste("`breaks` must be NULL or increasing class limits of 0 or more,",
                              "not %s."), describe_value(breaks)), call = call)
  }
  as.numeric(breaks)
}
