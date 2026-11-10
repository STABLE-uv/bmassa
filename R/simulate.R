# Totals of estimates by simulation with shared draws --------------------------

population_effect_levels <- c("shared", "independent")

#' @rdname summarize_biomass
#' @export
summarize_biomass.bmassa_estimate <- function(data, by = "site", level = 0.95, n_draws = 1000,
                                                population_effect = "shared", population_by = NULL,
                                                na_rm = FALSE, seed = NULL, ...) {
  check_dots_empty(...)
  check_by(by)
  check_level(level)
  check_n_draws(n_draws)
  check_choice(population_effect, population_effect_levels)
  check_flag(na_rm)
  check_seed(seed)
  e <- data$estimates
  pop_group <- population_groups(e, population_by, population_effect)
  by_label <- if (is.null(population_by)) NA_character_ else paste(population_by, collapse = ", ")
  check_columns(e, by, arg = "data$estimates")
  result_cols <- c("currency", "biomass_unit", "biomass_mean", "biomass_median", "biomass_lower",
                   "biomass_upper", "level", "interval", "n_rows", "n_estimated", "n_missing",
                   "status", "n_draws", "population_effect", "population_by")
  clash <- intersect(by, result_cols)
  if (length(clash) > 0L) {
    abort_input(sprintf("`by` cannot contain %s, which are columns of the result.",
                        format_items(quote_items(clash), max = Inf)))
  }
  if (nrow(e) == 0L) {
    return(data.frame(e[0, by, drop = FALSE], currency = character(), biomass_unit = character(),
                      biomass_mean = numeric(), biomass_median = numeric(),
                      biomass_lower = numeric(), biomass_upper = numeric(), level = numeric(),
                      interval = character(), n_rows = integer(), n_estimated = integer(),
                      n_missing = integer(), status = character(), n_draws = integer(),
                      population_effect = character(), population_by = character(),
                      stringsAsFactors = FALSE))
  }
  # Groups (C locale order) and, within them, currency and unit.
  key <- make_key(e, by)
  ukey <- unique(key)
  first <- match(ukey, key)
  ord <- do.call(order, c(unname(lapply(e[first, by, drop = FALSE], as_text)), list(method = "radix")))
  ukey <- ukey[ord]
  first <- first[ord]
  group <- match(key, ukey)
  known <- !is.na(e$biomass_mean)
  n_rows <- tabulate(group, nbins = length(ukey))
  n_missing <- tabulate(group[!known], nbins = length(ukey))
  est <- which(known)
  combo <- paste(group[est], e$currency[est], e$biomass_unit[est], sep = "\r")
  ucombo <- unique(combo)
  cfirst <- est[match(ucombo, combo)]
  cid <- match(combo, ucombo)
  # All the rows with a mass are drawn, as in simulate(), so that the same seed
  # gives the same draws; only those with a biomass enter the totals.
  drawn <- which(!is.na(e$ind_mean))
  draws <- if (length(est) == 0L) {
    matrix(numeric(0), nrow = n_draws, ncol = 0L)
  } else if (is.null(seed)) {
    draw_totals(data, est, cid, length(ucombo), n_draws, population_effect, pop_group, drawn)
  } else {
    withr::with_seed(seed, draw_totals(data, est, cid, length(ucombo), n_draws, population_effect,
                                       pop_group, drawn))
  }
  p <- c((1 - level) / 2, 1 - (1 - level) / 2)
  # The expectation of a sum is the sum of the expectations whatever the
  # correlations: the mean of a total is exact, only its quantiles are simulated.
  mean_total <- vapply(split(e$biomass_mean[est], factor(cid, levels = seq_along(ucombo))),
                       sum, numeric(1), USE.NAMES = FALSE)
  res <- data.frame(g = group[cfirst], currency = e$currency[cfirst],
                    biomass_unit = e$biomass_unit[cfirst],
                    biomass_mean = mean_total,
                    biomass_median = apply(draws, 2, stats::median),
                    biomass_lower = apply(draws, 2, stats::quantile, probs = p[1], names = FALSE),
                    biomass_upper = apply(draws, 2, stats::quantile, probs = p[2], names = FALSE),
                    n_estimated = tabulate(cid, nbins = length(ucombo)),
                    stringsAsFactors = FALSE)
  none <- setdiff(seq_along(ukey), res$g)
  if (length(none)) {
    res <- rbind(res, data.frame(g = none, currency = NA_character_, biomass_unit = NA_character_,
                                 biomass_mean = NA_real_, biomass_median = NA_real_,
                                 biomass_lower = NA_real_, biomass_upper = NA_real_,
                                 n_estimated = 0L, stringsAsFactors = FALSE))
  }
  res <- res[order(res$g, res$currency, res$biomass_unit, method = "radix"), ]
  incomplete <- n_missing[res$g] > 0L
  if (!na_rm) {
    res[incomplete, c("biomass_mean", "biomass_median", "biomass_lower", "biomass_upper")] <- NA_real_
  }
  out <- data.frame(e[first[res$g], by, drop = FALSE], currency = res$currency,
                    biomass_unit = res$biomass_unit, biomass_mean = res$biomass_mean,
                    biomass_median = res$biomass_median, biomass_lower = res$biomass_lower,
                    biomass_upper = res$biomass_upper, level = level,
                    interval = ifelse(is.na(res$currency), NA_character_, "prediction"),
                    n_rows = n_rows[res$g], n_estimated = res$n_estimated,
                    n_missing = n_missing[res$g],
                    status = ifelse(incomplete, "incomplete", "complete"),
                    n_draws = as.integer(n_draws), population_effect = population_effect,
                    population_by = by_label, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

# Draws of the totals: a matrix with one row per draw and one column per
# total. `est` are the rows with biomass and `cid` their total; `drawn`, the
# rows whose masses are drawn (`est` and possibly more, which change the draws
# of `est` but not their distribution); `group`, the group of the populations
# of every row (see population_groups()).
draw_totals <- function(x, est, cid, n_totals, n_draws, population_effect, group = NULL,
                        drawn = est) {
  pos <- match(est, drawn)
  if (anyNA(pos)) {
    abort_internal("rows with a biomass but without a drawn mass")
  }
  ld <- draw_log_mass(x, drawn, n_draws, population_effect, group)
  ld <- ld[, pos, drop = FALSE]
  a <- x$estimates$abundance[est]
  totals <- matrix(0, nrow = n_draws, ncol = n_totals)
  for (k in seq_len(n_totals)) {
    cols <- which(cid == k)
    totals[, k] <- exp(ld[, cols, drop = FALSE]) %*% a[cols]
  }
  totals
}

# Draws of the log mean individual mass of the rows `est` (matrix: draws by
# rows), with the errors shared as described in summarize_biomass(). `group`
# is NULL or the group of the populations of every row of the estimates: with
# "shared", the rows of a taxon and stage share their population only within
# a group.
draw_log_mass <- function(x, est, n_draws, population_effect, group = NULL) {
  sim <- x$simulation$rows[est, , drop = FALSE]
  keys <- x$simulation$missing_keys[est]
  ref <- x$reference
  nodes <- ref$nodes
  comps <- ref$components
  comp_value <- function(m, l, p) {
    v <- comps$value[comps$model == m & comps$level == l & comps$parameter == p]
    if (length(v)) v else 0
  }
  lwr <- sim$method %in% c("length", "typical_length")
  out <- matrix(sim$eta, nrow = n_draws, ncol = nrow(sim), byrow = TRUE)
  # 1. Errors of the nodes, from the root down.
  err <- draw_node_errors(unique(sim$node), nodes, n_draws)
  for (u in unique(sim$node)) {
    cols <- which(sim$node == u)
    eu <- err[[as.character(u)]]
    out[, cols] <- out[, cols] + eu[, 1]
    if (ncol(eu) == 2L) out[, cols] <- out[, cols] + eu[, 2] %o% sim$xc[cols]
  }
  # 1b. Typical lengths of tadpoles without lengths (version 0.6.0): the
  # errors of the node of the model of typical lengths and of its missing
  # levels, shared like those of the masses, enter the log mass times the
  # slope of the length-mass relation (their part between populations is in
  # `v_population`, step 4).
  tl <- if (is.null(sim$length_node)) integer(0) else which(!is.na(sim$length_node))
  if (length(tl)) {
    err_l <- draw_node_errors(unique(sim$length_node[tl]), nodes, n_draws)
    for (u in unique(sim$length_node[tl])) {
      cols <- tl[sim$length_node[tl] == u]
      out[, cols] <- out[, cols] + err_l[[as.character(u)]][, 1] %o% sim$length_slope[cols]
    }
    lkeys <- x$simulation$length_missing_keys[est]
    if (is.null(lkeys)) lkeys <- vector("list", length(est))
    rows_of_l <- rows_by_key(lkeys)
    for (k in names(rows_of_l)) {
      cols <- rows_of_l[[k]]
      parts <- strsplit(k, "\u001f", fixed = TRUE)[[1]]
      e <- stats::rnorm(n_draws, 0, sqrt(comp_value(parts[1], parts[3], "intercept")))
      out[, cols] <- out[, cols] + e %o% sim$length_slope[cols]
    }
  }
  # 2. Stage effects.
  for (s in unique(stats::na.omit(sim$stage_node))) {
    cols <- which(sim$stage_node == s)
    out[, cols] <- out[, cols] + stats::rnorm(n_draws, 0, sqrt(nodes$effect_var[s]))
  }
  # 2b. Fixed contrasts of stage (copepodites and nauplii of copepods): the
  # error of a contrast is shared by its rows (one standard normal per
  # contrast) and covaries with the error of the path of each node, so that
  # each row has the variance of its prediction (path + contrast).
  if (!is.null(sim$contrast)) {
    for (k in unique(stats::na.omit(sim$contrast))) {
      z <- stats::rnorm(n_draws)
      for (u in unique(sim$node[sim$contrast %in% k])) {
        cols <- which(sim$contrast %in% k & sim$node == u)
        m <- sim$model[cols[1]]
        vk <- ref$fixed$se[ref$fixed$model == m & ref$fixed$term == paste0("c_", k)]^2
        ck <- nodes[[paste0("contrast_cov_", k)]][u]
        vu <- nodes$path_var[u]
        b <- if (vu > 0) ck / vu else 0
        e <- b * err[[as.character(u)]][, 1] + sqrt(max(vk - b^2 * vu, 0)) * z
        out[, cols] <- out[, cols] + e
      }
    }
  }
  # 3. Missing levels: one deviation per missing node.
  rows_of <- rows_by_key(keys)
  for (k in names(rows_of)) {
    cols <- rows_of[[k]]
    parts <- strsplit(k, "\u001f", fixed = TRUE)[[1]]
    m <- parts[1]
    l <- parts[3]
    if (lwr[cols[1]]) {
      S <- matrix(c(comp_value(m, l, "intercept"), comp_value(m, l, "intercept_slope"),
                    comp_value(m, l, "intercept_slope"), comp_value(m, l, "slope")), 2L)
      d <- rmvnorm0(n_draws, S)
      out[, cols] <- out[, cols] + d[, 1] + d[, 2] %o% sim$xc[cols]
    } else {
      out[, cols] <- out[, cols] + stats::rnorm(n_draws, 0, sqrt(comp_value(m, l, "intercept")))
    }
  }
  # 4. Populations.
  pop <- if (population_effect == "shared") {
    paste(sim$model, sim$taxon_id, sim$stage_used, if (is.null(group)) "" else group[est],
          sep = "\u001f")
  } else {
    as.character(seq_len(nrow(sim)))
  }
  for (p in unique(pop)) {
    cols <- which(pop == p)
    out[, cols] <- out[, cols] + stats::rnorm(n_draws, 0, sqrt(sim$v_population[cols[1]]))
  }
  # 5. A new source or study, one per model.
  for (m in unique(sim$model)) {
    cols <- which(sim$model == m)
    out[, cols] <- out[, cols] + stats::rnorm(n_draws, 0, sqrt(sim$v_source[cols[1]]))
  }
  # 6. Currency conversion (convert_currency()): the linear map of the log mass
  # of each row, then the errors of the conversion factors, shared like those
  # of the masses.
  cv <- x$simulation$conversion
  if (!is.null(cv)) {
    out <- out * rep(cv$slope[est], each = n_draws) + rep(cv$shift[est], each = n_draws)
    out <- add_conversion_errors(out, x, est, cv, n_draws, population_effect, group)
  }
  out
}

# Errors of the conversion factors of the rows `est`: nodes (along the tree of
# each conversion model), missing levels (one per missing node), population
# (per taxon or per row), and source and habitat (one per conversion model).
# The steps of an inverse conversion enter with the opposite sign. Conventions
# with an uncertainty (amphibians) get the population error and, as the
# source, the error of the mean of their studies; the other conventions and
# the carbon equations add none.
add_conversion_errors <- function(out, x, est, cv, n_draws, population_effect, group = NULL) {
  st <- cv$steps
  keep <- which(st$row %in% est &
                  (st$type == "factor" |
                     (st$type == "convention" & (st$v_population > 0 | st$v_source > 0))))
  if (length(keep) == 0L) {
    return(out)
  }
  keys <- cv$missing_keys[keep]
  st <- st[keep, , drop = FALSE]
  cn <- x$reference$conversion_nodes
  cc <- x$reference$conversion_components
  # One column of errors per step, added to the log mass of its row at the end.
  e <- matrix(0, nrow = n_draws, ncol = nrow(st))
  err <- draw_node_errors(unique(st$node[!is.na(st$node)]), cn, n_draws)
  for (k in which(!is.na(st$node))) e[, k] <- err[[as.character(st$node[k])]][, 1]
  steps_of <- rows_by_key(keys)
  for (key in names(steps_of)) {
    parts <- strsplit(key, "\u001f", fixed = TRUE)[[1]]
    v <- cc$value[cc$model == parts[1] & cc$level == parts[2] & cc$parameter == "intercept"]
    ks <- steps_of[[key]]
    e[, ks] <- e[, ks] + stats::rnorm(n_draws, 0, sqrt(if (length(v)) v else 0))
  }
  sim <- x$simulation$rows
  pop <- if (population_effect == "shared") {
    paste(st$model, sim$taxon_id[st$row], sim$stage_used[st$row],
          if (is.null(group)) "" else group[st$row], sep = "\u001f")
  } else {
    paste(st$model, st$row, sep = "\u001f")
  }
  for (p in unique(pop)) {
    ks <- which(pop == p)
    e[, ks] <- e[, ks] + stats::rnorm(n_draws, 0, sqrt(st$v_population[ks[1]]))
  }
  for (m in unique(st$model)) {
    ks <- which(st$model == m)
    e[, ks] <- e[, ks] + stats::rnorm(n_draws, 0, sqrt(st$v_source[ks[1]])) +
      stats::rnorm(n_draws, 0, sqrt(st$v_habitat[ks[1]]))
  }
  col <- match(st$row, est)
  for (k in seq_len(nrow(st))) out[, col[k]] <- out[, col[k]] + st$sign[k] * e[, k]
  out
}

# For each key of missing levels (in the order of their first appearance), the
# elements of `keys` (a list of character vectors) that contain it.
rows_by_key <- function(keys) {
  all_keys <- unlist(keys)
  if (!length(all_keys)) {
    return(list())
  }
  idx <- rep(seq_along(keys), lengths(keys))
  lapply(split(idx, factor(all_keys, levels = unique(all_keys))), unique)
}

# Errors of the paths of the nodes `used` and of their ancestors, drawn from
# the root down: the error of a node given that of its parent has mean
# K V_p^-1 e_p and covariance V - K V_p^-1 K'. Returns a list of matrices
# (draws by 1 or 2 columns: intercept and slope) named by node row.
draw_node_errors <- function(used, nodes, n_draws) {
  chain <- integer(0)
  for (u in used) {
    cur <- u
    path <- integer(0)
    while (!is.na(cur)) {
      path <- c(cur, path)
      cur <- parent_row(cur, nodes)
    }
    chain <- c(chain, path)
  }
  chain <- unique(chain)
  depth <- vapply(chain, function(k) {
    if (nodes$level[k] == "root") 0L else length(strsplit(nodes$path[k], "|", fixed = TRUE)[[1]])
  }, integer(1))
  chain <- chain[order(depth)]
  err <- list()
  for (k in chain) {
    V <- node_cov(nodes, k)
    p <- parent_row(k, nodes)
    if (is.na(p)) {
      e <- rmvnorm0(n_draws, V)
    } else {
      Vp <- node_cov(nodes, p)
      K <- parent_cov_matrix(nodes, k, ncol(V))
      Vp_inv <- pinv_sym(Vp)
      R <- K %*% Vp_inv
      C <- V - R %*% t(K)
      e <- err[[as.character(p)]] %*% t(R) + rmvnorm0(n_draws, C)
    }
    err[[as.character(k)]] <- e
  }
  err
}

parent_row <- function(k, nodes) {
  if (nodes$level[k] == "root") return(NA_integer_)
  m <- nodes$model[k]
  parent_path <- nodes$parent[k]
  lt <- nodes$length_type[k]
  nxt <- if (is.na(parent_path) || identical(parent_path, "all")) {
    which(nodes$model == m & nodes$level == "root" &
            (if (is.na(lt)) TRUE else nodes$node == lt))
  } else {
    which(nodes$model == m & nodes$path == parent_path & nodes$level %in% main_ranks &
            (if (is.na(lt)) is.na(nodes$length_type) else nodes$length_type %in% lt))
  }
  if (length(nxt) != 1L) {
    abort_internal(sprintf("parent of node %s of model %s not found", nodes$node[k], m))
  }
  nxt
}

node_cov <- function(nodes, k) {
  if (is.na(nodes$path_var_slope[k])) {
    matrix(nodes$path_var[k], 1L)
  } else {
    matrix(c(nodes$path_var[k], nodes$path_cov[k], nodes$path_cov[k], nodes$path_var_slope[k]), 2L)
  }
}

parent_cov_matrix <- function(nodes, k, dim) {
  if (dim == 1L) {
    matrix(nodes$parent_cov[k], 1L)
  } else {
    matrix(c(nodes$parent_cov[k], nodes$parent_cov_si[k], nodes$parent_cov_is[k],
             nodes$parent_cov_ss[k]), 2L)
  }
}

# Pseudo-inverse of a small symmetric matrix.
pinv_sym <- function(S) {
  ev <- eigen(S, symmetric = TRUE)
  tol <- max(dim(S)) * max(abs(ev$values)) * .Machine$double.eps * 100
  inv <- ifelse(ev$values > tol, 1 / ev$values, 0)
  ev$vectors %*% (inv * t(ev$vectors))
}

# Draws from a centred normal distribution with covariance S (1 x 1 or
# 2 x 2), through its eigen decomposition; negative eigenvalues from
# rounding are set to zero.
rmvnorm0 <- function(n, S) {
  if (length(S) == 1L) {
    return(matrix(stats::rnorm(n, 0, sqrt(max(S[1], 0))), ncol = 1L))
  }
  ev <- eigen(S, symmetric = TRUE)
  A <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow = length(ev$values))
  matrix(stats::rnorm(n * ncol(S)), nrow = n) %*% t(A)
}

# Draws for other packages -------------------------------------------------------

#' Joint draws of individual masses and biomasses
#'
#' `simulate()` draws, for every row of the result of [estimate_biomass()] (or
#' of [convert_currency()] on it), the mean individual mass of the population
#' and its biomass, with the errors shared between rows as in the totals of
#' [summarize_biomass()]. The draws are meant for models that use the biomass
#' of every taxon, such as the energy fluxes of a food web: compute the
#' quantity of interest in each draw and summarise it over the draws, so that
#' the correlations between rows (the same taxon in several sites, taxa
#' estimated from the same family) are kept.
#'
#' @param object Result of [estimate_biomass()] or [convert_currency()].
#' @param nsim Number of draws (default 1000).
#' @param seed `NULL` (default: the current random state is used) or a number:
#'   the draws are then made with that seed and the random state of the session
#'   is restored afterwards (`withr::with_seed()`).
#' @param population_effect How the variation between populations is shared:
#'   `"shared"` (default; one draw per taxon and life stage for all its rows
#'   or, with `population_by`, for its rows of each group) or `"independent"`
#'   (one draw per row).
#' @param population_by `NULL` (default) or the names of columns of the
#'   abundance table (they are kept in the result) whose values separate
#'   populations, such as a habitat or a region: with
#'   `population_effect = "shared"`, the rows of a taxon and stage with the
#'   same values share their population, and those with different values have
#'   independent ones. Not allowed with `"independent"`; the columns must not
#'   have missing values.
#' @param ... Not used; an error if anything is given.
#' @param x A `bmassa_simulation` object, from `simulate()`.
#' @param row.names,optional Not used (for compatibility with
#'   [as.data.frame()]).
#' @param draws `NULL` (default: all) or the numbers of the draws to return.
#'
#' @return `simulate()` returns an object of class `bmassa_simulation`, a list
#'   with:
#'   * `rows`: a data frame with one row per row of `object`, in the same
#'     order: the columns of the abundance table and `broad_group`,
#'     `scientific_name` (if given), `life_stage`, `sex`, `count_unit`,
#'     `method`, `match_level`, `currency`,
#'     `unit`, `biomass_unit`, `ind_mean`, `biomass_mean`, `converted_from`
#'     (after a conversion), `status` and `reason`.
#'   * `ind_mass`: a numeric matrix with one row per row and one column per
#'     draw: the mean individual mass of the population, in `unit`; `NA` for
#'     rows without an estimate.
#'   * `settings`: `nsim`, `seed`, `population_effect`, `population_by` and the
#'     conversion, if any.
#'
#'   `as.data.frame()` returns the long table: one row per row and draw (draw
#'   after draw, the rows in their order) with `draw`, the columns of `rows`
#'   except `ind_mean`, `biomass_mean` and `reason`, `ind_mass` and `biomass`
#'   (`abundance` times `ind_mass`, in `biomass_unit`). Zero rows give zero
#'   rows.
#'
#' @details
#' The draws of each row are lognormal with the median and the spread of its
#' prediction interval: the mean, the median and the quantiles of many draws
#' approach `ind_mean`, `ind_median`, `ind_lower` and `ind_upper`. The errors
#' are shared as described in [summarize_biomass()]: nodes of the reference,
#' missing taxonomic levels, life stages, populations, the source and, after a
#' conversion, those of the conversion models. [summarize_biomass()] uses
#' these same draws: with the same `seed`, `n_draws = nsim` and the same
#' `population_effect` and `population_by`, its median and limits are
#' quantiles of the sums of these biomasses.
#'
#' Rows without an estimate get `NA` draws; rows with a mass and a missing
#' abundance get simulated masses and `NA` biomasses.
#'
#' Biomasses in different units (`biomass_unit`) must not be added up. To add
#' groups sampled per litre and per square metre, as in the nodes of a food
#' web, express all the abundances on the same basis (for instance per
#' mesocosm, with `per = "mesocosm"`) before [estimate_biomass()].
#'
#' **Size.** The long table has `nsim` times as many rows as `object`: use
#' `draws` to take a subset, or work with `ind_mass` directly.
#'
#' @section Assumptions:
#' The abundances are taken as known: the error of the counts is not included.
#' Each row has the mean mass of its population; the spread of sizes within a
#' population is not simulated (with different sizes, the metabolic demand of
#' a population is not its abundance times that of an individual of mean
#' mass). A taxon without lengths has the same distribution of mass in all its
#' rows, so that its differences of biomass between sites or dates come from
#' the abundances.
#'
#' @section Reproducibility:
#' With the same `seed`, the draws are the same on the same platform; other
#' versions of R or other BLAS and LAPACK libraries may give slightly
#' different values.
#'
#' @seealso [summarize_biomass()] for totals with prediction intervals.
#' @importFrom stats simulate
#' @export
#' @examples
#' est <- estimate_biomass(example_survey$abundance, example_survey$taxa,
#'                         example_survey$lengths, bias_correction = "lognormal")
#' wet <- convert_currency(est, to = "wet_mass", unit = "mg")
#' sim <- simulate(wet, nsim = 200, seed = 1)
#' sim
#' d <- as.data.frame(sim, draws = 1:2)
#' head(d[, c("draw", "site", "taxon_id", "ind_mass", "biomass", "biomass_unit")])
#' # Totals of each draw by site and group.
#' tot <- aggregate(biomass ~ draw + site + broad_group + biomass_unit, data = d, FUN = sum)
#' head(tot)
simulate.bmassa_estimate <- function(object, nsim = 1000, seed = NULL,
                                     population_effect = "shared", population_by = NULL, ...) {
  check_dots_empty(...)
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1 ||
      nsim != round(nsim)) {
    abort_input(sprintf("`nsim` must be a whole number of at least 1, not %s.",
                        describe_value(nsim)))
  }
  check_seed(seed)
  check_choice(population_effect, population_effect_levels)
  e <- object$estimates
  group <- population_groups(e, population_by, population_effect)
  nsim <- as.integer(nsim)
  drawn <- which(!is.na(e$ind_mean))
  ind_mass <- matrix(NA_real_, nrow = nrow(e), ncol = nsim)
  if (length(drawn)) {
    ld <- if (is.null(seed)) {
      draw_log_mass(object, drawn, nsim, population_effect, group)
    } else {
      withr::with_seed(seed, draw_log_mass(object, drawn, nsim, population_effect, group))
    }
    ind_mass[drawn, ] <- t(exp(ld))
  }
  structure(list(rows = simulation_rows(e), ind_mass = ind_mass,
                 settings = list(nsim = nsim, seed = seed, population_effect = population_effect,
                                 population_by = population_by,
                                 conversion = object$settings$conversion)),
            class = "bmassa_simulation")
}

#' @rdname simulate.bmassa_estimate
#' @export
as.data.frame.bmassa_simulation <- function(x, row.names = NULL, optional = FALSE, ...,
                                            draws = NULL) {
  check_dots_empty(...)
  nsim <- ncol(x$ind_mass)
  if (is.null(draws)) {
    draws <- seq_len(nsim)
  } else if (!is.numeric(draws) || length(draws) == 0L || anyNA(draws) ||
             any(draws != round(draws)) || any(draws < 1 | draws > nsim)) {
    abort_input(sprintf("`draws` must be NULL or whole numbers between 1 and %d, not %s.",
                        nsim, describe_value(draws)))
  }
  draws <- as.integer(draws)
  r <- x$rows
  keep <- setdiff(names(r), c("ind_mean", "biomass_mean", "reason"))
  clash <- intersect(c("draw", "ind_mass", "biomass"), keep)
  if (length(clash) > 0L) {
    abort_input(sprintf(
      "The abundance table has columns named %s, which the long table adds; rename them.",
      format_items(quote_items(clash), max = Inf)))
  }
  n <- nrow(r)
  out <- data.frame(draw = rep(draws, each = n), r[rep(seq_len(n), times = length(draws)), keep,
                                                     drop = FALSE],
                    stringsAsFactors = FALSE)
  out$ind_mass <- as.vector(x$ind_mass[, draws, drop = FALSE])
  out$biomass <- out$abundance * out$ind_mass
  rownames(out) <- NULL
  out
}

#' @export
print.bmassa_simulation <- function(x, ...) {
  r <- x$rows
  s <- x$settings
  has_mass <- !is.na(x$ind_mass[, 1L])
  cat(sprintf("<bmassa_simulation> %d draws of the mean individual mass of %d rows (%d with a mass, %d with a biomass)\n",
              s$nsim, nrow(r), sum(has_mass), sum(has_mass & !is.na(r$abundance))))
  if (any(has_mass)) {
    u <- unique(r[has_mass, c("currency", "unit")])
    cat("Currencies: ", paste(sprintf("%s (%s)", u$currency, u$unit), collapse = ", "), "\n",
        sep = "")
  }
  pop <- if (s$population_effect == "independent") {
    "one per row"
  } else if (is.null(s$population_by)) {
    "one per taxon and stage, shared by all its rows"
  } else {
    sprintf("one per taxon and stage within each value of %s",
            paste(sprintf("'%s'", s$population_by), collapse = " and "))
  }
  cat("Populations: ", pop, "\n", sep = "")
  cat("as.data.frame() gives one row per row and draw; x$ind_mass holds the masses\n")
  invisible(x)
}

# Columns of the rows of a simulation: those of the abundance table, then those
# that describe the estimate.
simulation_columns <- c("broad_group", "scientific_name", "life_stage", "sex", "count_unit",
                        "method", "match_level", "currency", "unit", "biomass_unit", "ind_mean",
                        "biomass_mean", "converted_from", "status", "reason")

simulation_rows <- function(e) {
  input <- input_columns(e)
  out <- e[, c(input, setdiff(intersect(simulation_columns, names(e)), input)), drop = FALSE]
  rownames(out) <- NULL
  out
}

# Columns of the estimates that came with the abundance table.
input_columns <- function(e) {
  setdiff(names(e), c(estimate_biomass_columns, conversion_columns, "scientific_name"))
}

# Group of the populations of every row (NULL without `population_by`).
population_groups <- function(e, population_by, population_effect, call = sys.call(-1)) {
  if (is.null(population_by)) {
    return(NULL)
  }
  if (!is.character(population_by) || length(population_by) == 0L || anyNA(population_by) ||
      anyDuplicated(population_by) > 0L) {
    abort_input("`population_by` must be NULL or a character vector of distinct column names.",
                call = call)
  }
  if (population_effect != "shared") {
    abort_input(paste("`population_by` sets which rows share a population, so it needs",
                      "`population_effect = \"shared\"`; with \"independent\" every row is",
                      "its own population."), call = call)
  }
  allowed <- setdiff(input_columns(e), "abundance")
  bad <- setdiff(population_by, allowed)
  if (length(bad) > 0L) {
    abort_input(sprintf("`population_by` must name columns of the abundance table (%s), not %s.",
                        format_items(quote_items(allowed), max = 10L),
                        format_items(quote_items(bad))), call = call)
  }
  for (col in population_by) {
    rows <- which(is.na(e[[col]]))
    if (length(rows) > 0L) {
      abort_input(sprintf(paste("Column '%s' of the abundance table (`population_by`) has missing",
                                "values in %s: give every row its population."),
                          col, format_rows(rows)), call = call)
    }
  }
  make_key(e, population_by)
}

check_seed <- function(seed, call = sys.call(-1)) {
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed))) {
    abort_input(sprintf("`seed` must be NULL or a single number, not %s.", describe_value(seed)),
                call = call)
  }
  invisible(seed)
}
