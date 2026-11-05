# Hierarchical models of version 0.2.0: fit with lme4 (REML) and export of the tables
#
# Reads data-raw/derived/model_typical.csv and model_lwr.csv (model_data.R) and
# writes, in data-raw/derived/:
#   model_nodes.csv       one row per model and node of the reference (and one
#                         root per model, or per length type): the estimate of
#                         the path from the root to the node (intercept, and
#                         slope in the length-mass models), its variance and its
#                         covariance with the path of the parent node; and the
#                         effect of the node itself
#   model_components.csv  variance components by model and level
#   model_fixed.csv       fixed effects, their covariance and the centre of
#                         the log lengths
#   model_fit_log.md      what was fitted, and diagnostics
# The functions are also used by dev/validacion/validacion_modelos.R (with
# BMASSA_FIT_ONLY_FUNCTIONS=true only the functions are defined).
#
# Model:
# - Typical masses, one model per domain and count unit: y = ln of the typical
#   mass (or biovolume) of a record; a fixed intercept, and in the model of
#   invertebrates fixed contrasts of the copepodites and nauplii of copepods
#   against their adults (column `contrast`); random intercepts for
#   each nested taxonomic level that splits the data; random intercepts of
#   stage group (invertebrates) or habitat (conversion factors), and of
#   source, when the model has two or more sources (a user's data are a new
#   source: what each source takes as the typical mass differs); the residual
#   is the variation between records of the same species and source
#   (populations, varieties, size classes).
# - Length-mass relations, one model per domain: three points per equation, at
#   the centres of the thirds of its log length range, and one point per
#   population mean (amphibians); x = ln length centred by length type; fixed
#   intercept and slope by length type; correlated random intercepts and slopes
#   for each nested taxonomic level within length type (an equation informs
#   only the nodes of its own length type; the variance components are shared),
#   with intercepts only at levels where most nodes have points at a single
#   length; random intercepts of equation (variation between populations) and,
#   with five or more studies, of study.
# - Prediction for a taxon: the path of its deepest node in the reference; the
#   variance of the path comes from the joint conditional covariance of fixed
#   and random effects (Henderson's mixed model equations, from the Cholesky
#   factors of lme4), so that the effects of a node and of its ancestors and the
#   fixed effects are not treated as independent; a contrast adds its
#   estimate, its variance and twice its covariance with the path.
# Natural logs throughout. Packages: lme4 (>= 1.1-35) and Matrix.
#
# Run from the package root, in a UTF-8 locale:  Rscript data-raw/fit_models.R

out_dir <- "data-raw/derived"
levels_all <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

# Path of each row down to a level: node names are unique within a model.
node_path <- function(d, level) {
  j <- match(level, levels_all)
  do.call(paste, c(unname(as.list(d[levels_all[seq_len(j)]])), sep = "|"))
}
# Lineage of each row from the hierarchy (ref_taxa.csv), with nodes of its own
# below its lowest rank ("Chironomus sp." of a source is a species of its own
# of Chironomus), so that every row has a full path from kingdom to species.
add_lineage <- function(d, rt) {
  k <- match(paste(d$source_id, d$taxon), paste(rt$source_id, rt$taxon))
  lost <- sum(is.na(k))
  d <- d[!is.na(k), , drop = FALSE]
  k <- k[!is.na(k)]
  d$col_id <- rt$col_id[k]
  d$domain <- rt$domain[k]
  d$lowest_rank <- rt$lowest_rank[k]
  for (r in levels_all) d[[r]] <- rt[[r]][k]
  low <- match(d$lowest_rank, levels_all)
  stopifnot(!anyNA(low))
  at_low <- as.matrix(d[levels_all])[cbind(seq_len(nrow(d)), low)]
  for (j in seq_along(levels_all)) {
    m <- low < j
    d[[levels_all[j]]][m] <- paste0(at_low[m], " [", levels_all[j], " not identified; ",
                                    d$source_id[m], ": ", d$taxon[m], "]")
  }
  stopifnot(!anyNA(as.matrix(d[levels_all])))
  attr(d, "lost") <- lost
  d
}
# Levels of a fitted model, from its variance components.
model_levels <- function(tabs) intersect(levels_all, tabs$components$level)
# Levels that split the data: more nodes than the level above, and more than one.
split_levels <- function(d) {
  n <- vapply(levels_all, function(l) length(unique(node_path(d, l))), integer(1))
  levels_all[n > 1 & n > c(1L, utils::head(n, -1))]
}
# Points of the length-mass relations.
lwr_points <- function(eq) {
  is_eq <- eq$kind == "equation"
  e <- eq[is_eq, ]
  lo <- log(e$length_min)
  w <- log(e$length_max) - lo
  p_eq <- e[rep(seq_len(nrow(e)), each = 3), ]
  p_eq$point <- rep(1:3, nrow(e))
  p_eq$x <- rep(lo, each = 3) + rep(w, each = 3) * (2 * p_eq$point - 1) / 6
  p_eq$y <- p_eq$ln_a + p_eq$b * p_eq$x
  m <- eq[!is_eq, ]
  m$point <- rep(0L, nrow(m))
  m$x <- log(m$length_mean)
  m$y <- log(m$mass_mean)
  p <- rbind(p_eq, m)
  rownames(p) <- NULL
  p
}
# Model frame: grouping factors (keys) and, for length-mass, the centred lengths.
model_frame <- function(d, type, levels, centre = NULL) {
  if (type == "lwr") {
    lt <- sort(unique(d$length_type))
    if (is.null(centre)) centre <- vapply(lt, function(t) mean(d$x[d$length_type == t]), numeric(1))
    d$xc <- d$x - centre[d$length_type]
    for (t in lt) {
      d[[paste0("i_", t)]] <- as.numeric(d$length_type == t)
      d[[paste0("s_", t)]] <- d[[paste0("i_", t)]] * d$xc
    }
    d$root <- d$length_type
    for (l in levels) d[[paste0("g_", l)]] <- factor(paste(d$length_type, node_path(d, l), sep = "|"))
  } else {
    d$root <- "all"
    for (l in levels) d[[paste0("g_", l)]] <- factor(node_path(d, l))
    # Fixed contrasts (stages of copepods): one indicator per contrast.
    if ("contrast" %in% names(d)) {
      for (k in sort(unique(stats::na.omit(d$contrast)))) {
        d[[paste0("c_", k)]] <- as.numeric(d$contrast %in% k)
      }
    }
  }
  attr(d, "centre") <- centre
  d
}
# Levels with random slopes: those where at least half of the nodes have
# points at two or more lengths.
slope_levels <- function(d, levels) {
  keep <- vapply(levels, function(l) {
    nx <- tapply(d$x, d[[paste0("g_", l)]], function(z) length(unique(round(z, 8))))
    mean(nx >= 2) >= 0.5
  }, logical(1))
  levels[keep]
}
# Crossed terms of a model.
model_options <- function(d, type) {
  if (type == "lwr") {
    return(list(study = length(unique(d$study)) >= 5, source = FALSE, context = character(0)))
  }
  # Context factors (crossed): stage groups, and habitat of the conversion
  # factors when two or more habitats have 10 or more rows each.
  context <- character(0)
  if ("stage" %in% names(d) && length(unique(d$stage)) > 1) context <- c(context, "stage")
  if ("habitat" %in% names(d) && sum(table(d$habitat) >= 10) >= 2) context <- c(context, "habitat")
  list(study = FALSE, source = length(unique(d$source_id)) >= 2, context = context)
}
model_formula <- function(d, type, levels, opt) {
  if (type == "lwr") {
    lt <- sort(unique(d$length_type))
    sl <- slope_levels(d, levels)
    fixed <- paste("0 +", paste(c(paste0("i_", lt), paste0("s_", lt)), collapse = " + "))
    rand <- c(sprintf(ifelse(levels %in% sl, "(1 + xc | g_%s)", "(1 | g_%s)"), levels),
              "(1 | equation_id)", if (opt$study) "(1 | study)")
  } else {
    fixed <- paste(c("1", grep("^c_", names(d), value = TRUE)), collapse = " + ")
    rand <- c(sprintf("(1 | g_%s)", levels), if (opt$source) "(1 | source_id)",
              sprintf("(1 | %s)", opt$context))
  }
  stats::as.formula(paste("y ~", fixed, "+", paste(rand, collapse = " + ")))
}
fit_model <- function(d, type, levels, opt) {
  lme4::lmer(model_formula(d, type, levels, opt), data = d, REML = TRUE,
             control = lme4::lmerControl(calc.derivs = FALSE, optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 1e5)))
}

# Tables of a fitted model ------------------------------------------------------------

model_tables <- function(fit, d, model, type, levels) {
  s2 <- stats::sigma(fit)^2
  beta <- lme4::fixef(fit)
  bhat <- lme4::getME(fit, "b")
  cn <- lme4::getME(fit, "cnms")
  gp <- lme4::getME(fit, "Gp")
  fl <- lme4::getME(fit, "flist")
  asg <- attr(fl, "assign")
  q <- length(bhat)
  p <- length(beta)
  # Position in b of a coefficient of a level of a grouping factor.
  bpos <- function(g, key, coef) {
    t <- match(g, names(cn))
    nc <- length(cn[[t]])
    gp[t] + (match(key, levels(fl[[asg[t]]])) - 1L) * nc + match(coef, cn[[t]])
  }
  has_slope <- vapply(levels, function(l) "xc" %in% cn[[paste0("g_", l)]], logical(1))
  # Tree of nodes: roots (one per length type, or one) and the nodes by level.
  keys <- vapply(levels, function(l) as.character(d[[paste0("g_", l)]]), character(nrow(d)))
  keys <- matrix(keys, nrow = nrow(d))
  tree <- list(data.frame(level = "root", key = sort(unique(d$root)), parent = NA_character_,
                          stringsAsFactors = FALSE))
  for (k in seq_along(levels)) {
    u <- unique(data.frame(level = levels[k], key = keys[, k],
                           parent = if (k == 1) d$root else keys[, k - 1], stringsAsFactors = FALSE))
    stopifnot(!anyDuplicated(u$key))
    tree[[k + 1]] <- u
  }
  tree <- do.call(rbind, tree)
  rownames(tree) <- NULL
  tree$root <- tree$key
  for (k in which(tree$level != "root")) tree$root[k] <- tree$root[match(tree$parent[k], tree$key)]
  n_nodes <- nrow(tree)
  slope <- type == "lwr"
  m <- if (slope) 2L * n_nodes else n_nodes
  # Coefficient vectors of each path: random part (in b) and fixed part.
  ci <- integer(0); cj <- integer(0)
  ai <- integer(0); aj <- integer(0)
  for (k in seq_len(n_nodes)) {
    a <- k
    while (tree$level[a] != "root") {
      g <- paste0("g_", tree$level[a])
      ci <- c(ci, bpos(g, tree$key[a], "(Intercept)"))
      cj <- c(cj, k)
      if (slope && has_slope[[tree$level[a]]]) {
        ci <- c(ci, bpos(g, tree$key[a], "xc"))
        cj <- c(cj, n_nodes + k)
      }
      a <- match(tree$parent[a], tree$key)
    }
    if (slope) {
      ai <- c(ai, match(paste0(c("i_", "s_"), tree$root[k]), names(beta)))
      aj <- c(aj, c(k, n_nodes + k))
    } else {
      ai <- c(ai, match("(Intercept)", names(beta)))
      aj <- c(aj, k)
    }
  }
  Cb <- Matrix::sparseMatrix(i = ci, j = cj, x = 1, dims = c(q, m))
  A <- Matrix::sparseMatrix(i = ai, j = aj, x = 1, dims = c(p, m))
  est <- as.vector(Matrix::crossprod(A, beta) + Matrix::crossprod(Cb, bhat))
  # Joint conditional covariance: Var(w'(theta_hat - theta)) = s2 * |M^-1 w|^2,
  # with M = [L 0; RZX' RX'] and w = [P Lambda' c; a] (lme4: L L' =
  # P (Lambda' Z'Z Lambda + I) P', L RZX = P Lambda' Z'X, RX'RX = X'X - RZX'RZX).
  L <- lme4::getME(fit, "L")
  Wu <- lme4::getME(fit, "Lambdat") %*% Cb
  Z1 <- as.matrix(Matrix::solve(L, Matrix::solve(L, Wu, system = "P"), system = "L"))
  RX <- as.matrix(lme4::getME(fit, "RX"))
  RZX <- as.matrix(lme4::getME(fit, "RZX"))
  Z2 <- forwardsolve(t(RX), as.matrix(A) - t(RZX) %*% Z1)
  dotc <- function(i, j) s2 * (colSums(Z1[, i, drop = FALSE] * Z1[, j, drop = FALSE]) +
                                 colSums(Z2[, i, drop = FALSE] * Z2[, j, drop = FALSE]))
  ii <- seq_len(n_nodes)
  par <- match(tree$parent, tree$key)
  pp <- ifelse(is.na(par), ii, par)
  nodes <- data.frame(model = model, level = tree$level, length_type = if (slope) tree$root else NA_character_,
                      node = ifelse(tree$level == "root", tree$key, sub("^.*\\|", "", tree$key)),
                      path = tree$key, parent = tree$parent, stringsAsFactors = FALSE)
  if (slope) nodes$path <- ifelse(tree$level == "root", NA_character_, sub("^[^|]*\\|", "", tree$key))
  if (slope) nodes$parent <- ifelse(tree$level %in% c("root", levels[1]), NA_character_,
                                    sub("^[^|]*\\|", "", tree$parent))
  unit <- if (slope) "equation_id" else "record_id"
  nodes$n_units <- vapply(seq_len(n_nodes), function(k) {
    if (tree$level[k] == "root") return(length(unique(d[[unit]][d$root == tree$key[k]])))
    length(unique(d[[unit]][as.character(d[[paste0("g_", tree$level[k])]]) == tree$key[k]]))
  }, integer(1))
  nodes$path_intercept <- est[ii]
  nodes$path_var <- dotc(ii, ii)
  nodes$parent_cov <- ifelse(is.na(par), NA_real_, dotc(ii, pp))
  # Fixed contrasts (typical models): covariance of each path with each
  # contrast, so that a prediction with a contrast has its exact variance
  # (path + contrast). A contrast has no random part (Z1 = 0).
  for (k in grep("^c_", names(beta), value = TRUE)) {
    z2k <- forwardsolve(t(RX), as.numeric(names(beta) == k))
    nodes[[paste0("contrast_cov_", sub("^c_", "", k))]] <- s2 * colSums(Z2[, ii, drop = FALSE] * z2k)
  }
  if (slope) {
    is_ <- n_nodes + ii
    ps <- n_nodes + pp
    nodes$path_slope <- est[is_]
    nodes$path_var_slope <- dotc(is_, is_)
    nodes$path_cov <- dotc(ii, is_)
    nodes$parent_cov_is <- ifelse(is.na(par), NA_real_, dotc(ii, ps))
    nodes$parent_cov_si <- ifelse(is.na(par), NA_real_, dotc(is_, pp))
    nodes$parent_cov_ss <- ifelse(is.na(par), NA_real_, dotc(is_, ps))
  }
  # Effect of each node itself (lme4 conditional modes and variances, given
  # the fixed effects).
  re <- lme4::ranef(fit, condVar = TRUE)
  own <- function(k, coef) {
    if (tree$level[k] == "root") {
      j <- if (coef == "(Intercept)") (if (slope) paste0("i_", tree$key[k]) else "(Intercept)") else
        paste0("s_", tree$key[k])
      return(unname(beta[j]))
    }
    g <- paste0("g_", tree$level[k])
    if (!coef %in% names(re[[g]])) return(0)
    re[[g]][tree$key[k], coef]
  }
  nodes$effect <- vapply(seq_len(n_nodes), own, numeric(1), coef = "(Intercept)")
  if (slope) nodes$effect_slope <- vapply(seq_len(n_nodes), own, numeric(1), coef = "xc")
  # Context factors (stage groups, habitat): their effects and conditional
  # variances, added to the path in a prediction.
  nodes$effect_var <- NA_real_
  for (cf in intersect(c("stage", "habitat"), names(re))) {
    r <- re[[cf]]
    extra <- data.frame(model = model, level = cf, length_type = NA_character_, node = rownames(r),
                        path = rownames(r), parent = NA_character_,
                        n_units = as.integer(table(d[[cf]])[rownames(r)]), path_intercept = NA_real_,
                        path_var = NA_real_, parent_cov = NA_real_, effect = r[["(Intercept)"]],
                        stringsAsFactors = FALSE)
    extra$effect_var <- attr(r, "postVar")[1, 1, ]
    for (j in setdiff(names(nodes), names(extra))) extra[[j]] <- NA_real_
    nodes <- rbind(nodes, extra[, names(nodes)])
  }
  vc <- as.data.frame(lme4::VarCorr(fit))
  comp <- data.frame(model = model, term = vc$grp, var1 = vc$var1, var2 = vc$var2, value = vc$vcov,
                     sd_cor = vc$sdcor, stringsAsFactors = FALSE)
  comp$level <- ifelse(startsWith(comp$term, "g_"), sub("^g_", "", comp$term),
                       ifelse(comp$term == "equation_id", "equation",
                              ifelse(comp$term == "source_id", "source", tolower(comp$term))))
  comp$quantity <- ifelse(!is.na(comp$var2), "covariance", "variance")
  comp$parameter <- ifelse(comp$quantity == "covariance", "intercept_slope",
                           ifelse(comp$var1 %in% "xc", "slope", "intercept"))
  comp <- comp[, c("model", "level", "parameter", "quantity", "value", "sd_cor")]
  V <- as.matrix(stats::vcov(fit))
  fixed <- data.frame(model = model, term = names(beta), estimate = unname(beta),
                      se = sqrt(diag(V)), stringsAsFactors = FALSE)
  fixed$vcov <- vapply(names(beta), function(j) paste(signif(V[j, ], 10), collapse = ";"), character(1))
  centre <- attr(d, "centre")
  fixed$centre <- if (is.null(centre)) NA_real_ else unname(centre[sub("^[is]_", "", names(beta))])
  list(nodes = nodes, components = comp, fixed = fixed)
}

# Prediction for new rows from the tables: deepest node of the reference on
# the path of the row; variance of the estimate of that path, of the levels
# below it that are missing, and of a new population (and of a new source or
# study, when the model has them).
predict_tables <- function(tabs, newdata, type, levels) {
  nodes <- tabs$nodes
  comp <- tabs$components
  getc <- function(level, parameter) {
    v <- comp$value[comp$level == level & comp$parameter == parameter]
    if (length(v)) v else 0
  }
  n <- nrow(newdata)
  slope <- type == "lwr"
  pre <- if (slope) newdata$length_type else rep("all", n)
  nkey <- ifelse(nodes$level == "root", nodes$node,
                 if (slope) paste(nodes$length_type, nodes$path, sep = "|") else nodes$path)
  nkey <- paste(nodes$level, nkey)
  k <- match(paste("root", pre), nkey)
  matched <- rep("root", n)
  for (l in levels) {
    key <- node_path(newdata, l)
    if (slope) key <- paste(pre, key, sep = "|")
    kk <- match(paste(l, key), nkey)
    hit <- !is.na(kk) & matched == c("root", levels)[match(l, levels)]
    k[hit] <- kk[hit]
    matched[hit] <- l
  }
  xc <- if (slope) {
    fx <- tabs$fixed
    newdata$x - fx$centre[match(paste0("i_", newdata$length_type), fx$term)]
  } else rep(0, n)
  eta <- nodes$path_intercept[k]
  v_param <- nodes$path_var[k]
  if (slope) {
    eta <- eta + nodes$path_slope[k] * xc
    v_param <- v_param + 2 * xc * nodes$path_cov[k] + xc^2 * nodes$path_var_slope[k]
  }
  for (cf in intersect(c("stage", "habitat"), nodes$level)) {
    ks <- match(paste(cf, newdata[[cf]]), paste(nodes$level, nodes$node))
    eta <- eta + nodes$effect[ks]
    v_param <- v_param + nodes$effect_var[ks]
  }
  if ("contrast" %in% names(newdata)) {
    for (ck in unique(stats::na.omit(newdata$contrast))) {
      fx <- tabs$fixed[tabs$fixed$term == paste0("c_", ck), ]
      stopifnot(nrow(fx) == 1L)
      w <- which(newdata$contrast %in% ck)
      eta[w] <- eta[w] + fx$estimate
      v_param[w] <- v_param[w] + fx$se^2 + 2 * nodes[[paste0("contrast_cov_", ck)]][k[w]]
    }
  }
  v_missing <- rep(0, n)
  depth <- match(matched, c("root", levels))
  for (j in seq_along(levels)) {
    miss <- depth <= j
    v_missing[miss] <- v_missing[miss] + getc(levels[j], "intercept") +
      if (slope) 2 * xc[miss] * getc(levels[j], "intercept_slope") + xc[miss]^2 * getc(levels[j], "slope") else 0
  }
  v_pop <- getc("residual", "intercept") + getc("equation", "intercept") + getc("study", "intercept") +
    getc("source", "intercept")
  data.frame(eta = eta, v_param = v_param, v_missing = v_missing, v_population = v_pop,
             v_total = v_param + v_missing + v_pop, matched_level = matched,
             stringsAsFactors = FALSE)
}

model_specs <- function(typ, lwr) {
  rbind(data.frame(model = sort(unique(typ$model)), type = "typical", stringsAsFactors = FALSE),
        data.frame(model = sort(unique(lwr$model)), type = "lwr", stringsAsFactors = FALSE))
}
read_model_data <- function() {
  typ <- read.csv(file.path(out_dir, "model_typical.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
  lwr <- read.csv(file.path(out_dir, "model_lwr.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
  # One length-mass model per domain, with all its length types.
  lwr$model <- paste0(lwr$domain, "_lwr")
  list(typ = typ, lwr = lwr)
}
model_rows <- function(dat, model, type) {
  if (type == "typical") dat$typ[dat$typ$model == model, ] else lwr_points(dat$lwr[dat$lwr$model == model, ])
}

if (identical(Sys.getenv("BMASSA_FIT_ONLY_FUNCTIONS"), "true")) {
  message("fit_models.R: functions loaded")
} else {
  dat <- read_model_data()
  log_lines <- c("# Fits of the hierarchical models", "",
                 paste0("Built by `data-raw/fit_models.R` on ", format(Sys.Date()), " (", R.version.string,
                        ", lme4 ", as.character(utils::packageVersion("lme4")), ", REML)."), "")
  say <- function(...) log_lines <<- c(log_lines, sprintf(...))
  out <- list(nodes = list(), components = list(), fixed = list())
  specs <- model_specs(dat$typ, dat$lwr)
  for (i in seq_len(nrow(specs))) {
    m <- specs$model[i]
    type <- specs$type[i]
    d <- model_rows(dat, m, type)
    levels <- split_levels(d)
    opt <- model_options(d, type)
    d <- model_frame(d, type, levels)
    warn <- character(0)
    t0 <- Sys.time()
    fit <- withCallingHandlers(fit_model(d, type, levels, opt),
                               warning = function(w) {
                                 warn <<- c(warn, conditionMessage(w))
                                 invokeRestart("muffleWarning")
                               },
                               message = function(w) {
                                 warn <<- c(warn, trimws(conditionMessage(w)))
                                 invokeRestart("muffleMessage")
                               })
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    tabs <- model_tables(fit, d, m, type, levels)
    for (k in names(out)) out[[k]][[m]] <- tabs[[k]]
    say("## %s", m)
    say("")
    say("- %d rows, %d %s; levels: %s%s%s%s; fitted in %.1f s; singular: %s%s.", nrow(d),
        length(unique(if (type == "lwr") d$equation_id else d$record_id)),
        if (type == "lwr") "equations and means" else "records", paste(levels, collapse = ", "),
        if (opt$study) "; plus study" else "", if (opt$source) "; plus source" else "",
        paste0(if (length(opt$context)) "; plus ", paste(opt$context, collapse = ", ")), secs,
        lme4::isSingular(fit),
        if (length(warn)) paste0(" (", paste(unique(warn), collapse = "; "), ")") else "")
    if (type == "lwr") say("- Levels with random slopes: %s.", paste(slope_levels(d, levels), collapse = ", "))
    cp <- tabs$components
    v <- cp[cp$quantity == "variance", ]
    say("- Standard deviations (ln): %s.",
        paste(sprintf("%s%s %.3f", v$level, ifelse(v$parameter == "slope", " (slope)", ""), v$sd_cor),
              collapse = "; "))
    # Spread for a new population by the lowest level known: population terms
    # plus the levels below (x/÷ at 95 %, at the centre of the lengths, without
    # the uncertainty of the estimates).
    vv <- stats::setNames(v$value[v$parameter == "intercept"], v$level[v$parameter == "intercept"])
    pop <- sum(vv[intersect(names(vv), c("residual", "equation", "study", "source"))])
    below <- rev(cumsum(rev(vv[levels])))
    fac <- exp(1.96 * sqrt(pop + c(below[-1], 0)))
    say("- 95 %% factor (x/÷) for a new population by the lowest level known%s: %s.",
        if (type == "lwr") " (at the centre of the lengths)" else "",
        paste(sprintf("%s %.1f", rev(levels), rev(fac)), collapse = ", "))
    say("- Fixed effects: %s.", paste(sprintf("%s %.3f (SE %.3f)", tabs$fixed$term, tabs$fixed$estimate,
                                               tabs$fixed$se), collapse = "; "))
    st <- tabs$nodes[tabs$nodes$level %in% c("stage", "habitat"), ]
    if (nrow(st)) say("- Effects of %s: %s.", paste(unique(st$level), collapse = " and "),
                      paste(sprintf("%s %.2f", st$node, st$effect), collapse = "; "))
    say("")
  }
  node_cols <- unique(unlist(lapply(out$nodes, names)))
  nodes <- do.call(rbind, lapply(out$nodes, function(x) {
    for (j in setdiff(node_cols, names(x))) x[[j]] <- NA_real_
    x[, node_cols]
  }))
  components <- do.call(rbind, out$components)
  fixed <- do.call(rbind, out$fixed)
  rownames(nodes) <- rownames(components) <- rownames(fixed) <- NULL
  write.csv(nodes, file.path(out_dir, "model_nodes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(components, file.path(out_dir, "model_components.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(fixed, file.path(out_dir, "model_fixed.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  writeLines(log_lines, file.path(out_dir, "model_fit_log.md"), useBytes = TRUE)
  cat(log_lines, sep = "\n")
}
