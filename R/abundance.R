# Abundances on a common basis ---------------------------------------------------

volume_unit_L <- c(mL = 1e-3, L = 1, m3 = 1e3)
area_unit_m2 <- c(cm2 = 1e-4, m2 = 1, ha = 1e4)
depth_unit_m <- c(cm = 0.01, m = 1)
abundance_added_columns <- c("abundance_original", "per_original", "per_factor", "per_assumption")

#' Abundances per a common volume, area or site
#'
#' `convert_abundance()` expresses the abundances of an abundance table on one
#' basis: per a volume (`"L"`, `"m3"`...), per an area (`"m2"`, `"100 m2"`...)
#' or per whole site (`"site"`), from the size of each sample in `per` and, when
#' needed, the volume, area or depth of each site. Use it before
#' [estimate_biomass()] when groups sampled per litre and per square metre are
#' to be added up, as in the nodes of a food web: the biomasses then share one
#' unit.
#'
#' @param abundance The abundance table of [estimate_biomass()]: one row per
#'   site and taxon, with `site`, `abundance` and `per`. `per` says what the
#'   abundance refers to: the size of the sample, a volume or an area, as a
#'   number and a unit or a unit alone (one of that unit), with the units
#'   `"mL"`, `"L"`, `"m3"` (volumes) and `"cm2"`, `"m2"`, `"ha"` (areas); or
#'   `"site"` for counts of the whole site (a mesocosm, say). Write what was
#'   sampled: the individuals in a trap of 30 x 10 cm are per `"300 cm2"`,
#'   those of two such traps counted together per `"600 cm2"`, those in 10 L
#'   of filtered water per `"10 L"`, a density per litre per `"L"`. Other
#'   columns are kept.
#' @param to The basis of the result, in the same form as `per`: a volume, an
#'   area or `"site"`.
#' @param sites `NULL` (default) or a data frame with one row per site:
#'   `site` and, as needed, `volume` and `volume_unit` (`"mL"`, `"L"`,
#'   `"m3"`), `area` and `area_unit` (`"cm2"`, `"m2"`, `"ha"`), and `depth` and
#'   `depth_unit` (`"cm"`, `"m"`). Needed to go between volumes and areas (the
#'   depth, or the volume and the area, whose ratio is the mean depth) and to
#'   or from `"site"` (its volume or area).
#'
#' @return The abundance table with the same rows in the same order,
#'   `abundance` and `per` on the new basis, and four more columns:
#'   `abundance_original`, `per_original`, `per_factor` (the abundance was
#'   multiplied by it) and `per_assumption` (what the conversion assumes, such
#'   as the share of the site that the sample is, or `NA` within one
#'   dimension). A zero-row table gives zero rows.
#'
#' @details
#' With the size of the sample \eqn{s} (in L or m\out{<sup>2</sup>}) the
#' abundance \eqn{a} is a density \eqn{a/s}. Within one dimension the factor
#' is the ratio of the sizes (`per = "mL"`, `to = "L"`: 1000). Between volume
#' and area it uses the depth \eqn{z} of the site (in m; the volume divided by
#' the area if no depth is given): a density per L times \eqn{1000 z} is a
#' density per m\out{<sup>2</sup>}, and the reverse. To `"site"`, a density
#' per L is multiplied by the volume of the site and a density per
#' m\out{<sup>2</sup>} by its area, that is, a count is divided by the share
#' of the site that the sample is (a trap of 300 cm\out{<sup>2</sup>} in a
#' mesocosm of 1.54 m\out{<sup>2</sup>}: 1.95 \%, so the count is multiplied
#' by 51.3); from `"site"`, the count is divided by them.
#'
#' @section Assumptions:
#' Going from volume to area (or to a site total) takes the density of the
#' sample as that of the whole water column; going from area to a site total
#' takes it as that of the whole bottom. Counts from semi-quantitative samples
#' (sweep nets, artificial substrates) have no area of their own: giving one in
#' `per` is an assumption of the user, and the result depends on it.
#'
#' @section Missing values and zeros:
#' A missing abundance stays missing, and a zero stays zero.
#'
#' @seealso [estimate_biomass()], which gives the biomass per `per`.
#' @export
#' @examples
#' ab <- example_survey$abundance
#' unique(ab[, c("taxon_id", "per")])[1:6, ]
#' # Everything per square metre of each site, with its mean depth.
#' sites <- data.frame(site = c("L1", "L2"), depth = c(1.5, 0.8), depth_unit = "m")
#' ab_m2 <- convert_abundance(ab, to = "m2", sites = sites)
#' head(ab_m2[, c("site", "taxon_id", "abundance_original", "per_original", "abundance",
#'                "per", "per_factor")])
#' # Or the total of each site, from its volume and area (e.g. mesocosms).
#' tanks <- data.frame(site = c("L1", "L2"), volume = 1.23, volume_unit = "m3",
#'                     area = 1.54, area_unit = "m2")
#' ab_site <- convert_abundance(ab, to = "site", sites = tanks)
#' unique(ab_site$per_assumption)
convert_abundance <- function(abundance, to, sites = NULL) {
  check_data_frame(abundance)
  check_columns(abundance, c("site", "abundance", "per"), arg = "abundance")
  check_numeric_column(abundance, "abundance", lower = 0, allow_na = TRUE, arg = "abundance")
  check_text_column(abundance, "site", "abundance")
  check_text_column(abundance, "per", "abundance")
  clash <- intersect(abundance_added_columns, names(abundance))
  if (length(clash) > 0L) {
    abort_input(sprintf(paste("`abundance` already has %s: it seems converted. Convert the",
                              "original table, once, to the basis you need."),
                        format_items(quote_items(clash), max = Inf)))
  }
  check_string(to)
  target <- parse_basis(to)
  if (is.na(target$dim)) {
    abort_input(sprintf(paste("`to` must be a volume, an area or \"site\" (%s), with or without",
                              "a number before the unit, not %s."),
                        basis_units_text(), describe_value(to)))
  }
  ab <- as_plain_data_frame(abundance)
  per <- as_text(ab$per)
  src <- parse_basis(per)
  bad <- which(is.na(src$dim))
  if (length(bad) > 0L) {
    abort_input(sprintf(paste("Column 'per' of `abundance` must give the size of the sample as a",
                              "unit (%s), with or without a number before it (\"5 L\", \"0.25",
                              "m2\"); check %s (%s)."),
                        basis_units_text(), format_rows(bad),
                        format_items(quote_items(unique(per[bad])))), rows = bad)
  }
  site <- as_text(ab$site)
  # The sizes of the sites are needed only between dimensions.
  size <- if (any(src$dim != target$dim)) site_sizes(sites) else NULL
  n <- nrow(ab)
  factor <- rep(NA_real_, n)
  assumption <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    conv <- basis_factor(src$dim[i], src$size[i], target$dim, target$size, size, site[i])
    factor[i] <- conv$factor
    assumption[i] <- conv$assumption
  }
  ab$abundance_original <- ab$abundance
  ab$per_original <- per
  ab$abundance <- ab$abundance * factor
  ab$per <- rep(target$label, n)
  ab$per_factor <- factor
  ab$per_assumption <- assumption
  ab
}

# Size and dimension of each basis: "volume" (size in L), "area" (in m2) or
# "site" (size 1); NA when it cannot be read. `label` is the normalised text.
parse_basis <- function(x) {
  x <- trimws(as_text(x))
  m <- regmatches(x, regexec("^([0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)?\\s*([A-Za-z][A-Za-z0-9]*)$",
                             x, perl = TRUE))
  ok <- lengths(m) == 3L
  has_num <- vapply(m, function(z) length(z) == 3L && nzchar(z[2]), logical(1))
  num <- rep(NA_real_, length(x))
  num[ok] <- 1
  num[has_num] <- as.numeric(vapply(m[has_num], `[`, character(1), 2L))
  unit <- rep(NA_character_, length(x))
  unit[ok] <- vapply(m[ok], `[`, character(1), 3L)
  known <- c(names(volume_unit_L), names(area_unit_m2), "site")
  canon <- known[match(tolower(unit), tolower(known))]
  dim <- ifelse(canon %in% names(volume_unit_L), "volume",
                ifelse(canon %in% names(area_unit_m2), "area",
                       ifelse(canon %in% "site", "site", NA_character_)))
  dim[is.na(num) | num <= 0 | (dim %in% "site" & has_num)] <- NA_character_
  size <- ifelse(dim %in% "volume", num * volume_unit_L[canon],
                 ifelse(dim %in% "area", num * area_unit_m2[canon], 1))
  num_text <- vapply(num, function(v) format(v, scientific = FALSE, trim = TRUE), character(1))
  label <- ifelse(has_num, paste(num_text, canon), canon)
  list(dim = unname(dim), size = unname(size), label = unname(label))
}

basis_units_text <- function() {
  sprintf("volumes %s; areas %s; or \"site\"",
          paste(sprintf("\"%s\"", names(volume_unit_L)), collapse = ", "),
          paste(sprintf("\"%s\"", names(area_unit_m2)), collapse = ", "))
}

# Volume (L), area (m2) and depth (m) of the sites.
site_sizes <- function(sites, call = sys.call(-1)) {
  force(call)
  if (is.null(sites)) {
    abort_input(paste("Going between volumes and areas, or to or from \"site\", needs `sites`",
                      "with the volume, area or depth of each site."), call = call)
  }
  check_data_frame(sites, call = call)
  check_columns(sites, "site", arg = "sites", call = call)
  s <- as_plain_data_frame(sites)
  check_text_column(s, "site", "sites", call = call)
  dup <- which(duplicated(as_text(s$site)))
  if (length(dup) > 0L) {
    abort_input(sprintf("Column 'site' of `sites` has repeated values in %s.", format_rows(dup)),
                call = call)
  }
  out <- data.frame(site = as_text(s$site), volume_L = NA_real_, area_m2 = NA_real_,
                    depth_m = NA_real_, stringsAsFactors = FALSE)
  for (q in list(c("volume", "volume_unit"), c("area", "area_unit"), c("depth", "depth_unit"))) {
    if (!q[1] %in% names(s)) next
    check_numeric_column(s, q[1], lower = 0, include_lower = FALSE, allow_na = TRUE,
                         arg = "sites", call = call)
    check_columns(s, q[2], arg = "sites", call = call)
    units <- switch(q[1], volume = volume_unit_L, area = area_unit_m2, depth = depth_unit_m)
    check_values_in(s, q[2], names(units), arg = "sites", call = call)
    u <- as_text(s[[q[2]]])
    bad <- which(!is.na(s[[q[1]]]) & is.na(u))
    if (length(bad) > 0L) {
      abort_input(sprintf("Column '%s' of `sites` is missing in %s, which have a %s.", q[2],
                          format_rows(bad), q[1]), call = call)
    }
    out[[switch(q[1], volume = "volume_L", area = "area_m2", depth = "depth_m")]] <-
      s[[q[1]]] * unname(units[u])
  }
  out
}

# Factor that turns an abundance per `from` into one per `to`, and its assumption.
basis_factor <- function(from_dim, from_size, to_dim, to_size, sizes, site,
                         call = sys.call(-1)) {
  force(call)
  if (from_dim == to_dim) {
    return(list(factor = if (to_dim == "site") 1 else to_size / from_size,
                assumption = NA_character_))
  }
  k <- match(site, sizes$site)
  if (is.na(k)) {
    abort_input(sprintf("Site '%s' of `abundance` is not in `sites`.", site), call = call)
  }
  get <- function(what, label) {
    v <- sizes[[what]][k]
    if (is.na(v)) {
      abort_input(sprintf("Site '%s' needs its %s in `sites`.", site, label), call = call)
    }
    v
  }
  depth <- function() {
    v <- sizes$depth_m[k]
    if (is.na(v)) {
      vol <- sizes$volume_L[k]
      area <- sizes$area_m2[k]
      if (is.na(vol) || is.na(area)) {
        abort_input(sprintf("Site '%s' needs its depth, or its volume and area, in `sites`.",
                            site), call = call)
      }
      v <- vol / 1000 / area
    }
    v
  }
  fmt <- function(x) format(signif(x, 3), scientific = FALSE, trim = TRUE)
  if (from_dim == "volume" && to_dim == "area") {
    z <- depth()
    return(list(factor = z * 1000 * to_size / from_size,
                assumption = sprintf("density uniform through the water column (depth %s m)",
                                     fmt(z))))
  }
  if (from_dim == "area" && to_dim == "volume") {
    z <- depth()
    return(list(factor = to_size / (from_size * z * 1000),
                assumption = sprintf("density uniform through the water column (depth %s m)",
                                     fmt(z))))
  }
  if (to_dim == "site") {
    if (from_dim == "volume") {
      v <- get("volume_L", "volume")
      return(list(factor = v / from_size,
                  assumption = sprintf(paste("sample of %s %% of the volume of the site (%s of %s",
                                             "L), its density taken for the whole site"),
                                       fmt(100 * from_size / v), fmt(from_size), fmt(v))))
    }
    a <- get("area_m2", "area")
    return(list(factor = a / from_size,
                assumption = sprintf(paste("sample of %s %% of the area of the site (%s of %s",
                                           "m2), its density taken for the whole site"),
                                     fmt(100 * from_size / a), fmt(from_size), fmt(a))))
  }
  # From "site" to a volume or an area.
  if (to_dim == "volume") {
    v <- get("volume_L", "volume")
    return(list(factor = to_size / v,
                assumption = sprintf("count of the site spread over its volume (%s L)", fmt(v))))
  }
  a <- get("area_m2", "area")
  list(factor = to_size / a,
       assumption = sprintf("count of the site spread over its area (%s m2)", fmt(a)))
}
