# Reader of the fish of NEON (National Ecological Observatory Network; Fish
# electrofishing, gill netting, and fyke netting counts, DP1.20107.001), used
# by data-raw/fit_lwr.R and data-raw/typical_records.R. It reads the download
# made by data-raw/neon_fish_download.R and returns the fish measured one by
# one (table fsh_perFish) that can be used: their mass for the means of
# populations, and which of them enter the length-mass equations.
#
# Decisions (2026-11-06; dev/decisions.md):
# - Only the data of RELEASE-2026, which has a DOI; the provisional data of
#   2016 and 2025-2026 are left out.
# - Fish identified to species: taxonRank species or subspecies (taken as its
#   species), without identification qualifier (cf., aff.) or morphospecies;
#   weight above zero. Total lengths below 5 mm are errors (taken as missing).
# - The fish of a bout are its sampling: site, year and season (eventID, with
#   ".netCompromised" removed: a compromised block net changes the counts, not
#   the masses of the fish caught); gears apart (electrofishing, with one or
#   two electrofishers; mini-fyke net; gill net).
# - Floor of the weights: no fish is recorded below 0.3 g or at 0.4 g, and
#   0.3 g (17 % of the fish) and 0.5 g are recorded for fish over a wide range
#   of small lengths (95 % of the blacknose dace of 15-35 mm weigh 0.3 g).
#   Weights of `floor_max` (0.5 g) or less are taken as the minimum of the
#   scale: the fish weighed less than `limit` (0.55 g) (censored).
# - Equation of a species, used to find errors and to give a mass to the
#   censored fish: least squares of ln(W) on ln(L), all sites pooled, after
#   one pass removing |studentized residual| > 4, with at least 10 fish over a
#   1.5-fold range of lengths. Only fish longer than the length at which the
#   equation predicts `2 * limit` (1.1 g) enter, so that the floor does not
#   bend it (at that length fewer than about 1 in 1000 fish weigh less than
#   the limit): first fit to the uncensored fish of 1.1 g or more, then twice
#   to the uncensored fish above the length given by the previous fit.
#   Equations with b outside 1.5-4.5 are not used. A genus equation (its fish
#   pooled) serves species without one.
# - Gross errors, left out: an uncensored fish whose weight is more than 3
#   times above or below the prediction from its length and beyond 4 residual
#   SD; a censored fish whose length predicts more than 3 times the limit and
#   4 SD above it.
# - Censored fish take the mean mass of the fish of their length lighter than
#   the limit (lognormal of the equation truncated at the limit); without an
#   equation or a valid length their mass is missing (NA), and the means of
#   populations that hold such a fish are not computed (typical_records.R): the
#   guppies and swordtails of Puerto Rico, mostly censored, have no equation.
# Packages: none beyond R.

neon_settings <- list(release = "RELEASE-2026", min_length = 5, floor_max = 0.5, limit = 0.55,
                      min_n = 10, min_range_ratio = 1.5, max_rstudent = 4, b_range = c(1.5, 4.5),
                      gross_ratio = 3, gross_sd = 4)

# Equation of ln(W) on ln(L) for the floor-free fish of one taxon (species or
# genus): list(ln_a, b, sigma, n, length_min, length_max, floor_length), or
# NULL.
neon_equation <- function(L, W, censored, s = neon_settings) {
  ok <- is.finite(L) & is.finite(W) & !censored
  L <- L[ok]
  W <- W[ok]
  fit <- function(keep) {
    if (sum(keep) < s$min_n || max(L[keep]) / min(L[keep]) < s$min_range_ratio) return(NULL)
    x <- log(L[keep])
    y <- log(W[keep])
    f <- stats::lm(y ~ x)
    rs <- stats::rstudent(f)
    drop <- is.finite(rs) & abs(rs) > s$max_rstudent
    if (any(drop)) {
      x <- x[!drop]
      y <- y[!drop]
      if (length(x) < s$min_n || exp(max(x) - min(x)) < s$min_range_ratio) return(NULL)
      f <- stats::lm(y ~ x)
    }
    list(ln_a = unname(stats::coef(f)[1]), b = unname(stats::coef(f)[2]),
         sigma = summary(f)$sigma, n = length(x), length_min = exp(min(x)),
         length_max = exp(max(x)))
  }
  e <- fit(W >= 2 * s$limit)
  floor_length <- NA_real_
  for (k in 1:2) {
    if (is.null(e) || !is.finite(e$b) || e$b <= 0) return(NULL)
    floor_length <- exp((log(2 * s$limit) - e$ln_a) / e$b)
    e <- fit(L >= floor_length)
  }
  if (is.null(e) || e$b < s$b_range[1] || e$b > s$b_range[2]) return(NULL)
  e$floor_length <- floor_length
  e
}

read_neon_fish <- function(path, s = neon_settings) {
  x <- readRDS(path)
  p <- x$fsh_perFish
  num <- function(v) suppressWarnings(as.numeric(v))
  blank <- function(v) is.na(v) | !nzchar(trimws(v))
  counts <- c(all = nrow(p))
  k <- p$release %in% s$release
  counts["release"] <- sum(k)
  k <- k & p$taxonRank %in% c("species", "subspecies") & blank(p$identificationQualifier) &
    blank(p$morphospeciesID)
  counts["species"] <- sum(k)
  W <- num(p$fishWeight)
  k <- k & is.finite(W) & W > 0
  counts["weighed"] <- sum(k)
  p <- p[k, ]
  W <- W[k]
  L <- num(p$fishTotalLength)
  counts["short_length"] <- sum(is.finite(L) & L < s$min_length)
  L[!is.finite(L) | L < s$min_length] <- NA_real_
  taxon <- sub("^(\\S+ \\S+).*$", "\\1", trimws(p$scientificName))
  genus <- sub(" .*$", "", taxon)
  gear <- ifelse(p$samplerType %in% c("electrofisher", "two electrofishers"), "electrofishing",
                 p$samplerType)
  stopifnot(all(gear %in% c("electrofishing", "mini-fyke net", "gill net")))
  bout <- sub("\\.netCompromised$", "", p$eventID)
  stopifnot(all(startsWith(bout, paste0(p$siteID, "."))))
  bout <- gsub(".", " ", substring(bout, nchar(p$siteID) + 2L), fixed = TRUE)
  st <- x$fsh_fieldData[!duplicated(x$fsh_fieldData$siteID), c("siteID", "aquaticSiteType")]
  censored <- W <= s$floor_max
  counts["censored"] <- sum(censored)
  counts["at_0.3"] <- sum(abs(W - 0.3) < 1e-9)
  # Equations of species and of genera.
  eq_of <- function(group) {
    idx <- split(seq_along(group), group)
    out <- lapply(idx, function(i) neon_equation(L[i], W[i], censored[i], s))
    out <- out[!vapply(out, is.null, logical(1))]
    if (!length(out)) return(NULL)
    data.frame(taxon = names(out), do.call(rbind, lapply(out, as.data.frame)),
               stringsAsFactors = FALSE, row.names = NULL)
  }
  eq_sp <- eq_of(taxon)
  eq_ge <- eq_of(genus)
  i_sp <- match(taxon, eq_sp$taxon)
  i_ge <- match(genus, eq_ge$taxon)
  use_ge <- is.na(i_sp) & !is.na(i_ge)
  pick <- function(col) ifelse(!is.na(i_sp), eq_sp[[col]][i_sp], eq_ge[[col]][i_ge])
  ln_a <- pick("ln_a")
  b <- pick("b")
  sigma <- pick("sigma")
  mu <- ln_a + b * log(L)
  z <- (log(W) - mu) / sigma
  has_eq <- is.finite(mu)
  gross <- has_eq & !censored & abs(log(W) - mu) > log(s$gross_ratio) & abs(z) > s$gross_sd
  above <- (mu - log(s$limit))
  gross <- gross | (has_eq & censored & above > log(s$gross_ratio) & above / sigma > s$gross_sd)
  counts["gross_uncensored"] <- sum(gross & !censored)
  counts["gross_censored"] <- sum(gross & censored)
  # Mass of the censored fish: mean of the lognormal truncated at the limit.
  s2 <- sigma^2
  tmean <- exp(mu + s2 / 2) * stats::pnorm((log(s$limit) - mu - s2) / sigma) /
    stats::pnorm((log(s$limit) - mu) / sigma)
  mass <- ifelse(censored, tmean, W)
  imputed <- censored & is.finite(mass) & !gross
  counts["imputed_species"] <- sum(imputed & !use_ge)
  counts["imputed_genus"] <- sum(imputed & use_ge)
  counts["censored_no_mass"] <- sum(censored & !gross & !is.finite(mass))
  mass[!is.finite(mass)] <- NA_real_
  keep <- !gross
  counts["kept"] <- sum(keep)
  floor_length <- ifelse(!is.na(i_sp), eq_sp$floor_length[i_sp], NA_real_)
  fish <- data.frame(taxon = taxon, name_as_given = p$scientificName, site = p$siteID,
                     site_type = st$aquaticSiteType[match(p$siteID, st$siteID)], bout = bout,
                     gear = gear, life_stage = p$fishLifeStage, length = L, weight = W,
                     censored = censored, mass = mass, imputed = imputed,
                     for_lwr = !censored & is.finite(L) & !is.na(floor_length) & L >= floor_length,
                     stringsAsFactors = FALSE)[keep, ]
  rownames(fish) <- NULL
  counts["for_lwr"] <- sum(fish$for_lwr)
  list(fish = fish, equations = eq_sp, genus_equations = eq_ge, counts = counts,
       citation = x$citation)
}
