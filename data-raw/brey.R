# Reader of a local copy of the data bank of Brey et al. (2010), used by
# data-raw/conversions.R (it was the function read_brey_databank() of bmassa
# 0.1.0). The data bank is not distributed: only the parameters of the
# conversion models fitted to it are.
# Packages used only here: readxl.

# Ratios read from the data bank (labels of rows 2 and 3 of the sheet "Data")
# and the conversions they give: value in `to` = factor x value in `from`,
# with masses in micrograms and energy in joules. Nitrogen and phosphorus per
# dry mass from version 0.6.0.
brey_ratios <- data.frame(
  label = c("DM / WM", "AFDM / DM", "AFDM / WM", "C / DM", "C / AFDM",
            "J / mgWM", "J / mgDM", "J / mgAFDM", "J / mgC", "N / DM", "P / DM"),
  from = c("wet_mass", "dry_mass", "wet_mass", "dry_mass", "afdm",
           "wet_mass", "dry_mass", "afdm", "carbon", "dry_mass", "dry_mass"),
  to = c("dry_mass", "afdm", "afdm", "carbon", "carbon",
         "energy", "energy", "energy", "energy", "nitrogen", "phosphorus"),
  scale = c(1, 1, 1, 1, 1, 1e-3, 1e-3, 1e-3, 1e-3, 1, 1), # J/mg -> J/ug
  invertible = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

# One row per record, ratio and direction (mass ratios also inverted), with
# columns taxon, from, to, factor, source, source_detail, habitat, genus,
# family, order, class and kingdom. `habitat`: "freshwater", "marine" or
# "both".
read_brey_databank <- function(path, habitat, quiet = FALSE) {
  stopifnot(is.character(path), length(path) == 1L,
            habitat %in% c("freshwater", "marine", "both"), length(habitat) == 1L)
  if (!file.exists(path)) stop("File '", path, "' does not exist.")
  if (!"Data" %in% readxl::excel_sheets(path)) {
    stop("'", path, "' has no sheet called 'Data'; is it the conversion data bank of ",
         "Brey et al. (2010)?")
  }
  raw <- readxl::read_excel(path, sheet = "Data", col_names = FALSE, col_types = "text",
                            .name_repair = "minimal")
  parse_brey_data(as.data.frame(raw, stringsAsFactors = FALSE), habitat, quiet = quiet)
}

# Numbers from text; cells with text that is not a number become NA and are
# counted (not silenced).
parse_numbers <- function(v) {
  v <- trimws(as.character(v))
  num <- grepl("^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$", v)
  out <- rep(NA_real_, length(v))
  out[num] <- as.numeric(v[num])
  attr(out, "n_unreadable") <- sum(!is.na(v) & nzchar(v) & !num)
  out
}

# Parse the sheet "Data" read as text without column names. Layout of version
# 3-2009: rows 1-5 are headers (row 2 and 3 hold the ratio labels), data start
# in row 6; columns 3-10 hold species, genus, family, superfamily, order,
# class, major taxon and kingdom, column 12 the habitat (1 marine, 0
# freshwater).
parse_brey_data <- function(raw, habitat, quiet = TRUE) {
  if (ncol(raw) < 12L || nrow(raw) < 6L) {
    stop("The sheet 'Data' does not have the layout of the Brey et al. (2010) data bank.")
  }
  labels <- gsub("\\s+", " ", trimws(paste(raw[2, ], raw[3, ])))
  col_of <- match(brey_ratios$label, labels)
  if (anyNA(col_of)) {
    stop("Columns not found in the sheet 'Data': ",
         paste(brey_ratios$label[is.na(col_of)], collapse = ", "),
         ". Is it version 3-2009 of the data bank?")
  }
  ref_col <- which(trimws(as.character(raw[5, ])) == "References")[1L]
  d <- raw[-(1:5), , drop = FALSE]
  hab <- parse_numbers(d[[12L]])
  unreadable <- attr(hab, "n_unreadable")
  keep <- switch(habitat,
    freshwater = hab %in% 0,
    marine = hab %in% 1,
    both = hab %in% c(0, 1)
  )
  d <- d[keep, , drop = FALSE]
  hab <- hab[keep]
  na_text <- function(v) {
    v <- trimws(as.character(v))
    v[!is.na(v) & !nzchar(v)] <- NA_character_
    v
  }
  species <- na_text(d[[3L]])
  genus <- na_text(d[[4L]])
  taxon <- ifelse(is.na(species), genus, species)
  lineage <- data.frame(genus = genus, family = na_text(d[[5L]]), order = na_text(d[[7L]]),
                        class = na_text(d[[8L]]), kingdom = na_text(d[[10L]]),
                        stringsAsFactors = FALSE)
  detail <- paste0("Brey et al. 2010, record ", na_text(d[[1L]]),
                   if (!is.na(ref_col)) paste0(": ", na_text(d[[ref_col]])) else "")
  pieces <- list()
  for (k in seq_len(nrow(brey_ratios))) {
    v <- parse_numbers(d[[col_of[k]]])
    unreadable <- unreadable + attr(v, "n_unreadable")
    ok <- !is.na(v) & is.finite(v) & v > 0 & !is.na(taxon)
    if (!any(ok)) next
    base <- data.frame(taxon = taxon[ok], stringsAsFactors = FALSE)
    one <- function(from, to, factor) {
      cbind(base, from = from, to = to, factor = factor, source = "brey2010",
            source_detail = detail[ok],
            habitat = ifelse(hab[ok] == 1, "marine", "freshwater"),
            lineage[ok, , drop = FALSE], stringsAsFactors = FALSE)
    }
    pieces[[length(pieces) + 1L]] <- one(brey_ratios$from[k], brey_ratios$to[k],
                                         v[ok] * brey_ratios$scale[k])
    if (brey_ratios$invertible[k]) {
      pieces[[length(pieces) + 1L]] <- one(brey_ratios$to[k], brey_ratios$from[k], 1 / v[ok])
    }
  }
  if (length(pieces) == 0L) {
    out <- data.frame(taxon = character(), from = character(), to = character(),
                      factor = numeric(), source = character(), source_detail = character(),
                      habitat = character(), genus = character(), family = character(),
                      order = character(), class = character(), kingdom = character(),
                      stringsAsFactors = FALSE)
  } else {
    out <- do.call(rbind, pieces)
  }
  rownames(out) <- NULL
  if (!quiet && unreadable > 0L) {
    message(sprintf("%d cell%s of the data bank could not be read as numbers and %s left out.",
                    unreadable, if (unreadable == 1L) "" else "s",
                    if (unreadable == 1L) "was" else "were"))
  }
  out
}
