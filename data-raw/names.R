# Taxon name normalisation used by data-raw/references.R (from R/names.R of
# version 0.1.0, which normalised the names of the user and of the reference
# tables in the same way). Every change is reported in `note`; nothing is
# matched approximately. Sourced by references.R; base R only.

# Ranks used by bmassa, from the lowest to the highest.
rank_levels <- c(
  "species", "genus", "subfamily", "family", "suborder", "order",
  "subclass", "class", "subphylum", "phylum", "kingdom"
)

# Qualifiers that make a species-level identification uncertain or partial.
# When one of them follows the genus, the name is reduced to the genus.
name_qualifiers <- c(
  "cf.", "cf", "aff.", "aff", "gr.", "gr", "group", "agg.", "agg",
  "s.l.", "s.str.", "sensu", "nr.", "nr", "near", "?"
)

# Words that mean "some species of this genus".
name_unspecific <- c(
  "sp.", "sp", "spp.", "spp", "sp.nov.", "indet.", "indet", "ssp.", "ssp",
  "species", "undet.", "undet"
)

# Markers of infraspecific ranks.
name_infraspecific <- c("var.", "var", "f.", "fo.", "forma", "subsp.", "subsp")

# Words for life stages that some data sets append to the taxon name. They
# are removed from the name (the stage belongs in its own column).
name_stage_words <- c(
  "larva", "larvae", "larval", "nymph", "nymphs", "pupa", "pupae", "adult",
  "adults", "imago", "imagines", "subimago", "juvenile", "juveniles",
  "nauplius", "nauplii", "copepodite", "copepodites", "copepodid", "copepodids",
  "female", "females", "male", "males", "ovigerous", "ephippial", "egg", "eggs",
  "instar"
)

# Suffixes of names above the genus (family, subfamily, tribe, superfamily,
# and algal classes, orders and divisions): words after such a name are never
# species epithets.
suprageneric_pattern <- "(idae|inae|ini|oidea|aceae|ales|phyceae|phyta|mycota)$"

# Latin-1 letters with diacritics (Unicode code points) and their ASCII
# equivalents. Working on code points (utf8ToInt()) makes the result
# independent of the locale, and keeps the R code ASCII-only.
diacritics_codes <- c(
  0xE1, 0xE0, 0xE2, 0xE4, 0xE3, 0xE5, 0xE9, 0xE8, 0xEA, 0xEB, 0xED, 0xEC, 0xEE, 0xEF,
  0xF3, 0xF2, 0xF4, 0xF6, 0xF5, 0xF8, 0xFA, 0xF9, 0xFB, 0xFC, 0xF1, 0xE7, 0xFD, 0xFF,
  0xC1, 0xC0, 0xC2, 0xC4, 0xC3, 0xC5, 0xC9, 0xC8, 0xCA, 0xCB, 0xCD, 0xCC, 0xCE, 0xCF,
  0xD3, 0xD2, 0xD4, 0xD6, 0xD5, 0xD8, 0xDA, 0xD9, 0xDB, 0xDC, 0xD1, 0xC7, 0xDD
)
diacritics_ascii <- utf8ToInt(paste0(
  "aaaaaaeeeeiiii", "oooooouuuuncyy", "AAAAAAEEEEIIII", "OOOOOOUUUUNCY"
))

remove_diacritics <- function(s) {
  s <- as.character(s)
  # Only strings with bytes outside printable ASCII need the replacements.
  todo <- which(!is.na(s) & grepl("[^ -~]", s, useBytes = TRUE))
  for (i in todo) {
    x <- s[i]
    if (Encoding(x) == "latin1") {
      x <- iconv(x, from = "latin1", to = "UTF-8")
    } else if (Encoding(x) == "unknown" && !validUTF8(x)) {
      # Native encoding other than UTF-8 (e.g. Latin-1 on older Windows).
      x <- iconv(x, from = "", to = "UTF-8")
    }
    if (is.na(x) || !validUTF8(x)) next
    cp <- utf8ToInt(x)
    hit <- match(cp, diacritics_codes)
    if (all(is.na(hit))) {
      s[i] <- x
      next
    }
    cp[!is.na(hit)] <- diacritics_ascii[hit[!is.na(hit)]]
    s[i] <- intToUtf8(cp)
  }
  s
}

# Normalise taxon names.
#
# Returns a data frame with one row per element of `x` (same order) and
# columns `taxon_clean` (character; NA for missing or empty names) and
# `name_note` (character; NA when the name was used as given, apart from
# white space).
normalize_taxon_names <- function(x) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  x <- as.character(x)
  # Each distinct name is normalised once.
  ux <- unique(x)
  clean <- rep(NA_character_, length(ux))
  note <- rep(NA_character_, length(ux))
  for (i in seq_along(ux)) {
    res <- normalize_one_name(ux[i])
    clean[i] <- res[["clean"]]
    note[i] <- res[["note"]]
  }
  idx <- match(x, ux)
  data.frame(taxon_clean = clean[idx], name_note = note[idx], stringsAsFactors = FALSE)
}

normalize_one_name <- function(name) {
  if (is.na(name)) {
    return(list(clean = NA_character_, note = NA_character_))
  }
  notes <- character()
  s <- gsub("[[:space:]]+", " ", name)
  s <- gsub("^ | $", "", s)
  if (!nzchar(s)) {
    return(list(clean = NA_character_, note = "empty name"))
  }
  ascii <- remove_diacritics(s)
  if (!identical(ascii, s)) {
    notes <- c(notes, "diacritics removed")
    s <- ascii
  }
  # Subgenus in parentheses right after the genus: "Daphnia (Ctenodaphnia) magna".
  if (grepl("^[A-Za-z-]+ \\([A-Z][a-z-]+\\)( |$)", s)) {
    s <- sub("^([A-Za-z-]+) \\([A-Z][a-z-]+\\)", "\\1", s)
    notes <- c(notes, "subgenus removed")
  }
  words <- strsplit(s, " ", fixed = TRUE)[[1]]
  first <- words[1]
  # A first word written in capitals ("CHIRONOMIDAE") is converted to the
  # usual form ("Chironomidae").
  if (nchar(first) > 1L && grepl("^[A-Z-]+$", first)) {
    first <- paste0(substr(first, 1L, 1L), tolower(substr(first, 2L, nchar(first))))
    notes <- c(notes, "capitals converted")
  }
  if (grepl("/", first, fixed = TRUE)) {
    return(list(clean = NA_character_,
                note = paste(c(notes, "alternative names separated by '/'"), collapse = "; ")))
  }
  if (!grepl("^[A-Z][a-z]+(-[a-z]+)?$", first)) {
    return(list(clean = NA_character_,
                note = paste(c(notes, sprintf("'%s' is not a scientific name", first)),
                             collapse = "; ")))
  }
  rest <- words[-1]
  if (length(rest) == 0L) {
    return(list(clean = first, note = collapse_notes(notes)))
  }
  second <- rest[1]
  second_low <- tolower(second)
  if (second_low %in% name_stage_words) {
    notes <- c(notes, sprintf("life-stage words removed ('%s'); give stages in their own column",
                              paste(rest, collapse = " ")))
    return(list(clean = first, note = collapse_notes(notes)))
  }
  if (grepl(suprageneric_pattern, first)) {
    notes <- c(notes, sprintf("text after '%s' removed ('%s')", first,
                              paste(rest, collapse = " ")))
    return(list(clean = first, note = collapse_notes(notes)))
  }
  if (second_low %in% name_unspecific) {
    if (length(rest) > 1L) {
      notes <- c(notes, sprintf("words after '%s' removed", second))
    }
    notes <- c(notes, sprintf("'%s' removed; genus-level name", second))
    return(list(clean = first, note = collapse_notes(notes)))
  }
  if (second_low %in% name_qualifiers) {
    notes <- c(notes, sprintf("qualifier '%s' makes the species uncertain; genus-level name", second))
    return(list(clean = first, note = collapse_notes(notes)))
  }
  if (grepl("/", second, fixed = TRUE)) {
    notes <- c(notes, sprintf("alternative epithets '%s'; genus-level name", second))
    return(list(clean = first, note = collapse_notes(notes)))
  }
  if (!grepl("^[a-z]+(-[a-z]+)?$", second)) {
    # Authority or other text after a uninomial ("Chironomidae Macquart").
    notes <- c(notes, sprintf("text after '%s' removed ('%s')", first,
                              paste(rest, collapse = " ")))
    return(list(clean = first, note = collapse_notes(notes)))
  }
  clean <- paste(first, second)
  extra <- rest[-1]
  if (length(extra) > 0L) {
    if (tolower(extra[1]) %in% name_stage_words) {
      notes <- c(notes, sprintf("life-stage words removed ('%s'); give stages in their own column",
                                paste(extra, collapse = " ")))
    } else if (tolower(extra[1]) %in% name_infraspecific || grepl("^[a-z]+(-[a-z]+)?$", extra[1])) {
      notes <- c(notes, sprintf("infraspecific name removed ('%s')", paste(extra, collapse = " ")))
    } else {
      notes <- c(notes, sprintf("authority removed ('%s')", paste(extra, collapse = " ")))
    }
  }
  list(clean = clean, note = collapse_notes(notes))
}

collapse_notes <- function(notes) {
  if (length(notes) == 0L) NA_character_ else paste(notes, collapse = "; ")
}

# Rank of a name from its form only: two words = species. Uninomials are
# not ranked here (they need the taxonomy table).
is_binomial <- function(x) {
  !is.na(x) & grepl("^[A-Z][a-z]+(-[a-z]+)? [a-z]+(-[a-z]+)?$", x)
}

genus_of_binomial <- function(x) {
  out <- rep(NA_character_, length(x))
  ok <- is_binomial(x)
  out[ok] <- sub(" .*$", "", x[ok])
  out
}
