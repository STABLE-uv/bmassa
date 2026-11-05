# Taxonomic standardization with the Catalogue of Life (COL XR) ----------------------

colxr_checklist_key <- "7ddf754f-d193-4cc9-b351-99906754a03b"
colxr_match_url <- "https://api.gbif.org/v2/species/match"
# Answers received in this session, by URL, so that a repeated call (or the
# rest of a call that failed half-way) does not send the same queries again.
colxr_memory <- new.env(parent = emptyenv())
hint_ranks <- c("kingdom", "phylum", "class", "order", "family")
name_suffix_rank <- c(idae = "FAMILY", aceae = "FAMILY", inae = "SUBFAMILY", oideae = "SUBFAMILY",
                      ini = "TRIBE", oidea = "SUPERFAMILY", ales = "ORDER", ineae = "SUBORDER",
                      phyceae = "CLASS", opsida = "CLASS", phyta = "PHYLUM")
stage_words <- c("larva", "larvae", "nymph", "nymphs", "pupa", "pupae", "adult", "adults",
                 "juvenile", "juveniles", "male", "males", "female", "females", "copepodite",
                 "copepodites", "nauplius", "nauplii", "imago", "imagines")
group_hints <- list(phytoplankton = character(0), zooplankton = c(kingdom = "Animalia"),
                    macroinvertebrate = c(kingdom = "Animalia"),
                    fish = c(kingdom = "Animalia", phylum = "Chordata"),
                    amphibian = c(kingdom = "Animalia", phylum = "Chordata", class = "Amphibia"))
synonym_status <- c("SYNONYM", "HETEROTYPIC_SYNONYM", "HOMOTYPIC_SYNONYM", "PROPARTE_SYNONYM")

#' Classification of a list of taxa in the Catalogue of Life
#'
#' `standardize_taxa()` looks up each name of `taxa` in the Catalogue of Life
#' eXtended Release (COL XR) through the name-matching service of GBIF, and
#' returns the table with the accepted name, the identifier of the COL
#' (`col_id`) and the classification by main ranks that [estimate_biomass()]
#' and [check_coverage()] use. It needs internet; check the result before
#' using it, in particular the rows with a `note`.
#'
#' @param taxa A data frame with one row per taxon: `taxon_id` (unique),
#'   `broad_group` (`"phytoplankton"`, `"zooplankton"`, `"macroinvertebrate"`,
#'   `"fish"` or `"amphibian"`) and `scientific_name` (the name as
#'   identified: a species, a genus, `"Genus sp."`, a family...). Optional
#'   columns `kingdom`, `phylum`, `class`, `order` and `family` are sent as
#'   hints that tell homonyms apart. Other columns are kept.
#' @param quiet If `FALSE` (default), a message says how many names are
#'   queried.
#'
#' @return `taxa` (as a data frame, same rows and order) with `scientific_name`
#'   unchanged and these columns set (replacing columns of the same name):
#'   * `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`:
#'     the classification of the accepted taxon in the COL XR, `NA` below its
#'     rank; for names not found, the values given in `taxa` (if any).
#'   * `col_id`, `col_name`, `col_rank`: identifier, accepted name and rank of
#'     the taxon in the COL XR.
#'   * `col_status`: status of the name given, e.g. `"accepted"` or
#'     `"synonym"` (then `col_name` is the accepted name).
#'   * `match_type`: `"exact"`, `"variant"` (the epithet differs only in its
#'     ending, e.g. *plumosa* and *plumosus*), `"higher_rank"` (the name was
#'     not accepted and the taxon is placed in its genus) or `"none"`.
#'   * `match_confidence`: confidence of the matching service (0-100).
#'   * `col_release`: version of the COL XR (e.g. `"COL26.6 XR"`).
#'   * `note`: what was changed or should be checked (qualifiers removed, a
#'     rank given in `taxa` that differs from the COL XR, the closest name to
#'     one that was not accepted...).
#'
#' @details
#' The name sent is the name given without subgenus, without the words of life
#' stages ("larvae", "males"...) and, for infraspecific names, reduced to the
#' species (the models of bmassa stop at species). "Genus sp.", "Genus
#' spp." and names with "cf." or "aff." are looked up as the genus: the
#' identification is uncertain beyond it. The broad group gives hints
#' (Animalia for animals; Chordata for fish and amphibians; Amphibia for
#' amphibians) and the answer must agree with it (phytoplankton outside
#' Animalia). Among homonyms, the taxa expected in the broad group come first:
#' in the zooplankton, rotifers, branchiopods, copepods and ostracods, then
#' other arthropods (the cladoceran *Latona*, not the bivalve); in the
#' macroinvertebrates, the phyla of freshwater invertebrates. The other taxa
#' with the same name are listed in `note`.
#'
#' An answer is accepted when it carries the same name (or the same epithet
#' with another ending, a matter of grammatical gender) at the rank the name
#' implies. Synonyms are replaced by the accepted name. Ambiguous synonyms and
#' misapplied names are followed only when they lead to the same species. A
#' name that is not accepted is retried with the rank it implies and with
#' fewer hints; a species still not accepted is placed in its genus
#' (`match_type = "higher_rank"`), and the note gives the closest name found,
#' if any: a misspelt name is not corrected automatically. Names that match
#' nothing keep the classification given in `taxa`, and a warning says how
#' many there are.
#'
#' @section Internet:
#' The names are sent to the GBIF API (<https://api.gbif.org>), four at a
#' time, with up to two retries. If the service cannot be reached, the
#' function stops with an error of class `bmassa_error_network`; answers
#' already received are kept for the rest of the session and are not asked
#' again.
#'
#' @section Reproducibility:
#' The COL XR is updated about monthly: a new release may change accepted
#' names and classifications, and `col_release` records the one used. Keep the
#' result with your data rather than querying again.
#'
#' @references
#' Catalogue of Life eXtended Release, through the GBIF species matching API
#' (version 2), <https://techdocs.gbif.org/en/openapi/v1/species>.
#'
#' @seealso [check_coverage()] to see from which level of the reference each
#'   taxon would then be estimated.
#' @export
#' @examplesIf requireNamespace("curl", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE) && curl::has_internet()
#' taxa <- data.frame(taxon_id = c("a", "b", "c"), broad_group = "macroinvertebrate",
#'                    scientific_name = c("Baetis rhodani", "Ephemerella ignita",
#'                                        "Chironomus sp."))
#' # If the GBIF service cannot be reached, the error says so (class
#' # "bmassa_error_network"); here it becomes a message.
#' std <- tryCatch(standardize_taxa(taxa), bmassa_error_network = function(e) {
#'   message(conditionMessage(e))
#'   NULL
#' })
#' if (!is.null(std)) std[, c("scientific_name", "col_name", "family", "match_type", "note")]
standardize_taxa <- function(taxa, quiet = FALSE) {
  check_data_frame(taxa)
  check_flag(quiet)
  check_columns(taxa, c("taxon_id", "broad_group", "scientific_name"), arg = "taxa")
  tx <- as_plain_data_frame(taxa)
  for (col in c("taxon_id", "broad_group", "scientific_name")) {
    check_text_column(tx, col, "taxa")
    tx[[col]] <- as_text(tx[[col]])
  }
  check_unique_rows(tx, "taxon_id", "taxa")
  check_values_in(tx, "broad_group", broad_group_levels, allow_na = FALSE, arg = "taxa")
  for (r in main_ranks) {
    if (!r %in% names(tx)) {
      tx[[r]] <- rep(NA_character_, nrow(tx))
      next
    }
    check_text_column(tx, r, "taxa", allow_na = TRUE)
    v <- trimws(gsub("\\s+", " ", as_text(tx[[r]])))
    v[!is.na(v) & !nzchar(v)] <- NA_character_
    tx[[r]] <- v
  }
  out <- tx
  n <- nrow(tx)
  chr <- rep(NA_character_, n)
  res <- data.frame(col_id = chr, col_name = chr, col_rank = chr, col_status = chr,
                    match_type = rep("none", n), match_confidence = rep(NA_integer_, n),
                    col_release = chr, note = chr, stringsAsFactors = FALSE)
  if (n == 0L) {
    out[names(res)] <- res
    return(out)
  }
  check_installed(c("curl", "jsonlite"), purpose = "query the Catalogue of Life")
  nm <- parse_taxon_names(tx$scientific_name)
  res$note <- nm$note
  if (!quiet) {
    message(sprintf("Querying the Catalogue of Life (COL XR) through the GBIF API for %d names...",
                    length(unique(paste(nm$target, tx$broad_group)))))
  }
  hints <- lapply(seq_len(n), function(i) {
    h <- unlist(tx[i, hint_ranks])
    h <- h[!is.na(h)]
    g <- group_hints[[tx$broad_group[i]]]
    c(h, g[setdiff(names(g), names(h))])
  })
  group_only <- lapply(tx$broad_group, function(g) group_hints[[g]])
  state <- list(answer = vector("list", n), first = vector("list", n))
  state <- colxr_rounds(state, seq_len(n), nm$query, nm$rank, nm$form, tx$broad_group, hints,
                        group_only, record_first = TRUE)
  # Species not accepted: their genus.
  to_genus <- which(vapply(state$answer, is.null, logical(1)) &
                      nm$form %in% c("species", "infraspecific"))
  if (length(to_genus)) {
    genus <- sub(" .*$", "", nm$target)
    for (i in to_genus) {
      if (!is.null(state$first[[i]])) {
        state$answer[i] <- list(genus_from_answer(state$first[[i]], genus[i], tx$broad_group[i]))
      }
    }
    state <- colxr_rounds(state, to_genus, genus, rep("GENUS", n), rep("genus_of_species", n),
                          tx$broad_group, hints, group_only, record_first = FALSE)
  }
  answer <- state$answer
  first <- state$first
  release <- colxr_release()
  for (i in seq_len(n)) {
    a <- answer[[i]]
    if (is.null(a)) {
      res$note[i] <- append_note(res$note[i], paste0("not found in the COL XR",
                                                     closest_name(first[[i]], nm$target[i])))
      next
    }
    cls <- a$classification
    for (r in main_ranks) out[[r]][i] <- unname(cls[toupper(r)])
    syn <- a$status %in% synonym_status || (a$status %in% c("AMBIGUOUS_SYNONYM", "MISAPPLIED"))
    res$col_id[i] <- if (syn && !is.na(a$accepted_key)) a$accepted_key else a$key
    res$col_name[i] <- if (syn && !is.na(a$accepted_name)) a$accepted_name else a$name
    res$col_rank[i] <- tolower(if (syn && !is.na(a$accepted_rank)) a$accepted_rank else a$rank)
    res$col_status[i] <- tolower(a$status)
    res$match_confidence[i] <- a$confidence
    res$col_release[i] <- release
    is_genus_of_species <- nm$form[i] %in% c("species", "infraspecific") && a$rank == "GENUS"
    res$match_type[i] <- if (is_genus_of_species) "higher_rank" else
      if (identical(tolower(a$name), tolower(nm$target[i]))) "exact" else "variant"
    if (is_genus_of_species) {
      res$note[i] <- append_note(res$note[i], paste0(
        "species not found in the COL XR: placed in its genus", closest_name(first[[i]], nm$target[i])))
    }
    if (res$match_type[i] == "variant") {
      res$note[i] <- append_note(res$note[i], sprintf("spelling: %s -> %s", nm$target[i], a$name))
    }
    if (!is.null(a$homonyms)) res$note[i] <- append_note(res$note[i], a$homonyms)
    if (syn) res$note[i] <- append_note(res$note[i], sprintf("%s is a %s of %s", a$name,
                                                             gsub("_", " ", tolower(a$status)),
                                                             res$col_name[i]))
    given <- unlist(tx[i, hint_ranks])
    differ <- hint_ranks[!is.na(given) & !is.na(cls[toupper(hint_ranks)]) &
                           given != cls[toupper(hint_ranks)]]
    if (length(differ)) {
      res$note[i] <- append_note(res$note[i], paste(sprintf(
        "%s given %s, COL XR %s", differ, given[differ], cls[toupper(differ)]), collapse = "; "))
    }
  }
  out[names(res)] <- res
  rownames(out) <- NULL
  n_none <- sum(res$match_type == "none")
  if (n_none > 0L) {
    warn_not_estimable(sprintf(
      "%d of %d names were not found in the COL XR and keep the classification given; see `note`.",
      n_none, n))
  }
  out
}

# Name to send, rank it implies and name the answer must carry -----------------------

parse_taxon_names <- function(x) {
  y <- gsub("\\s+", " ", trimws(x))
  note <- rep(NA_character_, length(y))
  # Subgenus in parentheses.
  sg <- grepl("^[A-Z][a-z-]+ \\([A-Z][a-z-]+\\) ", y)
  y[sg] <- sub("^([A-Z][a-z-]+) \\([A-Z][a-z-]+\\) ", "\\1 ", y[sg])
  # Authorship: words after the name that start with a capital or "(".
  auth <- grepl("^[A-Z][a-z-]+( [a-z-]+){0,2} [A-Z(].*$", y) & !grepl("^[A-Z][a-z-]+ \\(", y)
  y[auth] <- sub("^([A-Z][a-z-]+( [a-z-]+){0,2}) [A-Z(].*$", "\\1", y[auth])
  note[auth] <- append_note(note[auth], "authorship removed")
  # Words of life stages.
  words <- strsplit(y, " ", fixed = TRUE)
  st <- vapply(words, function(w) any(tolower(w[-1]) %in% stage_words), logical(1))
  y[st] <- vapply(words[st], function(w) {
    paste(w[!(tolower(w) %in% stage_words) | seq_along(w) == 1L], collapse = " ")
  }, character(1))
  note[st] <- append_note(note[st], "life-stage words removed")
  qual <- grepl("^[A-Z][a-z-]+ (cf|aff)\\.?( |$)", y)
  gsp <- grepl("^[A-Z][a-z-]+ (sp|spp)(\\.|[0-9]| |$)", y)
  note[qual] <- append_note(note[qual], "cf./aff.: identified to genus")
  genus_part <- sub(" .*$", "", y)
  binom <- grepl("^[A-Z][a-z-]+ [a-z-]+$", y) & !qual & !gsp
  infra <- grepl("^[A-Z][a-z-]+ [a-z-]+ ((var|subsp|ssp|f|fo|forma)\\.? )?[a-z-]+$", y) & !qual &
    !gsp & !binom
  species_part <- sub("^([A-Z][a-z-]+ [a-z-]+).*$", "\\1", y)
  note[infra] <- append_note(note[infra], "infraspecific name reduced to its species")
  single <- grepl("^[A-Z][a-z-]+$", y)
  suf <- vapply(y, function(w) {
    hit <- names(name_suffix_rank)[endsWith(w, names(name_suffix_rank))]
    if (length(hit)) name_suffix_rank[[hit[which.max(nchar(hit))]]] else NA_character_
  }, character(1), USE.NAMES = FALSE)
  form <- ifelse(qual | gsp, "genus_sp", ifelse(binom, "species", ifelse(infra, "infraspecific",
                 ifelse(single & !is.na(suf), "suffix", ifelse(single, "single", "other")))))
  target <- ifelse(form == "genus_sp", genus_part, ifelse(form == "infraspecific", species_part, y))
  rank <- ifelse(form %in% c("genus_sp", "single"), "GENUS",
                 ifelse(form %in% c("species", "infraspecific"), "SPECIES",
                        ifelse(form == "suffix", suf, NA_character_)))
  data.frame(query = target, target = target, form = form, rank = rank, note = note,
             stringsAsFactors = FALSE)
}

# Queries --------------------------------------------------------------------------

# Attempts, in order, for the rows still without an accepted answer: the name
# with all the hints; with the rank it implies; with the hints of the broad
# group only; with no hints. Each answer is checked with its alternatives.
colxr_rounds <- function(state, rows, name, rank, form, group, hints, group_only, record_first) {
  for (a in 1:4) {
    rows <- rows[vapply(state$answer[rows], is.null, logical(1))]
    if (a > 1L) rows <- rows[!is.na(rank[rows])]
    if (!length(rows)) break
    urls <- vapply(rows, function(i) {
      h <- switch(a, hints[[i]], hints[[i]], group_only[[i]], character(0))
      colxr_url(name[i], h, if (a == 1L) NA_character_ else rank[i])
    }, character(1))
    got <- lapply(colxr_get(urls), parse_colxr_answer)
    for (k in seq_along(rows)) {
      i <- rows[k]
      if (record_first && is.null(state$first[[i]])) state$first[[i]] <- got[[k]]
      ok <- accept_colxr(got[[k]], name[i], form[i], group[i])
      if (!is.null(ok)) state$answer[[i]] <- ok
    }
  }
  state
}


colxr_url <- function(name, hints, rank) {
  p <- c(scientificName = name, checklistKey = colxr_checklist_key, hints,
         taxonRank = if (!is.na(rank)) rank, verbose = "true")
  paste0(colxr_match_url, "?", paste(names(p), vapply(p, utils::URLencode, "", reserved = TRUE),
                                     sep = "=", collapse = "&"))
}

# Answers for a set of URLs, from this session's memory or from the API.
colxr_get <- function(urls) {
  new <- unique(urls[!vapply(urls, exists, logical(1), envir = colxr_memory, inherits = FALSE)])
  if (length(new)) {
    body <- colxr_fetch(new)
    for (k in seq_along(new)) assign(new[k], body[k], envir = colxr_memory)
  }
  vapply(urls, get, character(1), envir = colxr_memory, inherits = FALSE, USE.NAMES = FALSE)
}

# Fetch the URLs, four at a time, with up to two retries of those that fail
# (after 4 and 6 seconds).
colxr_fetch <- function(urls, tries = 3L, wait = 2) {
  st <- new.env(parent = emptyenv())
  st$body <- rep(NA_character_, length(urls))
  st$status <- rep(NA_integer_, length(urls))
  todo <- seq_along(urls)
  for (attempt in seq_len(tries)) {
    if (attempt > 1L) Sys.sleep(wait * attempt)
    pool <- curl::new_pool(total_con = 4, host_con = 4)
    for (k in todo) {
      local({
        kk <- k
        h <- curl::new_handle(useragent = "bmassa R package (standardize_taxa)", timeout = 30)
        curl::curl_fetch_multi(urls[kk], handle = h, pool = pool,
                               done = function(r) {
                                 st$status[kk] <- r$status_code
                                 st$body[kk] <- rawToChar(r$content)
                               },
                               fail = function(msg) st$status[kk] <- -1L)
      })
    }
    curl::multi_run(pool = pool)
    todo <- which(is.na(st$status) | st$status != 200L)
    if (!length(todo)) break
  }
  if (length(todo)) {
    codes <- unique(st$status[todo])
    abort_network(sprintf(paste(
      "Could not get answers from the GBIF name-matching service (%s) for %d of %d queries",
      "(%s). Check the internet connection and run the call again: the answers already",
      "received are kept for this session."), colxr_match_url, length(todo), length(urls),
      paste(ifelse(codes < 0, "no connection", paste("HTTP status", codes)), collapse = ", ")),
      call = NULL)
  }
  body <- st$body
  Encoding(body) <- "UTF-8"
  body
}

# Version of the COL XR used by the service (once per session).
colxr_release <- function() {
  if (!exists("release", envir = colxr_memory, inherits = FALSE)) {
    txt <- colxr_fetch(paste0(colxr_match_url, "/metadata?checklistKey=", colxr_checklist_key))
    meta <- colxr_json(txt)
    alias <- meta$mainIndex$datasetAlias
    assign("release", if (is.null(alias)) NA_character_ else as.character(alias), envir = colxr_memory)
  }
  get("release", envir = colxr_memory, inherits = FALSE)
}

# An answer that cannot be read as JSON (the service failed or changed) is a
# network error, not a crash.
colxr_json <- function(txt) {
  tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) {
    abort_network(sprintf(paste(
      "The GBIF name-matching service (%s) gave an answer that could not be read (%s).",
      "Run the call again later."), colxr_match_url, conditionMessage(e)), call = NULL)
  })
}

parse_colxr_answer <- function(txt) {
  r <- colxr_json(txt)
  list(main = parse_colxr_usage(r),
       alternatives = lapply(r$diagnostics$alternatives, parse_colxr_usage))
}

parse_colxr_usage <- function(r) {
  g <- function(x, f) if (is.null(x) || is.null(x[[f]])) NA_character_ else as.character(x[[f]])
  cls <- stats::setNames(rep(NA_character_, length(main_ranks)), toupper(main_ranks))
  keys <- cls
  for (z in r$classification) {
    if (!is.null(z$rank) && z$rank %in% names(cls)) {
      cls[[z$rank]] <- z$name
      keys[[z$rank]] <- if (is.null(z$key)) NA_character_ else as.character(z$key)
    }
  }
  # The confidence is a number from 0 to 100, or absent (NA).
  list(match_type = g(r$diagnostics, "matchType"),
       confidence = suppressWarnings(as.integer(g(r$diagnostics, "confidence"))),
       key = g(r$usage, "key"), name = g(r$usage, "canonicalName"), rank = g(r$usage, "rank"),
       status = g(r$usage, "status"), accepted_key = g(r$acceptedUsage, "key"),
       accepted_name = g(r$acceptedUsage, "canonicalName"),
       accepted_rank = g(r$acceptedUsage, "rank"), classification = cls, keys = keys)
}

# Acceptance of an answer ----------------------------------------------------------

# The acceptable usage among the answer and its alternatives (in the order of
# the service), or NULL. A genus is preferred to a homonymous subgenus
# (Chironomus). Other acceptable taxa with the same name in another class or
# order are reported in `homonyms`.
accept_colxr <- function(parsed, target, form, group) {
  cand <- Filter(function(a) usage_acceptable(a, target, form, group),
                 c(list(parsed$main), parsed$alternatives))
  if (!length(cand)) return(NULL)
  rk <- vapply(cand, function(a) a$rank, character(1))
  tier <- vapply(cand, preference_tier, integer(1), group = group)
  cand <- cand[order(tier, rk == "SUBGENUS")]
  best <- cand[[1]]
  where <- function(a) paste(a$classification[c("KINGDOM", "CLASS", "ORDER")], collapse = "|")
  other <- Filter(function(a) where(a) != where(best) && !identical(a$rank, "SUBGENUS"), cand[-1])
  best$homonyms <- if (length(other)) {
    paste0("other taxa with this name in the COL XR: ", paste(vapply(other, function(a) {
      sprintf("%s (%s, %s)", a$name, tolower(a$rank),
              paste(stats::na.omit(a$classification[c("CLASS", "ORDER", "FAMILY")]), collapse = ", "))
    }, character(1)), collapse = "; "))
  } else {
    NA_character_
  }
  best
}

# Order of preference among acceptable homonyms, by the taxa expected in the
# broad group (lower first; ties keep the order of the service): in the
# zooplankton, rotifers, branchiopods, copepods and ostracods, then other
# arthropods (the cladoceran Latona, not the bivalve); in the
# macroinvertebrates, the phyla of freshwater invertebrates.
preference_tier <- function(a, group) {
  cls <- a$classification
  phylum <- cls[["PHYLUM"]]
  switch(group,
    zooplankton = if (phylum %in% "Rotifera" || cls[["CLASS"]] %in%
                        c("Branchiopoda", "Copepoda", "Ostracoda", "Hexanauplia", "Maxillopoda")) {
      1L
    } else if (phylum %in% "Arthropoda") 2L else 3L,
    macroinvertebrate = if (phylum %in% c("Annelida", "Arthropoda", "Mollusca", "Platyhelminthes",
                                          "Nematoda", "Nematomorpha", "Nemertea", "Cnidaria",
                                          "Porifera", "Bryozoa")) 1L else 2L,
    1L)
}

# Genus of a species not accepted, from the classification of an answer to the
# species that is in that genus (a misspelt epithet): the genus as placed in
# the COL XR rather than a homonymous record of it.
genus_from_answer <- function(parsed, genus, group) {
  for (a in c(list(parsed$main), parsed$alternatives)) {
    if (is.na(a$name) || !identical(unname(a$classification["GENUS"]), genus) ||
        is.na(a$keys[["GENUS"]])) next
    g <- a
    g$name <- genus
    g$rank <- "GENUS"
    g$key <- a$keys[["GENUS"]]
    g$status <- "ACCEPTED"
    g$accepted_key <- g$accepted_name <- g$accepted_rank <- NA_character_
    g$classification[["SPECIES"]] <- NA_character_
    g$match_type <- "EXACT"
    if (usage_acceptable(g, genus, "genus_of_species", group)) return(g)
  }
  NULL
}

usage_acceptable <- function(a, target, form, group) {
  if (is.na(a$name) || !a$match_type %in% c("EXACT", "VARIANT", "FUZZY")) return(FALSE)
  if (!same_taxon_name(target, a$name)) return(FALSE)
  cls <- a$classification
  if (form %in% c("species", "infraspecific")) {
    if (!identical(a$rank, "SPECIES")) return(FALSE)
  } else if (form == "genus_of_species") {
    if (!identical(a$rank, "GENUS")) return(FALSE)
  } else if (identical(a$rank, "GENUS") && is.na(cls[["FAMILY"]]) && is.na(cls[["ORDER"]])) {
    # A one-word name is not taken as a genus without family or order: the COL
    # XR has placeholder genera named after higher taxa (Oligochaeta). Genera
    # with an order and no family (the cladocerans Latona and Sida) are kept.
    return(FALSE)
  }
  k <- cls[["KINGDOM"]]
  ok_group <- switch(group,
    phytoplankton = !is.na(k) && k != "Animalia",
    fish = identical(k, "Animalia") && identical(cls[["PHYLUM"]], "Chordata") &&
      !identical(cls[["CLASS"]], "Amphibia"),
    amphibian = identical(cls[["CLASS"]], "Amphibia"),
    identical(k, "Animalia"))
  if (!ok_group) return(FALSE)
  if (a$status %in% c("AMBIGUOUS_SYNONYM", "MISAPPLIED")) {
    return(!is.na(a$accepted_name) &&
             same_taxon_name(target, sub("^(\\S+ \\S+).*$", "\\1", a$accepted_name)))
  }
  TRUE
}

# Same name, or a binomial whose epithet differs only in its ending
# (grammatical gender: -us, -a, -um, -is, -e, -i, -y).
same_taxon_name <- function(a, b) {
  if (is.na(a) || is.na(b)) return(FALSE)
  if (identical(tolower(a), tolower(b))) return(TRUE)
  wa <- strsplit(a, " ", fixed = TRUE)[[1]]
  wb <- strsplit(b, " ", fixed = TRUE)[[1]]
  if (length(wa) != 2L || length(wb) != 2L || wa[1] != wb[1]) return(FALSE)
  stem <- function(w) sub("(us|um|is|a|e|i|y)$", "", w)
  stem(wa[2]) == stem(wb[2])
}

closest_name <- function(parsed, target) {
  a <- if (is.null(parsed)) NULL else parsed$main
  if (is.null(a) || is.na(a$name) || identical(a$name, target) ||
      !a$match_type %in% c("EXACT", "VARIANT", "FUZZY")) {
    return("")
  }
  sprintf("; closest name: %s (%s, not used)", a$name, tolower(a$rank))
}
