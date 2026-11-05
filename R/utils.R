# Shared helpers and vocabularies (internal) ---------------------------------

# Ranks of the taxonomy, from the lowest to the highest.
rank_levels <- c(
  "species", "genus", "subfamily", "family", "suborder", "order",
  "subclass", "class", "subphylum", "phylum", "kingdom"
)

count_unit_levels <- c("individual", "cell", "natural_unit")

# Append a note to existing notes (vectorised over `old`).
append_note <- function(old, new) {
  if (length(old) == 0L) {
    return(character())
  }
  ifelse(is.na(new) | !nzchar(new), old,
         ifelse(is.na(old) | !nzchar(old), new, paste(old, new, sep = "; ")))
}

# Stop when the new columns would overwrite columns of the input.
check_no_clash <- function(data, new_cols, arg = "data", call = sys.call(-1)) {
  clash <- intersect(new_cols, names(data))
  if (length(clash) > 0L) {
    abort_input(
      sprintf(
        paste(
          "`%s` already has %s %s, which this function adds. Rename or drop",
          "%s first (for example, results of an earlier call)."
        ),
        arg, if (length(clash) == 1L) "column" else "columns",
        format_items(quote_items(clash), max = Inf),
        if (length(clash) == 1L) "it" else "them"
      ),
      call = call
    )
  }
  invisible(data)
}

# Character vector of text values from a column (factors as text).
as_text <- function(x) {
  if (is.factor(x)) as.character(x) else as.character(x)
}

# A column of allowed text values; NA allowed or not.
check_values_in <- function(data, column, allowed, allow_na = TRUE, arg = "data",
                            call = sys.call(-1)) {
  x <- as_text(data[[column]])
  bad <- which(!is.na(x) & !x %in% allowed)
  if (length(bad) > 0L) {
    abort_input(
      sprintf(
        "Column '%s' of `%s` has values that are not allowed (%s in %s). Allowed values: %s.",
        column, arg, format_items(quote_items(unique(x[bad])), max = 5L),
        format_rows(bad), format_items(quote_items(allowed), max = Inf)
      ),
      column = column, rows = bad, call = call
    )
  }
  if (!allow_na && anyNA(x)) {
    rows <- which(is.na(x))
    abort_input(
      sprintf("Column '%s' of `%s` has missing values in %s.", column, arg, format_rows(rows)),
      column = column, rows = rows, call = call
    )
  }
  invisible(data)
}
