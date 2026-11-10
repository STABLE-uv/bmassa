# Reader of the tadpoles of Knapp et al. (2021), used by data-raw/fit_lwr.R and
# data-raw/typical_records.R. The tadpoles come from the compilation of the
# elemental composition of freshwater animals of 2026-11-09
# (Stoichiometry2026/ under BMASSA_DATABASES), which transcribes each row of
# the data set of Knapp et al. (Knapp_et_al_Stoichiometry_Dataset.xlsx, Sheet1)
# in `verbatim`, in the order of its cells: date, wetland, species, Gosner
# stage, simplified stage (1 to 3), total length (mm), a mass whose column
# name is not in the compilation (not used), dry mass (mg), %C, %N, %P and the
# molar ratios, and for pooled samples a comment.
#
# Decisions (2026-11-09; dev/decisions.md):
# - One row per tadpole measured alone (a row of the sheet; its C, N and P
#   are rows of their own in the compilation). The pooled homogenates are left
#   out: the compilation gives the size of their first member only.
# - Total length and dry mass are those of the body and tail after the gut
#   was dissected out (the compilation says so for every row), so the dry
#   mass is that of the animal without its gut contents.
# - The mass between total length and dry mass is not used until its column
#   is known: its ratio to the dry mass (median 5.1, 2.6-14.9) fits neither a
#   whole wet mass nor that of the dissected body (tadpoles hold about 90 %
#   water).
# - Stage 358 (one tadpole) is a typing error and is taken as unknown.

# Tadpoles of the compilation in `path` (its CSV file). Columns: row (the row
# of the sheet), taxon (accepted name in the COL XR, as in the compilation),
# name_as_given, date (Date, NA when missing), wetland, gosner (integer, NA
# when unknown), stage_simplified, total_length_mm, dry_mass_mg.
read_knapp_tadpoles <- function(path) {
  ba <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, encoding = "UTF-8",
                 na.strings = "")
  k <- ba[ba$doi %in% "10.6084/m9.figshare.13683910.v1", ]
  stopifnot(nrow(k) > 0, all(k$group == "Amphibians"), all(startsWith(k$life_stage, "tadpole")))
  k <- k[!grepl("pooled", k$life_stage, fixed = TRUE), ]
  stopifnot(all(grepl("^Knapp_et_al_Stoichiometry_Dataset\\.xlsx, Sheet1 row [0-9]+$",
                      k$table_or_figure)))
  k <- k[!duplicated(k$table_or_figure), ]
  v <- strsplit(k$verbatim, ",", fixed = TRUE)
  stopifnot(all(lengths(v) == 14L))
  v <- do.call(rbind, v)
  num <- function(x) suppressWarnings(as.numeric(x))
  date <- as.Date(substr(v[, 1], 1, 10), format = "%Y-%m-%d")
  gosner <- suppressWarnings(as.integer(v[, 4]))
  gosner[!gosner %in% 20:46] <- NA_integer_
  out <- data.frame(row = as.integer(sub("^.* row ", "", k$table_or_figure)),
                    taxon = k$accepted_name, name_as_given = k$name_as_given, date = date,
                    wetland = v[, 2], gosner = gosner, stage_simplified = num(v[, 5]),
                    total_length_mm = num(v[, 6]), dry_mass_mg = num(v[, 8]),
                    stringsAsFactors = FALSE)
  # The compilation repeats the length and the dry mass in `notes`.
  tl <- num(sub("^.*Total length ([0-9.]+) mm.*$", "\\1", k$notes))
  dm <- num(sub("^.*dry mass ([0-9.]+) mg.*$", "\\1", k$notes))
  stopifnot(isTRUE(all.equal(tl, out$total_length_mm)), isTRUE(all.equal(dm, out$dry_mass_mg)),
            !anyNA(out$taxon), !anyNA(out$wetland), all(out$total_length_mm > 0),
            all(out$dry_mass_mg > 0))
  out <- out[order(out$row), ]
  rownames(out) <- NULL
  out
}
