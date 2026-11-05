# Helpers for the source registry (data-raw/sources.csv), used by the build
# scripts of data-raw/. Run them from the package root.
#
# Original files are read from the folder given by the environment variable
# BMASSA_DATABASES (default "data-raw/databases"), with one subfolder per
# source as in the `files` column of the registry, e.g.
#   Sys.setenv(BMASSA_DATABASES = "path/to/databases")
#
# Permission states (column `permission`):
#   not_needed     the licence allows redistribution under CC BY 4.0 or MIT
#   pending        permission requested or to be requested from the authors
#   unverified     licence not yet checked
#   granted        permission received (note who and when in `note`)
#   refused        permission refused: the source must be removed
#   not_requested  source not used as data (status "excluded")
#
# The package and its public repository hold only the parameters of the
# models (by taxon, with the sources to cite), never the records of the
# sources; the records, the transcriptions and this registry stay in the
# development repository (decision of 2026-11-05; check_permissions.R).

databases_dir <- function() Sys.getenv("BMASSA_DATABASES", "data-raw/databases")

read_sources <- function(path = "data-raw/sources.csv") {
  s <- read.csv(path, stringsAsFactors = FALSE, encoding = "UTF-8", na.strings = "")
  stopifnot(!anyDuplicated(s$source_id),
            all(s$status %in% c("included", "planned", "excluded")),
            all(s$permission %in% c("not_needed", "pending", "unverified", "granted",
                                    "refused", "not_requested")))
  s
}

# TRUE when the records of a source can be distributed openly as things stand.
redistributable <- function(permission) permission %in% c("not_needed", "granted")

# Whether a source enters the tables built now. Sources without permission for
# open redistribution are included (BMASSA_INCLUDE_PENDING unset or "true":
# only the fitted parameters are distributed) and left out with
# BMASSA_INCLUDE_PENDING = "false". Refused sources never enter.
use_source <- function(source_id, sources = read_sources()) {
  include_pending <- tolower(Sys.getenv("BMASSA_INCLUDE_PENDING", "true")) == "true"
  p <- sources$permission[match(source_id, sources$source_id)]
  if (anyNA(p)) stop("Unknown source: ", paste(source_id[is.na(p)], collapse = ", "))
  ok <- redistributable(p) | (include_pending & p %in% c("pending", "unverified"))
  ok & p != "refused"
}

# Path of the files of a source, after checking their MD5 checksums.
source_files <- function(source_id, sources = read_sources()) {
  i <- match(source_id, sources$source_id)
  if (is.na(i)) stop("Unknown source: ", source_id)
  files <- strsplit(sources$files[i], ";", fixed = TRUE)[[1]]
  md5 <- strsplit(sources$md5[i], ";", fixed = TRUE)[[1]]
  path <- file.path(databases_dir(), files)
  found <- file.exists(path)
  if (!all(found)) {
    stop("Files of '", source_id, "' not found under ", databases_dir(), ":\n",
         paste(files[!found], collapse = "\n"),
         "\nSet BMASSA_DATABASES to the folder that holds them.")
  }
  got <- unname(tools::md5sum(path))
  if (!identical(got, md5)) {
    stop("Files of '", source_id, "' differ from the registered versions:\n",
         paste(files[got != md5], collapse = "\n"))
  }
  setNames(path, files)
}
