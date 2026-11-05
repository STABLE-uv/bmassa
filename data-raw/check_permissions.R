# Report of the permissions of the sources and of the MD5 checksums of the
# original files.
#
# Run from the package root:
#   Rscript data-raw/check_permissions.R            # report
#   Rscript data-raw/check_permissions.R --release  # error if a source was refused
#
# Criterion (decision of 2026-11-05): the records of the sources are used to
# fit the models and the sources are cited (cite_sources()), whatever their
# permission, unless it was refused. The package and its public repository
# hold only the parameters of the models (by taxon, with the sources to
# cite); the records, the transcriptions and the source registry stay in the
# development repository. Git keeps the history: data removed later stay in
# earlier commits.

source("data-raw/sources.R")

args <- commandArgs(trailingOnly = TRUE)
release <- "--release" %in% args

s <- read_sources()
used <- s[s$status %in% c("included", "planned"), ]
closed <- used[!redistributable(used$permission), ]

cat("Sources in use or planned:", nrow(used), "\n")
if (nrow(closed)) {
  cat("Without permission for open redistribution of their records (only their",
      "fitted parameters are distributed):\n")
  for (i in seq_len(nrow(closed))) {
    cat(sprintf("  - %s [%s, %s]: %s\n", closed$source_id[i], closed$status[i],
                closed$permission[i], closed$licence[i]))
  }
} else {
  cat("All can be redistributed.\n")
}

with_files <- used[!is.na(used$files), ]
for (id in with_files$source_id) {
  ok <- tryCatch({ source_files(id, s); TRUE }, error = function(e) {
    cat("  ! ", conditionMessage(e), "\n", sep = "")
    FALSE
  })
  if (ok) cat("  MD5 ok:", id, "\n")
}

refused <- used$source_id[used$status == "included" & used$permission == "refused"]
if (length(refused)) {
  msg <- paste("included sources with permission refused:", paste(refused, collapse = ", "))
  if (release) stop("Release blocked: ", msg)
  cat("Problem:", msg, "\n")
}
