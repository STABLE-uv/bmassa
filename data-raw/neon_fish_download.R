# Download of the fish of NEON (DP1.20107.001, Fish electrofishing, gill
# netting, and fyke netting counts) straight from the NEON API, without
# neonUtilities::loadByProduct() (version 4.0.2 gave "No response" with a
# valid token in 2026-11).
#
# Needs an API token in the environment variable NEON_TOKEN (.Renviron) and the
# packages httr and jsonlite. Saves, for every site and month, the basic tables
# fsh_perFish (fish measured one by one), fsh_perPass (passes and nets) and
# fsh_fieldData (reaches and bouts), with the release of each file, in
# out_file (the folder NEON2026 of BMASSA_DATABASES, data-raw/sources.R). It
# prints its progress and the requests that failed. The copy in the registry
# (data-raw/sources.csv) was downloaded on 2026-11-06 and holds RELEASE-2026
# and the provisional data of that day; data-raw/neon.R keeps only the release
# named in its settings. A new download differs in its provisional data (and
# MD5), and after a new release it holds that release instead.
#
# Run from the package root:  Rscript data-raw/neon_fish_download.R

out_file <- file.path(Sys.getenv("BMASSA_DATABASES", "data-raw/databases"), "NEON2026",
                      "neon_fish_DP1.20107.001.rds")
tables <- c("fsh_perFish", "fsh_perPass", "fsh_fieldData")

token <- Sys.getenv("NEON_TOKEN")
if (!nzchar(token)) stop("NEON_TOKEN is empty: add it to .Renviron and restart R.")
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
httr::set_config(httr::config(http_version = 2))  # HTTP/1.1

# One request with up to four attempts; NULL if all fail.
fetch <- function(url, api = TRUE, tries = 4) {
  for (k in seq_len(tries)) {
    r <- try(if (api) {
      httr::GET(url, httr::add_headers("X-API-Token" = token, accept = "application/json"),
                httr::timeout(120))
    } else {
      httr::GET(url, httr::timeout(300))
    }, silent = TRUE)
    if (inherits(r, "response") && httr::status_code(r) == 200) {
      return(httr::content(r, as = "text", encoding = "UTF-8"))
    }
    why <- if (inherits(r, "response")) paste("HTTP", httr::status_code(r)) else
      conditionMessage(attr(r, "condition"))
    message("  attempt ", k, " failed (", why, "): ", substr(url, 1, 110))
    Sys.sleep(5 * k)
  }
  NULL
}

prod_txt <- fetch("https://data.neonscience.org/api/v0/products/DP1.20107.001")
if (is.null(prod_txt)) stop("The NEON API did not answer the product request.")
prod <- jsonlite::fromJSON(prod_txt, simplifyVector = FALSE)
month_urls <- unlist(lapply(prod$data$siteCodes, function(s) unlist(s$availableDataUrls)))
message(length(month_urls), " site-months to download; this takes a while.")

res <- setNames(lapply(tables, function(t) list()), tables)
failed <- character()
for (i in seq_along(month_urls)) {
  txt <- fetch(paste0(month_urls[i], "?package=basic"))
  if (is.null(txt)) {
    failed <- c(failed, month_urls[i])
    next
  }
  m <- jsonlite::fromJSON(txt, simplifyVector = FALSE)$data
  release <- if (is.null(m$release)) NA_character_ else m$release
  for (t in tables) {
    files <- Filter(function(x) grepl(paste0(".", t, "."), x$name, fixed = TRUE) &&
                      grepl("\\.csv$", x$name), m$files)
    for (x in files) {
      csv <- fetch(x$url, api = FALSE)
      if (is.null(csv)) {
        failed <- c(failed, x$name)
        next
      }
      d <- utils::read.csv(text = csv, stringsAsFactors = FALSE, colClasses = "character")
      d$release <- rep(release, nrow(d))
      d$file_name <- rep(x$name, nrow(d))
      res[[t]][[length(res[[t]]) + 1L]] <- d
    }
  }
  if (i %% 20 == 0) message(i, " of ", length(month_urls), " site-months")
}

# Tables joined; columns missing in some months filled with NA.
bind <- function(l) {
  if (!length(l)) return(NULL)
  cols <- unique(unlist(lapply(l, names)))
  do.call(rbind, lapply(l, function(d) {
    for (cn in setdiff(cols, names(d))) d[[cn]] <- rep(NA_character_, nrow(d))
    d[cols]
  }))
}
fish <- lapply(res, bind)
fish$citation <- list(
  product = "DP1.20107.001",
  downloaded = format(Sys.time(), tz = "UTC", usetz = TRUE),
  releases = tryCatch(vapply(prod$data$releases, function(r) {
    paste(r$release, if (!is.null(r$productDoi$url)) r$productDoi$url else "")
  }, character(1)), error = function(e) character())
)
fish$failed <- failed
saveRDS(fish, out_file)
message("Saved ", out_file, ": ", if (is.null(fish$fsh_perFish)) 0 else nrow(fish$fsh_perFish),
        " fish measured; ", length(failed), " requests failed.")
