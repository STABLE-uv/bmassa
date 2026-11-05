# Taxonomic hierarchy of the reference taxa (version 0.2.0)
#
# Builds, in data-raw/derived/:
#   ref_taxa.csv        one row per source and name used, with the main ranks
#                       (kingdom, phylum, class, order, family, genus,
#                       species) of its COL taxon, where each rank comes from
#                       and the domain of the model
#   ref_taxa_log.md     how many ranks were filled, and how
# from data-raw/derived/taxa_colxr.csv (run taxonomy_colxr.R first).
#
# Rule (decided on 2026-11-02): the ranks of the Catalogue of
# Life (COL XR), and the ranks it lacks filled with the classification of the
# sources. For each taxon and missing main rank, in this order:
#   1. col           the rank in the COL classification of the taxon;
#   2. col_related   the rank that the COL gives to other taxa of the same
#                    node (other species of the genus, other genera with the
#                    same tribe or superfamily...), most frequent value;
#   2b. col_above    (orders and classes only; added on 2026-11-03)
#                    the COL rank closest to the missing one among those
#                    above it and below the next main rank (a superorder,
#                    subterclass, infraclass or subclass in place of an
#                    order; a subphylum in place of a class). The sources
#                    name the order of each family on their own, often with
#                    old or inconsistent names (Basommatophora, Pulmonata and
#                    Hygrophila for the families of Hygrophila; Gymnoplea for
#                    Oithonidae, a family of Podoplea), so that sister
#                    families ended under different nodes; a COL taxon above
#                    them is the same for all of them;
#   3. source        the value given by the sources for the taxa of the node
#                    just below (the family of a genus, the order of a
#                    family), most frequent value, ties broken by a fixed
#                    order of sources. A value is not used when it is the
#                    name of another main rank of the lineage, or when the
#                    COL places that name in another class (for an order or
#                    family) or phylum (for a class): Scenedesmaceae is not
#                    the family of Crucigenia, which the COL puts in
#                    Trebouxiophyceae;
#   4. col_intermediate  the COL rank closest to the missing one between its
#                    two neighbours (a superorder in place of an order, a
#                    subphylum in place of a class...);
#   5. placeholder   "<higher taxon> incertae sedis".
# Ranks below the rank of the taxon are left empty (a record of a family has
# no genus).
#
# Run from the package root, in a UTF-8 locale:  Rscript data-raw/hierarchy.R

out_dir <- "data-raw/derived"
log_lines <- c("# Taxonomic hierarchy of the reference taxa", "",
               paste0("Built by `data-raw/hierarchy.R` on ", format(Sys.Date()), " (",
                      R.version.string, ")."), "")
say <- function(...) log_lines <<- c(log_lines, sprintf(...))

x <- read.csv(file.path(out_dir, "taxa_colxr.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
x <- x[x$use, ]
main <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
rank_order <- c("domain", "realm", "kingdom", "subkingdom", "infrakingdom", "phylum", "subphylum",
                "infraphylum", "parvphylum", "gigaclass", "megaclass", "superclass", "class",
                "subclass", "infraclass", "subterclass", "superorder", "order", "suborder",
                "infraorder", "parvorder", "section_zoology", "subsection_zoology",
                "superfamily", "epifamily", "family", "subfamily", "supertribe", "tribe",
                "subtribe", "genus", "subgenus", "section_botany", "species_aggregate",
                "species", "subspecies", "variety", "subvariety", "form", "infraspecific_name")
pos <- function(r) match(tolower(r), rank_order)

# Classification of each taxon as a list of (rank, name).
cls <- lapply(strsplit(x$classification, "|", fixed = TRUE), function(p) {
  data.frame(rank = tolower(sub(":.*$", "", p)), name = sub("^[A-Z_]+:", "", p),
             stringsAsFactors = FALSE)
})
# Main rank of the taxon itself (subspecies, varieties and forms count as
# species; a subgenus as its genus; a tribe as its family...).
own_main <- vapply(x$col_rank, function(r) {
  p <- pos(r)
  m <- main[pos(main) <= p]
  m[length(m)]
}, "")
# Copepods given by their order (Calanoida, Cyclopoida): the COL XR has no
# orders of copepods, and these names stand for the superorder Gymnoplea
# (which holds only the Calanoida) and the suborder Cyclopida
# (data-raw/taxa_manual.csv); their records are records of the order.
copepod_order_of <- c(Gymnoplea = "Calanoida", Cyclopida = "Cyclopoida")
co <- x$col_name %in% names(copepod_order_of)
own_main[co] <- "order"
deepest <- match(own_main, main)

lin <- as.data.frame(matrix(NA_character_, nrow(x), length(main), dimnames = list(NULL, main)),
                     stringsAsFactors = FALSE)
from <- lin
for (r in main) {
  v <- x[[paste0("col_", r)]]
  ok <- !is.na(v) & match(r, main) <= deepest
  lin[[r]][ok] <- v[ok]
  from[[r]][ok] <- "col"
}
lin$order[co] <- unname(copepod_order_of[x$col_name[co]])
from$order[co] <- "col_label"
needs <- function(r) is.na(lin[[r]]) & match(r, main) <= deepest

# Step 2: the rank given by the COL to other taxa that share a node below it.
# Nodes are told apart by their phylum (or kingdom), so that homonyms of
# different phyla do not share a rank: the cladoceran genus Latona, without a
# family in the COL XR, took that of the bivalve genus Latona (Donacidae).
node_scope <- function(c1) {
  ph <- c1$name[c1$rank == "phylum"]
  if (length(ph)) paste0("phylum:", ph[1]) else paste0("kingdom:", c(c1$name[c1$rank == "kingdom"], "?")[1])
}
fill_related <- function(r) {
  k <- which(needs(r))
  if (!length(k)) return(0L)
  have <- which(!is.na(lin[[r]]) & from[[r]] == "col")
  # node (scope, rank:name) -> values of r among the taxa that have r
  idx <- new.env()
  for (i in have) {
    nodes <- cls[[i]][pos(cls[[i]]$rank) > pos(r), ]
    for (n in paste(node_scope(cls[[i]]), nodes$rank, nodes$name)) {
      assign(n, c(get0(n, idx, ifnotfound = NULL), lin[[r]][i]), idx)
    }
  }
  n_filled <- 0L
  for (i in k) {
    nodes <- cls[[i]][pos(cls[[i]]$rank) > pos(r), ]
    vals <- unlist(lapply(paste(node_scope(cls[[i]]), nodes$rank, nodes$name), get0, envir = idx,
                          ifnotfound = NULL))
    if (length(vals)) {
      lin[[r]][i] <<- names(sort(table(vals), decreasing = TRUE))[1]
      from[[r]][i] <<- "col_related"
      n_filled <- n_filled + 1L
    }
  }
  n_filled
}

# Step 2b (orders and classes): the COL rank closest to the missing one among
# those between it and the main rank above (the closest above wins: a
# superorder before a subclass). Some names are replaced:
#   - the infraclass of the lower heterobranchs (Valvatoidea...) comes from
#     the API as "Lower" and is written in full;
#   - the COL XR places the orders of free-living flatworms (Tricladida,
#     Rhabdocoela, Prolecithophora) in the subphylum Rhabditophora without a
#     class, while it keeps the class Turbellaria for the records identified
#     only to it: those orders take the class Turbellaria (as in the
#     sources), and a taxon given as Turbellaria shares their node;
#   - the COL XR has no orders of copepods, only the superorders Gymnoplea
#     and Podoplea, while copepods are counted by order (Calanoida,
#     Cyclopoida, Harpacticoida). Gymnoplea has only Calanoida; the families
#     of Podoplea take their order in WoRMS (`copepod_orders`; the former
#     Poecilostomatoida are the suborder Ergasilida of Cyclopoida). A family
#     of Podoplea missing from the list stops the script.
col_labels <- list(order = c(Lower = "Lower Heterobranchia", Gymnoplea = "Calanoida"),
                   class = c(Rhabditophora = "Turbellaria"))
copepod_orders <- c(Bomolochidae = "Cyclopoida", Corycaeidae = "Cyclopoida", Cyclopidae = "Cyclopoida",
                    Oithonidae = "Cyclopoida", Oncaeidae = "Cyclopoida", Sapphirinidae = "Cyclopoida",
                    Canthocamptidae = "Harpacticoida", Ectinosomatidae = "Harpacticoida",
                    Harpacticidae = "Harpacticoida", Peltidiidae = "Harpacticoida",
                    Tisbidae = "Harpacticoida")
fill_above <- function(r) {
  k <- which(needs(r))
  if (!length(k)) return(0L)
  lo <- pos(main[match(r, main) - 1])
  hi <- pos(r)
  n_filled <- 0L
  for (i in k) {
    c1 <- cls[[i]]
    p <- pos(c1$rank)
    cand <- c1[!is.na(p) & p > lo & p < hi, ]
    if (!nrow(cand)) next
    nm <- cand$name[which.max(pos(cand$rank))]
    lab <- col_labels[[r]]
    if (r == "order" && nm == "Podoplea") {
      fam <- lin$family[i]
      if (!fam %in% names(copepod_orders)) stop("Family of Podoplea without an order in copepod_orders: ", fam)
      nm <- copepod_orders[[fam]]
    }
    lin[[r]][i] <<- if (nm %in% names(lab)) unname(lab[nm]) else nm
    from[[r]][i] <<- "col_above"
    n_filled <- n_filled + 1L
  }
  n_filled
}

# Step 3: the value of the sources for the node just below, if the COL does
# not place that name elsewhere.
source_pref <- c("dumont1975", "vanleeuwen2025", "moreau2021", "cadses_boxcorer",
                 "cadses_leafpack", "kremer2014", "rimet2018", "santini2018", "trochet2014",
                 "arranz2022", "coughlan2021", "hrycik2024", "brey2010")
clean_hint <- function(v) {
  v <- sub("^\\[unassigned\\] ", "", v)
  v[!is.na(v) & !grepl("^[A-Z][a-z]+$", v)] <- NA
  v
}
# Higher rank used to check where the COL places a name.
check_rank <- c(class = "phylum", order = "class", family = "class")
placed_in <- function(name, rank) {
  # values of `rank` of the taxa whose COL classification contains `name`
  unique(unlist(lapply(seq_along(cls), function(i) {
    if (name %in% cls[[i]]$name) lin[[rank]][i] else NULL
  })))
}
fill_source <- function(r) {
  k <- which(needs(r))
  if (!length(k)) return(c(filled = 0L, rejected = 0L))
  below <- main[(match(r, main) + 1):length(main)]
  hint <- clean_hint(x[[paste0("hint_", r)]])
  # a hint equal to another main rank of the lineage is not a value of r
  for (m in setdiff(main, r)) hint[!is.na(hint) & !is.na(lin[[m]]) & hint == lin[[m]]] <- NA
  up <- check_rank[[r]]
  # node just below r: the first main rank below r that the taxon has
  key <- vapply(seq_len(nrow(x)), function(i) {
    b <- below[!is.na(unlist(lin[i, below]))]
    if (!length(b)) NA_character_ else paste(lin[[up]][i], b[1], lin[[b[1]]][i])
  }, "")
  cache <- list()
  filled <- rejected <- 0L
  for (i in k) {
    if (is.na(key[i])) next
    same <- which(key == key[i] & !is.na(hint))
    if (!length(same)) next
    tab <- table(hint[same])
    best <- names(tab)[tab == max(tab)]
    if (length(best) > 1) {
      first_src <- vapply(best, function(b) min(match(x$source_id[same][hint[same] == b], source_pref)), 0)
      best <- best[which.min(first_src)]
    }
    where <- cache[[best]]
    if (is.null(where)) where <- cache[[best]] <- placed_in(best, up)
    where <- where[!is.na(where)]
    if (length(where) && !is.na(lin[[up]][i]) && !(lin[[up]][i] %in% where)) {
      rejected <- rejected + 1L
      next
    }
    lin[[r]][i] <<- best
    from[[r]][i] <<- "source"
    filled <- filled + 1L
  }
  c(filled = filled, rejected = rejected)
}

# Step 4: the COL rank closest to the missing one, between its neighbours.
fill_intermediate <- function(r) {
  k <- which(needs(r))
  if (!length(k)) return(0L)
  j <- match(r, main)
  lo <- pos(main[j - 1])
  hi <- pos(main[j + 1])
  n_filled <- 0L
  for (i in k) {
    c1 <- cls[[i]]
    p <- pos(c1$rank)
    cand <- c1[!is.na(p) & p > lo & p < hi, ]
    if (!nrow(cand)) next
    d <- abs(pos(cand$rank) - pos(r)) + 0.1 * (pos(cand$rank) > pos(r))
    lin[[r]][i] <<- cand$name[which.min(d)]
    from[[r]][i] <<- "col_intermediate"
    n_filled <- n_filled + 1L
  }
  n_filled
}

# Families first, so that the order of a family is filled once for all its
# genera; placeholders last, from the top, once the parents are known.
report <- list()
for (r in c("family", "order", "class")) {
  missing0 <- sum(needs(r))
  a <- fill_related(r)
  ab <- if (r %in% c("order", "class")) fill_above(r) else 0L
  b <- fill_source(r)
  cc <- fill_intermediate(r)
  report[[r]] <- c(missing = missing0, col_related = a, col_above = ab, source = b[["filled"]],
                   source_rejected = b[["rejected"]], col_intermediate = cc)
}
for (r in c("class", "order", "family")) {
  k <- which(needs(r))
  for (i in k) {
    parent <- lin[[main[match(r, main) - 1]]][i]
    lin[[r]][i] <- paste(parent, "incertae sedis")
    from[[r]][i] <- "placeholder"
  }
  report[[r]] <- c(report[[r]], placeholder = length(k))
}

# One parent per node. Where the COL itself, or the filling, gives a genus
# two families (or a family two orders...), every taxon of the node takes the
# parent with the best origin (col, col_related, col_above, source,
# col_intermediate, placeholder) and, among those, the most frequent; the rows
# changed are marked "harmonized". From genus up, so that each change carries
# upwards. Genera, families and orders are nodes within their phylum:
# homonyms of different phyla stay apart (the cladoceran genus Latona and the
# bivalve genus Latona); classes within their kingdom.
origin_rank <- c(col = 1, col_related = 2, col_above = 2.5, source = 3, col_intermediate = 4,
                 placeholder = 5, harmonized = 6)
node_key <- function(d, r) {
  scope <- if (r == "class") d$kingdom else paste(d$kingdom, d$phylum)
  ifelse(is.na(d[[r]]), NA, paste(scope, d[[r]]))
}
n_harm <- c()
for (r in c("genus", "family", "order", "class")) {
  up <- main[match(r, main) - 1]
  node <- node_key(lin, r)
  changed <- 0L
  for (nd in unique(node[!is.na(node)])) {
    k <- which(node == nd)
    vals <- lin[[up]][k]
    if (length(unique(vals)) < 2) next
    o <- origin_rank[from[[up]][k]]
    best_o <- min(o)
    cand <- vals[o == best_o]
    pick <- names(sort(table(cand), decreasing = TRUE))[1]
    ch <- k[vals != pick]
    # ranks above the parent follow the row that already has the chosen parent
    model_row <- k[vals == pick][1]
    above <- main[seq_len(match(up, main))]
    for (m in above) {
      lin[[m]][ch] <- lin[[m]][model_row]
      from[[m]][ch] <- ifelse(lin[[m]][ch] == lin[[m]][model_row] & from[[m]][ch] == from[[m]][model_row],
                              from[[m]][ch], "harmonized")
    }
    from[[up]][ch] <- "harmonized"
    changed <- changed + length(ch)
  }
  n_harm[r] <- changed
}
# Genus and species: from the COL only (records above them stay above).
stopifnot(!anyNA(lin$kingdom), !anyNA(lin$phylum))

# Domain of the model.
vertebrate <- grepl("SUBPHYLUM:Vertebrata", x$classification, fixed = TRUE)
domain <- ifelse(lin$class %in% "Amphibia", "amphibian",
          ifelse(vertebrate & !(lin$class %in% c("Aves", "Mammalia", "Reptilia", "Squamata",
                                                   "Testudines", "Crocodylia")), "fish",
          ifelse(lin$kingdom == "Animalia", "invertebrate",
          ifelse(x$source_id %in% c("kremer2014", "rimet2018"), "phytoplankton", "other"))))

out <- data.frame(source_id = x$source_id, taxon = x$taxon, col_id = x$col_id,
                  col_name = x$col_name, col_rank = x$col_rank, lowest_rank = own_main,
                  domain = domain, lin, stringsAsFactors = FALSE)
names(from) <- paste0("from_", main)
out <- cbind(out, from)
out$col_release <- x$col_release
out <- out[order(out$source_id, out$taxon, method = "radix"), ]
write.csv(out, file.path(out_dir, "ref_taxa.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Log --------------------------------------------------------------------------------

say("`ref_taxa.csv`: %d names of %d sources; %d COL taxa.", nrow(out),
    length(unique(out$source_id)), length(unique(out$col_id)))
say("")
say("Domains: %s.", paste(sprintf("%s %d", names(table(out$domain)), as.integer(table(out$domain))),
                          collapse = ", "))
say("")
say("## Ranks filled")
say("")
say("| Rank | Missing in the COL | From related COL taxa | Nearest COL rank above | From the sources | Source value rejected | Nearest COL rank | Placeholder |")
say("|---|---|---|---|---|---|---|---|")
for (r in names(report)) {
  v <- report[[r]]
  say("| %s | %d | %d | %d | %d | %d | %d | %d |", r, v[["missing"]], v[["col_related"]], v[["col_above"]],
      v[["source"]], v[["source_rejected"]], v[["col_intermediate"]], v[["placeholder"]])
}
say("")
say("- Copepods given by their order (the COL XR taxa Gymnoplea and Cyclopida taken as the orders Calanoida and Cyclopoida): %d rows.",
    sum(co))
say("- Rows given another parent so that each node has one: %s.",
    paste(sprintf("%s of a %s %d", main[match(names(n_harm), main) - 1], names(n_harm), n_harm), collapse = ", "))
say("")
for (r in c("class", "order", "family")) {
  s <- out[out[[paste0("from_", r)]] %in% c("col_above", "source", "col_intermediate"), ]
  if (!nrow(s)) next
  tab <- sort(table(paste0(s[[r]], " (", s[[paste0("from_", r)]], ")")), decreasing = TRUE)
  say("- %s filled: %s.", r, paste(sprintf("%s %d", names(head(tab, 25)), as.integer(head(tab, 25))),
                                    collapse = ", "))
}
say("")
# Check: one parent per node (a family under one order, an order under one class).
for (r in c("class", "order", "family", "genus")) {
  up <- main[match(r, main) - 1]
  k <- !is.na(out[[r]])
  par <- tapply(out[[up]][k], node_key(out, r)[k], function(v) length(unique(v)))
  say("- %s with more than one %s: %d%s.", r, up, sum(par > 1),
      if (any(par > 1)) paste0(" (", paste(head(names(par)[par > 1], 15), collapse = ", "), ")") else "")
}
writeLines(log_lines, file.path(out_dir, "ref_taxa_log.md"), useBytes = TRUE)
cat(log_lines, sep = "\n")
