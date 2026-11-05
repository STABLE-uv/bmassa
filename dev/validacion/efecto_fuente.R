# Efecto de fuente en el modelo de invertebrados sin longitud (fuera del check)
#
# Compara cuatro formas de tratar las fuentes en invertebrate_individual:
#  - A: efecto aleatorio de fuente cruzado con toda la jerarquía (versión 0.2.0);
#  - C: fuente dentro de grupo de taxones (rotíferos, microcrustáceos y el
#    resto de invertebrados);
#  - N: sin efecto de fuente;
#  - M: sin efecto de fuente y con el contraste fijo de los registros de bentos
#    de Moreau et al. (2021) (`moreau_benthos`, model_data.R): versión 0.3.0,
#    decisiones del 2026-11-05.
# A, C y N se ajustan sin ese contraste.
# Para cada una: la mediana, por fuente, de la razón entre la predicción y la
# media de los registros de las especies medidas por una sola fuente, y una
# validación cruzada rápida (10 grupos de registros y cada fuente fuera por
# turno) con la cobertura del intervalo del 95 % y el factor medio
# datos/predicción por grupo de taxones.
# Escribe dev/validacion/efecto_fuente.md.
# Desde la raíz del paquete: Rscript dev/validacion/efecto_fuente.R  (unos 2 min)

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
dat <- read_model_data()
m <- "invertebrate_individual"
d0 <- model_rows(dat, m, "typical")
grupo <- ifelse(d0$phylum %in% "Rotifera", "rotíferos",
                ifelse(d0$class %in% c("Branchiopoda", "Copepoda", "Ostracoda"),
                       "microcrustáceos", "otros invertebrados"))
variantes <- c(A = "cruzada", C = "dentro de grupo", N = "sin fuente",
               M = "sin fuente, con contraste de Moreau en el bentos")
sin_moreau <- function(d, v) {
  if (v != "M") d$contrast[d$contrast %in% "moreau_benthos"] <- NA_character_
  d
}

ajusta <- function(d, g, v) {
  d <- sin_moreau(d, v)
  levels <- split_levels(d)
  opt <- model_options(d, "typical")
  opt$source <- v %in% c("A", "C")
  f <- model_frame(d, "typical", levels)
  if (v == "C") f$source_id <- paste(f$source_id, g, sep = "|")
  fit <- suppressMessages(suppressWarnings(fit_model(f, "typical", levels, opt)))
  list(tabs = model_tables(fit, f, m, "typical", levels), levels = levels, fit = fit)
}

set.seed(20261105)
u <- unique(d0$record_id)
fold <- stats::setNames(sample(rep(1:10, length.out = length(u))), u)
res <- list()
cv <- list()
param <- list()
for (v in names(variantes)) {
  a <- ajusta(d0, grupo, v)
  pr <- predict_tables(a$tabs, sin_moreau(d0, v), "typical", a$levels)
  res[[v]] <- data.frame(variante = v, species = d0$species, source = d0$source_id, y = d0$y,
                         eta = pr$eta, stringsAsFactors = FALSE)
  vc <- as.data.frame(lme4::VarCorr(a$fit))
  s <- vc$sdcor[vc$grp == "source_id"]
  fx <- lme4::fixef(a$fit)
  se <- sqrt(diag(as.matrix(stats::vcov(a$fit))))
  param[[v]] <- data.frame(variante = v, de_fuente = if (length(s)) round(s, 2) else NA,
                           de_residual = round(vc$sdcor[vc$grp == "Residual"], 2),
                           pesadas = round(exp(fx[["c_rotifer_weighing"]]), 1),
                           moreau_bentos = if ("c_moreau_benthos" %in% names(fx))
                             sprintf("%.1f (%.1f-%.1f)", exp(fx[["c_moreau_benthos"]]),
                                     exp(fx[["c_moreau_benthos"]] - 1.96 * se[["c_moreau_benthos"]]),
                                     exp(fx[["c_moreau_benthos"]] + 1.96 * se[["c_moreau_benthos"]]))
                           else "—")
  jobs <- c(lapply(1:10, function(f) list(esq = "unit", test = fold[d0$record_id] == f)),
            lapply(sort(unique(d0$source_id)), function(s) list(esq = "source", test = d0$source_id == s)))
  cv[[v]] <- do.call(rbind, lapply(jobs, function(j) {
    b <- tryCatch(ajusta(d0[!j$test, ], grupo[!j$test], v), error = function(e) NULL)
    if (is.null(b)) return(NULL)
    p2 <- tryCatch(predict_tables(b$tabs, sin_moreau(d0[j$test, ], v), "typical", b$levels),
                   error = function(e) NULL)
    if (is.null(p2)) return(NULL)
    data.frame(variante = v, esq = j$esq, grupo = grupo[j$test], y = d0$y[j$test],
               eta = p2$eta, vt = p2$v_total, stringsAsFactors = FALSE)
  }))
}
res <- do.call(rbind, res)
cv <- do.call(rbind, cv)
param <- do.call(rbind, param)

# Especies de una sola fuente (sin las no identificadas a especie).
res <- res[!grepl("not identified", res$species, fixed = TRUE), ]
nsrc <- tapply(res$source, res$species, function(s) length(unique(s)))
r1 <- res[nsrc[res$species] == 1, ]
sm <- stats::aggregate(cbind(y, eta) ~ variante + species + source, r1, mean)
sm$razon <- exp(sm$eta - sm$y)
t1 <- stats::aggregate(razon ~ source + variante, sm, stats::median)
t1 <- stats::reshape(t1, idvar = "source", timevar = "variante", direction = "wide")
n_sp <- table(sm$source[sm$variante == "A"])
t1 <- data.frame(fuente = t1$source, especies = as.integer(n_sp[t1$source]),
                 A = sprintf("%.2f", t1$razon.A), C = sprintf("%.2f", t1$razon.C),
                 N = sprintf("%.2f", t1$razon.N), M = sprintf("%.2f", t1$razon.M))

cv$z <- (cv$y - cv$eta) / sqrt(cv$vt)
cv <- cv[is.finite(cv$z), ]
t2 <- do.call(rbind, lapply(split(cv, list(cv$esq, cv$grupo, cv$variante), drop = TRUE), function(x) {
  data.frame(esquema = x$esq[1], grupo = x$grupo[1], variante = x$variante[1], n = nrow(x),
             cob95 = sprintf("%.0f %%", 100 * mean(abs(x$z) < 1.96)),
             datos_prediccion = sprintf("%.2f", exp(mean(x$y - x$eta))))
}))
t2 <- t2[order(t2$esquema, t2$grupo, t2$variante), ]
t3 <- do.call(rbind, lapply(split(cv, list(cv$esq, cv$variante), drop = TRUE), function(x) {
  data.frame(esquema = x$esq[1], variante = x$variante[1], n = nrow(x),
             cob95 = sprintf("%.0f %%", 100 * mean(abs(x$z) < 1.96)),
             factor95 = sprintf("%.1f", stats::median(exp(1.96 * sqrt(x$vt)))))
}))
t3 <- t3[order(t3$esquema, t3$variante), ]

md_table <- function(x) {
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep("---|", ncol(x)), collapse = "")),
    apply(x, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |")))
}
lines <- c("# Efecto de fuente en el modelo de invertebrados sin longitud", "",
           paste0("Generado por `dev/validacion/efecto_fuente.R` el ", format(Sys.Date()), " (",
                  R.version.string, ", lme4 ", as.character(utils::packageVersion("lme4")),
                  "; semilla 20261105)."), "",
           "Variantes: `A`, efecto de fuente cruzado (versión 0.2.0); `C`, fuente dentro de grupo de taxones (rotíferos, microcrustáceos y el resto); `N`, sin efecto de fuente; `M`, sin efecto de fuente y con un contraste fijo para los registros de bentos de Moreau et al. (2021) (versión 0.3.0). A, C y N se ajustan sin ese contraste. Las predicciones de las tablas llevan los contrastes de cada registro (pesadas de rotíferos, estadios de copépodos y, en M, el de Moreau).", "",
           "## Parámetros", "",
           "`pesadas`: factor del contraste de las pesadas de rotíferos frente a las masas por volumen. `moreau_bentos`: factor (IC 95 %) del contraste de los registros de bentos de Moreau frente a las medias de población de Hrycik; en M, una especie de bentos medida solo por Moreau se predice a su masa dividida por ese factor.", "",
           md_table(param), "",
           "## Predicción frente a los datos de cada especie", "",
           "Mediana, por fuente, de predicción/media de los registros en las especies medidas por una sola fuente (1 = la predicción coincide con sus datos).", "",
           md_table(t1), "",
           "## Validación cruzada por grupo", "",
           "`unit`: 10 grupos de registros; `source`: cada fuente fuera por turno. `datos_prediccion`: media geométrica de dato/predicción.", "",
           md_table(t2), "",
           "## Validación cruzada global", "",
           md_table(t3), "",
           "## Lectura", "",
           "- Con fuente (A y C), la predicción es para una fuente media, y las especies de una sola fuente se desplazan hacia ella: los rotíferos de Ejsmont-Karabin y Bottrell quedan 2-3 veces por encima de sus datos, y las especies de Moreau, 3-5 veces por debajo. C no lo arregla: las diferencias están dentro de cada grupo (pesadas frente a volúmenes en rotíferos; Moreau frente a Hrycik en el bentos), y en C los microcrustáceos pasan a quedar por encima de los suyos.",
           "- Sin fuente (N), cada especie queda en sus datos, el contraste de las pesadas recoge la diferencia de método de los rotíferos, la cobertura se mantiene y los intervalos son más estrechos. Queda sin corregir la diferencia entre Moreau y Hrycik en el bentos (cobertura del esquema `source` más baja en ese grupo).",
           "- Moreau et al. (2021) recopila masas de presas individuales de estudios de dieta de aves y peces y de calorimetría (sobre todo Driver et al. 1974 y Driver 1981), a menudo de animales grandes o adultos; Hrycik et al. (2024) da medias de todos los individuos de muestras de campo, que es lo que necesita una biomasa. M estima esa diferencia con las familias que comparten (pocas, de ahí su intervalo) y predice para medias de población. Con Hrycik fuera (esquema `source`), M no puede calibrar el contraste y sus predicciones de bentos no sirven: depende de esa fuente.")
writeLines(lines, "dev/validacion/efecto_fuente.md", useBytes = TRUE)
cat(lines, sep = "\n")
