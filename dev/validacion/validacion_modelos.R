# Validación de los modelos jerárquicos de la versión 0.2.0 (fuera del check)
#
# Validación cruzada de los intervalos de predicción:
#  - "unit": 10 grupos de registros (masas típicas) o de ecuaciones (talla-masa):
#    una población nueva de una especie de la referencia;
#  - "species": 10 grupos de nodos de especie: una especie que no está en la
#    referencia, predicha desde el nivel más bajo que queda;
#  - "source": cada fuente (masas típicas) o estudio (talla-masa), fuera por turno.
# Para cada unidad apartada (registro; centro del rango de cada ecuación; media
# de población): z = (y - eta) / sqrt(v_total); cobertura de los intervalos del
# 50, 90 y 95 %, sesgo y error absoluto mediano (ln), por modelo, esquema y
# nivel de la predicción.
# Escribe dev/validacion/validacion_modelos.csv (una fila por unidad apartada)
# y dev/validacion/validacion_modelos.md.
# Desde la raíz del paquete: Rscript dev/validacion/validacion_modelos.R  (unos 30 min)

Sys.setenv(BMASSA_FIT_ONLY_FUNCTIONS = "true")
source("data-raw/fit_models.R")
dat <- read_model_data()
specs <- model_specs(dat$typ, dat$lwr)
cores <- min(2L, parallel::detectCores())
k_folds <- 10L
seed <- 20261102L

one_fold <- function(d_all, type, test, model, scheme, fold) {
  train <- d_all[!test, ]
  levels <- split_levels(train)
  opt <- model_options(train, type)
  train <- model_frame(train, type, levels)
  fit <- suppressMessages(suppressWarnings(fit_model(train, type, levels, opt)))
  tabs <- model_tables(fit, train, model, type, levels)
  nd <- d_all[test, ]
  if (type == "lwr") nd <- nd[nd$point %in% c(0L, 2L), ]
  pr <- predict_tables(tabs, nd, type, levels)
  unit <- if (type == "lwr") nd$equation_id else nd$record_id
  data.frame(model = model, scheme = scheme, fold = fold, unit = unit, source_id = nd$source_id,
             taxon = nd$taxon, lowest_rank = nd$lowest_rank,
             length_type = if (type == "lwr") nd$length_type else NA_character_,
             y = nd$y, pr, stringsAsFactors = FALSE)
}

# Resultados parciales por modelo en dev/validacion/cv/ (se reutilizan si
# existen; BMASSA_CV_REDO=true los recalcula).
cv_dir <- "dev/validacion/cv"
dir.create(cv_dir, showWarnings = FALSE)
redo <- identical(Sys.getenv("BMASSA_CV_REDO"), "true")
res <- list()
for (i in seq_len(nrow(specs))) {
  m <- specs$model[i]
  type <- specs$type[i]
  part <- file.path(cv_dir, paste0(m, ".csv"))
  if (file.exists(part) && !redo) {
    res[[m]] <- read.csv(part, stringsAsFactors = FALSE, encoding = "UTF-8")
    next
  }
  d <- model_rows(dat, m, type)
  unit <- if (type == "lwr") d$equation_id else d$record_id
  sp <- node_path(d, "species")
  src <- if (type == "lwr") d$study else d$source_id
  set.seed(seed)
  u <- unique(unit)
  fold_unit <- stats::setNames(sample(rep(seq_len(k_folds), length.out = length(u))), u)
  s <- unique(sp)
  fold_sp <- stats::setNames(sample(rep(seq_len(k_folds), length.out = length(s))), s)
  jobs <- c(lapply(seq_len(k_folds), function(f) list(scheme = "unit", fold = f, test = fold_unit[unit] == f)),
            lapply(seq_len(k_folds), function(f) list(scheme = "species", fold = f, test = fold_sp[sp] == f)),
            if (length(unique(src)) > 1)
              lapply(sort(unique(src)), function(z) list(scheme = "source", fold = z, test = src == z)))
  t0 <- Sys.time()
  out <- parallel::mclapply(jobs, function(j) {
    tryCatch(one_fold(d, type, j$test, m, j$scheme, as.character(j$fold)),
             error = function(e) {
               data.frame(model = m, scheme = j$scheme, fold = as.character(j$fold), unit = NA_character_,
                          source_id = NA_character_, taxon = NA_character_, lowest_rank = NA_character_,
                          length_type = NA_character_, y = NA_real_, eta = NA_real_, v_param = NA_real_,
                          v_missing = NA_real_, v_population = NA_real_, v_total = NA_real_,
                          matched_level = paste("error:", conditionMessage(e)), stringsAsFactors = FALSE)
             })
  }, mc.cores = cores)
  res[[m]] <- do.call(rbind, out)
  write.csv(res[[m]], part, row.names = FALSE, fileEncoding = "UTF-8")
  message(sprintf("%s: %d fits in %.1f min", m, length(jobs),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, res)
rownames(res) <- NULL
res$z <- (res$y - res$eta) / sqrt(res$v_total)
write.csv(res, "dev/validacion/validacion_modelos.csv", row.names = FALSE, fileEncoding = "UTF-8")

# Resumen --------------------------------------------------------------------------
ok <- is.finite(res$z)
summ <- function(x) {
  data.frame(n = nrow(x), cob50 = mean(abs(x$z) < stats::qnorm(0.75)), cob90 = mean(abs(x$z) < stats::qnorm(0.95)),
             cob95 = mean(abs(x$z) < stats::qnorm(0.975)), z_media = mean(x$z), z_de = stats::sd(x$z),
             sesgo = mean(x$y - x$eta), eam = stats::median(abs(x$y - x$eta)),
             factor95 = stats::median(exp(1.96 * sqrt(x$v_total))))
}
r <- res[ok, ]
lv <- c("root", levels_all)
by_level <- do.call(rbind, lapply(split(r, list(r$model, r$scheme, r$matched_level), drop = TRUE), function(x)
  cbind(x[1, c("model", "scheme", "matched_level")], summ(x))))
by_level <- by_level[order(by_level$model, by_level$scheme, match(by_level$matched_level, lv)), ]
overall <- do.call(rbind, lapply(split(r, list(r$model, r$scheme), drop = TRUE), function(x)
  cbind(x[1, c("model", "scheme")], summ(x))))
overall <- overall[order(overall$model, overall$scheme), ]
fmt <- function(x) {
  x$cob50 <- sprintf("%.0f %%", 100 * x$cob50)
  x$cob90 <- sprintf("%.0f %%", 100 * x$cob90)
  x$cob95 <- sprintf("%.0f %%", 100 * x$cob95)
  for (j in c("z_media", "z_de", "sesgo", "eam")) x[[j]] <- sprintf("%.2f", x[[j]])
  x$factor95 <- sprintf("%.1f", x$factor95)
  x
}
md_table <- function(x) {
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep("---|", ncol(x)), collapse = "")),
    apply(x, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |")))
}
not_pred <- res[!ok, ]
not_pred$motivo <- ifelse(grepl("^error:", not_pred$matched_level),
                          paste("ajuste imposible:", sub("^error: ", "", not_pred$matched_level)),
                          ifelse(!is.na(not_pred$length_type) & not_pred$matched_level == "root",
                                 "tipo de longitud sin datos en el ajuste",
                                 "estadio sin datos en el ajuste"))
np <- if (nrow(not_pred)) {
  a <- stats::aggregate(list(n = rep(1L, nrow(not_pred))),
                        by = not_pred[, c("model", "scheme", "fold", "motivo")], FUN = sum)
  a[order(a$model, a$scheme, a$fold), ]
} else NULL
lines <- c("# Validación cruzada de los modelos jerárquicos (versión 0.2.0)", "",
           paste0("Generado por `dev/validacion/validacion_modelos.R` el ", format(Sys.Date()), " (",
                  R.version.string, ", lme4 ", as.character(utils::packageVersion("lme4")), "; semilla ", seed, ")."),
           "",
           "Esquemas: `unit`, población nueva de una especie de la referencia (10 grupos de registros o ecuaciones); `species`, especie nueva (10 grupos de especies); `source`, fuente o estudio nuevo (cada uno fuera por turno). Para las ecuaciones se evalúa el centro de su rango de tallas. `cob50`, `cob90`, `cob95`: proporción de valores observados dentro de los intervalos de predicción del 50, 90 y 95 %; `z_de`: desviación típica de z (1 si la varianza está bien calibrada); `sesgo` y `eam`: media de y − η y mediana de |y − η| (ln); `factor95`: mediana del factor ×/÷ del intervalo del 95 %.",
           "", "## Global", "", md_table(fmt(overall)), "", "## Por nivel de la predicción", "",
           md_table(fmt(by_level)), "",
           "## Lectura", "",
           "- `unit` y `species` miden lo que el modelo promete: una población nueva de una especie conocida y una especie nueva de un taxón superior conocido. La cobertura del 95 % debe quedar entre el 90 y el 97 %.",
           "- `source` es más exigente de lo que será el uso real: con dos fuentes (fitoplancton, anfibios, peces), al apartar una no queda varianza entre fuentes que estimar, de modo que los intervalos no la incluyen; en invertebrados las fuentes cubren grupos distintos (zooplancton de Dumont, bentos de Hrycik y Moreau) y apartar una obliga a predecir desde clases o filos con otra composición de tamaños. En el modelo completo la varianza entre fuentes sí entra en la predicción.",
           "",
           sprintf("## Unidades sin predicción (%d)", nrow(not_pred)), "",
           if (is.null(np)) "Ninguna." else md_table(np))
writeLines(lines, "dev/validacion/validacion_modelos.md", useBytes = TRUE)
cat(lines, sep = "\n")
