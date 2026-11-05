# Flujo completo de un usuario de bmassa 0.2.0 (ejecutable con Rscript --vanilla
# en un directorio vacío). Se convierte en el ejemplo del README, en la viñeta y
# en tests/testthat/test-workflow.R.
library(bmassa)

# 1. Dos tablas: las unidades de taxón (taxón, o taxón en un estadio y sexo, con
#    el grupo ecológico, la clasificación, la unidad de conteo y la talla media
#    si se midió) y sus abundancias por sitio (y fecha, si la hay). Una tercera,
#    opcional, con las tallas medidas en muestras concretas.
tx <- example_survey$taxa
ab <- example_survey$abundance
ln <- example_survey$lengths
head(tx[, c("taxon_id", "broad_group", "scientific_name", "life_stage", "count_unit",
            "length_mean", "length_unit", "length_type")])
head(ab)
ln

# 0. (Con internet) La clasificación y el col_id de una lista de nombres propia
#    salen de standardize_taxa(), que consulta el Catalogue of Life (COL XR):
#    mis_taxa <- data.frame(taxon_id = c("a", "b"), broad_group = "macroinvertebrate",
#                           scientific_name = c("Baetis rhodani", "Ephemerella ignita"))
#    standardize_taxa(mis_taxa)

# 2. Antes de estimar: de qué nivel taxonómico saldrá cada taxón y cuánto
#    ancho tendrá el intervalo (resultado intermedio que condiciona la estimación).
cov <- check_coverage(tx, ln)
cov
as.data.frame(cov)[, c("taxon_id", "method", "match_level", "match_name", "factor", "status")]

# 3. Masa media individual de una población nueva de cada taxón (media, mediana
#    e intervalo de predicción del 95 %) y biomasa por fila. Con tallas hay que
#    decidir la retransformación: "lognormal" da la media de los individuos de
#    cada talla; "none", su mediana.
est <- estimate_biomass(ab, tx, ln, bias_correction = "lognormal")
est
e <- as.data.frame(est)
e[, c("site", "taxon_id", "method", "match_level", "ind_median", "ind_lower", "ind_upper",
      "unit", "status")]
e[e$status != "estimated", c("site", "taxon_id", "life_stage", "reason")]

# 4. Otra moneda, con la incertidumbre de los factores: carbono en µg para
#    todos los grupos (fitoplancton por célula con Menden-Deuer y Lessard 2000;
#    anfibios sin conversión).
est_c <- convert_currency(est, to = "carbon", unit = "ug")
ec <- as.data.frame(est_c)
ec[ec$site == "L1", c("taxon_id", "converted_from", "ind_median", "ind_lower", "ind_upper",
                      "unit", "status")]

# 5. Totales por sitio y grupo con extracciones compartidas (dos especies
#    estimadas desde la misma familia tienen errores correlacionados); nunca se
#    suman monedas ni unidades distintas.
tot <- summarize_biomass(est_c, by = c("site", "broad_group"), na_rm = TRUE, seed = 1)
tot

# 6. Extracciones conjuntas para otros paquetes (p. ej., los flujos de energía
#    de una red trófica en bwebs): masa media individual y biomasa de cada fila
#    en cada extracción, con las mismas correlaciones que los totales. Con
#    population_by, cada hábitat tiene sus poblaciones. Para sumar grupos
#    muestreados por litro y por m², todas las abundancias en la misma base
#    antes de estimar: convert_abundance() con el tamaño de cada muestra en
#    `per` y la profundidad, el volumen o la superficie de cada sitio (por m²
#    aquí; por tanque con to = "site").
sites <- data.frame(site = c("L1", "L2"), depth = c(1.5, 0.8), depth_unit = "m")
ab_h <- convert_abundance(ab, to = "m2", sites = sites)
unique(ab_h[, c("per_original", "per", "per_factor", "per_assumption")])
ab_h$habitat <- ifelse(ab_h$site == "L1", "lake", "spring")
est_w <- convert_currency(estimate_biomass(ab_h, tx, ln, bias_correction = "lognormal"),
                          to = "wet_mass", unit = "mg")
sim <- simulate(est_w, nsim = 1000, seed = 1, population_by = "habitat")
sim
d <- as.data.frame(sim, draws = 1:2)
d[d$draw == 1, c("site", "taxon_id", "habitat", "ind_mass", "biomass", "biomass_unit")]

# 7. Qué fuentes citar.
cite_sources(est_c)[, c("source_id", "citation", "licence")]

# 8. Exportación (nunca al directorio de trabajo por defecto).
out <- tempfile(fileext = ".csv")
utils::write.csv(ec, out, row.names = FALSE)
