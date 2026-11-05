# Efecto de fuente en el modelo de invertebrados sin longitud

Generado por `dev/validacion/efecto_fuente.R` el 2026-11-05 (R version 4.3.3 (2024-02-29), lme4 1.1.35.1; semilla 20261105).

Variantes: `A`, efecto de fuente cruzado (versión 0.2.0); `C`, fuente dentro de grupo de taxones (rotíferos, microcrustáceos y el resto); `N`, sin efecto de fuente; `M`, sin efecto de fuente y con un contraste fijo para los registros de bentos de Moreau et al. (2021) (versión 0.3.0). A, C y N se ajustan sin ese contraste. Las predicciones de las tablas llevan los contrastes de cada registro (pesadas de rotíferos, estadios de copépodos y, en M, el de Moreau).

## Parámetros

`pesadas`: factor del contraste de las pesadas de rotíferos frente a las masas por volumen. `moreau_bentos`: factor (IC 95 %) del contraste de los registros de bentos de Moreau frente a las medias de población de Hrycik; en M, una especie de bentos medida solo por Moreau se predice a su masa dividida por ese factor.

| variante | de_fuente | de_residual | pesadas | moreau_bentos |
|---|---|---|---|---|
| A | 0.76 | 0.89 | 2.1 | — |
| C | 0.84 | 0.89 | 2.0 | — |
| N | NA | 0.92 | 4.0 | — |
| M | NA | 0.91 | 4.0 | 6.4 (3.6-11.4) |

## Predicción frente a los datos de cada especie

Mediana, por fuente, de predicción/media de los registros en las especies medidas por una sola fuente (1 = la predicción coincide con sus datos).

| fuente | especies | A | C | N | M |
|---|---|---|---|---|---|
| bottrell1976 |   9 | 1.94 | 1.67 | 0.93 | 0.93 |
| culver1985 |   2 | 0.80 | 1.60 | 1.02 | 1.03 |
| dumont1975 |  14 | 0.68 | 1.40 | 1.13 | 1.13 |
| ejsmontkarabin1998 | 206 | 2.32 | 2.00 | 1.12 | 1.11 |
| hebert2016 |  78 | 0.65 | 1.36 | 1.09 | 1.09 |
| hrycik2024 |  64 | 1.74 | 1.18 | 1.06 | 1.01 |
| moreau2021 |  62 | 0.31 | 0.22 | 0.83 | 0.85 |

## Validación cruzada por grupo

`unit`: 10 grupos de registros; `source`: cada fuente fuera por turno. `datos_prediccion`: media geométrica de dato/predicción.

| esquema | grupo | variante | n | cob95 | datos_prediccion |
|---|---|---|---|---|---|
| source | microcrustáceos | A | 240 | 97 % | 1.84 |
| source | microcrustáceos | C | 240 | 98 % | 0.65 |
| source | microcrustáceos | M | 211 | 100 % | 0.83 |
| source | microcrustáceos | N | 240 | 96 % | 1.08 |
| source | otros invertebrados | A | 317 | 77 % | 3.02 |
| source | otros invertebrados | C | 317 | 76 % | 2.56 |
| source | otros invertebrados | M | 191 | 94 % | 12.73 |
| source | otros invertebrados | N | 317 | 82 % | 0.42 |
| source | rotíferos | A | 347 | 93 % | 0.33 |
| source | rotíferos | C | 347 | 97 % | 0.09 |
| source | rotíferos | M | 347 | 96 % | 0.64 |
| source | rotíferos | N | 347 | 97 % | 0.64 |
| unit | microcrustáceos | A | 257 | 99 % | 1.73 |
| unit | microcrustáceos | C | 257 | 99 % | 0.78 |
| unit | microcrustáceos | M | 257 | 99 % | 1.03 |
| unit | microcrustáceos | N | 257 | 99 % | 1.03 |
| unit | otros invertebrados | A | 388 | 94 % | 1.31 |
| unit | otros invertebrados | C | 388 | 93 % | 2.26 |
| unit | otros invertebrados | M | 388 | 94 % | 1.00 |
| unit | otros invertebrados | N | 388 | 94 % | 1.00 |
| unit | rotíferos | A | 347 | 99 % | 0.55 |
| unit | rotíferos | C | 347 | 99 % | 0.61 |
| unit | rotíferos | M | 347 | 98 % | 1.00 |
| unit | rotíferos | N | 347 | 99 % | 1.00 |

## Validación cruzada global

| esquema | variante | n | cob95 | factor95 |
|---|---|---|---|---|
| source | A | 904 | 88 % | 35.0 |
| source | C | 904 | 90 % | 139.3 |
| source | M | 749 | 96 % | 20.3 |
| source | N | 904 | 91 % | 22.7 |
| unit | A | 992 | 97 % | 22.5 |
| unit | C | 992 | 97 % | 26.4 |
| unit | M | 992 | 96 % | 13.7 |
| unit | N | 992 | 97 % | 16.7 |

## Lectura

- Con fuente (A y C), la predicción es para una fuente media, y las especies de una sola fuente se desplazan hacia ella: los rotíferos de Ejsmont-Karabin y Bottrell quedan 2-3 veces por encima de sus datos, y las especies de Moreau, 3-5 veces por debajo. C no lo arregla: las diferencias están dentro de cada grupo (pesadas frente a volúmenes en rotíferos; Moreau frente a Hrycik en el bentos), y en C los microcrustáceos pasan a quedar por encima de los suyos.
- Sin fuente (N), cada especie queda en sus datos, el contraste de las pesadas recoge la diferencia de método de los rotíferos, la cobertura se mantiene y los intervalos son más estrechos. Queda sin corregir la diferencia entre Moreau y Hrycik en el bentos (cobertura del esquema `source` más baja en ese grupo).
- Moreau et al. (2021) recopila masas de presas individuales de estudios de dieta de aves y peces y de calorimetría (sobre todo Driver et al. 1974 y Driver 1981), a menudo de animales grandes o adultos; Hrycik et al. (2024) da medias de todos los individuos de muestras de campo, que es lo que necesita una biomasa. M estima esa diferencia con las familias que comparten (pocas, de ahí su intervalo) y predice para medias de población. Con Hrycik fuera (esquema `source`), M no puede calibrar el contraste y sus predicciones de bentos no sirven: depende de esa fuente.
