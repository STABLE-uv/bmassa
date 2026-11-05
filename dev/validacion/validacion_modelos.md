# Validación cruzada de los modelos jerárquicos (versión 0.2.0)

Generado por `dev/validacion/validacion_modelos.R` el 2026-11-05 (R version 4.3.3 (2024-02-29), lme4 1.1.35.1; semilla 20261102).

Esquemas: `unit`, población nueva de una especie de la referencia (10 grupos de registros o ecuaciones); `species`, especie nueva (10 grupos de especies); `source`, fuente o estudio nuevo (cada uno fuera por turno). Para las ecuaciones se evalúa el centro de su rango de tallas. `cob50`, `cob90`, `cob95`: proporción de valores observados dentro de los intervalos de predicción del 50, 90 y 95 %; `z_de`: desviación típica de z (1 si la varianza está bien calibrada); `sesgo` y `eam`: media de y − η y mediana de |y − η| (ln); `factor95`: mediana del factor ×/÷ del intervalo del 95 %.

## Global

| model | scheme | n | cob50 | cob90 | cob95 | z_media | z_de | sesgo | eam | factor95 |
|---|---|---|---|---|---|---|---|---|---|---|
| amphibian_adult | source |  290 | 53 % | 82 % | 86 % | -0.40 | 2.05 | -0.28 | 0.62 | 10.3 |
| amphibian_adult | species |  290 | 57 % | 91 % | 94 % | -0.01 | 0.95 | 0.00 | 0.56 | 7.3 |
| amphibian_adult | unit |  290 | 58 % | 90 % | 93 % | -0.11 | 1.06 | -0.10 | 0.42 | 6.4 |
| amphibian_lwr | source |   55 | 42 % | 76 % | 87 % | -0.11 | 1.55 | -0.07 | 0.25 | 1.8 |
| amphibian_lwr | species |  298 | 55 % | 89 % | 94 % | 0.00 | 1.04 | -0.01 | 0.20 | 1.9 |
| amphibian_lwr | unit |  298 | 57 % | 90 % | 95 % | -0.01 | 1.04 | -0.01 | 0.18 | 1.8 |
| fish_individual | source |   36 | 89 % | 100 % | 100 % | -0.16 | 0.46 | -0.34 | 0.76 | 49.5 |
| fish_individual | species |  421 | 48 % | 89 % | 91 % | 0.09 | 1.13 | 0.14 | 1.07 | 20.3 |
| fish_individual | unit |  421 | 60 % | 91 % | 94 % | 0.01 | 1.00 | 0.01 | 0.61 | 8.9 |
| fish_lwr | source |   44 | 66 % | 95 % | 98 % | 0.26 | 0.83 | 0.25 | 0.48 | 7.9 |
| fish_lwr | species |  355 | 72 % | 92 % | 96 % | -0.15 | 0.79 | -0.02 | 0.12 | 1.6 |
| fish_lwr | unit |  355 | 57 % | 93 % | 97 % | -0.00 | 0.91 | -0.00 | 0.07 | 1.2 |
| invertebrate_individual | source |  904 | 47 % | 85 % | 88 % | -0.02 | 1.21 | 0.13 | 1.45 | 35.0 |
| invertebrate_individual | species |  992 | 61 % | 94 % | 97 % | 0.03 | 0.84 | 0.08 | 0.88 | 26.0 |
| invertebrate_individual | unit |  992 | 62 % | 95 % | 97 % | 0.01 | 0.87 | 0.02 | 0.82 | 22.0 |
| invertebrate_lwr | source | 1107 | 58 % | 93 % | 96 % | -0.13 | 0.94 | -0.22 | 0.46 | 4.1 |
| invertebrate_lwr | species | 1107 | 64 % | 94 % | 96 % | -0.04 | 0.88 | -0.03 | 0.35 | 4.4 |
| invertebrate_lwr | unit | 1107 | 64 % | 93 % | 96 % | -0.06 | 0.90 | -0.06 | 0.34 | 3.7 |
| phytoplankton_cell | source | 1530 | 42 % | 82 % | 88 % | 0.24 | 1.27 | 0.39 | 1.11 | 16.4 |
| phytoplankton_cell | species | 1530 | 55 % | 92 % | 96 % | 0.11 | 0.94 | 0.17 | 0.94 | 21.7 |
| phytoplankton_cell | unit | 1530 | 57 % | 92 % | 96 % | 0.14 | 0.93 | 0.21 | 0.82 | 20.4 |
| phytoplankton_natural_unit | source | 1625 | 46 % | 86 % | 92 % | -0.16 | 1.12 | -0.31 | 1.39 | 51.7 |
| phytoplankton_natural_unit | species | 1625 | 58 % | 93 % | 97 % | -0.18 | 0.86 | -0.33 | 1.06 | 38.0 |
| phytoplankton_natural_unit | unit | 1625 | 58 % | 94 % | 97 % | -0.18 | 0.86 | -0.35 | 1.03 | 37.6 |

## Por nivel de la predicción

| model | scheme | matched_level | n | cob50 | cob90 | cob95 | z_media | z_de | sesgo | eam | factor95 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| amphibian_adult | source | order |  128 | 54 % | 84 % | 90 % | -0.25 | 1.15 | -0.34 | 0.76 | 13.4 |
| amphibian_adult | source | family |   75 | 56 % | 91 % | 92 % | -0.28 | 0.98 | -0.32 | 0.72 | 9.9 |
| amphibian_adult | source | genus |   52 | 62 % | 87 % | 87 % | -0.27 | 1.49 | -0.19 | 0.41 | 3.8 |
| amphibian_adult | source | species |   35 | 34 % | 51 % | 57 % | -1.42 | 4.91 | -0.13 | 0.44 | 1.4 |
| amphibian_adult | species | order |   11 | 45 % | 64 % | 73 % | 0.46 | 1.44 | 0.82 | 1.83 | 50.3 |
| amphibian_adult | species | family |   61 | 48 % | 92 % | 97 % | -0.10 | 0.97 | -0.11 | 0.86 | 9.6 |
| amphibian_adult | species | genus |  218 | 61 % | 92 % | 95 % | -0.00 | 0.91 | -0.01 | 0.51 | 6.9 |
| amphibian_adult | unit | order |    8 | 12 % | 50 % | 75 % | 0.04 | 1.85 | 0.05 | 3.15 | 39.4 |
| amphibian_adult | unit | family |   37 | 35 % | 86 % | 95 % | -0.12 | 1.11 | -0.13 | 1.06 | 9.5 |
| amphibian_adult | unit | genus |  115 | 59 % | 90 % | 92 % | -0.12 | 0.96 | -0.12 | 0.49 | 6.9 |
| amphibian_adult | unit | species |  130 | 66 % | 93 % | 94 % | -0.12 | 1.07 | -0.08 | 0.25 | 3.5 |
| amphibian_lwr | source | order |    8 | 38 % | 75 % | 75 % | -0.86 | 1.52 | -0.46 | 0.30 | 2.1 |
| amphibian_lwr | source | family |   12 | 50 % | 75 % | 92 % | -0.30 | 1.41 | -0.08 | 0.22 | 1.8 |
| amphibian_lwr | source | genus |   22 | 50 % | 86 % | 100 % | 0.04 | 1.02 | 0.00 | 0.19 | 1.8 |
| amphibian_lwr | source | species |   13 | 23 % | 62 % | 69 % | 0.27 | 2.28 | 0.08 | 0.31 | 1.6 |
| amphibian_lwr | species | order |   13 | 77 % | 85 % | 85 % | -0.11 | 1.19 | -0.17 | 0.14 | 3.8 |
| amphibian_lwr | species | family |   67 | 54 % | 88 % | 94 % | -0.13 | 1.04 | -0.03 | 0.26 | 2.0 |
| amphibian_lwr | species | genus |  218 | 54 % | 89 % | 94 % | 0.05 | 1.03 | 0.01 | 0.19 | 1.9 |
| amphibian_lwr | unit | order |    6 | 67 % | 83 % | 83 % | -0.68 | 1.21 | -0.59 | 0.22 | 4.9 |
| amphibian_lwr | unit | family |   40 | 48 % | 85 % | 95 % | -0.02 | 1.23 | -0.01 | 0.32 | 2.0 |
| amphibian_lwr | unit | genus |  109 | 58 % | 89 % | 95 % | 0.04 | 0.98 | 0.01 | 0.18 | 1.9 |
| amphibian_lwr | unit | species |  143 | 58 % | 92 % | 94 % | -0.01 | 1.02 | -0.00 | 0.16 | 1.7 |
| fish_individual | source | root |    8 | 62 % | 100 % | 100 % | -0.21 | 0.70 | -0.51 | 0.98 | 117.3 |
| fish_individual | source | order |   10 | 100 % | 100 % | 100 % | -0.19 | 0.40 | -0.38 | 0.91 | 66.6 |
| fish_individual | source | family |   18 | 94 % | 100 % | 100 % | -0.13 | 0.37 | -0.24 | 0.49 | 34.1 |
| fish_individual | species | root |   54 | 46 % | 96 % | 100 % | 0.54 | 0.68 | 1.15 | 1.51 | 65.0 |
| fish_individual | species | order |   73 | 73 % | 100 % | 100 % | -0.27 | 0.55 | -0.46 | 0.67 | 24.6 |
| fish_individual | species | family |  245 | 37 % | 82 % | 84 % | 0.07 | 1.36 | 0.08 | 1.22 | 14.2 |
| fish_individual | species | genus |   49 | 71 % | 98 % | 100 % | 0.21 | 0.60 | 0.26 | 0.59 | 10.7 |
| fish_individual | unit | root |    2 | 50 % | 100 % | 100 % | 0.82 | 0.22 | 1.78 | 1.78 | 69.6 |
| fish_individual | unit | order |    7 | 86 % | 100 % | 100 % | 0.14 | 0.61 | 0.24 | 0.84 | 32.4 |
| fish_individual | unit | family |   17 | 65 % | 100 % | 100 % | -0.07 | 0.56 | -0.12 | 0.53 | 19.4 |
| fish_individual | unit | genus |   14 | 71 % | 100 % | 100 % | -0.11 | 0.63 | -0.14 | 0.47 | 12.4 |
| fish_individual | unit | species |  381 | 59 % | 90 % | 93 % | 0.01 | 1.04 | 0.01 | 0.60 | 8.8 |
| fish_lwr | source | root |   34 | 65 % | 94 % | 97 % | 0.42 | 0.77 | 0.36 | 0.52 | 8.2 |
| fish_lwr | source | order |    7 | 86 % | 100 % | 100 % | -0.07 | 0.54 | -0.06 | 0.15 | 3.3 |
| fish_lwr | source | family |    3 | 33 % | 100 % | 100 % | -0.89 | 1.14 | -0.23 | 0.43 | 1.7 |
| fish_lwr | species | root |   47 | 74 % | 98 % | 98 % | 0.15 | 0.61 | 0.12 | 0.15 | 4.0 |
| fish_lwr | species | order |   68 | 76 % | 90 % | 97 % | -0.02 | 0.87 | -0.02 | 0.31 | 3.3 |
| fish_lwr | species | family |  188 | 81 % | 95 % | 97 % | -0.18 | 0.66 | -0.05 | 0.07 | 1.6 |
| fish_lwr | species | genus |   52 | 31 % | 81 % | 90 % | -0.50 | 1.08 | -0.08 | 0.18 | 1.4 |
| fish_lwr | unit | root |    4 | 0 % | 50 % | 75 % | -0.98 | 1.66 | -0.66 | 0.96 | 3.7 |
| fish_lwr | unit | order |   16 | 38 % | 88 % | 100 % | -0.05 | 1.11 | -0.02 | 0.62 | 3.5 |
| fish_lwr | unit | family |   27 | 52 % | 89 % | 93 % | 0.21 | 1.03 | 0.06 | 0.17 | 1.6 |
| fish_lwr | unit | genus |   27 | 63 % | 96 % | 96 % | -0.12 | 0.81 | -0.02 | 0.11 | 1.5 |
| fish_lwr | unit | species |  281 | 58 % | 94 % | 97 % | 0.00 | 0.88 | 0.00 | 0.05 | 1.2 |
| invertebrate_individual | source | root |    2 | 100 % | 100 % | 100 % | -0.09 | 0.55 | -0.45 | 1.91 | 16641.4 |
| invertebrate_individual | source | phylum |   25 | 52 % | 92 % | 92 % | 0.30 | 0.97 | 0.79 | 2.64 | 3253.5 |
| invertebrate_individual | source | class |  172 | 41 % | 78 % | 79 % | 0.25 | 1.47 | 0.55 | 2.34 | 237.5 |
| invertebrate_individual | source | order |  124 | 32 % | 73 % | 77 % | 0.49 | 1.47 | 1.39 | 2.36 | 88.2 |
| invertebrate_individual | source | family |  159 | 46 % | 89 % | 94 % | -0.23 | 1.05 | -0.27 | 1.36 | 38.6 |
| invertebrate_individual | source | genus |  235 | 47 % | 83 % | 89 % | -0.34 | 1.11 | -0.49 | 1.14 | 25.7 |
| invertebrate_individual | source | species |  187 | 61 % | 95 % | 97 % | -0.06 | 0.82 | -0.08 | 0.85 | 19.3 |
| invertebrate_individual | species | root |    2 | 100 % | 100 % | 100 % | 0.02 | 0.79 | 0.01 | 2.32 | 3491.6 |
| invertebrate_individual | species | phylum |    2 | 50 % | 100 % | 100 % | -0.86 | 0.50 | -2.77 | 2.77 | 770.9 |
| invertebrate_individual | species | class |   32 | 66 % | 91 % | 94 % | 0.10 | 0.97 | 0.25 | 0.85 | 145.5 |
| invertebrate_individual | species | order |   91 | 42 % | 84 % | 92 % | 0.55 | 1.02 | 1.22 | 1.87 | 69.2 |
| invertebrate_individual | species | family |  267 | 57 % | 91 % | 96 % | 0.08 | 0.92 | 0.14 | 1.06 | 31.5 |
| invertebrate_individual | species | genus |  598 | 65 % | 98 % | 99 % | -0.08 | 0.73 | -0.12 | 0.78 | 22.2 |
| invertebrate_individual | unit | root |    2 | 100 % | 100 % | 100 % | -0.00 | 0.71 | -0.00 | 2.08 | 3417.5 |
| invertebrate_individual | unit | phylum |    2 | 50 % | 100 % | 100 % | -0.85 | 0.51 | -2.75 | 2.75 | 780.6 |
| invertebrate_individual | unit | class |    7 | 57 % | 86 % | 86 % | 0.50 | 1.25 | 1.22 | 1.75 | 152.7 |
| invertebrate_individual | unit | order |   48 | 65 % | 96 % | 98 % | 0.02 | 0.78 | 0.02 | 0.97 | 69.6 |
| invertebrate_individual | unit | family |  128 | 49 % | 89 % | 96 % | 0.01 | 1.02 | 0.03 | 1.11 | 31.7 |
| invertebrate_individual | unit | genus |  359 | 63 % | 98 % | 99 % | -0.17 | 0.73 | -0.27 | 0.81 | 21.7 |
| invertebrate_individual | unit | species |  446 | 64 % | 94 % | 96 % | 0.16 | 0.91 | 0.24 | 0.78 | 21.1 |
| invertebrate_lwr | source | root |    4 | 25 % | 75 % | 100 % | 0.35 | 1.49 | 0.68 | 3.30 | 105.7 |
| invertebrate_lwr | source | phylum |    2 | 100 % | 100 % | 100 % | -0.52 | 0.09 | -0.73 | 0.73 | 15.8 |
| invertebrate_lwr | source | class |   41 | 49 % | 98 % | 100 % | -0.30 | 0.77 | -0.49 | 1.00 | 12.7 |
| invertebrate_lwr | source | order |  257 | 46 % | 95 % | 98 % | -0.55 | 0.79 | -0.86 | 1.06 | 18.2 |
| invertebrate_lwr | source | family |  239 | 59 % | 90 % | 94 % | -0.03 | 1.00 | -0.06 | 0.45 | 4.6 |
| invertebrate_lwr | source | genus |  175 | 65 % | 97 % | 97 % | -0.00 | 0.84 | -0.01 | 0.32 | 3.5 |
| invertebrate_lwr | source | species |  389 | 63 % | 91 % | 95 % | 0.05 | 0.96 | 0.03 | 0.33 | 3.3 |
| invertebrate_lwr | species | root |    4 | 25 % | 75 % | 100 % | 0.30 | 1.56 | 0.90 | 3.59 | 177.4 |
| invertebrate_lwr | species | phylum |    2 | 100 % | 100 % | 100 % | -0.42 | 0.04 | -0.58 | 0.58 | 14.9 |
| invertebrate_lwr | species | class |   28 | 64 % | 96 % | 100 % | -0.13 | 0.77 | -0.13 | 0.41 | 8.9 |
| invertebrate_lwr | species | order |  135 | 60 % | 92 % | 93 % | -0.04 | 1.05 | -0.04 | 0.50 | 5.7 |
| invertebrate_lwr | species | family |  365 | 62 % | 91 % | 94 % | 0.03 | 0.96 | 0.02 | 0.40 | 4.4 |
| invertebrate_lwr | species | genus |  573 | 66 % | 97 % | 98 % | -0.07 | 0.79 | -0.06 | 0.30 | 3.7 |
| invertebrate_lwr | unit | root |    2 | 50 % | 50 % | 50 % | -0.81 | 1.80 | -1.49 | 2.81 | 164.9 |
| invertebrate_lwr | unit | phylum |    2 | 100 % | 100 % | 100 % | -0.42 | 0.08 | -0.58 | 0.58 | 15.7 |
| invertebrate_lwr | unit | class |    9 | 78 % | 100 % | 100 % | -0.02 | 0.62 | -0.03 | 0.42 | 8.7 |
| invertebrate_lwr | unit | order |   60 | 55 % | 95 % | 95 % | -0.10 | 1.00 | -0.10 | 0.62 | 6.0 |
| invertebrate_lwr | unit | family |  209 | 59 % | 89 % | 93 % | 0.00 | 1.00 | -0.00 | 0.46 | 4.4 |
| invertebrate_lwr | unit | genus |  327 | 71 % | 96 % | 98 % | -0.18 | 0.75 | -0.14 | 0.28 | 4.8 |
| invertebrate_lwr | unit | species |  498 | 62 % | 93 % | 96 % | -0.01 | 0.93 | -0.01 | 0.32 | 3.3 |
| phytoplankton_cell | source | root |    1 | 100 % | 100 % | 100 % | -0.25 | NA | -0.76 | 0.76 | 386.8 |
| phytoplankton_cell | source | kingdom |   26 | 69 % | 100 % | 100 % | -0.43 | 0.51 | -1.00 | 0.95 | 101.7 |
| phytoplankton_cell | source | phylum |   39 | 59 % | 100 % | 100 % | 0.01 | 0.73 | 0.02 | 1.04 | 54.0 |
| phytoplankton_cell | source | class |   54 | 56 % | 81 % | 87 % | 0.33 | 1.07 | 0.61 | 0.87 | 32.4 |
| phytoplankton_cell | source | order |   69 | 46 % | 87 % | 93 % | 0.42 | 0.97 | 0.66 | 1.16 | 21.6 |
| phytoplankton_cell | source | family |  220 | 52 % | 86 % | 90 % | 0.25 | 1.06 | 0.37 | 1.01 | 18.3 |
| phytoplankton_cell | source | genus |  622 | 37 % | 82 % | 88 % | 0.48 | 1.23 | 0.67 | 1.24 | 16.0 |
| phytoplankton_cell | source | species |  499 | 40 % | 77 % | 85 % | -0.05 | 1.45 | 0.10 | 1.01 | 14.2 |
| phytoplankton_cell | species | root |    1 | 100 % | 100 % | 100 % | -0.43 | NA | -1.12 | 1.12 | 176.4 |
| phytoplankton_cell | species | kingdom |    2 | 100 % | 100 % | 100 % | 0.42 | 0.24 | 1.06 | 1.06 | 136.7 |
| phytoplankton_cell | species | phylum |    6 | 17 % | 100 % | 100 % | -0.27 | 0.76 | -0.59 | 1.47 | 59.5 |
| phytoplankton_cell | species | class |   18 | 56 % | 78 % | 89 % | -0.07 | 1.16 | -0.11 | 0.81 | 42.3 |
| phytoplankton_cell | species | order |   37 | 54 % | 92 % | 100 % | 0.14 | 0.98 | 0.23 | 1.10 | 34.6 |
| phytoplankton_cell | species | family |  155 | 58 % | 94 % | 96 % | -0.02 | 0.92 | -0.03 | 0.96 | 26.0 |
| phytoplankton_cell | species | genus | 1311 | 55 % | 92 % | 96 % | 0.12 | 0.94 | 0.20 | 0.93 | 21.3 |
| phytoplankton_cell | unit | root |    1 | 100 % | 100 % | 100 % | -0.44 | NA | -1.17 | 1.17 | 178.0 |
| phytoplankton_cell | unit | kingdom |    2 | 100 % | 100 % | 100 % | 0.44 | 0.27 | 1.11 | 1.11 | 147.8 |
| phytoplankton_cell | unit | phylum |    3 | 33 % | 100 % | 100 % | 0.22 | 0.93 | 0.45 | 1.49 | 68.1 |
| phytoplankton_cell | unit | class |   16 | 50 % | 75 % | 88 % | -0.12 | 1.15 | -0.23 | 1.29 | 42.8 |
| phytoplankton_cell | unit | order |   26 | 46 % | 88 % | 100 % | 0.22 | 0.96 | 0.37 | 1.33 | 33.3 |
| phytoplankton_cell | unit | family |  103 | 58 % | 93 % | 97 % | 0.05 | 0.90 | 0.09 | 0.89 | 26.2 |
| phytoplankton_cell | unit | genus |  806 | 56 % | 92 % | 95 % | 0.19 | 0.93 | 0.30 | 0.93 | 21.6 |
| phytoplankton_cell | unit | species |  573 | 60 % | 93 % | 96 % | 0.08 | 0.93 | 0.10 | 0.68 | 13.4 |
| phytoplankton_natural_unit | source | kingdom |   87 | 69 % | 100 % | 100 % | 0.07 | 0.68 | 0.17 | 1.16 | 116.3 |
| phytoplankton_natural_unit | source | phylum |    6 | 17 % | 83 % | 100 % | -1.08 | 0.48 | -2.63 | 2.95 | 116.3 |
| phytoplankton_natural_unit | source | class |  833 | 52 % | 91 % | 96 % | -0.27 | 0.92 | -0.55 | 1.38 | 67.6 |
| phytoplankton_natural_unit | source | order |  105 | 51 % | 86 % | 92 % | -0.40 | 1.06 | -0.66 | 1.26 | 40.4 |
| phytoplankton_natural_unit | source | family |  100 | 33 % | 81 % | 88 % | -0.20 | 1.33 | -0.35 | 1.54 | 25.3 |
| phytoplankton_natural_unit | source | genus |  247 | 34 % | 74 % | 84 % | -0.05 | 1.39 | -0.06 | 1.53 | 19.9 |
| phytoplankton_natural_unit | source | species |  247 | 35 % | 78 % | 85 % | 0.15 | 1.37 | 0.28 | 1.40 | 17.0 |
| phytoplankton_natural_unit | species | kingdom |    2 | 100 % | 100 % | 100 % | -0.01 | 0.27 | -0.03 | 0.49 | 150.7 |
| phytoplankton_natural_unit | species | phylum |    3 | 67 % | 100 % | 100 % | -0.60 | 0.38 | -1.53 | 1.28 | 145.3 |
| phytoplankton_natural_unit | species | class |   19 | 42 % | 89 % | 95 % | -0.21 | 0.96 | -0.47 | 1.63 | 80.1 |
| phytoplankton_natural_unit | species | order |   39 | 51 % | 90 % | 95 % | -0.08 | 0.97 | -0.17 | 1.40 | 70.6 |
| phytoplankton_natural_unit | species | family |  128 | 56 % | 91 % | 98 % | -0.14 | 0.90 | -0.28 | 1.18 | 48.4 |
| phytoplankton_natural_unit | species | genus | 1434 | 58 % | 94 % | 97 % | -0.18 | 0.85 | -0.34 | 1.04 | 37.3 |
| phytoplankton_natural_unit | unit | kingdom |    2 | 100 % | 100 % | 100 % | -0.03 | 0.28 | -0.07 | 0.50 | 149.2 |
| phytoplankton_natural_unit | unit | phylum |    3 | 67 % | 100 % | 100 % | -0.60 | 0.37 | -1.53 | 1.30 | 146.8 |
| phytoplankton_natural_unit | unit | class |   19 | 47 % | 89 % | 100 % | -0.19 | 0.98 | -0.42 | 1.58 | 81.7 |
| phytoplankton_natural_unit | unit | order |   33 | 52 % | 91 % | 94 % | -0.22 | 0.96 | -0.48 | 1.42 | 72.4 |
| phytoplankton_natural_unit | unit | family |  116 | 59 % | 92 % | 97 % | -0.15 | 0.91 | -0.30 | 1.17 | 49.3 |
| phytoplankton_natural_unit | unit | genus | 1191 | 57 % | 94 % | 98 % | -0.26 | 0.82 | -0.47 | 1.06 | 37.8 |
| phytoplankton_natural_unit | unit | species |  261 | 63 % | 92 % | 95 % | 0.14 | 0.91 | 0.25 | 0.89 | 33.2 |

## Lectura

- `unit` y `species` miden lo que el modelo promete: una población nueva de una especie conocida y una especie nueva de un taxón superior conocido. La cobertura del 95 % debe quedar entre el 90 y el 97 %.
- `source` es más exigente de lo que será el uso real: con dos fuentes (fitoplancton, anfibios, peces), al apartar una no queda varianza entre fuentes que estimar, de modo que los intervalos no la incluyen; en invertebrados las fuentes cubren grupos distintos (zooplancton de Dumont, bentos de Hrycik y Moreau) y apartar una obliga a predecir desde clases o filos con otra composición de tamaños. En el modelo completo la varianza entre fuentes sí entra en la predicción.

## Unidades sin predicción (385)

| model | scheme | fold | motivo | n |
|---|---|---|---|---|
| amphibian_lwr | source | santini2018 | ajuste imposible: number of levels of each grouping factor must be < number of observations (problems: equation_id) |   1 |
| fish_individual | source | arranz2022 | ajuste imposible: number of levels of each grouping factor must be < number of observations (problems: g_species) |   1 |
| fish_lwr | source | arranz2022 | tipo de longitud sin datos en el ajuste | 276 |
| fish_lwr | source | koning2024 | tipo de longitud sin datos en el ajuste |  35 |
| invertebrate_individual | source | culver1985 | ajuste imposible: nrow(fx) == 1L is not TRUE |   1 |
| invertebrate_individual | source | moreau2021 | estadio sin datos en el ajuste |  71 |
