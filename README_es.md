# MorphoPlus

🌐 [English](README.md) · **Español** · [Português](README_pt.md)

---

**MorphoPlus** es un pipeline semi-automático para el ajuste morfológico multi-banda de galaxias en campos de cúmulos de S-PLUS usando [GALFITM](https://www.nottingham.ac.uk/astronomy/megamorph/). Descarga imágenes de S-PLUS, construye PSFs, crea máscaras de segmentación y genera archivos de entrada para GALFITM — todo organizado en torno a una grilla espacial que subdivide cada campo de S-PLUS en sub-imágenes manejables.

---

## Tabla de contenidos

1. [Dependencias](#1-dependencias)
2. [Instalación](#2-instalación)
3. [Estructura de directorios](#3-estructura-de-directorios)
4. [Catálogo de entrada](#4-catálogo-de-entrada)
5. [Configuración](#5-configuración)
6. [Ejecución del pipeline](#6-ejecución-del-pipeline)
7. [Descripción del pipeline](#7-descripción-del-pipeline)
8. [Archivos de salida](#8-archivos-de-salida)
9. [Muestra de galaxias de campo (control)](#9-muestra-de-galaxias-de-campo-control)
10. [Query SQL de S-PLUS](#10-query-sql-de-s-plus)
11. [Cita](#11-cita)
12. [Solución de problemas](#12-solución-de-problemas)

---

## 1. Dependencias

### Paquetes de Python

```bash
pip install astropy numpy matplotlib splusdata
```

### Herramientas externas

| Herramienta | Versión | Función |
|-------------|---------|---------|
| [SExtractor](https://github.com/astromatic/sextractor) | ≥ 2.25 | Detección de fuentes y mapas de segmentación |
| [PSFEx](https://github.com/astromatic/psfex) | ≥ 3.21 | Modelado de la PSF |
| [GALFITM](https://www.nottingham.ac.uk/astronomy/megamorph/) | 1.4.4 | Ajuste multi-banda de perfiles de Sérsic |

---

## 2. Instalación

### Paso 1 — Bibliotecas del sistema

```bash
sudo apt-get update
sudo apt-get install libatlas-base-dev libblas-dev liblapack-dev   # ATLAS 3.6
sudo apt-get install libfftw3-dev                                  # FFTW 3.3.8
sudo apt-get install libplplot-dev                                 # PLPlot 5.9
```

### Paso 2 — SExtractor

```bash
git clone https://github.com/astromatic/sextractor.git
cd sextractor
sh autogen.sh
./configure
make -j
sudo make install
```

Verificar la instalación:
```bash
sex --version
```

> **Alternativa (conda):**
> ```bash
> conda install -c conda-forge astromatic-source-extractor
> ```
> Más información: https://anaconda.org/conda-forge/astromatic-source-extractor

### Paso 3 — PSFEx

```bash
git clone https://github.com/astromatic/psfex.git
cd psfex
sh autogen.sh
./configure
make -j
sudo make install
```

Verificar la instalación:
```bash
psfex --version
```

### Paso 4 — GALFITM

Descargar el binario para tu plataforma desde:
https://www.nottingham.ac.uk/astronomy/megamorph/

Colocar el binario (`galfitm-1.4.4-linux-x86_64` o equivalente) en el **directorio raíz de MorphoPlus**.

---

## 3. Estructura de directorios

Antes de correr el pipeline, crear los siguientes directorios vacíos:

```
MorphoPlus/
├── Catalogos/                  # Catálogo de entrada y tablas de salida
├── Field_Img/                  # Imágenes descargadas y procesadas
│   ├── det/                    # Imágenes de detección y mapas de segmentación
│   ├── mask/                   # Máscaras binarias en formato FITS
│   └── psf/                    # PSFs Moffat en formato FITS por filtro
├── Out_img/                    # Imágenes de visualización de salida (SVG)
├── ejecutable.py
├── Recortar.py
├── mascara.py
├── segmetation.py
├── table_generation.py
├── leer_header_output.py
├── Img_galxgal.py
├── psf_new.py
├── gauss_5_0_9x9.conv          # Filtro de convolución para SExtractor
├── sextopsfex.param            # Archivo de parámetros de SExtractor
└── galfitm-1.4.4-linux-x86_64 # Binario de GALFITM
```

Comando rápido de configuración:
```bash
mkdir -p Catalogos Field_Img/det Field_Img/mask Field_Img/psf Out_img
```

---

## 4. Catálogo de entrada

Colocar la tabla de entrada en:
```
Catalogos/SPLUS_Table.csv
```

El catálogo debe estar **corregido por extinción** y contener las siguientes columnas:

| Columna | Descripción |
|---------|-------------|
| `ID` | Identificador único de la fuente |
| `Field` | Nombre del campo S-PLUS (ej. `SPLUS-S14-S01`) |
| `ra` / `dec` | Ascensión recta y declinación (grados) |
| `RA` / `DEC` | Igual que lo anterior (usado para recortes) |
| `X` / `Y` | Coordenadas en píxeles dentro del campo S-PLUS |
| `ELONGATION` | Relación entre ejes semi-mayor y semi-menor (A/B) |
| `THETA` | Ángulo de posición de SExtractor (grados) |
| `FLUX_RADIUS_50` | Radio de semi-luz en píxeles |
| `FLUX_RADIUS_90` | Radio que contiene el 90% de la luz en píxeles |
| `r_auto` | Magnitud AUTO en banda r |
| `u_auto`, `g_auto`, `i_auto`, `z_auto` | Magnitudes AUTO en filtros de banda ancha |
| `J0378_auto` … `J0861_auto` | Magnitudes AUTO en filtros angostos de S-PLUS |

### Generación del catálogo con el notebook de Jupyter

El notebook `splus-table-morfoplus.ipynb` automatiza la preparación del catálogo a partir de un archivo parquet o CSV bruto de S-PLUS. Maneja el renombramiento de columnas, calcula `ELONGATION` a partir de A/B y genera el `SPLUS_Table.csv` con el formato correcto.

Si el catálogo no contiene las columnas `X`, `Y` o `ID`, `ejecutable.py` las calculará automáticamente usando el WCS de los frames descargados.

### Query SQL para descargar fotometría de S-PLUS (DR4/DR5)

La siguiente query recupera las columnas morfológicas y fotométricas necesarias. Reemplazar `ra`, `dec` y `r_g` con las coordenadas del centro del cúmulo y el radio de búsqueda (en grados):

```sql
SELECT
  det.ID, det.Field, det.RA, det.DEC, det.X, det.Y,
  det.ELONGATION, det.ELLIPTICITY, det.THETA, det.A, det.B,
  det.FLUX_RADIUS_50, det.FLUX_RADIUS_90,
  u.u_petro, u.e_u_petro,
  J0378.J0378_petro, J0378.e_J0378_petro,
  J0395.J0395_petro, J0395.e_J0395_petro,
  J0410.J0410_petro, J0410.e_J0410_petro,
  J0430.J0430_petro, J0430.e_J0430_petro,
  g.g_petro, g.e_g_petro,
  J0515.J0515_petro, J0515.e_J0515_petro,
  r.r_petro, r.e_r_petro,
  J0660.J0660_petro, J0660.e_J0660_petro,
  i.i_petro, i.e_i_petro,
  J0861.J0861_petro, J0861.e_J0861_petro,
  z.z_petro, z.e_z_petro,
  pz.zml, pz.odds
FROM idr4_dual.idr4_detection_image AS det
  JOIN idr4_dual.idr4_dual_u      AS u      ON u.ID      = det.ID
  JOIN idr4_dual.idr4_dual_J0378  AS J0378  ON J0378.ID  = det.ID
  JOIN idr4_dual.idr4_dual_J0395  AS J0395  ON J0395.ID  = det.ID
  JOIN idr4_dual.idr4_dual_J0410  AS J0410  ON J0410.ID  = det.ID
  JOIN idr4_dual.idr4_dual_J0430  AS J0430  ON J0430.ID  = det.ID
  JOIN idr4_dual.idr4_dual_g      AS g      ON g.ID      = det.ID
  JOIN idr4_dual.idr4_dual_J0515  AS J0515  ON J0515.ID  = det.ID
  JOIN idr4_dual.idr4_dual_r      AS r      ON r.ID      = det.ID
  JOIN idr4_dual.idr4_dual_J0660  AS J0660  ON J0660.ID  = det.ID
  JOIN idr4_dual.idr4_dual_i      AS i      ON i.ID      = det.ID
  JOIN idr4_dual.idr4_dual_J0861  AS J0861  ON J0861.ID  = det.ID
  JOIN idr4_dual.idr4_dual_z      AS z      ON z.ID      = det.ID
  JOIN idr4_vacs.idr4_photoz      AS pz     ON pz.ID     = det.ID
WHERE 1 = CONTAINS(
  POINT('ICRS', det.RA, det.DEC),
  CIRCLE('ICRS', ra, dec, r_g)
)
```

---

## 5. Configuración

Todos los parámetros de la grilla se definen **únicamente en `ejecutable.py`**. No es necesario editar ningún otro archivo para una corrida estándar.

| Parámetro | Archivo | Descripción |
|-----------|---------|-------------|
| `size` | `ejecutable.py` | Lado de cada sub-imagen en píxeles. Por defecto: `550`. |
| `c` | `ejecutable.py` | Lista de centros de la grilla (píxeles). Por defecto: 20 posiciones de 275 a 10725 en pasos de 550. |
| `MAX_WORKERS` | `Recortar.py` | Hilos paralelos para la descarga de filtros. Por defecto: `6`. Reducir a 3–4 si la API de S-PLUS devuelve errores de rate limit. |
| `DR` | `Recortar.py` | Data release preferido de S-PLUS. Por defecto: `"dr6"`. Los fallbacks se prueban automáticamente. |
| `MAX_RETRIES` | `Recortar.py` | Intentos máximos de descarga por filtro. Por defecto: `5`. |

---

## 6. Ejecución del pipeline

### Paso 1 — Generar el script de ejecución

```bash
python ejecutable.py
```

Este script:
- Lee `Catalogos/SPLUS_Table.csv`
- Calcula las coordenadas en píxeles (X, Y) via WCS si no están presentes
- Genera `ejecutable.sh` con todos los comandos para la corrida completa

### Paso 2 — Ejecutar el pipeline completo

```bash
./ejecutable.sh
```

`ejecutable.sh` corre la siguiente secuencia de forma automática:

1. Descarga las 12 imágenes por campo en paralelo
2. Recorta sub-imágenes en la grilla espacial
3. Construye PSFs Moffat a partir del header FITS
4. Corre SExtractor para generar mapas de segmentación
5. Crea máscaras binarias (`mascara.py`)
6. Genera los archivos de entrada de GALFITM (`.input`)
7. Corre GALFITM en cada sub-imagen

### Paso 3 — Leer los resultados de GALFITM

Este paso se ejecuta automáticamente al final de `ejecutable.sh`, pero también puede correrse de forma independiente:

```bash
python leer_header_output.py
```

Genera:
- `Catalogos/GalfitM_output.csv` — parámetros ajustados (posición, Re, n, b/a, AP, magnitud) para las 12 bandas con incertidumbres
- `Out_img/*.svg` — imágenes de tres paneles (entrada / modelo / residuo)

---

## 7. Descripción del pipeline

```
SPLUS_Table.csv
      │
      ▼
ejecutable.py ──────────────────────────── genera ejecutable.sh
      │
      ▼
./ejecutable.sh
  │
  ├─ [Stage 1] Recortar.py
  │     Descarga 12 filtros (paralelo), recorta grilla, construye PSFs
  │     Genera: dopsfex_mask_{campo}.sh
  │
  ├─ [Stage 2] dopsfex_mask_*.sh
  │     Corre SExtractor → mapas de segmentación (.seg.fits)
  │
  ├─ [Stage 3] mascara.py
  │     Máscaras binarias + inputs de GALFITM + ejecutable_gal.sh
  │
  ├─ [Stage 4] galfitm galfit_*_*_*.input
  │     Ajuste multi-banda de Sérsic para sub-imágenes de grupos
  │
  ├─ [Stage 5] ejecutable_gal.sh  (solo si existe)
  │     GALFITM galaxia por galaxia
  │
  └─ [Stage 6] leer_header_output.py
        GalfitM_output.csv + imágenes SVG
```

### Grilla espacial

Dado que GALFITM no puede ajustar un campo completo de S-PLUS (≈ 10800 × 10800 px), el campo se divide en sub-imágenes centradas en una grilla regular. La configuración por defecto usa 20 posiciones por eje (400 sub-imágenes por campo), cada una de 550 × 550 px. Las galaxias dentro de 25 px del borde de una sub-imagen son excluidas de ese ajuste para evitar efectos de borde.

### Convención de máscaras

| Valor del píxel | Significado |
|-----------------|-------------|
| `0` | Galaxia del grupo — **no** enmascarada, GALFITM la ajustará |
| `1` | Otra fuente — enmascarada, excluida del ajuste |

### Puntos cero fotométricos

Los puntos cero se recuperan espacialmente desde la API de S-PLUS (`conn.get_zp`) para cada campo, banda y posición en el cielo. Usa DR6 por defecto, que provee correcciones de ZP por posición que tienen en cuenta los gradientes de iluminación del CCD.

### Modo galaxia individual

Para galaxias en la lista `Catalogos/g_S.csv` (detectadas como aisladas o en grupos dispersos), `Img_galxgal.py` descarga stamps de 200 × 200 px directamente via `conn.stamp()` y ajusta cada galaxia individualmente en lugar de usar la grilla.

---

## 8. Archivos de salida

| Archivo | Descripción |
|---------|-------------|
| `Catalogos/SPLUS_Table.csv` | Catálogo de entrada (actualizado con X, Y, ID si faltaban) |
| `Catalogos/g_S.csv` | Sub-catálogo de galaxias procesadas individualmente |
| `Catalogos/procesados.txt` | Log de campos ya procesados (para retomar corridas) |
| `Catalogos/procesados_gal.txt` | Log de galaxias individuales ya procesadas |
| `Catalogos/GalfitM_output.csv` | Parámetros ajustados por GALFITM para todas las galaxias y bandas |
| `Field_Img/*.fits` | Recortes de sub-imágenes (un archivo por posición y filtro) |
| `Field_Img/det/*.fits` | Imágenes de detección (suma G+R+Z) y mapas de segmentación |
| `Field_Img/mask/*.fits` | Máscaras binarias |
| `Field_Img/psf/*.fits` | PSFs Moffat |
| `galfit_*.input` | Archivos de entrada de GALFITM |
| `Galfitm_*.fits` | FITS de salida de GALFITM (entrada / modelo / residuo) |
| `Out_img/*.svg` | Imágenes de visualización de tres paneles |
| `ejecutable.sh` | Script de ejecución generado automáticamente |
| `dopsfex_mask_*.sh` | Scripts de SExtractor generados automáticamente |

---

## 9. Muestra de galaxias de campo (control)

Para procesar una **muestra de galaxias de campo** (muestra de control, no en cúmulos):

1. Usar `imagenes_field.py` para descargar los recortes multi-banda de cada objeto.
2. Usar `Img_galxgal.py` (función `galporgal`) para generar los archivos de entrada de GALFITM.

Estos scripts siguen la misma estructura de 12 filtros pero trabajan sobre objetos individuales en lugar de sobre la grilla.

---

## 10. Query SQL de S-PLUS

Ver la [Sección 4](#4-catálogo-de-entrada) para la query SQL completa para descargar catálogos fotométricos de la base de datos de S-PLUS (DR4/DR5). La query devuelve magnitudes petrosianas en las 12 bandas junto con parámetros morfológicos y redshifts fotométricos.

---


## 11. Cita

Si usas MorphoPlus en tu investigación, por favor cita el siguiente artículo:

> Montaguth, G. P., O'Mill, A. L., Mendes de Oliveira, C., et al. (2026).
> *Galaxy Evolution in Compact Groups. III. Structural Analysis of Galaxies and Dynamical State of Non-isolated Compact Groups.*
> The Astrophysical Journal, 998(1), 91.
> https://doi.org/10.3847/1538-4357/ae2bdd

Entrada BibTeX:

```bibtex
@article{Montaguth_2026,
  doi       = {10.3847/1538-4357/ae2bdd},
  url       = {https://doi.org/10.3847/1538-4357/ae2bdd},
  year      = {2026},
  month     = {feb},
  publisher = {The American Astronomical Society},
  volume    = {998},
  number    = {1},
  pages     = {91},
  author    = {Montaguth, Gissel P. and O'Mill, Ana Laura and
               Mendes de Oliveira, Claudia and Lima-Dias, Ciria and
               Torres-Flores, Sergio and Monachesi, Antonela and
               Olave-Rojas, D. E. and Pallero, Diego and
               Humire, Pedro K. and Demarco, Ricardo and
               Telles, Eduardo and Lopes, Paulo A. A. and
               Panda, Swayamtrupta and Haack, Rodrigo F. and
               Lopes, Amanda R. and Alvarez-Candal, Alvaro and
               Smith Castelli, Analia V. and Kanaan, Antonio and
               Ribeiro, Tiago and Schoenell, William},
  title     = {Galaxy Evolution in Compact Groups. {III}. Structural Analysis
               of Galaxies and Dynamical State of Non-isolated Compact Groups},
  journal   = {The Astrophysical Journal}
}
```

---

## 12. Solución de problemas

| Síntoma | Causa probable | Solución |
|---------|---------------|----------|
| `mascara.py` colapsa / la memoria se llena | Archivos FITS no se cierran en el loop | Usar `mascara_fixed.py` que emplea `with fits.open()` y `del` explícito tras cada lectura de `CCDData` |
| Las descargas de S-PLUS se cuelgan o fallan | Timeout de red o rate limit de la API | Reducir `MAX_WORKERS` en `Recortar.py` a 3–4; el pipeline reintenta automáticamente con backoff exponencial |
| `sex: command not found` | SExtractor no instalado o no en el PATH | Seguir los pasos de instalación de la Sección 2; probar la alternativa conda |
| El input de GALFITM tiene magnitudes > 30 | Problema con flags de fotometría | El pipeline reemplaza automáticamente las magnitudes malas (> 30) con el valor de la banda r |
| El pipeline se interrumpe a mitad de la corrida | Cualquier error | Volver a correr `python ejecutable.py` y `./ejecutable.sh`; los campos y galaxias ya procesados se saltean via los logs en `Catalogos/` |
| `KeyError: 'X'` o `KeyError: 'Y'` | Columnas de píxeles faltantes en el catálogo | Correr `ejecutable.py` primero; calcula X e Y automáticamente via WCS |
| PSF no creada para un filtro | FWHM/BETA no encontrados en el header FITS | Se imprime una advertencia; GALFITM correrá sin esa PSF — verificar la versión del DR |
