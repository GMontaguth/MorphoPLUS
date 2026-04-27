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
9. [Scripts de recuperación y utilidades](#9-scripts-de-recuperación-y-utilidades)
10. [Muestra de galaxias de campo (control)](#10-muestra-de-galaxias-de-campo-control)
11. [Datos de S-PLUS](#11-datos-de-s-plus)
12. [Cita](#12-cita)
13. [Solución de problemas](#13-solución-de-problemas)

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

La estructura del proyecto se muestra a continuación. `ejecutable.py` crea todos los directorios de salida automáticamente — **solo es necesario crear `Catalogos/` y colocar el catálogo de entrada antes de la primera corrida**.

```
MorphoPlus/
├── Catalogos/                  # ← crear este y colocar SPLUS_Table.csv adentro
│   └── SPLUS_Table.csv         # tu catálogo de entrada
├── Field_Img/                  # creado automáticamente por ejecutable.py
│   ├── det/                    # imágenes de detección y mapas de segmentación
│   ├── mask/                   # máscaras binarias en formato FITS
│   └── psf/                    # PSFs Moffat en formato FITS por filtro
├── Out_img/                    # creado automáticamente por ejecutable.py
├── config.py                   # ← TUS credenciales (nunca subir a git)
├── ejecutable.py               # driver principal — genera ejecutable.sh
├── Recortar.py
├── mascara.py
├── segmetation.py
├── table_generation.py
├── leer_header_output.py
├── Img_galxgal.py
├── psf_new.py
├── count_fits.py               # chequeo de progreso (galaxias listas vs. faltantes)
├── generate_missing_fits.py    # recuperación: regenera .sh con los GALFITM pendientes
├── generate_psf_from_field.py  # recuperación: rearma una PSF faltante desde una imagen de campo
├── gauss_5_0_9x9.conv          # Filtro de convolución para SExtractor
├── sextopsfex.param            # Archivo de parámetros de SExtractor
└── galfitm-1.4.4-linux-x86_64  # Binario de GALFITM
```

```bash
mkdir Catalogos   # solo este directorio es necesario antes de correr
```

---

## 4. Catálogo de entrada

Colocar la tabla de entrada en:
```
Catalogos/SPLUS_Table.csv
```

El catálogo debe contener las siguientes columnas:

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

### Descarga de datos de S-PLUS

Los catálogos fotométricos e imágenes de S-PLUS están disponibles en **https://splus.cloud**.

---

## 5. Configuración

### Paso 0 — Configurar tus credenciales de S-PLUS

Abre `config.py` y reemplaza los valores de ejemplo con los tuyos:

```python
SPLUS_USERNAME = "tu_usuario_aqui"     # ← tu nombre de usuario en S-PLUS
SPLUS_PASSWORD = "tu_contraseña_aqui"  # ← tu contraseña de S-PLUS
```

Registrate en **https://splus.cloud** si todavía no tenés cuenta.

> **Seguridad:** `config.py` está incluido en `.gitignore` y **nunca** se sube al repositorio. No lo compartas ni pegues su contenido en issues o pull requests.

### Parámetros de grilla y descarga

Todos los parámetros de la grilla se definen **únicamente en `ejecutable.py`**. No es necesario editar ningún otro archivo para una corrida estándar.

| Parámetro | Archivo | Descripción |
|-----------|---------|-------------|
| `size` | `ejecutable.py` | Lado de cada sub-imagen en píxeles. Por defecto: `550`. |
| `c` | `ejecutable.py` | Lista de centros de la grilla (píxeles). Por defecto: 20 posiciones de 275 a 10725 en pasos de 550. |
| `GALFITM_BIN` | `ejecutable.py` | Ruta/nombre del binario de GALFITM. Por defecto: `./galfitm-1.4.4-linux-x86_64`. |
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
- Calcula las coordenadas en píxeles (X, Y) via WCS si no están presentes (usando un stamp pequeño de 15×15 px en banda R — mucho más rápido que descargar el campo completo)
- Precalcula la lista explícita de comandos GALFITM iterando sobre la grilla espacial
- Genera `ejecutable.sh` con todos los comandos para la corrida completa

### Paso 2 — Ejecutar el pipeline completo

```bash
chmod +wrx ejecutable.sh
./ejecutable.sh
```

`ejecutable.sh` corre las siguientes etapas de forma automática, con logging a `morphoplus_run.log`:

1. **Etapa 1** — Descarga las 12 imágenes por campo en paralelo, recorta sub-imágenes en la grilla espacial y construye PSFs Moffat a partir del header FITS (`Recortar.py`)
2. **Etapa 2** — Corre los scripts de SExtractor (`dopsfex_mask_*.sh`) para generar mapas de segmentación
3. **Etapa 3** — Crea máscaras binarias y archivos de entrada de GALFITM (`mascara.py`)
4. **Etapa 4** — Corre GALFITM **secuencialmente** sobre cada sub-imagen de grupo usando una lista de comandos explícita y precalculada. Las fallas de ajustes individuales se loguean pero no abortan el pipeline (`set +e` está activo en esta etapa)
5. **Etapa 5** — Corre GALFITM galaxia por galaxia via `ejecutable_gal.sh` (solo si existe)
6. **Etapa 6** — Lee los headers de salida de GALFITM y escribe el catálogo final de resultados (`leer_header_output.py`)

> **¿Por qué GALFITM secuencial?** Versiones anteriores del pipeline intentaban paralelizar las llamadas a GALFITM, pero esto causaba que el proceso se colgara. La versión actual usa una lista de comandos explícita y secuencial — más lenta, pero robusta.

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
      │                                    (con lista explícita de comandos GALFITM)
      ▼
./ejecutable.sh
  │
  ├─ [Etapa 1] Recortar.py
  │     Descarga 12 filtros (paralelo), recorta grilla, construye PSFs
  │     Genera: dopsfex_mask_{campo}.sh
  │
  ├─ [Etapa 2] dopsfex_mask_*.sh
  │     Corre SExtractor → mapas de segmentación (.seg.fits)
  │
  ├─ [Etapa 3] mascara.py
  │     Máscaras binarias + inputs de GALFITM + ejecutable_gal.sh
  │
  ├─ [Etapa 4] galfitm galfit_*_*_*.input   (secuencial, set +e)
  │     Ajuste multi-banda de Sérsic para sub-imágenes de grupos
  │
  ├─ [Etapa 5] ejecutable_gal.sh   (solo si existe)
  │     GALFITM galaxia por galaxia
  │
  └─ [Etapa 6] leer_header_output.py
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
| `galfit_*.input` | Archivos de entrada de GALFITM (sub-imágenes de grupo) |
| `Gal_*.input` | Archivos de entrada de GALFITM (galaxias individuales) |
| `Galfitm_*.fits` | FITS de salida de GALFITM para grupos (entrada / modelo / residuo) |
| `Gal_*.fits` | FITS de salida de GALFITM para galaxias individuales |
| `Out_img/*.svg` | Imágenes de visualización de tres paneles |
| `ejecutable.sh` | Script de ejecución generado automáticamente |
| `dopsfex_mask_*.sh` | Scripts de SExtractor generados automáticamente |
| `morphoplus_run.log` | Log completo de la corrida (creado por `ejecutable.sh`) |

---

## 9. Scripts de recuperación y utilidades

GALFITM corre sobre cientos de sub-imágenes y galaxias individuales, y una corrida completa puede tomar varias horas. Para soportar corridas largas y recuperarse de interrupciones (cortes de luz, procesos matados, PSFs faltantes, ...), el pipeline incluye tres utilidades pequeñas.

### 9.1 `count_fits.py` — chequeo de progreso

Cuenta cuántos archivos `.fits` de salida de GALFITM existen en disco, comparado con cuántos se esperan según el catálogo. Útil en cualquier momento para chequear el avance del pipeline.

```bash
python count_fits.py
```

Salida típica:

```
GALFITM .fits: 312/400
  Grupos: 290/378
  Indiv:  22/22
⏳ 78.0%
```

Cuando todo está completo:

```
GALFITM .fits: 400/400
  Grupos: 378/378
  Indiv:  22/22
✅ COMPLETO
```

### 9.2 `generate_missing_fits.py` — retomar tras un crash

Si GALFITM se interrumpe a mitad de la corrida (corte de luz, proceso matado, reinicio del sistema, ...), este script revisa el directorio de trabajo buscando archivos `.input` **sin** su correspondiente salida `.fits` y escribe un script de recuperación `missing_fits.sh` que contiene **solo** los comandos GALFITM faltantes.

```bash
python generate_missing_fits.py
chmod +x missing_fits.sh
./missing_fits.sh
```

El script de recuperación:
- Usa `set +e` para que las fallas individuales no aborten la corrida
- Imprime una línea de progreso `[N/total] Completado` después de cada ajuste terminado
- Puede correrse las veces que sea necesario — los `.fits` ya generados se saltean

Si no falta nada, el script lo informa y termina sin crear `missing_fits.sh`.

### 9.3 `generate_psf_from_field.py` — reconstruir una PSF faltante

Si `Recortar.py` no logró construir una PSF para un filtro en particular (por ejemplo porque `FWHM` o `BETA` faltaban en el header FITS en la primera descarga, o el archivo fue borrado), este script busca cualquier imagen S-PLUS de campo de ese filtro en disco, lee `FWHM` y `BETA` del header y regenera la PSF Moffat faltante en `Field_Img/psf/`.

Uso:

```bash
python generate_psf_from_field.py CAMPO FILTRO
```

Ejemplos:

```bash
python generate_psf_from_field.py STRIPE82-0001 F378
python generate_psf_from_field.py HYDRA-0001 R
python generate_psf_from_field.py SPLUS-n03s27 F395
```

El script:
- Busca archivos que coincidan con el patrón `*_FILTRO.fits` en el directorio de trabajo y bajo `Field_Img/`
- Prueba varios keywords comunes del header (`FWHM`, `PSF_FWHM`, `SEEING`, ...; `BETA`, `MOFFAT_BETA`, ...)
- Escribe la nueva PSF en `Field_Img/psf/psf_{CAMPO}_{FILTRO}.fits`

Después de regenerar la PSF, podés volver a correr los ajustes de GALFITM afectados con `generate_missing_fits.py`.

---

## 10. Muestra de galaxias de campo (control)

Para procesar una **muestra de galaxias de campo** (muestra de control, no en cúmulos):

1. Usar `imagenes_field.py` para descargar los recortes multi-banda de cada objeto.
2. Usar `Img_galxgal.py` (función `galporgal`) para generar los archivos de entrada de GALFITM.

Estos scripts siguen la misma estructura de 12 filtros pero trabajan sobre objetos individuales en lugar de sobre la grilla.

---

## 11. Datos de S-PLUS

Los catálogos fotométricos e imágenes de S-PLUS están disponibles en **https://splus.cloud**.

---

## 12. Cita

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

## 13. Solución de problemas

| Síntoma | Causa probable | Solución |
|---------|---------------|----------|
| Pipeline interrumpido a mitad de la corrida (corte de luz, kill, reboot) | Cualquier error | Correr `python generate_missing_fits.py` y luego `./missing_fits.sh` para retomar **solo** los ajustes pendientes. Ver [Sección 9.2](#92-generate_missing_fitspy--retomar-tras-un-crash). |
| Querés saber cuánto avanzó el pipeline | — | Correr `python count_fits.py` en cualquier momento. Ver [Sección 9.1](#91-count_fitspy--chequeo-de-progreso). |
| PSF no creada para un filtro | FWHM/BETA faltantes en el header FITS en la descarga inicial, o archivo borrado | Correr `python generate_psf_from_field.py CAMPO FILTRO` para reconstruirla desde otra imagen del mismo filtro. Ver [Sección 9.3](#93-generate_psf_from_fieldpy--reconstruir-una-psf-faltante). |
| `mascara.py` colapsa / la memoria se llena | Archivos FITS no se cierran en el loop | Usar `mascara_fixed.py` que emplea `with fits.open()` y `del` explícito tras cada lectura de `CCDData` |
| Las descargas de S-PLUS se cuelgan o fallan | Timeout de red o rate limit de la API | Reducir `MAX_WORKERS` en `Recortar.py` a 3–4; el pipeline reintenta automáticamente con backoff exponencial |
| `sex: command not found` | SExtractor no instalado o no en el PATH | Seguir los pasos de instalación de la Sección 2; probar la alternativa conda |
| El input de GALFITM tiene magnitudes > 30 | Problema con flags de fotometría | El pipeline reemplaza automáticamente las magnitudes malas (> 30) con el valor de la banda r |
| `KeyError: 'X'` o `KeyError: 'Y'` | Columnas de píxeles faltantes en el catálogo | Correr `ejecutable.py` primero; calcula X e Y automáticamente via WCS |
| Un ajuste de GALFITM falla pero el pipeline sigue | Problema de convergencia en una sub-imagen | Esperado — la Etapa 4 corre con `set +e` para que fallas individuales no aborten la corrida. Usar `count_fits.py` y `generate_missing_fits.py` para revisar y reintentar. |
