# MorphoPlus

🌐 [Español](README_es.md) · **English** · [Português](README_pt.md)

---

**MorphoPlus** is a semi-automatic pipeline for multi-band morphological fitting of galaxies in S-PLUS cluster fields using [GALFITM](https://www.nottingham.ac.uk/astronomy/megamorph/). It downloads S-PLUS images, builds PSFs, creates segmentation masks, and generates GALFITM input files — all organized around a spatial grid that subdivides each S-PLUS field into manageable sub-images.

---

## Table of Contents

1. [Dependencies](#1-dependencies)
2. [Installation](#2-installation)
3. [Directory Structure](#3-directory-structure)
4. [Input Catalog](#4-input-catalog)
5. [Configuration](#5-configuration)
6. [Running the Pipeline](#6-running-the-pipeline)
7. [Pipeline Description](#7-pipeline-description)
8. [Output Files](#8-output-files)
9. [Recovery & Utility Scripts](#9-recovery--utility-scripts)
10. [Field Galaxy Sample (control)](#10-field-galaxy-sample-control)
11. [S-PLUS Data](#11-s-plus-data)
12. [Citation](#12-citation)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. Dependencies

### Python packages

```bash
pip install astropy numpy matplotlib splusdata
```

### External tools

| Tool | Version | Purpose |
|------|---------|---------|
| [SExtractor](https://github.com/astromatic/sextractor) | ≥ 2.25 | Source detection and segmentation maps |
| [PSFEx](https://github.com/astromatic/psfex) | ≥ 3.21 | PSF modelling |
| [GALFITM](https://www.nottingham.ac.uk/astronomy/megamorph/) | 1.4.4 | Multi-band Sérsic profile fitting |

---

## 2. Installation

### Step 1 — System libraries

```bash
sudo apt-get update
sudo apt-get install libatlas-base-dev libblas-dev liblapack-dev   # ATLAS 3.6
sudo apt-get install libfftw3-dev                                  # FFTW 3.3.8
sudo apt-get install libplplot-dev                                 # PLPlot 5.9
```

### Step 2 — SExtractor

```bash
git clone https://github.com/astromatic/sextractor.git
cd sextractor
sh autogen.sh
./configure
make -j
sudo make install
```

Verify the installation:
```bash
sex --version
```

> **Alternative (conda):**
> ```bash
> conda install -c conda-forge astromatic-source-extractor
> ```
> See: https://anaconda.org/conda-forge/astromatic-source-extractor

### Step 3 — PSFEx

```bash
git clone https://github.com/astromatic/psfex.git
cd psfex
sh autogen.sh
./configure
make -j
sudo make install
```

Verify the installation:
```bash
psfex --version
```

### Step 4 — GALFITM

Download the binary for your platform from:
https://www.nottingham.ac.uk/astronomy/megamorph/

Place the binary (`galfitm-1.4.4-linux-x86_64` or equivalent) in the **root MorphoPlus directory**.

---

## 3. Directory Structure

The expected layout after cloning is shown below. `ejecutable.py` creates all output directories automatically — **you only need to create `Catalogos/` and place your input catalog there** before the first run.

```
MorphoPlus/
├── Catalogos/                  # ← create this and place SPLUS_Table.csv inside
│   └── SPLUS_Table.csv         # your input catalog
├── Field_Img/                  # created automatically by ejecutable.py
│   ├── det/                    # detection images and segmentation maps
│   ├── mask/                   # binary mask FITS files
│   └── psf/                    # Moffat PSF FITS files per filter
├── Out_img/                    # created automatically by ejecutable.py
├── config.py                   # ← YOUR credentials (never commit this file)
├── ejecutable.py               # main pipeline driver — generates ejecutable.sh
├── Recortar.py
├── mascara.py
├── segmetation.py
├── table_generation.py
├── leer_header_output.py
├── Img_galxgal.py
├── psf_new.py
├── count_fits.py               # progress checker (galaxies done vs. missing)
├── generate_missing_fits.py    # recovery: regenerates .sh for unfinished GALFITM fits
├── generate_psf_from_field.py  # recovery: rebuilds a missing PSF from a field image
├── gauss_5_0_9x9.conv          # SExtractor convolution filter
├── sextopsfex.param            # SExtractor parameter file
└── galfitm-1.4.4-linux-x86_64  # GALFITM binary
```

```bash
mkdir Catalogos   # only this one is needed before running
```

---

## 4. Input Catalog

Place your input table at:
```
Catalogos/SPLUS_Table.csv
```

The catalog must contain the following columns:

| Column | Description |
|--------|-------------|
| `ID` | Unique source identifier |
| `Field` | S-PLUS field name (e.g. `SPLUS-S14-S01`) |
| `ra` / `dec` | Right ascension and declination (degrees) |
| `RA` / `DEC` | Same as above (used for cutout requests) |
| `X` / `Y` | Pixel coordinates within the S-PLUS field frame |
| `ELONGATION` | Semi-major / semi-minor axis ratio (A/B) |
| `THETA` | Position angle from SExtractor (degrees) |
| `FLUX_RADIUS_50` | Half-light radius in pixels |
| `FLUX_RADIUS_90` | 90% light radius in pixels |
| `r_auto` | AUTO magnitude in r band |
| `u_auto`, `g_auto`, `i_auto`, `z_auto` | AUTO magnitudes in broadband filters |
| `J0378_auto` … `J0861_auto` | AUTO magnitudes in S-PLUS narrow bands |

### Generating the catalog with the Jupyter notebook

The notebook `splus-table-morfoplus.ipynb` automates the catalog preparation from a raw S-PLUS parquet or CSV file. It handles column renaming, computes `ELONGATION` from A/B, and outputs the correctly formatted `SPLUS_Table.csv`.

If the catalog does not already contain `X`, `Y`, or `ID` columns, `ejecutable.py` will compute them automatically using the WCS from downloaded field frames.

### Downloading S-PLUS data

S-PLUS photometric catalogs and images can be accessed at **https://splus.cloud**.

---

## 5. Configuration

### Step 0 — Set your S-PLUS credentials

Open `config.py` and replace the placeholder values with your own:

```python
SPLUS_USERNAME = "your_username_here"   # ← your S-PLUS username
SPLUS_PASSWORD = "your_password_here"   # ← your S-PLUS password
```

Register at **https://splus.cloud** if you do not have an account yet.

> **Security note:** `config.py` is listed in `.gitignore` and will never be committed to version control. Never share it or paste its contents in issues or pull requests.

### Grid and download parameters

All grid parameters are set in **`ejecutable.py`** only. You should not need to edit any other file for a standard run.

| Parameter | Location | Description |
|-----------|----------|-------------|
| `size` | `ejecutable.py` | Side length of each sub-image cutout in pixels. Default: `550`. |
| `c` | `ejecutable.py` | List of grid centre positions (pixels). Default: 20 positions from 275 to 10725 in steps of 550, covering the full S-PLUS field. |
| `GALFITM_BIN` | `ejecutable.py` | Path/name of the GALFITM binary. Default: `./galfitm-1.4.4-linux-x86_64`. |
| `MAX_WORKERS` | `Recortar.py` | Number of parallel threads for downloading filters. Default: `6`. Reduce to 3–4 if the S-PLUS API returns rate-limit errors. |
| `DR` | `Recortar.py` | Preferred S-PLUS data release. Default: `"dr6"`. Fallbacks are tried automatically. |
| `MAX_RETRIES` | `Recortar.py` | Maximum download attempts per filter before giving up. Default: `5`. |

---

## 6. Running the Pipeline

### Step 1 — Prepare the execution script

```bash
python ejecutable.py
```

This script:
- Reads `Catalogos/SPLUS_Table.csv`
- Computes pixel coordinates (X, Y) via WCS if they are missing (using a small 15×15 px R-band stamp — much faster than downloading the full field)
- Precomputes the explicit list of GALFITM commands by iterating the spatial grid
- Generates `ejecutable.sh` with all the commands for the full pipeline run

### Step 2 — Give execution permissions and run

```bash
chmod +wrx ejecutable.sh
./ejecutable.sh
```

`ejecutable.sh` runs the following stages automatically, with logging to `morphoplus_run.log`:

1. **Stage 1** — Downloads all 12 S-PLUS filter images per field (parallel per filter), cuts sub-images on the spatial grid, and builds Moffat PSFs from FITS header keywords (`Recortar.py`)
2. **Stage 2** — Runs SExtractor scripts (`dopsfex_mask_*.sh`) to produce segmentation maps
3. **Stage 3** — Creates binary masks and GALFITM input files (`mascara.py`)
4. **Stage 4** — Runs GALFITM **sequentially** on each group sub-image using an explicit, precomputed command list. Failures of individual fits are logged but do not abort the pipeline (`set +e` is active around this stage)
5. **Stage 5** — Runs GALFITM on individual galaxies via `ejecutable_gal.sh` (only if it exists)
6. **Stage 6** — Reads GALFITM output headers and writes the final results catalog (`leer_header_output.py`)

> **Why sequential GALFITM?** Earlier versions of the pipeline tried to parallelize GALFITM calls, but this caused hangs. The current version uses an explicit sequential command list — slower, but robust.

### Step 3 — Read the GALFITM output

This step runs automatically at the end of `ejecutable.sh`, but it can also be run on its own:

```bash
python leer_header_output.py
```

This reads the FITS headers of all GALFITM output files and produces:
- `Catalogos/GalfitM_output.csv` — fitted parameters (position, Re, n, b/a, PA, magnitude) for all 12 bands, with uncertainties
- `Out_img/*.svg` — three-panel visualization images (input / model / residual) for each fitted sub-image

---

## 7. Pipeline Description

```
SPLUS_Table.csv
      │
      ▼
ejecutable.py ──────────────────────────── generates ejecutable.sh
      │                                    (with explicit GALFITM command list)
      ▼
./ejecutable.sh
  │
  ├─ [Stage 1] Recortar.py
  │     Downloads 12 filters (parallel), cuts grid, builds PSFs
  │     Generates: dopsfex_mask_{field}.sh
  │
  ├─ [Stage 2] dopsfex_mask_*.sh
  │     Runs SExtractor → segmentation maps (.seg.fits)
  │
  ├─ [Stage 3] mascara.py
  │     Binary masks + GALFITM input files + ejecutable_gal.sh
  │
  ├─ [Stage 4] galfitm galfit_*_*_*.input   (sequential, set +e)
  │     Multi-band Sérsic fit for group sub-images
  │
  ├─ [Stage 5] ejecutable_gal.sh   (only if it exists)
  │     GALFITM galaxy by galaxy
  │
  └─ [Stage 6] leer_header_output.py
        GalfitM_output.csv + SVG images
```

### Spatial grid

Because GALFITM cannot fit an entire S-PLUS field (≈ 10800 × 10800 px), the field is divided into overlapping sub-images centred on a regular grid. The default configuration uses 20 positions per axis (400 sub-images per field), each 550 × 550 px. Galaxies within 25 px of a sub-image border are excluded from that sub-image's fit to avoid edge effects.

### Mask convention

| Pixel value | Meaning |
|-------------|---------|
| `0` | Group galaxy — **not** masked, will be fitted by GALFITM |
| `1` | Other source — masked, excluded from the fit |

### Zero-points

Zero-points are retrieved spatially from the S-PLUS API (`conn.get_zp`) for each field, band, and sky position. This uses DR6 by default, which provides per-position ZP corrections that account for CCD illumination gradients.

### Individual galaxy mode

For galaxies in the `Catalogos/g_S.csv` list (detected as isolated or in sparse groups), `Img_galxgal.py` downloads 200 × 200 px stamps directly via `conn.stamp()` and fits each galaxy individually instead of using the grid approach.

---

## 8. Output Files

| File | Description |
|------|-------------|
| `Catalogos/SPLUS_Table.csv` | Input catalog (updated with X, Y, ID columns if missing) |
| `Catalogos/g_S.csv` | Sub-catalog of galaxies processed individually |
| `Catalogos/procesados.txt` | Log of fields already processed (used for resuming) |
| `Catalogos/procesados_gal.txt` | Log of individual galaxies already processed |
| `Catalogos/GalfitM_output.csv` | GALFITM fitted parameters for all galaxies and bands |
| `Field_Img/*.fits` | Sub-image cutouts (one file per position per filter) |
| `Field_Img/det/*.fits` | Detection images (G+R+Z sum) and segmentation maps |
| `Field_Img/mask/*.fits` | Binary mask images |
| `Field_Img/psf/*.fits` | Moffat PSF images |
| `galfit_*.input` | GALFITM input files (group sub-images) |
| `Gal_*.input` | GALFITM input files (individual galaxies) |
| `Galfitm_*.fits` | GALFITM output FITS blocks for groups (input / model / residual) |
| `Gal_*.fits` | GALFITM output FITS blocks for individual galaxies |
| `Out_img/*.svg` | Three-panel visualization images |
| `ejecutable.sh` | Auto-generated execution script |
| `dopsfex_mask_*.sh` | Auto-generated SExtractor scripts |
| `morphoplus_run.log` | Full log of the pipeline run (created by `ejecutable.sh`) |

---

## 9. Recovery & Utility Scripts

GALFITM runs on hundreds of sub-images and individual galaxies, and a single full run can take many hours. To support long runs and to recover gracefully from interruptions (power outages, killed processes, missing PSFs, ...), the pipeline ships three small utility scripts.

### 9.1 `count_fits.py` — progress checker

Counts how many GALFITM output `.fits` files exist on disk, compared to how many are expected from the catalog. Useful at any time to check the progress of the pipeline.

```bash
python count_fits.py
```

Sample output:

```
GALFITM .fits: 312/400
  Groups: 290/378
  Indiv:  22/22
⏳ 78.0%
```

If everything is complete:

```
GALFITM .fits: 400/400
  Groups: 378/378
  Indiv:  22/22
✅ COMPLETO
```

### 9.2 `generate_missing_fits.py` — resume after a crash

If GALFITM is interrupted mid-run (power outage, killed process, system reboot, ...), this script scans the working directory for `.input` files **without** the corresponding `.fits` output and writes a recovery script `missing_fits.sh` containing only the missing GALFITM commands.

```bash
python generate_missing_fits.py
chmod +x missing_fits.sh
./missing_fits.sh
```

The recovery script:
- Uses `set +e` so individual failures do not abort the run
- Prints a `[N/total] Done` progress line after each completed fit
- Can be re-run as many times as needed — already-generated `.fits` are skipped

If nothing is missing, the script reports it and exits without creating `missing_fits.sh`.

### 9.3 `generate_psf_from_field.py` — rebuild a missing PSF

If `Recortar.py` failed to build a PSF for a particular filter (e.g. because `FWHM` or `BETA` were missing from the FITS header on first download, or the file was deleted), this script searches for any S-PLUS field image of that filter on disk, reads `FWHM` and `BETA` from its header, and regenerates the missing Moffat PSF in `Field_Img/psf/`.

Usage:

```bash
python generate_psf_from_field.py FIELD FILTER
```

Examples:

```bash
python generate_psf_from_field.py STRIPE82-0001 F378
python generate_psf_from_field.py HYDRA-0001 R
python generate_psf_from_field.py SPLUS-n03s27 F395
```

The script:
- Looks for files matching `*_FILTER.fits` under the working directory and `Field_Img/`
- Tries several common header keywords (`FWHM`, `PSF_FWHM`, `SEEING`, ...; `BETA`, `MOFFAT_BETA`, ...)
- Writes the new PSF to `Field_Img/psf/psf_{FIELD}_{FILTER}.fits`

After regenerating the PSF, you can re-run the affected GALFITM fits with `generate_missing_fits.py`.

---

## 10. Field Galaxy Sample (control)

To process a **field galaxy sample** (control sample, not in clusters):

1. Use `imagenes_field.py` to download multi-band cutouts for each object.
2. Use `Img_galxgal.py` (`galporgal` function) to generate the GALFITM input files.

These scripts follow the same 12-filter structure but work on individual objects rather than on a grid.

---

## 11. S-PLUS Data

S-PLUS photometric catalogs and images can be accessed at **https://splus.cloud**.

---

## 12. Citation

If you use MorphoPlus in your research, please cite the following paper:

> Montaguth, G. P., O'Mill, A. L., Mendes de Oliveira, C., et al. (2026).
> *Galaxy Evolution in Compact Groups. III. Structural Analysis of Galaxies and Dynamical State of Non-isolated Compact Groups.*
> The Astrophysical Journal, 998(1), 91.
> https://doi.org/10.3847/1538-4357/ae2bdd

BibTeX entry:

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

## 13. Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| Pipeline interrupted mid-run (power cut, kill, reboot) | Any crash | Run `python generate_missing_fits.py` and then `./missing_fits.sh` to resume only the unfinished fits. See [Section 9.2](#92-generate_missing_fitspy--resume-after-a-crash). |
| Want to know how far the pipeline has progressed | — | Run `python count_fits.py` at any time. See [Section 9.1](#91-count_fitspy--progress-checker). |
| PSF not created for a filter | FWHM/BETA missing from FITS header on first download, or file deleted | Run `python generate_psf_from_field.py FIELD FILTER` to rebuild it from any other field image of the same filter. See [Section 9.3](#93-generate_psf_from_fieldpy--rebuild-a-missing-psf). |
| `mascara.py` crashes / memory collapses | FITS files not being closed in loop | Use the fixed `mascara_fixed.py` which uses `with fits.open()` and explicit `del` after each `CCDData` read |
| S-PLUS download hangs or fails | Network timeout or API rate limit | Reduce `MAX_WORKERS` in `Recortar.py` to 3–4; the pipeline retries automatically with exponential backoff |
| `sex: command not found` | SExtractor not installed or not in PATH | Follow installation steps in Section 2; try the conda alternative |
| GALFITM input has magnitudes > 30 | Photometry flag issue | The pipeline automatically replaces bad magnitudes (> 30) with the r-band value |
| `KeyError: 'X'` or `KeyError: 'Y'` | Pixel columns missing from catalog | Run `ejecutable.py` first; it computes X and Y automatically via WCS |
| Single GALFITM fit fails but pipeline continues | Convergence problem on one sub-image | Expected — Stage 4 runs with `set +e` so single failures don't abort the run. Use `count_fits.py` and `generate_missing_fits.py` to inspect and retry. |
