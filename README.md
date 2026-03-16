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
9. [Field Galaxy Sample (control)](#9-field-galaxy-sample-control)
10. [S-PLUS SQL Query](#10-s-plus-sql-query)
11. [Citation](#11-citation)
12. [Troubleshooting](#12-troubleshooting)

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

Before running the pipeline, create the following empty directories:

```
MorphoPlus/
├── Catalogos/                  # Input catalog and output tables
├── Field_Img/                  # All downloaded and processed images
│   ├── det/                    # Detection images and segmentation maps
│   ├── mask/                   # Binary mask FITS files
│   └── psf/                    # Moffat PSF FITS files per filter
├── Out_img/                    # Output visualization images (SVG)
├── ejecutable.py
├── Recortar.py
├── mascara.py
├── segmetation.py
├── table_generation.py
├── leer_header_output.py
├── Img_galxgal.py
├── psf_new.py
├── gauss_5_0_9x9.conv          # SExtractor convolution filter
├── sextopsfex.param            # SExtractor parameter file
└── galfitm-1.4.4-linux-x86_64 # GALFITM binary
```

Quick setup command:
```bash
mkdir -p Catalogos Field_Img/det Field_Img/mask Field_Img/psf Out_img
```

---

## 4. Input Catalog

Place your input table at:
```
Catalogos/SPLUS_Table.csv
```

The catalog must be **extinction-corrected** and contain the following columns:

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

### SQL query to download S-PLUS photometry (DR4/DR5)

The following query retrieves the morphological and photometric columns needed. Replace `ra`, `dec`, and `r_g` with your cluster centre coordinates and search radius (in degrees):

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

## 5. Configuration

All grid parameters are set in **`ejecutable.py`** only. You should not need to edit any other file for a standard run.

| Parameter | Location | Description |
|-----------|----------|-------------|
| `size` | `ejecutable.py` | Side length of each sub-image cutout in pixels. Default: `550`. |
| `c` | `ejecutable.py` | List of grid centre positions (pixels). Default: 20 positions from 275 to 10725 in steps of 550, covering the full S-PLUS field. |
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
- Computes pixel coordinates (X, Y) via WCS if they are missing
- Generates `ejecutable.sh` with all the commands for the full pipeline run

### Step 2 — Give execution permissions and run

```bash
chmod +wrx ejecutable.sh
./ejecutable.sh
```

`ejecutable.sh` runs the following sequence automatically:

1. Downloads all 12 S-PLUS filter images per field (parallel per filter)
2. Cuts sub-images on the spatial grid
3. Builds Moffat PSFs from FITS header keywords
4. Runs SExtractor to produce segmentation maps
5. Creates binary masks (`mascara.py`)
6. Generates GALFITM input files (`.input`) for each sub-image
7. Runs GALFITM on each sub-image

### Step 3 — Read the GALFITM output

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
ejecutable.py ──────────────────────────────────── generates ejecutable.sh
      │
      ▼
Recortar.py
  ├─ Downloads 12-filter field frames (parallel, with retry + backoff)
  ├─ Builds Moffat PSFs from FITS header FWHM/BETA keywords
  ├─ Cuts 550×550 px sub-images on a 20×20 spatial grid
  └─ Calls segmetation.py ──► SExtractor segmentation maps
      │
      ▼
mascara.py
  ├─ Reads segmentation maps
  ├─ Zeros out group-galaxy segments → binary mask (0=galaxy, 1=mask)
  ├─ Estimates per-band sky background (sigma-clipped mean)
  └─ Writes GALFITM input files (.input) with spatially-varying ZPs (DR6)
      │
      ▼
galfitm-1.4.4-linux-x86_64 galfit_*.input
  └─ Multi-band Sérsic fit (12 bands simultaneously)
      │
      ▼
leer_header_output.py
  ├─ Reads GALFITM output FITS headers
  └─ Produces GalfitM_output.csv + SVG visualization images
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
| `galfit_*.input` | GALFITM input files |
| `Galfitm_*.fits` | GALFITM output FITS blocks (input / model / residual) |
| `Out_img/*.svg` | Three-panel visualization images |
| `ejecutable.sh` | Auto-generated execution script |
| `dopsfex_mask_*.sh` | Auto-generated SExtractor scripts |

---

## 9. Field Galaxy Sample (control)

To process a **field galaxy sample** (control sample, not in clusters):

1. Use `imagenes_field.py` to download multi-band cutouts for each object.
2. Use `Img_galxgal.py` (`galporgal` function) to generate the GALFITM input files.

These scripts follow the same 12-filter structure but work on individual objects rather than on a grid.

---

## 10. S-PLUS SQL Query

See [Section 4](#4-input-catalog) for the full SQL query to download photometric catalogs from the S-PLUS database (DR4/DR5). The query returns petrosian magnitudes in all 12 bands along with morphological parameters and photometric redshifts.

---


## 11. Citation

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

## 12. Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| `mascara.py` crashes / memory collapses | FITS files not being closed in loop | Use the fixed `mascara_fixed.py` which uses `with fits.open()` and explicit `del` after each `CCDData` read |
| S-PLUS download hangs or fails | Network timeout or API rate limit | Reduce `MAX_WORKERS` in `Recortar.py` to 3–4; the pipeline retries automatically with exponential backoff |
| `sex: command not found` | SExtractor not installed or not in PATH | Follow installation steps in Section 2; try the conda alternative |
| GALFITM input has magnitudes > 30 | Photometry flag issue | The pipeline automatically replaces bad magnitudes (> 30) with the r-band value |
| Pipeline interrupted mid-run | Any crash | Re-run `ejecutable.py` and `./ejecutable.sh`; already-processed fields and galaxies are skipped via the log files in `Catalogos/` |
| `KeyError: 'X'` or `KeyError: 'Y'` | Pixel columns missing from catalog | Run `ejecutable.py` first; it computes X and Y automatically via WCS |
| PSF not created for a filter | FWHM/BETA not found in FITS header | A warning is printed; GALFITM will run without that PSF — check your DR version |
