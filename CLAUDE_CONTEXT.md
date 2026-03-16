# MorphoPlus — Context for Claude

## What is this project?
MorphoPlus is a Python pipeline for semi-automatic multi-band morphological fitting
of galaxies in S-PLUS galaxy cluster fields using GALFITM (MegaMorph).
It was written by Gissel Pardo. The code is shared with collaborators, so all
comments and print statements must be in English.

## Language rules
- ALL code comments, docstrings, and print() messages must be in English.
- Variable names and function names can remain as-is (many are in Spanish/Portuguese
  for historical reasons — do not rename them unless asked).

## File map
| File | Role |
|------|------|
| `ejecutable.py` | Orchestrator: reads catalog, computes WCS pixel coords, generates ejecutable.sh |
| `Recortar.py` | Downloads 12-filter field frames (parallel), cuts spatial grid, builds PSFs |
| `segmetation.py` | Creates G+R+Z detection images, generates SExtractor scripts |
| `mascara.py` | Builds binary masks from segmentation maps; generates GALFITM .input files |
| `Img_galxgal.py` | Downloads stamps and generates GALFITM inputs for individual galaxies |
| `psf_new.py` | Builds a Moffat PSF FITS from FWHM and beta |
| `table_generation.py` | `tables()` function: filters catalog by field and grid position |
| `leer_header_output.py` | Reads GALFITM output headers, writes result CSV and SVG images |
| `imagenes_field.py` | Downloads cutouts for field (control) galaxy sample |
| `splus-table-morfoplus.ipynb` | Prepares SPLUS_Table.csv from raw parquet/CSV |

## Key global parameters (defined in ejecutable.py)
- `size = 550` — sub-image cutout side length in pixels
- `c` — list of 20 grid center positions (275 to 10725, step 550)
- These are imported by mascara.py, Recortar.py, segmetation.py via `from ejecutable import c, size`

## S-PLUS specifics
- 12 filters: U, F378, F395, F410, F430, G, F515, R, F660, I, F861, Z
- GALFITM band labels: u, J0378, J0395, J0410, J0430, g, J0515, r, J0660, i, J0861, z
- Plate scale: 0.55 arcsec/pixel
- Field size: ~10800 × 10800 px
- Data release in use: DR6 (fallback to DR4/DR5 automatically)
- API: splusdata.Core('gpardo', 'gNGC5054')
- Zero-points: retrieved spatially via conn.get_zp(field, band, ra, dec)

## Mask convention
- 0 = group galaxy pixel → NOT masked → GALFITM will fit it
- 1 = other source pixel → masked → excluded from fit

## Known fixed bugs (already in the optimized versions)
1. Recortar.py: 12 filters now downloaded in parallel (ThreadPoolExecutor, MAX_WORKERS=6)
   with skip-if-exists and exponential backoff retry.
2. Img_galxgal.py: 12 stamps per galaxy downloaded in parallel.
3. mascara.py: all fits.open() calls use 'with' context managers + .copy() to prevent
   memory leaks. CCDData objects are deleted after each read. gc.collect() called
   after every grid position.

## Optimized file names
- Recortar_optimized.py  (replaces Recortar.py)
- Img_galxgal_optimized.py  (replaces Img_galxgal.py)
- mascara_fixed.py  (replaces mascara.py)

## Important patterns to preserve
- `tables(S, field, position, size)` returns (Tablef, Tabled):
    - Tablef: galaxies in the interior of the cutout (to be fitted)
    - Tabled: galaxies near the border (excluded from fit, used for mask only)
- GALFITM input parameter line format: `N) val,val,...,val flag` (12 comma-separated values)
- Sky component flag: 0 for edge grid positions, 1 for interior positions
- Bad magnitude fallback: if mag > 30 in any band, replace with r_auto

## Things NOT to change without asking
- Credentials in splusdata.Core() calls
- The `c` and `size` values (pipeline-wide constants)
- GALFITM binary name: galfitm-1.4.4-linux-x86_64
- Output directory names: Field_Img/, Field_Img/det/, Field_Img/mask/, Field_Img/psf/, Out_img/
- Catalog paths: Catalogos/SPLUS_Table.csv, Catalogos/g_S.csv
