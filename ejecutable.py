# -*- coding: utf-8 -*-
"""
ejecutable.py
=============
Entry point for MorphoPlus.

What this script does:
  1. Creates all required output directories.
  2. Reads Catalogos/SPLUS_Table.csv and adds ID, X, Y columns if missing
     (X/Y are computed via WCS from the S-PLUS R-band frame).
  3. Generates ejecutable.sh — a fully autonomous bash script that runs
     the complete pipeline from image download to GALFITM output reading.

Usage:
    python ejecutable.py          # generates ejecutable.sh
    chmod +x ejecutable.sh
    ./ejecutable.sh               # runs everything end-to-end

Pipeline stages inside ejecutable.sh:
  [1] python Recortar.py
        Downloads all 12-filter field frames (parallel),
        cuts the spatial grid, builds PSFs,
        and generates dopsfex_mask_{field}.sh (SExtractor scripts).
  [2] SExtractor scripts (dopsfex_mask_*.sh)
        Run SExtractor on every detection image to produce segmentation maps.
  [3] python mascara.py
        Builds binary masks from segmentation maps and generates
        GALFITM input files (galfit_*.input) for every group position.
        Also generates ejecutable_gal.sh for individual galaxies.
  [4] GALFITM — group sub-images
        Runs galfitm on every galfit_*_*_*.input file found on disk.
  [5] GALFITM — individual galaxies  (only if ejecutable_gal.sh exists)
        Runs galfitm galaxy-by-galaxy.
  [6] python leer_header_output.py
        Reads GALFITM output headers and writes the final CSV + SVG images.

Note: stages [2] and [4] use shell 'find' loops so they discover files
dynamically — no file names are hard-coded in the generated script.
"""

import os
import numpy as np
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import splusdata

# ===================== CONFIGURATION =====================
# These two parameters control the spatial grid used by the whole pipeline.
# Edit ONLY here — all other scripts import c and size from this file.
c    = [275, 825, 1375, 1925, 2475, 3025, 3575, 4125, 4675, 5225,
        5775, 6325, 6875, 7425, 7975, 8525, 9075, 9625, 10175, 10725]
size = 550

# GALFITM binary name — update if you use a different version or platform
GALFITM_BIN = "./galfitm-1.4.4-linux-x86_64"

# ===================== API CONNECTION =====================
conn = splusdata.Core()

# ===================== CATALOG =====================
S       = Table.read('Catalogos/SPLUS_Table.csv')
Dados_S = S.group_by('Field')
Fields  = Dados_S.groups.keys


# ===================== HELPERS =====================

def _download_r_frame(field):
    """
    Download the R-band science frame for WCS computation.
    Tries DR5 first, then DR4 as fallback.
    Returns an HDUList or raises on total failure.
    """
    for dr in ("dr5", "dr4", "dr6"):
        try:
            return conn.field_frame(field, 'R', data_release=dr)
        except Exception:
            continue
    raise RuntimeError(f"Could not download R frame for field {field}")


def iauname(ra_deg, dec_deg):
    """
    Generate IAU-style source IDs of the form JHHMMSS.ss+DDMMSS.s
    from arrays of RA and DEC in degrees.
    """
    coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
    return np.array([
        f"J{c_obj.to_string('hmsdms', precision=2, sep='', pad=True).replace(' ', '')}"
        for c_obj in coords
    ])


def add_missing_columns(S):
    """
    Add ID, X, and Y columns to the catalog if they are not already present.

    ID  : IAU-style name from RA/DEC.
    X/Y : Pixel coordinates computed from the S-PLUS R-band WCS for each field.

    The updated catalog is written back to Catalogos/SPLUS_Table.csv.
    """
    updated = False

    # --- ID column ---
    if 'ID' not in S.colnames:
        print("[INFO] Creating ID column from RA/DEC...")
        S['ID'] = iauname(np.array(S['ra']), np.array(S['dec']))
        updated = True
    else:
        print("[INFO] ID column already present.")

    # --- X and Y columns ---
    if 'X' not in S.colnames or 'Y' not in S.colnames:
        print("[INFO] Computing pixel coordinates X, Y via WCS...")

        grouped = S.group_by('Field')
        field_names = grouped.groups.keys['Field']

        X_all = np.full(len(S), np.nan)
        Y_all = np.full(len(S), np.nan)

        for field in field_names:
            mask_field = S['Field'] == field
            subtable   = S[mask_field]

            try:
                hdul   = _download_r_frame(field)
                header = hdul[1].header
                wcs    = WCS(header)
                coords = SkyCoord(
                    ra=np.array(subtable['ra'])  * u.deg,
                    dec=np.array(subtable['dec']) * u.deg
                )
                x, y = wcs.world_to_pixel(coords)
                X_all[mask_field] = x
                Y_all[mask_field] = y
                print(f"  [OK] WCS computed for field {field}")
            except Exception as e:
                print(f"  [ERROR] Could not compute WCS for field {field}: {e}")

        S['X'] = X_all
        S['Y'] = Y_all
        updated = True
    else:
        print("[INFO] X and Y columns already present.")

    if updated:
        print("[INFO] Saving updated catalog to Catalogos/SPLUS_Table.csv ...")
        S.write('Catalogos/SPLUS_Table.csv', overwrite=True)
        print("[INFO] Catalog saved.")
    else:
        print("[INFO] No changes made to the catalog.")


# ===================== SCRIPT GENERATOR =====================

def generate_ejecutable_sh():
    """
    Write ejecutable.sh — a fully autonomous bash script that runs the
    complete MorphoPlus pipeline without any manual intervention.

    Design principles:
      - Uses 'set -euo pipefail' so the script stops on any error.
      - Each stage is wrapped in a timestamped log block.
      - SExtractor scripts and GALFITM inputs are discovered dynamically
        with 'find' loops — no filenames are hard-coded.
      - A log file (morphoplus_run.log) records stdout + stderr.
      - Exit codes are checked after each critical stage.
    """

    lines = [
        "#!/bin/bash",
        "# ==========================================================================",
        "# ejecutable.sh — MorphoPlus autonomous pipeline",
        "# Generated automatically by ejecutable.py — do not edit by hand.",
        "# ==========================================================================",
        "",
        "# Stop immediately if any command fails, an unset variable is used,",
        "# or any command in a pipeline fails.",
        "set -euo pipefail",
        "",
        "# Redirect all output (stdout + stderr) to a log file AND the terminal.",
        "exec > >(tee -a morphoplus_run.log) 2>&1",
        "",
        "# Helper: print a timestamped stage banner.",
        'log_stage() { echo ""; echo "========================================'
        '================================"; '
        'echo "[$(date \'+%Y-%m-%d %H:%M:%S\')] STAGE $1: $2"; '
        'echo "================================================================'
        '================"; }',
        "",
        'log_stage 0 "Pipeline start"',
        'echo "Working directory: $(pwd)"',
        'echo "Python: $(python --version)"',
        "",

        # ── STAGE 1: Recortar.py ──────────────────────────────────────────────
        'log_stage 1 "Download images, cut grid, build PSFs, run segmentation"',
        "python Recortar.py",
        'echo "[OK] Recortar.py finished."',
        "",

        # ── STAGE 2: SExtractor — group positions ────────────────────────────
        'log_stage 2 "Run SExtractor on group sub-images (dopsfex_mask_*.sh)"',
        "# dopsfex_mask_{field}.sh scripts are generated dynamically by",
        "# segmetation.py (called from Recortar.py). We discover them here.",
        "NSEX=0",
        'for sh_file in dopsfex_mask_*.sh; do',
        '    [ -f "$sh_file" ] || continue   # skip if glob matched nothing',
        '    echo "[SEX] Running $sh_file ..."',
        '    chmod +x "$sh_file"',
        '    bash "$sh_file"',
        '    NSEX=$((NSEX + 1))',
        "done",
        'echo "[OK] SExtractor done — ran $NSEX script(s)."',
        "",

        # ── STAGE 3: mascara.py ───────────────────────────────────────────────
        'log_stage 3 "Build masks and generate GALFITM input files (mascara.py)"',
        "python mascara.py",
        'echo "[OK] mascara.py finished."',
        "",

        # ── STAGE 4: GALFITM — group sub-images ──────────────────────────────
        'log_stage 4 "Run GALFITM on group sub-images (galfit_*_*_*.input)"',
        "# galfit_*.input files are generated by mascara.py.",
        "# We discover them dynamically so no position is hard-coded.",
        "NGALFIT=0",
        'for inp in galfit_*_*_*.input; do',
        '    [ -f "$inp" ] || continue',
        '    echo "[GALFITM] $inp ..."',
        f'    {GALFITM_BIN} "$inp"',
        '    NGALFIT=$((NGALFIT + 1))',
        "done",
        'echo "[OK] GALFITM group stage done — ran $NGALFIT fit(s)."',
        "",

        # ── STAGE 5: GALFITM — individual galaxies (optional) ────────────────
        'log_stage 5 "Run GALFITM on individual galaxies (ejecutable_gal.sh)"',
        "# ejecutable_gal.sh is generated by mascara.py / Img_galxgal.py",
        "# only when isolated galaxies are present. Skip if it doesn't exist.",
        'if [ -f "ejecutable_gal.sh" ]; then',
        '    chmod +x ejecutable_gal.sh',
        '    bash ejecutable_gal.sh',
        '    echo "[OK] Individual galaxy GALFITM stage done."',
        "else",
        '    echo "[SKIP] ejecutable_gal.sh not found — no individual galaxies to fit."',
        "fi",
        "",

        # ── STAGE 6: leer_header_output.py ───────────────────────────────────
        'log_stage 6 "Read GALFITM outputs and write result catalog"',
        "python leer_header_output.py",
        'echo "[OK] leer_header_output.py finished."',
        "",

        # ── Done ──────────────────────────────────────────────────────────────
        'log_stage 7 "Pipeline complete"',
        'echo "Results: Catalogos/GalfitM_output.csv"',
        'echo "Images:  Out_img/"',
        'echo "Log:     morphoplus_run.log"',
    ]

    with open("ejecutable.sh", "w") as fic:
        for line in lines:
            fic.write(line + "\n")

    # Make it executable right away so the user doesn't need chmod
    os.chmod("ejecutable.sh", 0o755)
    print("[INFO] ejecutable.sh written and marked as executable.")
    print("[INFO] Run with:  ./ejecutable.sh")


# ===================== MAIN =====================

if __name__ == "__main__":

    # --- Create all required directories ---
    for d in ("Field_Img", "Field_Img/det", "Field_Img/mask",
              "Field_Img/psf", "Out_img", "Catalogos"):
        os.makedirs(d, exist_ok=True)
    print("[INFO] Output directories ready.")

    # --- Ensure catalog has ID, X, Y ---
    add_missing_columns(S)

    # --- Generate the autonomous pipeline script ---
    generate_ejecutable_sh()
