# -*- coding: utf-8 -*-
"""
ejecutable.py
=============
Entry point for MorphoPlus.

This version keeps all the improvements of the parallel rewrite
(stamp-based WCS, config.py, logging, safe-exec bash header, ...)
but runs GALFITM **sequentially** via an explicit command list
— same style as the original ejecutable.py — because parallel
GALFITM was hanging.

What this script does:
  1. Creates all required output directories.
  2. Reads Catalogos/SPLUS_Table.csv and adds ID, X, Y columns if missing
     (X/Y are computed via WCS from a tiny R-band stamp — fast).
  3. Generates ejecutable.sh — a fully autonomous bash script that runs
     the complete pipeline from image download to GALFITM output reading.

Pipeline stages inside ejecutable.sh:
  [1] python Recortar.py
  [2] dopsfex_mask_*.sh   (SExtractor scripts, discovered dynamically)
  [3] python mascara.py
  [4] GALFITM — group sub-images      <-- explicit sequential list
  [5] GALFITM — individual galaxies   <-- chmod + ./ejecutable_gal.sh
  [6] python leer_header_output.py
"""

import os
import numpy as np
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import splusdata

# Needed to precompute which galfit_*.input files mascara.py will create
from table_generation import tables

# ===================== CONFIGURATION =====================
# Spatial grid — edit ONLY here. All scripts import c and size from this file.
c    = [275, 825, 1375, 1925, 2475, 3025, 3575, 4125, 4675, 5225,
        5775, 6325, 6875, 7425, 7975, 8525, 9075, 9625, 10175, 10725]
size = 550

# GALFITM binary name — update if you use a different version or platform
GALFITM_BIN = "./galfitm-1.4.4-linux-x86_64"

# ===================== API CONNECTION =====================
from config import SPLUS_USERNAME, SPLUS_PASSWORD
conn = splusdata.Core(SPLUS_USERNAME, SPLUS_PASSWORD)

# ===================== CATALOG =====================
S       = Table.read('Catalogos/SPLUS_Table.csv')
Dados_S = S.group_by('Field')
Fields  = Dados_S.groups.keys


# ===================== HELPERS =====================

def _get_wcs_header(field, ra, dec):
    """
    Get the WCS header for a field by downloading a tiny 15x15 px R-band stamp.
    Orders of magnitude faster than downloading the full field frame
    just to extract a header. Tries DR6, DR5, DR4 in order.
    """
    for dr in ("dr6", "dr5", "dr4"):
        try:
            hdul = conn.stamp(ra, dec, size=15, band='R', data_release=dr)
            return hdul[1].header
        except Exception:
            continue
    raise RuntimeError(f"Could not get WCS header for field {field}")


def iauname(ra_deg, dec_deg):
    """Generate IAU-style IDs JHHMMSS.ss+DDMMSS.s from RA/DEC arrays."""
    coords = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
    return np.array([
        f"J{c_obj.to_string('hmsdms', precision=2, sep='', pad=True).replace(' ', '')}"
        for c_obj in coords
    ])


def add_missing_columns(S):
    """Add ID, X, Y columns to the catalog if they are not already present."""
    updated = False

    if 'ID' not in S.colnames:
        print("[INFO] Creating ID column from RA/DEC...")
        S['ID'] = iauname(np.array(S['ra']), np.array(S['dec']))
        updated = True
    else:
        print("[INFO] ID column already present.")

    if 'X' not in S.colnames or 'Y' not in S.colnames:
        print("[INFO] Computing pixel coordinates X, Y via WCS...")

        grouped     = S.group_by('Field')
        field_names = grouped.groups.keys['Field']

        X_all = np.full(len(S), np.nan)
        Y_all = np.full(len(S), np.nan)

        for field in field_names:
            mask_field = S['Field'] == field
            subtable   = S[mask_field]

            try:
                ra0    = float(subtable['ra'][0])
                dec0   = float(subtable['dec'][0])
                header = _get_wcs_header(field, ra0, dec0)
                wcs    = WCS(header)
                coords = SkyCoord(
                    ra  = np.array(subtable['ra'])  * u.deg,
                    dec = np.array(subtable['dec']) * u.deg
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

def _build_galfit_command_list():
    """
    Old-style: precompute the list of GALFITM commands by iterating the
    spatial grid and calling tables() for every (field, position).
    One explicit command per input file — no find loops, no parallelism.

    Returns
    -------
    group_cmds : list of str
        Pairs of ('chmod 777 <input>', '<galfitm> <input>') lines.
    has_individual : bool
        True if at least one position contains individual galaxies,
        i.e. ejecutable_gal.sh will be generated by mascara.py.
    """
    group_cmds     = []
    has_individual = False
    n_positions    = 0

    for f in Fields:
        field = f[0]
        for j in range(len(c)):
            for k in range(len(c)):
                position = (c[j], c[k])
                Tablef, Tabled = tables(S, field, position, size)

                if len(Tablef) > 0:
                    inp = f"galfit_{position[0]}_{position[1]}_{field}.input"
                    group_cmds.append(f"chmod 777 {inp}")
                    group_cmds.append(f"{GALFITM_BIN} {inp}")
                    n_positions += 1

                if len(Tabled) > 0:
                    has_individual = True

    print(f"[INFO] Predicted {n_positions} GALFITM group fits.")
    print(f"[INFO] Individual galaxies present: {has_individual}")
    return group_cmds, has_individual


def generate_ejecutable_sh():
    """
    Write ejecutable.sh — runs the full pipeline end-to-end.

    Stages 1-3 and 6 use dynamic file discovery (same as the parallel
    version). Stages 4 and 5 use an explicit command list built at
    ejecutable.py time (old-style, sequential).

    `set -euo pipefail` is active for housekeeping stages, but is
    temporarily relaxed around GALFITM calls so that a single failed
    fit does not kill the whole pipeline.
    """
    galfit_group_cmds, has_individual_galaxies = _build_galfit_command_list()

    lines = [
        "#!/bin/bash",
        "# ==========================================================================",
        "# ejecutable.sh - MorphoPlus autonomous pipeline",
        "# Generated automatically by ejecutable.py - do not edit by hand.",
        "#",
        "# GALFITM stages run SEQUENTIALLY via an explicit command list.",
        "# ==========================================================================",
        "",
        "set -euo pipefail",
        "",
        "exec > >(tee -a morphoplus_run.log) 2>&1",
        "",
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

        # --- STAGE 1 -------------------------------------------------
        'log_stage 1 "Download images, cut grid, build PSFs, run segmentation"',
        "python Recortar.py",
        'echo "[OK] Recortar.py finished."',
        "",

        # --- STAGE 2 -------------------------------------------------
        'log_stage 2 "Run SExtractor on group sub-images (dopsfex_mask_*.sh)"',
        "NSEX=0",
        'for sh_file in dopsfex_mask_*.sh; do',
        '    [ -f "$sh_file" ] || continue',
        '    echo "[SEX] Running $sh_file ..."',
        '    chmod +x "$sh_file"',
        '    bash "$sh_file"',
        '    NSEX=$((NSEX + 1))',
        "done",
        'echo "[OK] SExtractor done - ran $NSEX script(s)."',
        "",

        # --- STAGE 3 -------------------------------------------------
        'log_stage 3 "Build masks and generate GALFITM input files (mascara.py)"',
        "python mascara.py",
        'echo "[OK] mascara.py finished."',
        "",

        # --- STAGE 4 : EXPLICIT SEQUENTIAL GALFITM LIST --------------
        'log_stage 4 "Run GALFITM on group sub-images (sequential, explicit list)"',
        '# One explicit chmod + galfitm call per predicted input file.',
        '# Any single failure is logged but does not abort the pipeline.',
        "set +e",
    ]

    # Inject the precomputed old-style command list
    lines.extend(galfit_group_cmds)

    lines.extend([
        "set -e",
        'echo "[OK] GALFITM group stage done."',
        "",

        # --- STAGE 5 : individual galaxies, old style ---------------
        'log_stage 5 "Run GALFITM on individual galaxies (ejecutable_gal.sh)"',
    ])

    if has_individual_galaxies:
        lines.extend([
            'if [ -f "ejecutable_gal.sh" ]; then',
            '    chmod 777 ejecutable_gal.sh',
            '    set +e',
            '    ./ejecutable_gal.sh',
            '    set -e',
            '    echo "[OK] Individual galaxies done."',
            'else',
            '    echo "[WARN] ejecutable_gal.sh was predicted but not found."',
            'fi',
            "",
        ])
    else:
        lines.extend([
            'echo "[SKIP] No individual galaxies predicted by tables()."',
            "",
        ])

    # --- STAGE 6 ----------------------------------------------------
    lines.extend([
        'log_stage 6 "Read GALFITM outputs and write result catalog"',
        "python leer_header_output.py",
        'echo "[OK] leer_header_output.py finished."',
        "",
        'log_stage 7 "Pipeline complete"',
        'echo "Results: Catalogos/GalfitM_output.csv"',
        'echo "Images:  Out_img/"',
        'echo "Log:     morphoplus_run.log"',
    ])

    with open("ejecutable.sh", "w") as fic:
        for line in lines:
            fic.write(line + "\n")

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
