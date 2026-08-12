# -*- coding: utf-8 -*-
"""
leer_header_output.py
=====================
Reads GALFITM output files and writes the final result catalog + SVG images.

Two types of output files are processed:
  1. Group sub-images  : Galfitm_{cx}_{cy}_{field}.fits
  2. Individual galaxies: Gal_{ID}.fits

RESUME SUPPORT
--------------
If `Catalogos/GalfitM_output.csv` already exists, the galaxies listed there
(by column `ID`) are skipped — both header reading and SVG generation.
The script also skips any individual SVG that already exists in `Out_img/`,
so it can be safely re-run after a crash.

The CSV is re-written after every FITS file is processed, so partial progress
is preserved if the run is interrupted.

Usage:
    python leer_header_output.py
"""

from astropy.table import Table
from astropy.io import ascii, fits
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import colors

#from ejecutable import c, size, S
from table_generation import tables
c    = [275, 825, 1375, 1925, 2475, 3025, 3575, 4125, 4675, 5225,
        5775, 6325, 6875, 7425, 7975, 8525, 9075, 9625, 10175, 10725]
size = 550
S       = Table.read('Catalogos/SPLUS_Table.csv')
Dados_S = S.group_by('Field')
Fields  = Dados_S.groups.keys


# =============================================================================
# CONFIG
# =============================================================================
OUT_CSV     = "Catalogos/GalfitM_output.csv"
OUT_IMG_DIR = "Out_img"
MAKE_IMAGES = True

os.makedirs(OUT_IMG_DIR, exist_ok=True)
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

# Filter arrays — order must match GALFITM band order
Fil_name = np.array(['R','J0378','J0395','J0410','J0430','J0515',
                     'J0660','J0861','G','I','Z','U'])
Bands    = np.array(['U','F378','F395','F410','F430','G',
                     'F515','R','F660','I','F861','Z'])

# =============================================================================
# HELPERS
# =============================================================================

def correct_image_orientation(image):
    """Flip image vertically to match display convention."""
    return np.flipud(image)


def get_chi2nu(band_path):
    """Read Chi^2/nu from a .galfit.01.band file. Returns nan if not found."""
    try:
        with open(band_path, "r") as f:
            for line in f:
                if "Chi^2/nu =" in line:
                    return float(line.split("=")[1].split(",")[0].strip())
    except (FileNotFoundError, OSError):
        pass
    return np.nan


def leer_header(NGAL, HEADER, ID, Field_ID, data):
    """
    Append GALFITM Sersic parameters for all 12 bands from FITS header.
    NGAL: galaxy index within the group (1-based).
    """
    data.append(ID)
    data.append(Field_ID)
    for i in range(12):
        xc = HEADER[f"{NGAL}_XC_{Fil_name[i]}"].split()
        yc = HEADER[f"{NGAL}_YC_{Fil_name[i]}"].split()
        R  = HEADER[f"{NGAL}_RE_{Fil_name[i]}"].split()
        m  = HEADER[f"{NGAL}_MAG_{Fil_name[i]}"].split()
        N  = HEADER[f"{NGAL}_n_{Fil_name[i]}"].split()
        ar = HEADER[f"{NGAL}_AR_{Fil_name[i]}"].split()
        PA = HEADER[f"{NGAL}_PA_{Fil_name[i]}"].split()
        data.extend([xc[0], xc[2],
                     yc[0], yc[2],
                     R[0],  R[2],
                     m[0],  m[2],
                     N[0],  N[2],
                     ar[0], ar[2],
                     PA[0], PA[2]])
    return data


def build_header_names():
    """Column names for the output CSV."""
    names = ["CHI2NU", "ID", "Field_ID"]
    for i in range(12):
        names += [
            f"XC_{Fil_name[i]}",  f"e_XC_{Fil_name[i]}",
            f"YC_{Fil_name[i]}",  f"e_YC_{Fil_name[i]}",
            f"RE_{Fil_name[i]}",  f"e_RE_{Fil_name[i]}",
            f"MAG_{Fil_name[i]}", f"e_MAG_{Fil_name[i]}",
            f"n_{Fil_name[i]}",   f"e_n_{Fil_name[i]}",
            f"AR_{Fil_name[i]}",  f"e_AR_{Fil_name[i]}",
            f"PA_{Fil_name[i]}",  f"e_PA_{Fil_name[i]}",
        ]
    return names


def save_csv(rows, names, path):
    """Write the full results table to disk. Safe to call repeatedly."""
    if not rows:
        return
    tbl = Table(rows=rows, names=names)
    ascii.write(tbl, path, format='csv', fast_writer=False, overwrite=True)


def _imshow_band(ax, data2d, j):
    """Display one band image with SymLogNorm."""
    lt = [0.003,0.003,0.003,0.003,0.003,0.04,0.04,0.01,0.004,0.01,0.08,0.05]
    ls = [0.003,0.003,0.003,0.003,0.003,0.04,0.04,0.01,0.004,0.01,0.08,0.05]
    ax.imshow(data2d, cmap='gray',
              norm=colors.SymLogNorm(linthresh=lt[j], linscale=ls[j],
                                     vmin=-1.0, vmax=5.0))
    ax.set_xticks([]); ax.set_yticks([])


def grafico_g(fi, position, field, Tablef, k):
    """
    Plot for a galaxy inside a group sub-image.
    Centres the cutout on the galaxy pixel coordinates from the catalog.
    """
    npix = 50
    fig, axs = plt.subplots(3, 12, figsize=(20, 8), sharex=True, sharey=True)

    for j in range(12):
        x = Tablef[k]['X']
        y = len(fi[j].data) - Tablef[k]['Y']
        y1 = max(y - npix, 0)
        x1 = max(x - npix, 0)

        fic   = correct_image_orientation(fi[j].data)
        fic2  = correct_image_orientation(fi[12 + j].data)
        fic3  = correct_image_orientation(fi[24 + j].data)

        for row_ax, img in zip([0, 1, 2], [fic, fic2, fic3]):
            _imshow_band(axs[row_ax, j],
                         img[int(y1):int(y+npix), int(x1):int(x+npix)], j)

        axs[0, j].set_title(f'Filter {Bands[j]}', fontsize=8)

    # Scale bar on first panel
    axs[0, 0].text(1, 11, '5.5"', fontsize=10, color='blue')
    axs[0, 0].plot([1, 11], [2, 2], 'b-', lw=3)

    out_svg = os.path.join(
        OUT_IMG_DIR,
        f"{int(position[0])}_{int(position[1])}_{field}_{Tablef['ID'][k]}_{k}.svg"
    )
    plt.tight_layout()
    plt.savefig(out_svg, format='svg', dpi=1200)
    plt.close()
    return out_svg


def grafico_i(fi, position, field, Tabled, r):
    """
    Plot for an individual (isolated) galaxy.
    Centres at pixel (100, 100) — the stamp centre.
    """
    npix = 50
    x, y = 100, 100
    fig, axs = plt.subplots(3, 12, figsize=(20, 8), sharex=True, sharey=True)

    for j in range(12):
        fic  = correct_image_orientation(fi[j].data)
        fic2 = correct_image_orientation(fi[12 + j].data)
        fic3 = correct_image_orientation(fi[24 + j].data)

        for row_ax, img in zip([0, 1, 2], [fic, fic2, fic3]):
            _imshow_band(axs[row_ax, j],
                         img[int(y-npix):int(y+npix), int(x-npix):int(x+npix)], j)

        axs[0, j].set_title(f'Filter {Bands[j]}', fontsize=8)

    out_svg = os.path.join(
        OUT_IMG_DIR,
        f"{int(position[0])}_{int(position[1])}_{field}_{Tabled['ID'][r]}.svg"
    )
    plt.tight_layout()
    plt.savefig(out_svg, format='svg', dpi=1200)
    plt.close()
    return out_svg


# =============================================================================
# MAIN
# =============================================================================

header_names = build_header_names()

# -------------------------------------------------------------------
# RESUME: load galaxies already processed so we can skip them.
# -------------------------------------------------------------------
processed_ids = set()
Tabla = []
if os.path.exists(OUT_CSV):
    try:
        existing = Table.read(OUT_CSV, format='csv')
        processed_ids = {str(x) for x in existing['ID']}
        # Keep the existing rows so the final CSV preserves them
        Tabla = [list(row) for row in existing]
        print(f"[RESUME] {len(processed_ids)} galaxies already in '{OUT_CSV}' "
              f"will be skipped.")
    except Exception as e:
        print(f"[WARN] Could not read existing CSV ({e}). Starting from scratch.")
        processed_ids = set()
        Tabla = []
else:
    print(f"[INFO] No existing '{OUT_CSV}' found — processing everything.")


S_cat   = Table.read('Catalogos/SPLUS_Table.csv')
Datos_S = S_cat.group_by('Field')
Fields  = Datos_S.groups.keys

for f in Fields:
    field = f[0]

    for j in range(len(c)):
        for k in range(len(c)):
            position  = (c[j], c[k])
            Tablef, Tabled = tables(S_cat, field, position, size)

            # ------------------------------------------------------------------
            # Type 1 — Group sub-image: Galfitm_{cx}_{cy}_{field}.fits
            # ------------------------------------------------------------------
            if len(Tablef) > 0:
                # Which galaxies in this group still need processing?
                pending_idx = [i for i in range(len(Tablef))
                               if str(Tablef['ID'][i]) not in processed_ids]

                if not pending_idx:
                    # Nothing new here — don't even open the FITS.
                    pass
                else:
                    fits_path = f"Galfitm_{position[0]}_{position[1]}_{field}.fits"
                    band_path = f"Galfitm_{position[0]}_{position[1]}_{field}.galfit.01.band"
                    fi = None
                    new_in_group = 0
                    try:
                        fi = fits.open(fits_path)
                        chi2_nu = get_chi2nu(band_path)

                        for i in pending_idx:
                            gid = Tablef['ID'][i]
                            data = [chi2_nu]
                            Tabla.append(
                                leer_header(
                                    i + 1,
                                    fi[12].header,
                                    gid,
                                    f"{position[0]}_{position[1]}_{field}",
                                    data
                                )
                            )
                            processed_ids.add(str(gid))
                            new_in_group += 1

                            if MAKE_IMAGES:
                                if len(fi) < 36:
                                    print(f"[SKIP IMG] {fits_path} — only {len(fi)} HDUs, need >= 36.")
                                else:
                                    expected_svg = os.path.join(
                                        OUT_IMG_DIR,
                                        f"{int(position[0])}_{int(position[1])}_{field}_{gid}_{i}.svg"
                                    )
                                    if os.path.exists(expected_svg):
                                        print(f"[SKIP IMG] {expected_svg} already exists")
                                    else:
                                        out_svg = grafico_g(fi, position, field, Tablef, i)
                                        print(f"[IMG] {out_svg}")

                    except (FileNotFoundError, OSError) as e:
                        print(f"[SKIP] {fits_path}: {e}")
                    finally:
                        if fi is not None:
                            fi.close()

                    # Incremental save after every group file
                    if new_in_group > 0:
                        save_csv(Tabla, header_names, OUT_CSV)

            # ------------------------------------------------------------------
            # Type 2 — Individual galaxies: Gal_{ID}.fits
            # ------------------------------------------------------------------
            if len(Tabled) > 0:
                for r in range(len(Tabled)):
                    gid = Tabled['ID'][r]

                    # Skip galaxies already in the catalog
                    if str(gid) in processed_ids:
                        continue

                    fits_path = f"Gal_{gid}.fits"
                    band_path = f"Gal_{gid}.galfit.01.band"
                    fi = None
                    added_row = False
                    print(f"[GAL] {fits_path}")
                    try:
                        fi = fits.open(fits_path)
                        chi2_nu = get_chi2nu(band_path)
                        data = [chi2_nu]
                        Tabla.append(
                            leer_header(
                                1,
                                fi[12].header,
                                gid,
                                f"{position[0]}_{position[1]}_{field}",
                                data
                            )
                        )
                        processed_ids.add(str(gid))
                        added_row = True

                        if MAKE_IMAGES:
                            if len(fi) < 36:
                                print(f"[SKIP IMG] {fits_path} — only {len(fi)} HDUs.")
                            else:
                                expected_svg = os.path.join(
                                    OUT_IMG_DIR,
                                    f"{int(position[0])}_{int(position[1])}_{field}_{gid}.svg"
                                )
                                if os.path.exists(expected_svg):
                                    print(f"[SKIP IMG] {expected_svg} already exists")
                                else:
                                    out_svg = grafico_i(fi, position, field, Tabled, r)
                                    print(f"[IMG] {out_svg}")

                    except (FileNotFoundError, OSError) as e:
                        print(f"[SKIP] {fits_path}: {e}")
                    finally:
                        if fi is not None:
                            fi.close()

                    # Incremental save after every individual galaxy
                    if added_row:
                        save_csv(Tabla, header_names, OUT_CSV)

# =============================================================================
# FINAL WRITE
# =============================================================================
save_csv(Tabla, header_names, OUT_CSV)
print(f"\n[Done] {len(Tabla)} total rows saved to '{OUT_CSV}'.")
