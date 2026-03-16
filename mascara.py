# -*- coding: utf-8 -*-
"""
mascara.py
==========
Builds binary masks from SExtractor segmentation maps and generates
GALFITM multi-band input files for each group position on the spatial grid.

Mask convention:
  0 = group galaxy pixel  (not masked — will be fitted by GALFITM)
  1 = other source pixel  (masked — excluded from the fit)

Memory management notes:
  All FITS files are opened with 'with' context managers so they are
  guaranteed to be closed (and their memory freed) after each use.
  CCDData objects in Median_sky() are explicitly deleted after reading
  to prevent accumulation across the 400-position loop per field.
  gc.collect() is called at the end of each grid position to release
  any cyclic references that the garbage collector might miss.

Pipeline position:
  segmetation.py (SExtractor) -> mascara.py -> GALFITM input files
"""

import gc
import csv
from math import isnan

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.nddata import CCDData

import splusdata

from ejecutable import c, size
from table_generation import tables

# ===================== API CONNECTION =====================
conn = splusdata.Core()

# ===================== INPUT CATALOG =====================
S = Table.read('Catalogos/SPLUS_Table.csv')

# S-PLUS filter names (order used by GALFITM internally)
FILTROS = np.array(['U', 'F378', 'F395', 'F410', 'F430',
                    'G', 'F515', 'R', 'F660', 'I', 'F861', 'Z'])


# ===================== MASK FUNCTIONS =====================

def mask_gal(galaxy):
    """
    Build a binary mask for an individual galaxy from its SExtractor
    segmentation map (det_{galaxy}.seg.fits).

    The pixel at the center of the stamp (100, 100) is assumed to belong
    to the target galaxy.  All pixels with the same segment ID are set to 0
    (unmasked); all other non-zero pixels are set to 1 (masked).

    Output: Field_Img/mask/mask_{galaxy}.fits
    """
    seg_path  = f'Field_Img/det/det_{galaxy}.seg.fits'
    mask_path = f'Field_Img/mask/mask_{galaxy}.fits'

    # Use a context manager so the file is closed as soon as we are done.
    # We copy the data array out of the HDU before leaving the 'with' block.
    with fits.open(seg_path) as f:
        seg_data = f[0].data.copy()   # copy: lets the HDUList be freed immediately
        hdr      = f[0].header.copy()

    # Zero out the target galaxy's segment (center pixel defines the segment ID)
    target_id              = seg_data[100, 100]
    seg_data[seg_data == target_id] = 0

    # Convert all remaining segments to 1 (masked)
    seg_data[seg_data > 0] = 1

    fits.PrimaryHDU(seg_data, header=hdr).writeto(mask_path, overwrite=True)

    # Explicitly release the arrays; not strictly necessary but good practice
    del seg_data, hdr


def mask(field, Table, position):
    """
    Build a binary mask for a group cutout from its SExtractor segmentation map
    (det_{x}_{y}_{field}.seg.fits).

    For each galaxy in Table, the segment that overlaps pixel (X[i], Y[i]) is
    zeroed out (not masked).  All remaining non-zero segments are set to 1.

    Parameters
    ----------
    field    : S-PLUS field name string
    Table    : astropy Table — group galaxies with pixel columns X, Y
    position : (x, y) tuple — grid position used in the filename
    """
    seg_path  = ('Field_Img/det/det_%i_%i_%s.seg.fits'
                 % (position[0], position[1], field))
    mask_path = ('Field_Img/mask/mask_%i_%i_%s.fits'
                 % (position[0], position[1], field))

    with fits.open(seg_path) as f:
        seg_data = f[0].data.copy()   # copy data out before closing
        hdr      = f[0].header.copy()

    # Zero out every segment that contains a group galaxy
    X = np.array(Table['X'], dtype=int)
    Y = np.array(Table['Y'], dtype=int)
    for i in range(len(X)):
        gal_id = seg_data[Y[i], X[i]]   # segment ID at galaxy center
        seg_data[seg_data == gal_id] = 0

    # Mask everything else
    seg_data[seg_data > 0] = 1

    fits.PrimaryHDU(seg_data, header=hdr).writeto(mask_path, overwrite=True)

    del seg_data, hdr


# ===================== SKY BACKGROUND ESTIMATION =====================

def Median_sky(position, field):
    """
    Estimate the sky background level (sigma-clipped mean) in each S-PLUS
    band for a given group cutout.

    Each CCDData object is explicitly deleted after its statistics are
    extracted to avoid accumulating 12 full images in memory per call.

    Sky-fitting flag:
      - Edge positions (first or last element of the 'c' grid) get flag=0
        (sky fixed during GALFITM fitting).
      - Interior positions get flag=1 (sky free to vary).

    Returns a GALFITM-formatted parameter line (parameter 1, sky level).
    """
    sk = []
    for filtro in FILTROS:
        img_path = ('Field_Img/%i_%i_%s_%s.fits'
                    % (position[0], position[1], field, filtro))
        image = CCDData.read(img_path, unit="adu")
        mean, _median, _std = sigma_clipped_stats(image.data, sigma=3.0)
        sk.append(mean)
        # Explicitly free the CCDData object — critical in a 400-position loop
        del image

    # Edge positions: fix the sky (flag=0); interior positions: free sky (flag=1)
    is_edge = (position[0] in (c[0], c[-1])) or (position[1] in (c[0], c[-1]))
    flag = 0 if is_edge else 1

    vals = ','.join([f'{v:.6f}' for v in sk])
    return f'1) {vals} {flag}'


# ===================== COORDINATE HELPER =====================

def _get_ra_dec_for_stamp(Table):
    """
    Return (RA_deg, DEC_deg) for the first galaxy in Table.
    Looks for RA/DEC columns first, then RAdeg/DECdeg.
    """
    if 'RA' in Table.colnames and 'DEC' in Table.colnames:
        return float(Table['RA'][0]), float(Table['DEC'][0])
    if 'RAdeg' in Table.colnames and 'DECdeg' in Table.colnames:
        return float(Table['RAdeg'][0]), float(Table['DECdeg'][0])
    raise KeyError("Table must have RA/DEC or RAdeg/DECdeg columns.")


# ===================== GALFITM INPUT FILE GENERATOR =====================

def arh_galfit(position, Table, field, size):
    """
    Generate a GALFITM multi-band input file for a group cutout.

    Zero-points are retrieved spatially from S-PLUS DR6 via conn.get_zp().
    Initial structural parameters come from SExtractor catalog columns.

    Parameters
    ----------
    position : (x, y) tuple — grid center in pixels
    Table    : astropy Table — group galaxies (one row per galaxy)
    field    : S-PLUS field name string
    size     : cutout side length in pixels

    Output file: galfit_{x}_{y}_{field}.input
    """
    # Filter suffixes used in FITS filenames
    Filtros = np.array(['U', 'F378', 'F395', 'F410', 'F430',
                        'G', 'F515', 'R', 'F660', 'I', 'F861', 'Z'])

    # ---------- GALFITM header blocks ----------
    A1 = 'A1) u,J0378,J0395,J0410,J0430,g,J0515,r,J0660,i,J0861,z'
    A2 = 'A2) 3536,3770,3940,4094,4292,4751,5133,6258,6614,7690,8611,8831'
    B  = 'B) Galfitm_%i_%i_%s.fits'  % (position[0], position[1], field)
    C  = 'C) none'
    E  = 'E) 1'
    F  = 'F) ' + ','.join(
             ['Field_Img/mask/mask_%i_%i_%s.fits' % (position[0], position[1], field)]
             * len(Filtros))
    G  = 'G) none'
    H  = 'H) 0 %i 0 %i' % (size, size)
    I  = 'I) 150 150'
    K  = 'K) 0.55 0.55'
    O  = 'O) regular'
    P  = 'P) 0'
    U  = 'U) 0'
    W  = 'W) input,model,residual'

    # ---------- Spatially-varying zero-points from DR6 ----------
    ra0, dec0 = _get_ra_dec_for_stamp(Table)
    zps = [conn.get_zp(field=field, band=b, ra=ra0, dec=dec0) for b in Filtros]
    J = 'J) ' + ','.join([str(z) for z in zps])

    # ---------- Input images and PSFs ----------
    a = ['Field_Img/%i_%i_%s_%s.fits'    % (position[0], position[1], field, b)
         for b in Filtros]
    d = ['Field_Img/psf/psf_%s_%s.fits'  % (field, b) for b in Filtros]
    A = 'A) ' + ','.join(a)
    D = 'D) ' + ','.join(d)

    Data = [
        '=' * 80,
        '# GALFITM multi-band input — generated by MorphoPlus (mascara.py)',
        '=' * 80,
        A, A1, A2, B, C, D, E, F, G, H, I, J, K, O, P, U, W,
    ]

    # ---------- One Sersic component per group galaxy ----------
    for i in range(len(Table)):
        Data.append(f'# ---- Galaxy {i} ----')
        Data.append('0) sersic')

        x1 = float(Table['X'][i])
        y1 = float(Table['Y'][i])
        Data.append('1) ' + ','.join([f'{x1:.6f}'] * 12) + ' 1')  # x, free
        Data.append('2) ' + ','.join([f'{y1:.6f}'] * 12) + ' 1')  # y, free

        # Magnitudes — replace unphysical values (>30 mag) with r_auto
        r_auto = float(Table['r_auto'][i])
        band_cols = ['u_auto', 'J0378_auto', 'J0395_auto', 'J0410_auto',
                     'J0430_auto', 'g_auto', 'J0515_auto', 'r_auto',
                     'J0660_auto', 'i_auto', 'J0861_auto', 'z_auto']
        mags = [float(Table[col][i]) for col in band_cols]
        mags = [m if m <= 30 else r_auto for m in mags]
        Data.append('3) ' + ','.join([f'{m:.6f}' for m in mags]) + ' 12')

        # Effective radius and Sersic index (Andrae et al. 2010 approximation)
        R  = float(Table['FLUX_RADIUS_50'][i])
        Cn = 5.0 * np.log10(float(Table['FLUX_RADIUS_90'][i]) / R)
        n  = (Cn / 2.77) ** (1.0 / 0.466)
        if n < 1 or isnan(n):
            n = 2.5
            R = 7.0
        Data.append('4) ' + ','.join([f'{R:.6f}'] * 12) + ' 2')  # Re, coupled
        Data.append('5) ' + ','.join([f'{n:.6f}'] * 12) + ' 2')  # n,  coupled

        # Axis ratio b/a
        el = 1.0 / float(Table['ELONGATION'][i])
        Data.append('9) ' + ','.join([f'{el:.6f}'] * 12) + ' 2')

        # Position angle (SExtractor -> GALFITM convention)
        theta = float(Table['THETA'][i])
        Th = theta - 90.0 if theta >= 0 else theta + 90.0
        Data.append('10) ' + ','.join([f'{Th:.6f}'] * 12) + ' 2')
        Data.append('Z) 0')

    # ---------- Sky component ----------
    Data += [
        '# ---- Sky ----',
        '0) sky',
        Median_sky(position, field),  # sigma-clipped mean per band
        '2) 0.000  0',                # dsky/dx, fixed
        '3) 0.000  0',                # dsky/dy, fixed
        'Z) 0',
    ]

    out_path = 'galfit_%i_%i_%s.input' % (position[0], position[1], field)
    with open(out_path, 'w') as fic:
        for line in Data:
            print(line, file=fic)
    # 'with' guarantees the file is closed even if an exception occurs


# ===================== MAIN LOOP — GROUP MASKS & GALFITM INPUTS =====================

Datos_S = S.group_by('Field')
Fields  = Datos_S.groups.keys
print(f"[INFO] Total fields to process: {len(Fields)}")

conta = 0
for f in Fields:
    field_name = f[0]
    print(f"[FIELD] {field_name}  ({conta + 1}/{len(Fields)})")

    for j in range(len(c)):
        for k in range(len(c)):
            p = (c[j], c[k])
            Tablef, Tabled = tables(S, field_name, p, size)

            if len(Tablef) > 0:
                print(f"  [MASK] position={p}  n_gal={len(Tablef)}")
                try:
                    # Build the binary mask from the segmentation map
                    mask(field_name, Tablef, p)
                except Exception as e:
                    print(f"  [ERROR] mask() failed at {p}: {e}")

                try:
                    # Generate the GALFITM input file for this group position
                    arh_galfit(p, Tablef, field_name, size)
                except Exception as e:
                    print(f"  [ERROR] arh_galfit() failed at {p}: {e}")

            # Force garbage collection after every grid position to free
            # any arrays or file handles that are no longer referenced.
            gc.collect()

    conta += 1
    print(f"[OK] {field_name} done.")


# ===================== INDIVIDUAL GALAXY MASKS =====================
# If Catalogos/g_S.csv exists, build per-galaxy masks from their
# individual segmentation maps (produced by segmetation.main_gal).

print("\n[INFO] Checking for individual galaxy catalog (Catalogos/g_S.csv)...")

archivo_csv  = "Catalogos/g_S.csv"
has_gal_catalog = False

try:
    with open(archivo_csv, "r") as archivo:
        reader = csv.reader(archivo)
        # Check that the file has at least one data row
        has_gal_catalog = any(True for _ in reader)
except FileNotFoundError:
    pass

if has_gal_catalog:
    g_S = Table.read(archivo_csv)
    print(f"[INFO] Individual galaxy catalog found: {len(g_S)} galaxies.")

    for j in range(len(g_S)):
        gid = g_S['ID'][j]
        print(f"  [GAL MASK] ID={gid}  ({j + 1}/{len(g_S)})")
        try:
            mask_gal(gid)
        except Exception as e:
            print(f"  [ERROR] mask_gal() failed for ID={gid}: {e}")

        # Release memory after each galaxy mask
        gc.collect()

    print("[INFO] Individual galaxy masks completed.")
else:
    print("[INFO] No individual galaxy catalog found; skipping per-galaxy masks.")
