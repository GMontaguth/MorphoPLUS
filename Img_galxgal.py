# -*- coding: utf-8 -*-
"""
Img_galxgal.py — optimized version
=====================================
Downloads multi-band S-PLUS stamps for individual galaxies and generates
GALFITM input files using spatially-varying zero-points (DR6).

Key optimizations over the original version:
  1. All 12 stamps for one galaxy are downloaded IN PARALLEL (ThreadPoolExecutor)
  2. Skip-if-exists: existing FITS files are not re-downloaded
  3. Exponential backoff retry on each individual stamp download
  4. MAX_WORKERS is configurable to match your network bandwidth

Pipeline position:
  Recortar.py -> Img_galxgal.py (galporgal) -> mascara.py -> GALFITM
"""

import os
import time
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.nddata import CCDData
from math import isnan
import concurrent.futures

import splusdata

conn = splusdata.Core('gpardo', 'gNGC5054')


# ===================== CONFIGURATION =====================
MAX_WORKERS  = 6       # Parallel threads for downloading the 12 bands of one galaxy.
                       # Reduce to 3-4 if the API returns rate-limit errors (HTTP 429).
MAX_RETRIES  = 4       # Maximum download attempts per stamp
BACKOFF_BASE = 2.0     # Exponential backoff base in seconds: wait = BACKOFF_BASE^attempt

# S-PLUS filter names used as filename suffixes
Filters = np.array(['R', 'F378', 'F395', 'F410', 'F430',
                    'F515', 'F660', 'F861', 'G', 'I', 'Z', 'U'])


# ===================== SINGLE-STAMP DOWNLOAD WORKER =====================

def _download_stamp(ra_deg, dec_deg, src_id, field_name, banda,
                    size_pix=200, dr='dr6'):
    """
    Download and save ONE stamp (one band, one galaxy) to disk.

    Features:
      - Skip-if-exists: returns immediately if the file is already on disk.
      - Retry with exponential backoff on any download failure.
      - Writes useful metadata (ID, RA, DEC, band, field) into the FITS header.

    Designed to run inside a ThreadPoolExecutor thread.

    Returns: (banda, success: bool, status: str)
             status is 'skip', 'ok', or an error message string.
    """
    out = f'Field_Img/Gal_{src_id}_{banda}.fits'

    # --- Skip if the file already exists on disk ---
    if os.path.exists(out):
        return banda, True, "skip"

    last_err = None
    for attempt in range(MAX_RETRIES):
        try:
            hdulist = conn.stamp(
                ra=ra_deg, dec=dec_deg,
                size=int(size_pix), band=banda, weight=False
            )

            # Science data is usually in extension [1]; fall back to [0]
            try:
                hdu_data = hdulist[1].data
                hdr = hdulist[1].header
            except Exception:
                hdu_data = hdulist[0].data
                hdr = hdulist[0].header

            # Add useful metadata to the header
            hdr['SRC_ID']  = (str(src_id),    'Galaxy ID')
            hdr['RA_DEG']  = (ra_deg,          'RA (deg)')
            hdr['DEC_DEG'] = (dec_deg,         'Dec (deg)')
            hdr['CUTSIZE'] = (int(size_pix),   'Cutout size (pixels)')
            hdr['BAND']    = (banda,            'S-PLUS band')
            hdr['DRUSE']   = (dr,               'Data release used')
            hdr['FIELD']   = (field_name,       'S-PLUS field')

            fits.PrimaryHDU(hdu_data, header=hdr).writeto(out, overwrite=True)
            return banda, True, "ok"

        except Exception as e:
            last_err = e
            wait = BACKOFF_BASE ** attempt
            print(f"[RETRY] Galaxy {src_id}/{banda} — attempt {attempt + 1}/{MAX_RETRIES} "
                  f"— waiting {wait:.0f}s. Error: {e}")
            time.sleep(wait)

    print(f"[FAIL] Galaxy {src_id}/{banda} — exhausted {MAX_RETRIES} attempts. "
          f"Last error: {last_err}")
    return banda, False, str(last_err)


# ===================== PARALLEL STAMP DOWNLOAD FOR ONE GALAXY =====================

def _download_galaxy_stamps(ra_deg, dec_deg, src_id, field_name,
                             size_pix=200, dr='dr6'):
    """
    Download all 12 S-PLUS stamps for a single galaxy IN PARALLEL.

    Submits one thread per filter to a ThreadPoolExecutor and collects results.

    Returns:
        results : dict mapping each band name to its status string
                  ('skip', 'ok', or an error message)
    """
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_band = {
            executor.submit(
                _download_stamp,
                ra_deg, dec_deg, src_id, field_name, banda, size_pix, dr
            ): banda
            for banda in Filters
        }
        for fut in concurrent.futures.as_completed(future_to_band):
            banda = future_to_band[fut]
            try:
                band_done, ok, status = fut.result()
                results[band_done] = status
            except Exception as e:
                results[banda] = f"exception: {e}"
                print(f"[ERROR] Galaxy {src_id}/{banda}: {e}")

    return results


# ===================== img() — PUBLIC INTERFACE =====================

def img(Fi, size_pix=200, dr='dr6'):
    """
    Download multi-band cutouts for all galaxies in table Fi.

    OPTIMIZED: the 12 filters for each galaxy are downloaded in parallel
    (see _download_galaxy_stamps). Files already on disk are skipped.

    Parameters
    ----------
    Fi       : astropy Table with columns ID, Field, and RA/DEC (or RA_1/DEC_1)
    size_pix : cutout size in pixels (default 200)
    dr       : S-PLUS data release string (default 'dr6')
    """
    # Detect which RA/DEC column names are present
    if 'RA' in Fi.colnames and 'DEC' in Fi.colnames:
        ra_col, dec_col = 'RA', 'DEC'
    elif 'RA_1' in Fi.colnames and 'DEC_1' in Fi.colnames:
        ra_col, dec_col = 'RA_1', 'DEC_1'
    else:
        raise KeyError("Table must have RA/DEC or RA_1/DEC_1 columns.")

    for i in range(len(Fi)):
        ra_deg     = float(Fi[ra_col][i])
        dec_deg    = float(Fi[dec_col][i])
        src_id     = Fi['ID'][i]
        field_name = str(Fi['Field'][i])

        print(f"[GALAXY {i + 1}/{len(Fi)}] ID={src_id} — "
              f"downloading {len(Filters)} bands in parallel...")

        results = _download_galaxy_stamps(
            ra_deg, dec_deg, src_id, field_name, size_pix, dr)

        # Summary log for this galaxy
        ok_count   = sum(1 for s in results.values() if s in ('ok', 'skip'))
        skip_count = sum(1 for s in results.values() if s == 'skip')
        fail_count = sum(1 for s in results.values() if s not in ('ok', 'skip'))
        print(f"  -> {ok_count}/{len(Filters)} OK "
              f"({skip_count} already on disk, {fail_count} failed)")


# ===================== SKY BACKGROUND ESTIMATION =====================

def Median_sky(galaxy):
    """
    Estimate the sky background (sigma-clipped mean) for each S-PLUS band
    from the stamp of a given galaxy.

    Returns a GALFITM-formatted sky parameter string (parameter line 1).
    """
    # Ordered filter list used internally by GALFITM
    Filtros = np.array(['U', 'F378', 'F395', 'F410', 'F430',
                        'G', 'F515', 'R', 'F660', 'I', 'F861', 'Z'])
    sk = []
    for band in Filtros:
        image = CCDData.read(f'Field_Img/Gal_{galaxy}_{band}.fits', unit="adu")
        mean, median, std = sigma_clipped_stats(image.data, sigma=3.0)
        sk.append(mean)

    return '1) ' + ','.join([f'{v:.6f}' for v in sk]) + ' 1'


# ===================== HELPERS =====================

def _get_ra_dec_from_row(row):
    """
    Extract (RA_deg, DEC_deg) from a table row.
    Looks for RA/DEC first, then RA_1/DEC_1.
    Raises KeyError if neither pair is found.
    """
    colnames = row.colnames
    if 'RA' in colnames and 'DEC' in colnames:
        return float(row['RA']), float(row['DEC'])
    elif 'RA_1' in colnames and 'DEC_1' in colnames:
        return float(row['RA_1']), float(row['DEC_1'])
    raise KeyError("Could not find RA/DEC or RA_1/DEC_1 columns in row.")


# ===================== GALFITM INPUT FILE GENERATOR =====================

def arh_galfit(GRf):
    """
    Generate a GALFITM multi-band input file for one galaxy.

    Uses spatially-varying zero-points retrieved from S-PLUS DR6 via
    conn.get_zp(field, band, ra, dec).

    Initial structural parameters are estimated from SExtractor columns:
      - Position  : center of the 200x200 stamp (x=100, y=100)
      - Magnitude : *_auto columns (bad values > 30 are replaced by r_auto)
      - Re        : FLUX_RADIUS_50
      - n (Sersic): estimated from FLUX_RADIUS_90/FLUX_RADIUS_50 ratio
                    using the Andrae et al. (2010) approximation
      - b/a       : 1 / ELONGATION
      - PA        : THETA (converted to GALFITM convention)

    Output file: galfit_{ID}.input
    """
    galaxy = GRf['ID']
    Field  = str(GRf['Field'])

    # Filter names: used as FITS filename suffixes and for get_zp calls
    Filtros = np.array(['U', 'F378', 'F395', 'F410', 'F430',
                        'G', 'F515', 'R', 'F660', 'I', 'F861', 'Z'])

    # ---------- GALFITM header blocks ----------
    A1 = 'A1) u,J0378,J0395,J0410,J0430,g,J0515,r,J0660,i,J0861,z'  # band labels
    A2 = 'A2) 3536,3770,3940,4094,4292,4751,5133,6258,6614,7690,8611,8831'  # pivot wavelengths (A)
    B  = f'B) Gal_{galaxy}.fits'         # Output model FITS
    C  = 'C) none'                        # No sigma image
    E  = 'E) 1'                           # PSF fine-sampling factor
    F  = 'F) ' + ','.join(               # Bad-pixel mask (one per band)
             [f'Field_Img/mask/mask_{galaxy}.fits'] * 12)
    G  = 'G) none'                        # No constraints file
    H  = 'H) 0 200 0 200'                 # Fitting region (full 200x200 stamp)
    I  = 'I) 150 150'                     # Convolution box size
    K  = 'K) 0.55 0.55'                   # Plate scale [arcsec/pixel]
    O  = 'O) regular'
    P  = 'P) 0'
    U  = 'U) 0'                           # Standard parametric fitting
    W  = 'W) input,model,residual'        # Output image types

    # ---------- Spatially-varying zero-points from DR6 ----------
    ra_deg, dec_deg = _get_ra_dec_from_row(GRf)
    zps = [conn.get_zp(field=Field, band=b, ra=ra_deg, dec=dec_deg)
           for b in Filtros]
    J = 'J) ' + ','.join([f'{z}' for z in zps])

    # ---------- Input image and PSF file lists ----------
    a = [f'Field_Img/Gal_{galaxy}_{b}.fits'          for b in Filtros]
    d = [f'Field_Img/psf/psf_{Field}_{b}.fits'        for b in Filtros]
    A = 'A) ' + ','.join(a)
    D = 'D) ' + ','.join(d)

    # ---------- Assemble the input file ----------
    Data = [
        '=' * 80,
        '# GALFITM multi-band input file — generated by MorphoPlus',
        '=' * 80,
        A, A1, A2, B, C, D, E, F, G, H, I, J, K, O, P, U, W,
        '# ---- Galaxy component ----',
        '0) sersic',
    ]

    # Initial centroid (center of the 200x200 stamp)
    x1 = y1 = 100.0
    Data.append('1) ' + ','.join([f'{x1:.6f}'] * 12) + ' 1')  # x position, free
    Data.append('2) ' + ','.join([f'{y1:.6f}'] * 12) + ' 1')  # y position, free

    # Initial magnitudes from SExtractor AUTO photometry
    bands_auto = ['u', 'J0378', 'J0395', 'J0410', 'J0430',
                  'g', 'J0515', 'r', 'J0660', 'i', 'J0861', 'z']
    mags  = [float(GRf[f'{b}_auto']) for b in bands_auto]
    r_mag = float(GRf['r_auto'])
    # Replace unphysical magnitudes (>30) with r_auto as fallback
    mags  = [m if m <= 30 else r_mag for m in mags]
    Data.append('3) ' + ','.join([f'{m:.6f}' for m in mags]) + ' 12')  # mags, free in all bands

    # Initial effective radius and Sersic index
    R  = float(GRf['FLUX_RADIUS_50'])
    Cn = 5.0 * np.log10(float(GRf['FLUX_RADIUS_90']) / R)
    n  = (Cn / 2.77) ** (1.0 / 0.466)  # Andrae et al. (2010) approximation
    if n < 1 or isnan(n):
        # Fallback to a generic early-type galaxy profile
        n = 2.5
        R = 7.0
    Data.append('4) ' + ','.join([f'{R:.6f}'] * 12) + ' 2')   # Re, free (wavelength-coupled)
    Data.append('5) ' + ','.join([f'{n:.6f}'] * 12) + ' 2')   # n,  free (wavelength-coupled)

    # Axis ratio b/a = 1 / elongation
    el = 1.0 / float(GRf['ELONGATION'])
    Data.append('9) ' + ','.join([f'{el:.6f}'] * 12) + ' 2')  # b/a, free (wavelength-coupled)

    # Position angle: convert from SExtractor convention (N->E) to GALFITM (E of N)
    theta = float(GRf['THETA'])
    Th = theta - 90.0 if theta >= 0 else theta + 90.0
    Data.append('10) ' + ','.join([f'{Th:.6f}'] * 12) + ' 2') # PA, free (wavelength-coupled)
    Data.append('Z) 0')

    # ---------- Sky component ----------
    Data += [
        '# ---- Sky component ----',
        '0) sky',
        Median_sky(galaxy),   # Background level (sigma-clipped mean per band)
        '2) 0.000  0',        # dsky/dx, fixed
        '3) 0.000  0',        # dsky/dy, fixed
        'Z) 0',
    ]

    with open(f"galfit_{galaxy}.input", "w") as fic:
        for line in Data:
            print(line, file=fic)


# ===================== MAIN ENTRY POINT =====================

def galporgal(Fi):
    """
    Full per-galaxy pipeline:
      1. Download multi-band stamps in parallel (img)
      2. Generate a GALFITM input file for each galaxy (arh_galfit)
      3. Append GALFITM run commands to ejecutable_gal.sh

    Parameters
    ----------
    Fi : astropy Table — one or more galaxies to process
    """
    # Step 1: Download stamps (parallel internally)
    img(Fi)

    script_lines = []
    for i in range(len(Fi)):
        # Step 2: Generate GALFITM input file
        arh_galfit(Fi[i])
        galaxy = Fi['ID'][i]
        # Collect the run commands for the shell script
        script_lines.append(f"chmod 777 galfit_{galaxy}.input")
        script_lines.append(f"./galfitm-1.4.4-linux-x86_64 galfit_{galaxy}.input")

    # Step 3: Append commands to the execution script
    with open('ejecutable_gal.sh', 'a') as fic:
        for line in script_lines:
            print(line, file=fic)
