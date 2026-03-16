# -*- coding: utf-8 -*-
"""
Recortar.py — optimized version
================================
Downloads S-PLUS field frames, builds PSFs, and cuts sub-images on a
spatial grid defined by the (c, size) parameters from ejecutable.py.

Key optimizations over the original version:
  1. Parallel filter downloads per field (ThreadPoolExecutor)
  2. Skip-if-exists: individual cutout files are not re-downloaded
  3. Exponential backoff retry on failed downloads
  4. MAX_WORKERS is configurable to match your network bandwidth

Pipeline position:
  ejecutable.py -> Recortar.py -> segmetation.py -> mascara.py -> GALFITM
"""

import os
import time
import numpy as np
import numpy.ma as ma
import concurrent.futures
from functools import partial

from astropy.io import fits, ascii
from astropy.nddata.utils import Cutout2D
from astropy.table import Table, join

import splusdata

# Local modules
from segmetation import main, main_gal
from psf_new import make_psf
from table_generation import tables
from ejecutable import c, size, S
from Img_galxgal import galporgal


# ===================== CONFIGURATION =====================
DR           = "dr6"       # Preferred data release; fallbacks are tried automatically
FILTROS      = ['R', 'F378', 'F395', 'F410', 'F430',
                'F515', 'F660', 'F861', 'G', 'I', 'Z', 'U']
MAX_RETRIES  = 5           # Maximum download attempts per filter/field combination
TIMEOUT      = 300         # Seconds allowed per download attempt
BACKOFF_BASE = 2.0         # Base for exponential backoff: wait = BACKOFF_BASE^attempt (s)
MAX_WORKERS  = 6           # Parallel threads for filter downloads.
                           # Recommended: 4-8. Reduce if you hit API rate limits (HTTP 429).
DIR_OUT            = "Field_Img"
DIR_PSF            = os.path.join(DIR_OUT, "psf")
LOG_PROCESSADOS    = "Catalogos/procesados.txt"        # Fields already processed
LOG_PROCESADOS_GAL = "Catalogos/procesados_gal.txt"    # Individual galaxies already processed
TABELA_SPLUS       = "Catalogos/SPLUS_Table.csv"
EDGE_FILTER        = "R"   # Filter used to probe frame dimensions (cheap single download)

# Ensure output directories exist
os.makedirs(DIR_OUT, exist_ok=True)
os.makedirs(DIR_PSF, exist_ok=True)
os.makedirs("Catalogos", exist_ok=True)


# ===================== API CONNECTION =====================
conn = splusdata.Core('gpardo', 'gNGC5054')


# ===================== ROBUST EXCEPTION ALIAS =====================
# splusdata changed its exception class across versions; this finds whichever exists.
try:
    SPLUS_ERR = splusdata.SplusdataError
except AttributeError:
    try:
        from splusdata.core import SplusdataError as SPLUS_ERR
    except Exception:
        try:
            from splusdata.exceptions import SplusdataError as SPLUS_ERR
        except Exception:
            SPLUS_ERR = Exception  # Generic fallback


# ===================== HELPERS =====================

def normalize_field_name(field: str):
    """
    Generate alternative field-name candidates by swapping '_' <-> '-'
    and handling SPLUS_/SPLUS- prefixes.
    Returns an ordered list without duplicates.
    """
    cands, seen = [], set()
    def add(x):
        if x not in seen:
            seen.add(x)
            cands.append(x)
    add(field)
    add(field.replace('_', '-'))
    add(field.replace('-', '_'))
    if field.startswith('SPLUS_'):
        add('SPLUS-' + field.split('_', 1)[1])
    if field.startswith('SPLUS-'):
        add('SPLUS_' + field.split('-', 1)[1])
    return cands


def tentar_download_field_frame(field: str, band: str):
    """
    Download a science frame for (field, band) with retry + exponential backoff.

    Strategy:
      - Try every field-name variant from normalize_field_name()
      - Try every available data release (preferred DR first)
      - On failure, wait BACKOFF_BASE^attempt seconds before the next round
      - Give up after MAX_RETRIES total rounds

    Returns an HDUList on success, or None on total failure.
    """
    # Build release list: preferred DR first, then all available ones
    releases = [DR]
    try:
        avail = conn.check_available_images_releases()
        for r in avail:
            if r not in releases:
                releases.append(r)
    except Exception:
        pass  # If we cannot check, just use the preferred DR

    for attempt in range(MAX_RETRIES):
        last_err = None
        for f_try in normalize_field_name(field):
            for dr_try in releases:
                def _call(f=f_try, dr=dr_try):
                    return conn.field_frame(
                        field=f, weight=False,
                        band=band, data_release=dr,
                        verbose=False
                    )
                try:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as ex:
                        fut = ex.submit(_call)
                        return fut.result(timeout=TIMEOUT)
                except concurrent.futures.TimeoutError as e:
                    last_err = e
                    print(f"[TIMEOUT] field={f_try}, band={band}, DR={dr_try} "
                          f"(attempt {attempt + 1})")
                except SPLUS_ERR as e:
                    last_err = e
                except Exception as e:
                    last_err = e

        # Exponential backoff before the next attempt
        wait = BACKOFF_BASE ** attempt
        print(f"[RETRY] field={field}, band={band} — attempt {attempt + 1}/{MAX_RETRIES} "
              f"— waiting {wait:.0f}s. Last error: {last_err}")
        time.sleep(wait)

    print(f"[FAIL] field={field}, band={band} — exhausted {MAX_RETRIES} attempts.")
    return None


def extrair_fwhm_beta(header):
    """
    Extract FWHM (arcsec) and Moffat beta from a FITS header.
    Tries several keyword variants that appear across S-PLUS data releases.
    Returns (fwhm, beta); either value may be None if not found.
    """
    cand_fwhm = [
        "HIERARCH MAR PRO FWHMMEAN", "HIERARCH OAJ PRO FWHMMEAN",
        "HIERARCH OAJ PRO FWHM_MEAN", "HIERARCH OAJ PRO FWHM",
        "SEEING", "FWHMMEAN",
    ]
    cand_beta = [
        "HIERARCH OAJ PRO FWHMBETA", "HIERARCH MAR PRO FWHMBETA",
        "HIERARCH OAJ PRO BETA", "BETA", "FWHMBETA",
    ]
    fwhm = next((header[k] for k in cand_fwhm if k in header), None)
    beta = next((header[k] for k in cand_beta if k in header), None)
    return fwhm, beta


# ===================== PARALLEL DOWNLOAD WORKER =====================

def _download_and_cut_filtro(filtro, field, positions_seen_per_field):
    """
    Worker function: download ONE filter for a given field, build its PSF,
    and save all cutouts on the spatial grid.

    Designed to run inside a ThreadPoolExecutor thread.

    Skip logic:
      - All cutouts AND the PSF already exist  -> skip the download entirely.
      - Only some cutouts are missing          -> download and fill the gaps.

    Returns: (filtro, data_array, header)
             data/header are None if the filter was skipped or the download failed.
    """
    psf_path = os.path.join(DIR_PSF, f"psf_{field}_{filtro}.fits")

    # --- Check whether anything needs to be done for this filter ---
    any_cutout_missing = False
    for j in range(len(c)):
        for k in range(len(c)):
            position = (c[j], c[k])
            outpath = os.path.join(
                DIR_OUT, f'{position[0]}_{position[1]}_{field}_{filtro}.fits')
            if not os.path.exists(outpath):
                any_cutout_missing = True
                break
        if any_cutout_missing:
            break

    psf_missing = not os.path.exists(psf_path)

    if not any_cutout_missing and not psf_missing:
        print(f"[SKIP] {field}/{filtro} — all cutouts and PSF already on disk.")
        return filtro, None, None

    # --- Download the full science frame ---
    hdu_list = tentar_download_field_frame(field, filtro)
    if hdu_list is None:
        return filtro, None, None

    try:
        data = hdu_list[1].data
        hdr  = hdu_list[1].header
    except Exception:
        data = hdu_list[0].data
        hdr  = hdu_list[0].header

    # --- Build PSF if it does not exist yet ---
    if psf_missing:
        fwhm, beta = extrair_fwhm_beta(hdr)
        if fwhm is not None and beta is not None:
            try:
                make_psf(fwhm, beta, psf_path)
                print(f"[PSF] Created: {field}/{filtro}")
            except Exception as e:
                print(f"[WARNING] PSF creation failed for {field}/{filtro}: {e}")
        else:
            print(f"[WARNING] FWHM/BETA not found in header for "
                  f"{field}/{filtro}; PSF not created.")

    # --- Save cutouts, skipping positions that already exist on disk ---
    for j in range(len(c)):
        for k in range(len(c)):
            position = (c[j], c[k])
            outpath = os.path.join(
                DIR_OUT, f'{position[0]}_{position[1]}_{field}_{filtro}.fits')

            if os.path.exists(outpath):
                continue  # Skip-if-exists: avoid re-writing the file

            try:
                Tablef, _ = tables(S, field, position, size)
            except Exception as e:
                print(f"[WARNING] tables() failed for {field} at {position}: {e}")
                continue

            if len(Tablef) > 0:
                # Record this grid position as containing group members
                positions_seen_per_field.setdefault(field, set()).add(position)
                hdr_cut = hdr.copy()
                hdr_cut['CUTINFO'] = f"{position} / group center and cutout size"
                try:
                    # mode='strict' raises if the cutout touches the frame border
                    cut = Cutout2D(data, position, size, mode='strict')
                    fits.PrimaryHDU(cut.data, header=hdr_cut).writeto(
                        outpath, overwrite=True)
                except Exception as e:
                    print(f"[ERROR] Cutout failed for {field}/{filtro} "
                          f"at {position}: {e}")

    return filtro, data, hdr


# ===================== FIELD PROCESSING — PARALLEL ENTRY POINT =====================

def recortar_subcampos(field: str, positions_seen_per_field: dict):
    """
    Process one S-PLUS field:
      1. Download all 12 filters IN PARALLEL using ThreadPoolExecutor.
      2. For each filter, build the PSF and save all grid cutouts.
      3. Run segmentation/mask post-processing (sequential, called once per field).

    The degree of parallelism is set by MAX_WORKERS.
    """
    print(f"[DOWNLOAD] {field} — starting parallel download "
          f"({MAX_WORKERS} workers, {len(FILTROS)} filters)")

    # Bind fixed arguments so only 'filtro' varies across submitted tasks
    worker = partial(_download_and_cut_filtro,
                     field=field,
                     positions_seen_per_field=positions_seen_per_field)

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(worker, filtro): filtro for filtro in FILTROS}
        for fut in concurrent.futures.as_completed(futures):
            filtro = futures[fut]
            try:
                filtro_done, _, _ = fut.result()
                print(f"[DONE] {field}/{filtro_done}")
            except Exception as e:
                print(f"[ERROR] {field}/{filtro}: {e}")

    # Post-processing: build detection image and SExtractor scripts for this field
    try:
        main(field, c, size, S)
    except Exception as e:
        print(f"[WARNING] main() failed for field {field}: {e}")


# ===================== GALAXY CATALOG HELPERS =====================

def gal(field: str):
    """
    Return (ids, positions) for all galaxies detected by tables()
    in a given field. Used to build the individual-galaxy catalog.
    """
    ids, positions = [], []
    for j in range(len(c)):
        for k in range(len(c)):
            position = (c[j], c[k])
            try:
                Tablef, Tabled = tables(S, field, position, size)
            except Exception as e:
                print(f"[WARNING] tables() failed in gal() for "
                      f"{field} at {position}: {e}")
                continue
            if len(Tabled) > 0:
                for h in range(len(Tabled)):
                    ids.append(Tabled['ID'][h])
                    positions.append(position)
    return ids, positions


def coletar_de_todos_os_campos(Fields):
    """
    Rebuild the full galaxy ID/position lists by iterating over ALL fields,
    regardless of the processed-fields log. No image downloads are performed.

    Returns:
        IDS_all        : flat list of galaxy IDs
        P_all          : flat list of corresponding grid positions
        positions_seen : dict mapping field_name -> set of (x, y) positions
                         that contain at least one group member
    """
    IDS_all, P_all = [], []
    positions_seen = {}
    for frow in Fields:
        raw = frow['Field'] if 'Field' in frow.colnames else frow[0]
        if ma.is_masked(raw):
            print("[WARNING] Skipping row with masked 'Field' value.")
            continue
        field_name = str(raw)
        for j in range(len(c)):
            for k in range(len(c)):
                position = (c[j], c[k])
                try:
                    Tablef, Tabled = tables(S, field_name, position, size)
                except Exception:
                    continue
                if len(Tabled) > 0:
                    positions_seen.setdefault(field_name, set()).add(position)
                    for h in range(len(Tabled)):
                        IDS_all.append(Tabled['ID'][h])
                        P_all.append(position)
    return IDS_all, P_all, positions_seen


def inside_frame(nx: int, ny: int, position, size):
    """
    Return True if a cutout of 'size' centered at 'position'
    fits entirely inside the frame with dimensions (nx, ny).
    'size' can be a scalar int or a (sx, sy) sequence.
    """
    if isinstance(size, (list, tuple, np.ndarray)):
        sx = float(size[0])
        sy = float(size[1]) if len(size) > 1 else float(size[0])
    else:
        sx = sy = float(size)
    hx, hy = sx / 2.0, sy / 2.0
    x, y = float(position[0]), float(position[1])
    return (x - hx >= 0) and (x + hx <= nx - 1) and \
           (y - hy >= 0) and (y + hy <= ny - 1)


# ===================== INPUT TABLES =====================
S = Table.read(TABELA_SPLUS)
Dados_S = S.group_by('Field')
Fields  = Dados_S.groups.keys
print(f"[INFO] Total fields in catalog: {len(Fields)}")


# ===================== LOAD PROCESSED LOGS =====================
# Fields already processed (cutouts exist on disk)
campos_processados = set()
if os.path.exists(LOG_PROCESSADOS):
    with open(LOG_PROCESSADOS) as f:
        campos_processados = {line.strip() for line in f if line.strip()}
print(f"[INFO] Fields already processed: {len(campos_processados)}")

# Individual galaxies already fitted by GALFITM
gals_procesadas = set()
if os.path.exists(LOG_PROCESADOS_GAL):
    with open(LOG_PROCESADOS_GAL) as f:
        gals_procesadas = {line.strip() for line in f if line.strip()}
print(f"[INFO] Individual galaxies already processed: {len(gals_procesadas)}")


def marcar_gal_procesada(gal_id: str):
    """
    Append a galaxy ID to the processed-galaxies log and in-memory set.
    Uses fsync to guarantee the write survives a crash or keyboard interrupt.
    """
    gal_id = str(gal_id).strip()
    if not gal_id or gal_id in gals_procesadas:
        return
    with open(LOG_PROCESADOS_GAL, "a") as f:
        f.write(f"{gal_id}\n")
        f.flush()
        os.fsync(f.fileno())
    gals_procesadas.add(gal_id)


# ===================== MAIN LOOP — FIELDS =====================
IDS_g, p_g = [], []
positions_seen_per_field = {}

for frow in Fields:
    raw = frow['Field'] if 'Field' in frow.colnames else frow[0]
    if ma.is_masked(raw):
        print("[WARNING] Skipping row with masked 'Field' value in main loop.")
        continue
    field_name = str(raw)

    if field_name in campos_processados:
        print(f"[SKIP] Field {field_name} already processed.")
        continue

    print(f"\n{'=' * 60}")
    print(f"[START] Field: {field_name}  |  DR: {DR}")
    print(f"{'=' * 60}")

    try:
        # Download all 12 filters in parallel and produce the spatial grid cutouts
        recortar_subcampos(field_name, positions_seen_per_field)

        # Collect IDs of galaxies found in this field
        I, P = gal(field_name)
        if len(I) > 0:
            IDS_g.append(I)
            p_g.append(P)

        # Mark field as successfully processed
        with open(LOG_PROCESSADOS, "a") as flog:
            flog.write(f"{field_name}\n")
        print(f"[OK] Field {field_name} completed.")

    except Exception as e:
        print(f"[ERROR] Failure processing field {field_name}: {e}")


# ===================== FINAL REBUILD — ALL FIELDS =====================
print("\n[REBUILD] Reconstructing galaxy catalog from all fields "
      "(ignoring the processed-fields log)...")
IDS_g, p_g, positions_seen_all = coletar_de_todos_os_campos(Fields)

if len(IDS_g) > 0:
    # Map galaxy IDs to their grid positions
    g_S = Table([IDS_g, p_g], names=['ID', 'position'])
    try:
        gal_S = join(S, g_S, keys='ID', join_type='inner')
    except Exception as e:
        print(f"[ERROR] Table join failed: {e}")
        gal_S = g_S

    # Remove the helper column before saving
    if 'position' in gal_S.colnames:
        gal_S.remove_column('position')

    try:
        ascii.write(gal_S, 'Catalogos/g_S.csv', format='csv',
                    fast_writer=False, overwrite=True)
        print(f"[INFO] Saved galaxy catalog: Catalogos/g_S.csv ({len(gal_S)} rows)")
    except Exception as e:
        print(f"[WARNING] Could not save Catalogos/g_S.csv: {e}")

    # ===================== INDIVIDUAL GALAXY LOOP (RESUMABLE) =====================
    print(f"\n[GALAXIES] Total in catalog : {len(gal_S)}")
    print(f"[GALAXIES] Already processed: {len(gals_procesadas)}")

    if 'ID' not in gal_S.colnames:
        raise RuntimeError("[FATAL] gal_S has no 'ID' column.")

    # Build a boolean mask to select only pending galaxies
    pending_mask = np.array(
        [str(gid).strip() not in gals_procesadas for gid in gal_S['ID']],
        dtype=bool
    )
    gal_pending = gal_S[pending_mask]
    print(f"[GALAXIES] Pending: {len(gal_pending)}")

    for i, row in enumerate(gal_pending, start=1):
        gid = str(row['ID']).strip()
        if not gid or gid in gals_procesadas:
            continue

        print(f"[GALAXY START] {i}/{len(gal_pending)}  ID={gid}")

        # Single-row sub-table for this galaxy
        gal_1 = gal_S[gal_S['ID'] == row['ID']]

        try:
            # Step 1: Download multi-band stamps & generate GALFITM input file
            galporgal(gal_1)
            # Step 2: Build detection image and SExtractor script for this galaxy
            main_gal(gal_1)
            # Mark as done only if both steps succeeded
            marcar_gal_procesada(gid)
            print(f"[GALAXY OK] ID={gid}")
        except Exception as e:
            print(f"[GALAXY ERROR] ID={gid}: {e}")
            # Do NOT mark as processed; it will be retried on the next run

else:
    print("[INFO] No galaxies collected in this run (IDS_g is empty).")
