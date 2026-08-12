# -*- coding: utf-8 -*-
"""
verify_downloads.py
===================
Run AFTER Recortar.py.

1) Checks which fields were actually downloaded (each field must have all
   12 SPLUS-*_FILTER.fits bands on disk).
2) Removes from Catalogos/procesados.txt the fields that were NOT fully
   downloaded, so the next Recortar.py run retries them.

Usage:
    python verify_downloads.py            # check and clean procesados.txt
    python verify_downloads.py --dry-run  # check only, modify nothing
"""

import os
import sys
import glob
from astropy.table import Table

# The 12 bands each fully-downloaded field must have
BANDS = ['U', 'F378', 'F395', 'F410', 'F430', 'G',
         'F515', 'R', 'F660', 'I', 'F861', 'Z']

# Where to look for the field images
SEARCH_DIRS = [".", "Field_Img", "Field_Img/det"]


def get_all_fields():
    """Get the list of all expected fields from SPLUS_Table.csv."""
    splus_file = "Catalogos/SPLUS_Table.csv"
    if not os.path.exists(splus_file):
        print(f"ERROR: {splus_file} not found")
        return []
    S = Table.read(splus_file, format='csv')
    if 'Field' not in S.colnames:
        print(f"ERROR: No 'Field' column in {splus_file}")
        return []
    return sorted(set(str(f).strip() for f in S['Field']))


def find_field_bands(field):
    """Return the set of bands found on disk for a field."""
    found = set()
    for band in BANDS:
        for d in SEARCH_DIRS:
            pattern = os.path.join(d, f"*{field}*_{band}.fits")
            if glob.glob(pattern):
                found.add(band)
                break
    return found


def find_procesados_file():
    """Locate procesados.txt in a few possible places."""
    for path in ("Catalogos/procesados.txt", "procesados.txt"):
        if os.path.exists(path):
            return path
    matches = glob.glob("**/procesados*.txt", recursive=True)
    return matches[0] if matches else None


def main():
    dry_run = "--dry-run" in sys.argv

    print("Verifying field downloads (12 bands)...")
    print("=" * 60)
    if dry_run:
        print("DRY-RUN MODE: procesados.txt will not be modified\n")

    # 1. Expected fields
    all_fields = get_all_fields()
    if not all_fields:
        return
    print(f"Total fields in catalog: {len(all_fields)}\n")

    # 2. Check each field
    complete = []
    incomplete = []   # (field, n_found, missing_bands)
    missing = []

    for field in all_fields:
        found = find_field_bands(field)
        n = len(found)
        if n == len(BANDS):
            complete.append(field)
        elif n == 0:
            missing.append(field)
        else:
            incomplete.append((field, n, sorted(set(BANDS) - found)))

    # 3. Report
    print(f"COMPLETE ({len(complete)}/{len(all_fields)})")
    if incomplete:
        print(f"\nINCOMPLETE ({len(incomplete)}):")
        for field, n, missing_bands in incomplete:
            print(f"  {field}: {n}/{len(BANDS)} bands (missing: {', '.join(missing_bands)})")
    if missing:
        print(f"\nNOT DOWNLOADED ({len(missing)}):")
        for field in missing:
            print(f"  {field}: 0 bands")

    # Fields to retry = incomplete + not downloaded
    fields_to_retry = [f for f, _, _ in incomplete] + missing

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  Complete:       {len(complete)}")
    print(f"  Incomplete:     {len(incomplete)}")
    print(f"  Not downloaded: {len(missing)}")
    print(f"  TO RETRY:       {len(fields_to_retry)}")

    if not fields_to_retry:
        print("\nAll fields downloaded completely. Nothing to retry.")
        return

    # 4. Clean procesados.txt
    log_file = find_procesados_file()
    if not log_file:
        print("\nprocesados.txt not found. Cannot clean automatically.")
        print("Fields to retry:")
        for f in sorted(fields_to_retry):
            print(f"  {f}")
        return

    print(f"\nLog file: {log_file}")

    with open(log_file, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    retry_set = set(fields_to_retry)
    cleaned = []
    removed = []
    for line in lines:
        if any(field in line for field in retry_set):
            removed.append(line)
        else:
            cleaned.append(line)

    print(f"procesados.txt: {len(lines)} fields -> remove {len(removed)}")

    if dry_run:
        print("\n(DRY-RUN) Fields that would be removed:")
        for line in removed:
            print(f"  {line}")
        print("\nRun without --dry-run to apply the changes.")
        return

    if removed:
        # Backup
        backup = f"{log_file}.backup"
        with open(backup, "w") as f:
            for line in lines:
                f.write(f"{line}\n")
        print(f"OK Backup: {backup}")

        # Rewrite cleaned
        with open(log_file, "w") as f:
            for line in cleaned:
                f.write(f"{line}\n")
        print(f"OK {log_file} updated: {len(lines)} -> {len(cleaned)} fields")
    else:
        print("None of the fields to retry were in procesados.txt.")

    print("\nNEXT STEP:")
    print("  python Recortar.py   # will retry downloading the removed fields")


if __name__ == "__main__":
    main()
