# -*- coding: utf-8 -*-
"""
count_fits.py
=============
Cuenta SOLO archivos .fits de GALFITM - versión mínima.

Uso:
    python count_fits.py
"""

import os
from ejecutable import c, size, Fields, S
from table_generation import tables

def main():
    # Contar archivos esperados y existentes
    expected_groups = 0
    expected_indiv = 0
    existing_groups = 0
    existing_indiv = 0
    
    for f in Fields:
        field = f[0]
        for j in range(len(c)):
            for k in range(len(c)):
                position = (c[j], c[k])
                Tablef, Tabled = tables(S, field, position, size)
                
                # Grupos
                if len(Tablef) > 0:
                    expected_groups += 1
                    fits_file = f"Galfitm_{position[0]}_{position[1]}_{field}.fits"
                    if os.path.exists(fits_file):
                        existing_groups += 1
                
                # Individuales
                for r in range(len(Tabled)):
                    expected_indiv += 1
                    gid = str(Tabled['ID'][r])
                    fits_file = f"Gal_{gid}.fits"
                    if os.path.exists(fits_file):
                        existing_indiv += 1
    
    total_expected = expected_groups + expected_indiv
    total_existing = existing_groups + existing_indiv
    
    print(f"GALFITM .fits: {total_existing}/{total_expected}")
    print(f"  Grupos: {existing_groups}/{expected_groups}")
    print(f"  Indiv:  {existing_indiv}/{expected_indiv}")
    
    if total_existing == total_expected:
        print("✅ COMPLETO")
    else:
        progress = (total_existing / total_expected) * 100 if total_expected > 0 else 0
        print(f"⏳ {progress:.1f}%")

if __name__ == "__main__":
    main()
