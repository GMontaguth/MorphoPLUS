# -*- coding: utf-8 -*-
"""
generate_missing_individual_fits.py
===================================
Genera script bash con comandos GALFITM solo para GALAXIAS INDIVIDUALES.
Lee directamente el catálogo g_S.csv (generado por Recortar.py).

Uso:
    python generate_missing_individual_fits.py
    # Genera missing_individual_fits.sh
"""

import os
from astropy.table import Table

GALFITM_BIN = "./galfitm-1.4.4-linux-x86_64"

def main():
    print("Identificando archivos .fits faltantes para GALAXIAS INDIVIDUALES...")
    
    # Cargar catálogo de galaxias individuales
    catalog_file = "Catalogos/g_S.csv"
    if not os.path.exists(catalog_file):
        print(f"❌ ERROR: No se encuentra {catalog_file}")
        print("   Este archivo se genera al ejecutar Recortar.py")
        print("   Contiene el catálogo de galaxias individuales detectadas.")
        return
    
    try:
        gal_catalog = Table.read(catalog_file, format='csv')
    except Exception as e:
        print(f"❌ ERROR leyendo {catalog_file}: {e}")
        return
    
    print(f"  Catálogo de galaxias individuales: {len(gal_catalog)} objetos")
    
    if 'ID' not in gal_catalog.colnames:
        print(f"❌ ERROR: No se encuentra columna 'ID' en {catalog_file}")
        print(f"   Columnas disponibles: {gal_catalog.colnames}")
        return
    
    missing_commands = []
    missing_count = 0
    total_checked = 0
    
    # Verificar cada galaxia individual
    for row in gal_catalog:
        gid = str(row['ID']).strip()
        if not gid:
            continue
            
        total_checked += 1
        input_file = f"Gal_{gid}.input"
        fits_file = f"Gal_{gid}.fits"
        
        # Verificar si falta el .fits pero existe el .input
        if os.path.exists(input_file) and not os.path.exists(fits_file):
            missing_commands.append(f"chmod 777 {input_file}")
            missing_commands.append(f"{GALFITM_BIN} {input_file}")
            missing_count += 1
            print(f"  FALTA: {fits_file}")
        elif not os.path.exists(input_file):
            print(f"  [INFO] Sin .input: {input_file} (saltado)")
        else:
            print(f"  [OK] Completo: {fits_file}")
    
    print(f"\nResumen:")
    print(f"  Galaxias verificadas: {total_checked}")
    print(f"  Archivos .fits faltantes: {missing_count}")
    
    # Generar script
    if not missing_commands:
        print("\n✅ ¡No faltan archivos .fits individuales! Todos están completos.")
        # Eliminar script anterior si existe
        script_name = "missing_individual_fits.sh"
        if os.path.exists(script_name):
            os.remove(script_name)
            print(f"   (Removido {script_name} anterior)")
        return
    
    script_name = "missing_individual_fits.sh"
    with open(script_name, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"# Comandos para generar {missing_count} archivos .fits de galaxias individuales\n")
        f.write(f"# Generado por generate_missing_individual_fits.py desde {catalog_file}\n\n")
        f.write("set +e  # No parar si falla un fit individual\n")
        f.write("echo 'Ejecutando GALFITM para galaxias individuales pendientes...'\n\n")
        
        for i, cmd in enumerate(missing_commands):
            f.write(f"{cmd}\n")
            # Mensaje de progreso cada par de comandos
            if i % 2 == 1:  # Después del comando galfitm
                file_num = (i + 1) // 2
                f.write(f"echo '[{file_num}/{missing_count}] Galaxia completada'\n\n")
        
        f.write("echo '¡Todas las galaxias individuales pendientes ejecutadas!'\n")
        f.write("echo 'Verificar con: python count_individual_fits.py'\n")
    
    # Hacer ejecutable
    os.chmod(script_name, 0o755)
    
    print(f"\n✅ {script_name} creado")
    print(f"   Galaxias individuales faltantes: {missing_count}")
    print(f"   Comandos generados: {len(missing_commands)}")
    print(f"   Ejecutar con: ./{script_name}")

if __name__ == "__main__":
    main()
