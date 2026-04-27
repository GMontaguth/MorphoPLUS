# -*- coding: utf-8 -*-
"""
generate_missing_fits.py
========================
Genera script bash con los comandos GALFITM para archivos .fits que faltan.
Solo mira qué .fits no existen (independiente del CSV).

Uso:
    python generate_missing_fits.py
    # Genera missing_fits.sh
"""

import os
from ejecutable import c, size, Fields, S
from table_generation import tables

GALFITM_BIN = "./galfitm-1.4.4-linux-x86_64"

def main():
    print("Identificando archivos .fits faltantes...")
    
    missing_commands = []
    missing_count = 0
    
    # Recorrer posiciones y verificar qué .fits faltan
    for f in Fields:
        field = f[0]
        for j in range(len(c)):
            for k in range(len(c)):
                position = (c[j], c[k])
                Tablef, Tabled = tables(S, field, position, size)
                
                # Grupos: galfit_*.input -> Galfitm_*.fits
                if len(Tablef) > 0:
                    input_file = f"galfit_{position[0]}_{position[1]}_{field}.input"
                    fits_file = f"Galfitm_{position[0]}_{position[1]}_{field}.fits"
                    
                    if os.path.exists(input_file) and not os.path.exists(fits_file):
                        missing_commands.append(f"chmod 777 {input_file}")
                        missing_commands.append(f"{GALFITM_BIN} {input_file}")
                        missing_count += 1
                        print(f"  FALTA: {fits_file}")
                
                # Individuales: Gal_{ID}.input -> Gal_{ID}.fits
                if len(Tabled) > 0:
                    for r in range(len(Tabled)):
                        gid = str(Tabled['ID'][r])
                        input_file = f"Gal_{gid}.input"
                        fits_file = f"Gal_{gid}.fits"
                        
                        if os.path.exists(input_file) and not os.path.exists(fits_file):
                            missing_commands.append(f"chmod 777 {input_file}")
                            missing_commands.append(f"{GALFITM_BIN} {input_file}")
                            missing_count += 1
                            print(f"  FALTA: {fits_file}")
    
    # Generar script
    if not missing_commands:
        print("\n✅ ¡No faltan archivos .fits! Todos están completos.")
        # Eliminar script anterior si existe
        if os.path.exists("missing_fits.sh"):
            os.remove("missing_fits.sh")
            print("   (Removido missing_fits.sh anterior)")
        return
    
    script_name = "missing_fits.sh"
    with open(script_name, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"# Comandos para generar {missing_count} archivos .fits faltantes\n")
        f.write("# Generado por generate_missing_fits.py\n\n")
        f.write("set +e  # No parar si falla un fit individual\n")
        f.write("echo 'Ejecutando comandos GALFITM pendientes...'\n\n")
        
        for i, cmd in enumerate(missing_commands):
            f.write(f"{cmd}\n")
            # Agregar mensaje de progreso cada par de comandos (cada archivo)
            if i % 2 == 1:  # Después del comando galfitm (cada segundo comando)
                file_num = (i + 1) // 2
                f.write(f"echo '[{file_num}/{missing_count}] Completado'\n\n")
        
        f.write("echo '¡Todos los comandos pendientes ejecutados!'\n")
        f.write("echo 'Verificar con: python count_fits.py'\n")
    
    # Hacer ejecutable
    os.chmod(script_name, 0o755)
    
    print(f"\n✅ {script_name} creado")
    print(f"   Archivos .fits faltantes: {missing_count}")
    print(f"   Comandos generados: {len(missing_commands)}")
    print(f"   Ejecutar con: ./{script_name}")
    print(f"   Verificar después con: python count_fits.py")

if __name__ == "__main__":
    main()
