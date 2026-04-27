# -*- coding: utf-8 -*-
"""
generate_psf_from_field.py
==========================
Genera PSF faltante buscando imágenes SPLUS de campo completo.

Busca archivos con patrón SPLUS-*_FILTRO.fits, lee los parámetros PSF
del header y genera la PSF faltante.

Uso:
    python generate_psf_from_field.py CAMPO FILTRO
    
Ejemplos:
    python generate_psf_from_field.py STRIPE82-0001 F378
    python generate_psf_from_field.py HYDRA-0001 R
"""

import sys
import os
import glob
import numpy as np
from astropy.io import fits

# Función make_psf copiada del script original
radius = 22

def make_psf(fwhm, beta, outfile):
    """Genera PSF con perfil Moffat usando parámetros FWHM y beta."""
    fwhm = fwhm/0.5
    alpha = fwhm / (2 * np.sqrt(np.power(2., 1/beta) - 1.))
    r = np.linspace(-radius, radius, 2 * radius + 1)
    X, Y = np.meshgrid(r, r)
    R = np.sqrt(X**2 + Y**2)
    I = (beta - 1.) / (np.pi * alpha**2) * \
        np.power(1. + np.power(R / alpha, 2), -beta)
    hdu = fits.PrimaryHDU(I)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(outfile, overwrite=True)
    return


def find_splus_images(filter_name):
    """
    Busca imágenes SPLUS que terminen con _FILTRO.fits
    """
    print(f"Buscando imágenes SPLUS con filtro {filter_name}...")
    
    # Patrones de búsqueda en diferentes ubicaciones
    search_patterns = [
        f"*_{filter_name}.fits",
        f"SPLUS-*_{filter_name}.fits", 
        f"Field_Img/*_{filter_name}.fits",
        f"Field_Img/det/*_{filter_name}.fits",
        f"**/*_{filter_name}.fits",  # Búsqueda recursiva
    ]
    
    found_images = []
    
    for pattern in search_patterns:
        matches = glob.glob(pattern, recursive=True)
        for match in matches:
            if match not in found_images:
                found_images.append(match)
    
    # Filtrar solo las que parecen imágenes SPLUS de campo completo
    splus_images = []
    for img in found_images:
        basename = os.path.basename(img)
        if ('SPLUS' in basename.upper() or 
            'STRIPE' in basename.upper() or 
            'HYDRA' in basename.upper() or
            len(basename.split('_')) >= 2):  # Formato field_filter.fits
            splus_images.append(img)
    
    print(f"  Encontradas {len(splus_images)} imágenes:")
    for img in splus_images:
        print(f"    {img}")
    
    return splus_images


def extract_psf_from_field_image(image_file):
    """
    Extrae parámetros PSF del header de una imagen de campo completo.
    """
    print(f"\nExtrayendo PSF de: {image_file}")
    
    try:
        with fits.open(image_file) as hdul:
            # Probar diferentes HDUs
            for i, hdu in enumerate(hdul):
                if hdu.header is None:
                    continue
                    
                header = hdu.header
                print(f"  Examinando HDU {i}...")
                
                # Keywords posibles para FWHM
                fwhm_keys = [
                    'FWHM', 'PSF_FWHM', 'SEEING', 'PSF_SEEING', 
                    'FWHM_PSF', 'SEESMIC', 'HIERARCH FWHM'
                ]
                
                # Keywords posibles para beta
                beta_keys = [
                    'BETA', 'PSF_BETA', 'MOFFAT_BETA', 'PSF_MOFFAT_BETA',
                    'HIERARCH BETA', 'HIERARCH PSF_BETA'
                ]
                
                fwhm = None
                beta = None
                
                # Buscar FWHM
                for key in fwhm_keys:
                    if key in header:
                        try:
                            fwhm = float(header[key])
                            print(f"    ✓ FWHM = {fwhm} (keyword: {key})")
                            break
                        except ValueError:
                            continue
                
                # Buscar beta
                for key in beta_keys:
                    if key in header:
                        try:
                            beta = float(header[key])
                            print(f"    ✓ BETA = {beta} (keyword: {key})")
                            break
                        except ValueError:
                            continue
                
                # Si encontramos ambos, listo
                if fwhm is not None and beta is not None:
                    return fwhm, beta
                
                # Si no encontramos con nombres obvios, buscar valores numéricos típicos
                if fwhm is None or beta is None:
                    print(f"    Buscando valores numéricos típicos...")
                    for key in header.keys():
                        try:
                            val = float(header[key])
                            # FWHM típico: 0.5 - 5 arcsec
                            if fwhm is None and 0.5 <= val <= 5.0:
                                if 'FWHM' in key.upper() or 'SEEING' in key.upper():
                                    fwhm = val
                                    print(f"    ? FWHM = {fwhm} (keyword: {key})")
                            # Beta típico: 1.5 - 8
                            if beta is None and 1.5 <= val <= 8.0:
                                if 'BETA' in key.upper() or 'MOFFAT' in key.upper():
                                    beta = val
                                    print(f"    ? BETA = {beta} (keyword: {key})")
                        except (ValueError, TypeError):
                            continue
                
                if fwhm is not None and beta is not None:
                    return fwhm, beta
        
        print("    ❌ No se encontraron parámetros PSF completos")
        return None, None
        
    except Exception as e:
        print(f"    ❌ Error leyendo {image_file}: {e}")
        return None, None


def main():
    if len(sys.argv) != 3:
        print("Uso: python generate_psf_from_field.py CAMPO FILTRO")
        print("")
        print("Ejemplos:")
        print("  python generate_psf_from_field.py STRIPE82-0001 F378")
        print("  python generate_psf_from_field.py HYDRA-0001 R")
        print("  python generate_psf_from_field.py SPLUS-n03s27 F395")
        print("")
        print("El script buscará imágenes con patrón: *_FILTRO.fits")
        return
    
    field = sys.argv[1]
    filter_name = sys.argv[2]
    
    print(f"Generando PSF para campo {field}, filtro {filter_name}")
    print("="*60)
    
    # 1. Buscar imágenes SPLUS del filtro especificado
    splus_images = find_splus_images(filter_name)
    
    if not splus_images:
        print(f"❌ ERROR: No se encontraron imágenes con filtro {filter_name}")
        print("Patrones buscados:")
        print(f"  - *_{filter_name}.fits")
        print(f"  - SPLUS-*_{filter_name}.fits")
        print("  - Field_Img/**/*_FILTRO.fits")
        return
    
    # 2. Intentar extraer parámetros PSF de cualquier imagen encontrada
    fwhm, beta = None, None
    used_image = None
    
    for image in splus_images:
        fwhm, beta = extract_psf_from_field_image(image)
        if fwhm is not None and beta is not None:
            used_image = image
            break
    
    if fwhm is None or beta is None:
        print(f"\n❌ ERROR: No se pudieron extraer parámetros PSF")
        print(f"Se probaron {len(splus_images)} imágenes")
        print(f"Valores encontrados: FWHM={fwhm}, BETA={beta}")
        
        # Mostrar headers para debug
        print(f"\n🔍 Para debug, keywords en la primera imagen:")
        if splus_images:
            try:
                with fits.open(splus_images[0]) as hdul:
                    for i, hdu in enumerate(hdul):
                        if hdu.header:
                            keys = [k for k in hdu.header.keys() if k.strip()]
                            print(f"  HDU {i}: {len(keys)} keywords")
                            # Mostrar keywords que podrían contener PSF info
                            psf_keys = [k for k in keys if any(x in k.upper() for x in ['PSF', 'FWHM', 'SEEING', 'BETA', 'MOFFAT'])]
                            if psf_keys:
                                print(f"    PSF-related: {psf_keys}")
            except:
                pass
        return
    
    print(f"\n✅ Parámetros PSF extraídos de: {used_image}")
    print(f"   FWHM: {fwhm}")
    print(f"   BETA: {beta}")
    
    # 3. Generar PSF
    psf_dir = "Field_Img/psf"
    os.makedirs(psf_dir, exist_ok=True)
    psf_filename = f"{psf_dir}/psf_{field}_{filter_name}.fits"
    
    print(f"\n📐 Generando PSF: {psf_filename}")
    
    try:
        make_psf(fwhm, beta, psf_filename)
        
        # Verificar que se creó
        if os.path.exists(psf_filename):
            file_size = os.path.getsize(psf_filename)
            print(f"✅ PSF generada exitosamente!")
            print(f"   Archivo: {psf_filename}")
            print(f"   Tamaño:  {file_size} bytes")
            print(f"   Dimensiones: {2*radius + 1} x {2*radius + 1} pixels")
        else:
            print(f"❌ ERROR: El archivo PSF no se creó")
            
    except Exception as e:
        print(f"❌ ERROR generando PSF: {e}")


if __name__ == "__main__":
    main()
