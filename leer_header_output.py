from astropy.table import Table
from astropy.io import ascii, fits
import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib import colors

# -----------------------------
# Config
# -----------------------------
CAT_IN   = "Catalogos/g_S.csv"
OUT_CSV  = "Catalogos/GalfitM_output.csv"

FITS_FMT = "Gal_{id}.fits"
BAND_FMT = "Gal_{id}.galfit.01.band"

Fil_name = np.array(['R','J0378','J0395','J0410','J0430','J0515','J0660','J0861','G','I','Z','U'])
Bands = Fil_name  # tu grafico_i usa Bands[j]

OUT_IMG_DIR = "Out_img"
MAKE_IMAGES = True
os.makedirs(OUT_IMG_DIR, exist_ok=True)

# -----------------------------
# Helpers
# -----------------------------
def correct_image_orientation(img2d):
    # Poné tu versión real si corresponde
    return img2d


def get_chi2nu(band_path):
    """
    Lee Chi^2/nu desde el archivo .band.
    Devuelve float o np.nan si no lo encuentra.
    """
    try:
        with open(band_path, "r") as f:
            for line in f:
                if "Chi^2/nu" in line:
                    return float(line.split("=")[1].split(",")[0].strip())
    except (FileNotFoundError, OSError):
        return np.nan
    return np.nan


def leer_header(NGAL, HEADER, ID, Field_ID, data):
    """
    Agrega parámetros por cada filtro desde el header.
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

        data.append(xc[0]); data.append(xc[2])
        data.append(yc[0]); data.append(yc[2])
        data.append(R[0]);  data.append(R[2])
        data.append(m[0]);  data.append(m[2])
        data.append(N[0]);  data.append(N[2])
        data.append(ar[0]); data.append(ar[2])
        data.append(PA[0]); data.append(PA[2])

    return data


def build_header_names():
    header_names = ["CHI2NU", "ID", "Field_ID"]
    for i in range(12):
        header_names += [
            f"XC_{Fil_name[i]}",  f"e_XC_{Fil_name[i]}",
            f"YC_{Fil_name[i]}",  f"e_YC_{Fil_name[i]}",
            f"RE_{Fil_name[i]}",  f"e_RE_{Fil_name[i]}",
            f"MAG_{Fil_name[i]}", f"e_MAG_{Fil_name[i]}",
            f"n_{Fil_name[i]}",   f"e_n_{Fil_name[i]}",
            f"AR_{Fil_name[i]}",  f"e_AR_{Fil_name[i]}",
            f"PA_{Fil_name[i]}",  f"e_PA_{Fil_name[i]}",
        ]
    return header_names


# -----------------------------
# Tu función (tal cual, solo cambié el path a OUT_IMG_DIR)
# -----------------------------
def grafico_i(f, position, field, Tabled, r):
    fig, axs = plt.subplots(3, 12, figsize=(20, 8), sharex=True, sharey=True)

    vmax = [5.0]*12
    vmin = [-1.0]*12
    lt   = [0.003,0.003,0.003,0.003,0.003,0.04,0.04,0.01,0.004,0.01,0.08,0.05]
    ls   = [0.003,0.003,0.003,0.003,0.003,0.04,0.04,0.01,0.004,0.01,0.08,0.05]
    npix = 50

    for j in range(12):
        fin = f[j].data
        fic_n = correct_image_orientation(fin)

        x = 100
        y = 100

        axs[0, j].imshow(
            fic_n[int(y-npix):int(y+npix), int(x-npix):int(x+npix)],
            cmap='gray',
            norm=colors.SymLogNorm(linthresh=lt[j], linscale=ls[j], vmin=vmin[j], vmax=vmax[j])
        )
        axs[0, j].set_title(f'Filter {Bands[j]}', fontsize=8)

        fin_2 = f[12 + j].data
        fic_2n = correct_image_orientation(fin_2)
        axs[1, j].imshow(
            fic_2n[int(y-npix):int(y+npix), int(x-npix):int(x+npix)],
            cmap='gray',
            norm=colors.SymLogNorm(linthresh=lt[j], linscale=ls[j], vmin=vmin[j], vmax=vmax[j])
        )

        fin_3 = f[24 + j].data
        fic_3n = correct_image_orientation(fin_3)
        axs[2, j].imshow(
            fic_3n[int(y-npix):int(y+npix), int(x-npix):int(x+npix)],
            cmap='gray',
            norm=colors.SymLogNorm(linthresh=lt[j], linscale=ls[j], vmin=vmin[j], vmax=vmax[j])
        )

        for rr in range(3):
            axs[rr, j].set_xticks([])
            axs[rr, j].set_yticks([])

    out_svg = os.path.join(OUT_IMG_DIR, f"{int(position[0])}_{int(position[1])}_{field}.svg")
    plt.tight_layout()
    plt.savefig(out_svg, format='svg', dpi=1200)
    plt.close()
    return out_svg


# -----------------------------
# Main
# -----------------------------
galcat = Table.read(CAT_IN)

ID_COL = "ID"
FIELD_COL = "Field_ID" if "Field_ID" in galcat.colnames else None

Tabla = []

for k, row in enumerate(galcat):
    gid = str(row[ID_COL])
    fits_path = FITS_FMT.format(id=gid)
    band_path = BAND_FMT.format(id=gid)

    field_id = str(row[FIELD_COL]) if FIELD_COL is not None else gid

    # solo para el nombre del archivo (tu grafico_i no usa position para centrar)
    position = (k, k)

    fi = None
    try:
        fi = fits.open(fits_path)

        # CSV: header como venías
        hdr = fi[12].header
        chi2_nu = get_chi2nu(band_path)
        data = [chi2_nu]
        Tabla.append(leer_header(1, hdr, gid, field_id, data))

        # IMGs: llamar tu función
        if MAKE_IMAGES:
            if len(fi) < 36:
                print(f"[SKIP IMG] ID={gid} | FITS tiene {len(fi)} HDUs, pero grafico_i necesita >= 36.")
            else:
                out_svg = grafico_i(fi, position, field_id, galcat, k)
                print(f"[IMG] guardada: {out_svg}")

        fi.close()

    except (FileNotFoundError, OSError, IndexError, KeyError) as e:
        print(f"[SKIP] ID={gid} | problema leyendo {fits_path} o header: {e}")
        try:
            if fi is not None:
                fi.close()
        except Exception:
            pass
        continue


header_names = build_header_names()
table_out = Table(rows=Tabla, names=header_names)

ascii.write(table_out, OUT_CSV, format="csv", fast_writer=False, overwrite=True)
print(f"Listo: guardé {len(table_out)} filas en {OUT_CSV}")

