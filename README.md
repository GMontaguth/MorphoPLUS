# MorphoPLUS

Este programa requiere tener Python, PSFEx, SExtractor y GalfitM.

Para instalar SExtractor y PSFEx, sigue los siguientes pasos:

Primero, actualiza todo:

Instala estas tres bibliotecas:
Atlas 3.6 en Ubuntu : sudo apt-get install libatlas-base-dev libblas-dev liblapack-dev
FFTw3.3.8 : sudo apt-get install libfftw3-dev
PLPlot 5.9 : sudo apt-get update, sudo apt-get install libplplot-dev libplplot_data

Luego descargar Psfex y sextractor de los repositorios:

Psfex: git clone https://github.com/astromatic/psfex.git
Sextractor: git clone https://github.com/astromatic/sextractor.git 

Para instalar psf:

- cd psfex
- sh autogen.sh
- ./configure
- make -j
- sudo make install

Verifica la instalación escribiendo en la terminal (psfex).

Haz lo mismo para SExtractor y verifica que esté instalado escribiendo sex en la terminal.

Si al intentar instalar SExtractor no funciona, existe la alternativa de hacer la instalación vía conda:
> conda install -c conda-forge astromatic-source-extractor

La información está en: https://anaconda.org/conda-forge/astromatic-source-extractor

Adicionalmente, descarga GalfitM desde esta página: https://www.nottingham.ac.uk/astronomy/megamorph/ y colócalo dentro de esta carpeta.

Se deben crear dos carpetas vacias, la primera es Field_img, que debe contener adentro tres carpetas vacias: det, mask, psf. La segunda carpeta  vacia dentro de MorphPLUS es Out_img.

Luego, en la carpeta "Catalogos", ingresa la lista de objetos que deseas ajustar. Esta tabla debe tener el nombre "SPLUS_Table.csv" y ya debe estar corregida por la extinción. Además, la tabla que contiene los zero points del survey es del DR4. Luego, debes modificar en los archivos los parámetros size y c, que es el tamaño que tendrá cada subimagen de cada campo, ya que GalfitM no es capaz de ajustar sobre el campo completo de SPLUS; c es el centro de cada subimagen. Esto se debe modificar solo en: ejecutable.py.

Para ejecutar GalfitM_Modificado, corre el programa ejecutable.py. Este crea un archivo .sh, al cual debes darle los permisos con:

> chmod +wrx ejecutable.sh

Luego, ejecútalo en bash:

> ./ejecutable.sh

Este script se encarga de descargar cada imagen en los doce filtros, sus PSFs, generar las máscaras y el punto de impunt para GalfitM para cada subcampo. Y por último, corre GalfitM.

Para generar las tablas con la información del ajuste, debes ejecutar el programa leer_outputs.py, que generará tres tablas. En la primera se encuentran los datos de SPLUS y del ajuste para cada galaxia donde el ajuste y el chi2_red es mayor a dos, en la segunda las galaxias donde el chi es mayor, y en la tercera se encuentran los Chi2 por subcampo.

Diagrama de flujo de GalfitM_modificado

Para las galaxias de campo, es decir, para la muestra de control, para generar y recortar las imágenes, debes usar el script imagenes_field.py, y para generar los inputs de GalfitM, usa galfit_inputs_field.py.

Por otra parte, la manera en que generé y recorté las imágenes de los grupos compactos también es diferente, porque es la primera versión del código.

Para descargar las tablas de SPLUS en el dr4: 
SELECT det.ID, det.Field ,det.RA, det.DEC, det.X, det.Y, det.ELONGATION, det.ELLIPTICITY, det.THETA, det.A, det.B, det.FLUX_RADIUS_50, det.FLUX_RADIUS_90, u.u_petro, u.e_u_petro, u.SEX_FLAGS_u, J0378.J0378_petro, J0378.e_J0378_petro, J0378.SEX_FLAGS_J0378, J0395.J0395_petro, J0395.e_J0395_petro, J0395.SEX_FLAGS_J0395, J0410.J0410_petro, J0410.e_J0410_petro, J0410.SEX_FLAGS_J0410, J0430.J0430_petro, J0430.e_J0430_petro, J0430.SEX_FLAGS_J0430, g.g_petro, g.e_g_petro, g.SEX_FLAGS_g, J0515.J0515_petro, J0515.e_J0515_petro, J0515.SEX_FLAGS_J0515, r.r_petro, r.e_r_petro, r.SEX_FLAGS_r, J0660.J0660_petro, J0660.e_J0660_petro, J0660.SEX_FLAGS_J0660, i.i_petro, i.e_i_petro, i.SEX_FLAGS_i, J0861.J0861_petro, J0861.e_J0861_petro, J0861.SEX_FLAGS_J0861, z.z_petro, z.e_z_petro, z.SEX_FLAGS_z, pz.zml, pz.odds FROM idr4_dual.idr4_detection_image as det JOIN idr4_dual.idr4_dual_u as u ON (u.ID = det.ID) JOIN idr4_dual.idr4_dual_J0378 as J0378 ON (J0378.ID = det.ID) JOIN idr4_dual.idr4_dual_J0395 as J0395 ON (J0395.ID = det.ID) JOIN idr4_dual.idr4_dual_J0410 as J0410 ON (J0410.ID = det.ID) JOIN idr4_dual.idr4_dual_J0430 as J0430 ON (J0430.ID = det.ID) JOIN idr4_dual.idr4_dual_g as g ON (g.ID = det.ID) JOIN idr4_dual.idr4_dual_J0515 as J0515 ON (J0515.ID = det.ID) JOIN idr4_dual.idr4_dual_r as r ON (r.ID = det.ID) JOIN idr4_dual.idr4_dual_J0660 as J0660 ON (J0660.ID = det.ID) JOIN idr4_dual.idr4_dual_i as i ON (i.ID = det.ID) JOIN idr4_dual.idr4_dual_J0861 as J0861 ON (J0861.ID = det.ID) JOIN idr4_dual.idr4_dual_z as z ON (z.ID = det.ID) JOIN idr4_vacs.idr4_photoz as pz ON (pz.ID = det.ID) WHERE 1 = CONTAINS( POINT('ICRS', det.RA, det.DEC), CIRCLE('ICRS', ra, dec, r_g))

Donde ra, dec, r_g son los parametros a cambiar, ra y dec son las posiciones del centro y r_g es es el radio en grados donde se buscaran todas las fuentes

