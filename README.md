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
