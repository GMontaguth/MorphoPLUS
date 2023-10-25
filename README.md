# MonrphPLUS

Este programa requiere tener python, psfex, sextractor y galfitm.

Para instalar sextrctor y psfex debo seguir los siguientes pasos:

Primero actualaizar todo:
sudo apt-get update

Intalar estas tres librerias:
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

Y escribiendo en la terminal (psfex) me aseguro que este instalado 

Y lo mismo se hace para sextractor y para verificar que esta instalado se de dar sex

Si al intentar instalar SEXtractor y no funciona existe la alternativa de hacer la instalacion via conda:
> conda install -c conda-forge astromatic-source-extractor

La informacion esta en: https://anaconda.org/conda-forge/astromatic-source-extractor

Adicionalmente hay que descargar galafitm en este pagina: https://www.nottingham.ac.uk/astronomy/megamorph/
y Colocarlo dentro de esta carpeta

Luego en la carpeta Catalogos se ingresa la lista de objetos que se quiere ajustar, esta tabla debe tener el nombre "SPLUS_Table.csv". Esta tabla ya bebe estar corregida por la extinci. Aidionalmente la tabla que contiene los zero points de survey es del DR4. Luego hay que modificar en los archivos los parametros size y c, que es el tamaño que tengra cada sub imgen de cada campo ya que GalfitM no es capaz de ajustar sobre el campo completo de splus; c es el centro de cada sub imagen. Esto se debe modificar solo en: ejecutable.py.

Para correr galfitM_Modificado, se corre el program ejecutable.py, este crea un archivo .sh, al cual se le debe dar los permisos con:

> chmod +wrx ejecutable.sh

Luego se ejecuta en bash

> ./ejecutable.sh

El cual se encarga de descargar cada imgen en los dosce filtros, sus PSFs, generar las mascaras y el impunt para GalfitM para cada subcampo. Y por ultimo corre galfitM.

Para generar las tablas con la infomación del ajuste se debe correr el programa leer_outputs.py, que generara tres tablas. En la primera se encuentra lso datos de splus y del ajuste para cada galaxia donde el ajuste el chi2_red es mor a dos, en la seunda las galaxias donde el chi es mayor, y en la tercera se encuentran los Chi2 por subcampo.

[Diagrama de flujo de GalfitM_modificado](https://github.com/GMontaguth/GalfitM_modificado/blob/main/diagrama_flujo.png)

Para las galaxias de campo, es decir para la muestra de control para generar y recortar las imagenes se debe usar el imagenes_field.py y para generar los inputs de galfitm es galfit_inputs_field.py.

Por otra parte la menera como genere y recorte las imagenes de los grupos compactos tambien es diferente, porque es la primera version del codigo.
