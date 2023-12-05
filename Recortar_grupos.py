### Este cogigo va en la carpeta de imagenes para crear las imagnes de los grupos en cada filtro


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import *
from astropy.io import fits
from astropy.nddata.utils import Cutout2D



Grupos_sub=Table.read('Catalogos/grupos_completos_splus_12-02-2021.csv') #grupos que tienen todas la galaxias reconocidas en el catalogo de splus
S=Table.read('Catalogos/SPLUS_DR3_Grupos_18072021_cor.csv')
Datos_S= S.group_by('Groups')
Grupos=Datos_S.groups.keys #numero de grupos
GS=Grupos_sub['Groups'][195:257]


def recortar(position,size,grupo,field): #le ingreso el nombre de la región del stripe82 y la posicion central de cada grupo
    Filtros=np.array(['R','F378','F395','F410','F430','F515','F660','F861','G','I','Z','U']) #filtros de splus
    #print(type(F[0]))
    for i in range(len(Filtros)):
        #print(Filtros[i])
        hdulist = fits.open('Splus_Fields/%s_%s.fits'%(field,Filtros[i]))
        hdu = hdulist[1].data
        hdr = hdulist[1].header #header
        #hdulist.close()
        hdr['Cut'] ="%s %s / position center group and size"%(position,size)
        Recortada=Cutout2D(hdu, position, size)
        #Convertir a fits la imagen how to convertrecortada
        im_Re = fits.PrimaryHDU(Recortada.data, header=hdr)

        im_Re.writeto('%s_%s_%s.fits'%(grupo,field,Filtros[i]))
    return

for i in range(len(GS)):#(295,len(Grupos_sub)):#len(G['Groups'])):
    #Para selecionar los datos que corresponden a cada grupo
    #print(Grupos_sub['Groups'][i])
    mask = Datos_S.groups.keys['Groups'] == GS[i] 
    GR=Datos_S.groups[mask]
    Da_f= GR.group_by('Field')
    Field=Da_f.groups.keys
    if len(Field)==1:
        position=(np.mean(GR['X']),np.mean(GR['Y'])) #Calcula pla posicion media del grupo
        size=(500,500)
        recortar(position,size,GR['Groups'][0],Field['Field'][0])
    else:
        print(1)
        GR1=GR[GR['Field']==Field['Field'][0]]
        print(2)
        GR2=GR[GR['Field']==Field['Field'][1]]
        print(Field['Field'][0],Field['Field'][1])

        position1=(np.mean(GR1['X']),np.mean(GR1['Y'])) #Calcula pla posicion media del grupo  
        position2=(np.mean(GR2['X']),np.mean(GR2['Y'])) #Calcula pla posicion media del grupo  
        size=(500,500)
        recortar(position1,size,GR['Groups'][0],Field['Field'][0])
        recortar(position2,size,GR['Groups'][0],Field['Field'][1])

