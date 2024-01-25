#Recorta las imagenes es todos lo grupos, y genera las psf de cada campo, ademas crea las imagenes de deteccion y el archivo para correr sextrctor
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import *
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
import splusdata
conn = splusdata.connect('gpardo','gNGC5054') #(usuario,contraseña) ## from splus.cloud
from psf_new import make_psf
from ejecutable import S,size


Datos_S= S.group_by('Groups')
GS=Datos_S.groups.keys #grupos

Datos_S_field= S.group_by('Field')
Fields=Datos_S_field.groups.keys #numero de grupos
Bands=np.array(['R','F378','F395','F410','F430','F515','F660','F861','G','I','Z','U']) #filtros de splus



def recortar(position,size2,grupo,field,hdu,hdr,B): #le ingreso el nombre de la región del stripe82 y la posicion central de cada grupo
	print('recortar')
	hdr['Cut'] ="%s %s / position center group and size"%(position,size2)
	Recortada=Cutout2D(hdu, position, size2)
	#Convertir a fits la imagen how to convertrecortada
	im_Re = fits.PrimaryHDU(Recortada.data, header=hdr)
	im_Re.writeto('Field_Img/%s_%s_%s.fits'%(grupo,field,B))
	im_Re.closed()
	return



def por_campo(field):
	SF=S[S['Field']==field]
	Datos_SF= SF.group_by('Groups')
	GS=Datos_SF.groups.keys #grupos
	for b in Bands:
		hdulist = conn.get_field(field, b)
		hdu = hdulist[1].data
		print("2")
		hdr = hdulist[1].header
		fwhm, beta = (hdr["HIERARCH OAJ PRO FWHMMEAN"],hdr["HIERARCH OAJ PRO FWHMBETA"])
		make_psf(fwhm, beta, "Field_Img/psf/psf_"  + field + "_" + b + ".fits")
		for g in GS['Groups']:#(295,len(Grupos_sub)):#len(G['Groups'])):
			mask = Datos_SF.groups.keys['Groups'] == g 
			GR=Datos_SF.groups[mask]
			position=(np.mean(GR['X']),np.mean(GR['Y'])) #Calcula pla posicion media del grupo
			size2=(size,size)
			recortar(position,size2,g,field,hdu,hdr,b)
		hdulist.close()	
	return
	
def img_det(field,grupo):
	hdulist = fits.open('Field_Img/%s_%s_R.fits'%(grupo,field))
	hdulist1 = fits.open('Field_Img/%s_%s_G.fits'%(grupo,field))
	hdulist2 =fits.open('Field_Img/%s_%s_Z.fits'%(grupo,field))
	hdu = hdulist[0].data + hdulist1[0].data + hdulist2[0].data
	hdr = hdulist[0].header 
	#print(j,k)
	hdr['FILTER'] ="Sum(G,R,Z) / combinacion de filtros det imag"
	#Recortada=Cutout2D(hdu, position, size)
	im_Re = fits.PrimaryHDU(hdu, header=hdr) ##Convertir a fits la imagen how to convertrecortada
	im_Re.writeto('Field_Img/det/det_%s_%s.fits'%(grupo,field),overwrite=True)
	hdulist.close()
	hdulist1.close()
	hdulist2.close()
	return
    

Data=[]
for f in Fields["Field"]:	
	print(f)
	por_campo(f)
	SF=S[S['Field']==f]
	Datos_SF= SF.group_by('Groups')
	GS=Datos_SF.groups.keys #grupos
	for g in GS["Groups"]:
		img_det(f,g)
		Data.append('sex  Field_Img/det/det_%s_%s.fits -c sextopsfex.conf -CATALOG_NAME R.fits -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME ./sextopsfex.param -DETECT_MINAREA 10 -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -FILTER_NAME gauss_5.0_9x9.conv -PHOT_APERTURES 20.55 -SATUR_LEVEL 25000  -MAG_ZEROPOINT 20.833 -GAIN_KEY GAIN -GAIN 652.7072652846 -PIXEL_SCALE 0.55 -SEEING_FWHM 1.38 -BACK_SIZE 900 -BACKPHOTO_TYPE GLOBAL -BACKPHOTO_THICK 24 -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME Field_Img/det/det_%s_%s.seg.fits'%(g,f,g,f))

fic = open("dopsfex_mask.sh", "w")
for line in Data:
	print(line, file=fic)
fic.close()
	#subprocess.	


