
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.io import ascii
import subprocess
from table_generation import tables




def img_det(field,c,size):
	Filtros=np.array(['R','I','Z'])#,'F378','F395','F410','F430','F515','F660','F861','G','I','Z','U']) #filtros de splus
	#print('Filtros')
	#print(type(F[0]))
	for j in range(len(c)):
		for k in range(len(c)):
			position=(c[j],c[k])
			Tablef,Tabled=tables(S,field,position,size)
			print(len(Tablef),position,field)
			if len(Tablef)>0:
				hdulist = fits.open('Field_Img/%i_%i_%s_R.fits'%(position[0],position[1],field))
				hdulist1 = fits.open('Field_Img/%i_%i_%s_G.fits'%(position[0],position[1],field))
				hdulist2 =fits.open('Field_Img/%i_%i_%s_Z.fits'%(position[0],position[1],field))
				hdu = hdulist[0].data + hdulist1[0].data + hdulist2[0].data
				hdr = hdulist[0].header 
				#print(j,k)
				hdr['FILTER'] ="Sum(G,R,Z) / combinacion de filtros det imag"
				#Recortada=Cutout2D(hdu, position, size)
				im_Re = fits.PrimaryHDU(hdu, header=hdr) ##Convertir a fits la imagen how to convertrecortada
				im_Re.writeto('Field_Img/det/det_%i_%i_%s.fits'%(position[0],position[1],field),overwrite=True)
			else:
				a=0
	return
	

def sextractor(field,c,size):
	Data=[]
	for j in range(len(c)):
		for k in range(len(c)):
			position=(c[j],c[k])
			Tablef,Tabled=tables(S,field,position,size)
			print(len(Tablef),position,field)
			if len(Tablef)>0:
  				#a=('sex  %s -c sextopsfex.conf -CATALOG_NAME R.fits -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME ./sextopsfex.param -DETECT_MINAREA 10 -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -FILTER_NAME gauss_5.0_9x9.conv -PHOT_APERTURES 20.55 -SATUR_LEVEL 25000  -MAG_ZEROPOINT 20.833 -GAIN_KEY GAIN -GAIN 652.7072652846 -PIXEL_SCALE 0.55 -SEEING_FWHM 1.38 -BACK_SIZE 900 -BACKPHOTO_TYPE GLOBAL -BACKPHOTO_THICK 24 -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME %s_%s.seg.fits'%(a[0],b[0],b[1])) 
				Data.append('sex  Field_Img/det/det_%i_%i_%s.fits -c sextopsfex.conf -CATALOG_NAME R.fits -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME ./sextopsfex.param -DETECT_MINAREA 10 -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -FILTER_NAME gauss_5.0_9x9.conv -PHOT_APERTURES 20.55 -SATUR_LEVEL 25000  -MAG_ZEROPOINT 20.833 -GAIN_KEY GAIN -GAIN 652.7072652846 -PIXEL_SCALE 0.55 -SEEING_FWHM 1.38 -BACK_SIZE 900 -BACKPHOTO_TYPE GLOBAL -BACKPHOTO_THICK 24 -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME Field_Img/det/det_%i_%i_%s.seg.fits'%(position[0],position[1],field,position[0],position[1],field))
			else:
				a=0
	fic = open("dopsfex_mask_%s.sh"%(field), "w")
	for line in Data:
    		print(line, file=fic)
	fic.close()
	#subprocess.run(['ls'])
	return
	
def main(field,c,size,S):
	Datos_S= S.group_by('Field')
	Fields=Datos_S.groups.keys 
	Fields= Fields[0:1]
	img_det(field,c,size)
	sextractor(field,c,size)
	return
S=Table.read('Catalogos/SPLUS_Table.csv')	
if __name__ == "__main__":
    main()
############################### crear los inputs de sextrator para analizar galaxi por galaxia##############################


def img_det_gal(Tabla):

	for j in range(len(Tabla)):
		hdulist = fits.open('Field_Img/Gal_%s_R.fits'%(Tabla['ID'][j]))
		hdulist1 = fits.open('Field_Img/Gal_%s_G.fits'%(Tabla['ID'][j]))
		hdulist2 =fits.open('Field_Img/Gal_%s_Z.fits'%(Tabla['ID'][j]))
		hdu = hdulist[0].data + hdulist1[0].data + hdulist2[0].data
		hdr = hdulist[0].header 
		#print(j,k)
		hdr['FILTER'] ="Sum(G,R,Z) / combinacion de filtros det imag"
		#Recortada=Cutout2D(hdu, position, size)
		im_Re = fits.PrimaryHDU(hdu, header=hdr) ##Convertir a fits la imagen how to convertrecortada
		im_Re.writeto('Field_Img/det/det_%s.fits'%(Tabla['ID'][j]),overwrite=True)
	return
	

def sextractor_gal(Tabla):
	Data=[]
	for j in range(len(Tabla)):
		Data.append('sex  Field_Img/det/det_%s.fits -c sextopsfex.conf -CATALOG_NAME R.fits -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME ./sextopsfex.param -DETECT_MINAREA 10 -DETECT_THRESH 1.5 -ANALYSIS_THRESH 1.5 -FILTER_NAME gauss_5.0_9x9.conv -PHOT_APERTURES 20.55 -SATUR_LEVEL 25000  -MAG_ZEROPOINT 20.833 -GAIN_KEY GAIN -GAIN 652.7072652846 -PIXEL_SCALE 0.55 -SEEING_FWHM 1.38 -BACK_SIZE 900 -BACKPHOTO_TYPE GLOBAL -BACKPHOTO_THICK 24 -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME Field_Img/det/det_%s.seg.fits'%(Tabla['ID'][j],Tabla['ID'][j]))
	fic = open("dopsfex_mask_gal.sh", "w")
	for line in Data:
    		print(line, file=fic)
	fic.close()
	#subprocess.run(['ls'])
	return
	
def main_gal(Tabla):
	img_det_gal(Tabla)
	sextractor_gal(Tabla)
	return
	

