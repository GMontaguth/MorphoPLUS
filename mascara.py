
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.table import *
from ejecutable import c,size
from table_generation import tables

#Esta parte del codigo me crea la mascara apartir de la imagen segmentation generada por sextractor, a la 
#region donde se encuentra las galaxias del grupo les asigna un 0 para no mascararlas y al resto un 1 macarar
#las funtes que no pertenecen al grupo

S=Table.read('Catalogos/SPLUS_Table.csv')

def mask_gal(galaxy):
	f = fits.open('Field_Img/det/det_%s.seg.fits'%(galaxy))
	print(f[0])
	arr_mask = f[0].data == f[0].data[100][100]
	print(arr_mask)
	f[0].data[arr_mask] = 0
	arr_mask = f[0].data >0
	f[0].data[arr_mask] = 1
	hdr=f[0].header
	imgF = fits.PrimaryHDU(f[0].data, header=hdr)
	imgF.writeto('Field_Img/mask/mask_%s.fits'%(galaxy))
	return
	

def mask(field,Table,position):
	f = fits.open('Field_Img/det/det_%i_%i_%s.seg.fits'%(position[0],position[1],field))
	X=np.array(Table['X'])
	Y=np.array(Table['Y'])
	for i in range(len(X)):
		arr_mask = f[0].data == f[0].data[int(Y[i])][int(X[i])]
		f[0].data[arr_mask] = 0
	arr_mask = f[0].data >0
	f[0].data[arr_mask] = 1
	hdr=f[0].header
	imgF = fits.PrimaryHDU(f[0].data, header=hdr)
	imgF.writeto('Field_Img/mask/mask_%i_%i_%s.fits'%(position[0],position[1],field))
	return


def Median_sky(position,field):
    '''
    Funcion que medie el valor de la mediana del cielo para los grupos en cada filtro
    ---------------------------------------------------------------------------------
    Inputs:
    grupo: nombre del grupo(C), 
    Field: el Field 
    
    Output:
    
    Un arreglo (I) de la mediana del cielo en cada filtro
    
    '''

    Filtros=np.array(['U','F378','F395','F410','F430','G','F515','R','F660','I','F861','Z',])
    sk=[]
    for i in range(len(Filtros)):
        #print(F[i])
        image = CCDData.read('Field_Img/%i_%i_%s_%s.fits'%(position[0],position[1],field,Filtros[i]),unit="adu") #Abre la imagen donde esta el grupo
        mean, median, std = sigma_clipped_stats(image.data, sigma=3.0)
        sk.append(mean)
    I='1) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(sk[0],sk[1],sk[2],sk[3],sk[4],sk[5],sk[6],sk[7],sk[8],sk[9],sk[10],sk[11])
    #print(I)
    return I

def arh_galfit(position,Table,Z,field,size):
	#print(field)
	z=Z[Z['Field']==field]
	#print(z)
	Filtros=np.array(['U','F378','F395','F410','F430','G','F515','R','F660','I','F861','Z']) #filtros de splus
	Filtros2=np.array(['u','J0378','J0395','J0410','J0430','g','J0515','r','J0660','i','J0861','z'])
	# Band labels (CSL of <nbands> labels containing no whitespace)
	# (these must be unique in a case-insensitive manner)
	# (can be omitted if fitting a single band)
	A1='A1) u,J0378,J0395,J0410,J0430,g,J0515,r,J0660,i,J0861,z' ## filtros
	# Band wavelengths (CSL of values)
	# (choice of wavelength units is arbitrary, as long as consistent,
	#  but affects the resulting wavelength-dependence parameters)
	A2='A2) 3536,3770,3940,4094,4292,4751,5133,6258,6614,7690,8611,8831' #longitd de onda de los filtros
	# Output data image block (FITS filename)
	B='B) Galfitm_%i_%i_%s.fits'%(position[0],position[1],field)
	# Sigma image name (CSL of <nbands> FITS filenames or "none")
	# (if an individual filename is specified as "none", then that sigma
	#  image will be made from data; if the whole entry consists of just a
	#  single "none", then all sigma images will be made from data.)
	# One can also add a minimum sigma value, such that any galfit-created
	# sigma image will have a minimum of that value times the sky-subtracted
	# input data.
	C='C) none'
	# PSF fine sampling factor relative to data 
	E='E) 1'
	# Bad pixel mask (CSL of <nbands> FITS image or ASCII coord list)
	# (if an individual filename is specified as "none", then a blank
	#  mask will be used; if the whole entry consists of just a single
	#  "none", then all masks will be blank.)
	F='F) Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits,Field_Img/mask/mask_%i_%i_%s.fits'%(position[0],position[1],field,position[0],position[1],field,position[0],position[1],field,position[0],position[1],field,position[0],position[1],field,position[0],position[1],field,position[0],position[1],field,position[0],position[1],field,position[0],position[1],field,position[0],position[1],field,position[0],position[1],field,position[0],position[1],field)
	# File with parameter constraints (ASCII file)
	G='G) none'
	xh1=np.mean(np.array(Table['X']))-(max(np.array(Table['X']))-min(np.array(Table['X'])))-100
	xh2=np.mean(np.array(Table['X']))-(max(np.array(Table['X']))-min(np.array(Table['X'])))+100
	yh1=np.mean(np.array(Table['Y']))-(max(np.array(Table['Y']))-min(np.array(Table['Y'])))-100
	yh2=np.mean(np.array(Table['Y']))-(max(np.array(Table['Y']))-min(np.array(Table['Y'])))+100
	# Image region to fit (xmin xmax ymin ymax)
	H='H) 0 %i 0 %i'%(size,size)
	# Size of the convolution box (x y)
	I='I) %i %i'%(int(size/2),int(size/2))
	# Magnitude photometric zeropoint (CSL of <nbands> values)
	J='J) %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(z["ZP_u"][0],z["ZP_J0378"][0],z["ZP_J0395"][0],z["ZP_J0410"][0],z["ZP_J0430"][0],z["ZP_g"][0],z["ZP_J0515"][0],z["ZP_r"][0],z["ZP_J0660"][0],z["ZP_i"][0],z["ZP_J0861"][0],z["ZP_z"][0])
	# Plate scale (dx dy)   [arcsec per pixel]
	K='K) 0.55 0.55'
	# Display type (regular, curses, both)
	O='O) regular'
	# Options: 0=normal run; 1,2=make model/imgblock & quit
	P='P) 0  '
	# Non-parametric component
	U='U) 0 ' # Standard parametric fitting
	W='W) input,model,residual '
	a=[]
	d=[]
	# Input data images (CSL of FITS filenames)
	# the number of input data images defines <nbands>
	# the order of the bands must be maintained in all multi-band options
	# the first band in the list is the 'reference band'
	for j in range(len(Filtros)):
		a.append('Field_Img/%i_%i_%s_%s.fits'%(position[0],position[1],field,Filtros[j]))# 'CGs_FITS/R_%i_%s_%s_g_%s.fits'%(frame,Field,Filters[j],grupo))
		d.append('Field_Img/psf/psf_%s_%s.fits'%(field,Filtros[j]))
	A='A) %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11])
	# Input PSF image (CSL of <nbands> FITS filenames) 
	# and a single diffusion kernel (FITS filename, # or omitted)
	D='D) %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10],d[11])
	Data=['================================================================================',
	'# GUIDE TO INPUT FILE FOR GALFITM (a product of the MegaMorph project)',
	'Including multi-band fitting, non-parametric component and MultiNest sampling.',
	'CSL = comma separated list (must not contain any whitespace)',
	'Where several lines are given for the same letter code, these are alternative',
	'examples. The behaviour for multiple lines with the same letter is undefined.',
	'================================================================================',
	A,A1,A2,B,C,D,E,F,G,H,I,J,K,O,P,U,W]
	# INITIAL FITTING PARAMETERS
	for i in range(len(Table)):
		Data.append('#----------Galaxia %i-----------'%(i))
		Data.append('0) sersic') # Object type
		x1=(Table['X'][i])
		y1=(Table['Y'][i])
		Data.append('1) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(x1,x1,x1,x1,x1,x1,x1,x1,x1,x1,x1,x1))  # position x
		Data.append('2) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(y1,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1))  # position y
		Data.append('3) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 12'%(Table['u_auto'][i],Table['J0378_auto'][i],Table['J0395_auto'][i],Table['J0410_auto'][i],Table['J0430_auto'][i],Table['g_auto'][i],
		Table['J0515_auto'][i],Table['J0515_auto'][i],Table['J0660_auto'][i],Table['i_auto'][i],Table['J0861_auto'][i],Table['z_auto'][i]))# total magnitude in each band
		R=Table['FLUX_RADIUS_50'][i]
		C=5*np.log10(Table['FLUX_RADIUS_90'][i]/R)
		n=(C/2.77)**(1/0.466) #aproximado (4R. Andrae, K. Jahnke & P. Melchior (2010))
		Data.append('4) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 2'%(R,R,R,R,R,R,R,R,R,R,R,R))# R_e in each band
		Data.append('5) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 2'%(n,n,n,n,n,n,n,n,n,n,n,n))  # Sersic exponent in each band
		el=1/Table['ELONGATION'][i]
		Data.append('9) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 2'%(el,el,el,el,el,el,el,el,el,el,el,el))  # axis ratio (b/a) in each band
		if Table['THETA'][i]>=0:
			Th=(Table['THETA'][i]-90)
		else:
			Th=(Table['THETA'][i]+90)

		Data.append('10) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 2'%(Th,Th,Th,Th,Th,Th,Th,Th,Th,Th,Th,Th)) # position angle (PA), same value in each band
		Data.append('Z) 0')   #  Skip this model in output image?  (yes=1, no=0)
	Data.append('#-------sky----------')
	Data.append('0) sky')
	Data.append(Median_sky(position,field)) # sky background       [ADU counts]
	Data.append('2) 0.000      0 ') # dsky/dx (sky gradient in x) 
	Data.append('3) 0.000      0 ') # dsky/dy (sky gradient in y)
	Data.append('Z) 0')   # Skip this model in output image?  (yes=1, no=0)
	# Guarda cada linea de data en un archivo
	fic = open("galfit_%i_%i_%s.input"%(position[0],position[1],field,), "w")
	for line in Data:
		print(line, file=fic)
	fic.close()
	return
    
    
Z= Table.read('Catalogos/ZPs.csv')

Datos_S= S.group_by('Field')
Fields=Datos_S.groups.keys 
#Fields=["SPLUS-n14s32"]
Fields= Fields[0:1]
print(Fields)
size=1100 #Fue el ultimo que use con el cluster abell 1644, este es el tamano de la imagen
for f in Fields:
	for j in range(len(c)):
		for k in range(len(c)):
			p=(c[j],c[k])
			Tablef, Tabled=tables(S,f[0],p,size)
			if len(Tablef)>0:
				print(f[0])
				mask(f[0],Tablef,p)
				arh_galfit(p,Tablef,Z,f[0],size)
			else:
				no=0

###################################### PARA QUE GENERE LA MASCARA DE LAS GALAXIAS INDIVIDUALES###########################

archivo_csv = "Catalogos/g_S.csv"

try:
    with open(archivo_csv, "r") as archivo:
        lector_csv = csv.reader(archivo)
        # Procesar los datos del archivo CSV aquí
        for fila in lector_csv:
            # Haz algo con cada fila del archivo CSV
            archivo_g_S=1
except FileNotFoundError:
    # Archivo no encontrado, continuar con el código
     archivo_g_S=0
 
if archivo_g_S==1:
	g_S=Table.read('Catalogos/g_S.csv')
	for j in range(len(g_S)):
		mask_gal(g_S['ID'][j])
else:
	nada=0



