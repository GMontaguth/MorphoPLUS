#Genera la mascara y el archvo de imput para galfitm
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import *
from astropy.table import Table
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.table import *
from ejecutable import S,size


Z=Table.read('Catalogos/ZPs.csv')

def mask_ge(field,grupo,X,Y):
	fig = fits.open('Field_Img/det/det_%s_%s.seg.fits'%(grupo,field))
	for i in range(len(X)):
		arr_mask = fig[0].data == fig[0].data[int(Y[i])][int(X[i])]
		fig[0].data[arr_mask] = 0
	arr_mask = fig[0].data >0
	fig[0].data[arr_mask] = 1
	hdr=fig[0].header
	imgF = fits.PrimaryHDU(fig[0].data, header=hdr)
	imgF.writeto('Field_Img/mask/mask_%s_%s.fits'%(grupo,field),overwrite=True)
	return
	
def Median_sky(grupo,Field):
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
        image = CCDData.read('Field_Img/%s_%s_%s.fits'%(grupo,Field,Filtros[i]),unit="adu") #Abre la imagen donde esta el grupo
        mean, median, std = sigma_clipped_stats(image.data, sigma=3.0)
        sk.append(median)
    I='1) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(sk[0],sk[1],sk[2],sk[3],sk[4],sk[5],sk[6],sk[7],sk[8],sk[9],sk[10],sk[11])
    return I
    
    
def arh_galfit(GRf,z,field,size):
	grupo=GRf['Groups'][0] #nombre del grupo
	Filtros=np.array(['U','F378','F395','F410','F430','G','F515','R','F660','I','F861','Z']) #filtros de splus
	# Band labels (CSL of <nbands> labels containing no whitespace)
	# (these must be unique in a case-insensitive manner)
	# (can be omitted if fitting a single band)
	A1='A1) u,J0378,J0395,J0410,J0430,g,J0515,r,J0660,i,J0861,z' ## filtros
	# Band wavelengths (CSL of values)
	# (choice of wavelength units is arbitrary, as long as consistent,
	#  but affects the resulting wavelength-dependence parameters)
	A2='A2) 3536,3770,3940,4094,4292,4751,5133,6258,6614,7690,8611,8831' #longitd de onda de los filtros
	# Output data image block (FITS filename)
	B='B) galfitm_%s_%s.fits'%(grupo,field) 
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
	maskima='Field_Img/mask/mask_%s_%s.fits'%(grupo,field)
	F='F) %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(maskima,maskima,maskima,maskima,maskima,maskima,maskima,maskima,maskima,maskima,maskima,maskima)
	# File with parameter constraints (ASCII file)
	G='G) none'
	xh1=min(np.array(GRf['X'])-np.mean(GRf['X'])+size/2)-30
	xh2=max(np.array(GRf['X'])-np.mean(GRf['X'])+size/2)+30
	yh1=min(np.array(GRf['Y'])-np.mean(GRf['Y'])+size/2)-30
	yh2=max(np.array(GRf['Y'])-np.mean(GRf['Y'])+size/2)+30
	# Image region to fit (xmin xmax ymin ymax)
	H='H) %f %f %f %f'%(0,size,0,size)
	# Size of the convolution box (x y)
	I='I) 150 150' #%(2*(int(max(np.array(GRf['X'])-np.mean(GRf['X'])+size/2)+20)-int(min(np.array(GRf['X'])-np.mean(GRf['X'])+size/2)-20)),2*(int(max(np.array(GRf['Y'])-np.mean(GRf['Y'])+size/2)+20)-int(min(np.array(GRf['Y'])-np.mean(GRf['Y'])+size/2)-20)))
	# Magnitude photometric zeropoint (CSL of <nbands> values)
	J='J) %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s'%(z['ZP_u'][0],z['ZP_J0378'][0],z['ZP_J0395'][0],z['ZP_J0410'][0],z['ZP_J0430'][0],z['ZP_g'][0],z['ZP_J0515'][0],z['ZP_r'][0],z['ZP_J0660'][0],z['ZP_i'][0],z['ZP_J0861'][0],z['ZP_z'][0])
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
		a.append('Field_Img/%s_%s_%s.fits'%(grupo,field,Filtros[j]))# 'CGs_FITS/R_%i_%s_%s_g_%s.fits'%(frame,Field,Filters[j],grupo))
		d.append('Field_Img/PSF/psf_%s_%s.fits'%(field,Filtros[j]))
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
	X=[]
	Y=[]
	for i in range(len(GRf)):
		Data.append('#----------Galaxia %i-----------'%(i))
		Data.append('0) sersic') # Object type
		x1=(GRf['X'][i]-np.mean(GRf['X'])+size/2)
		y1=(GRf['Y'][i]-np.mean(GRf['Y'])+size/2)
		X.append(x1)
		Y.append(y1)
		Data.append('1) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(x1,x1,x1,x1,x1,x1,x1,x1,x1,x1,x1,x1))  # position x
		Data.append('2) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(y1,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1))  # position y
		Data.append('3) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 12'%(GRf['u_auto'][i],GRf['J0378_auto'][i],GRf['J0395_auto'][i],GRf['J0410_auto'][i],GRf['J0430_auto'][i],GRf['g_auto'][i],GRf['J0515_auto'][i],GRf['J0515_auto'][i],GRf['J0660_auto'][i],GRf['i_auto'][i],GRf['J0861_auto'][i],GRf['z_auto'][i])) # total magnitude in each band
		R=GRf['FLUX_RADIUS_50'][i]
		C=5*np.log10(GRf['FLUX_RADIUS_90'][i]/R)
		n=(C/2.77)**(1/0.466) #aproximado (4R. Andrae, K. Jahnke & P. Melchior (2010))
		Data.append('4) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 2'%(R,R,R,R,R,R,R,R,R,R,R,R))# R_e in each band
		Data.append('5) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 2'%(n,n,n,n,n,n,n,n,n,n,n,n)) 
		el=1/GRf['ELONGATION'][i]
		Data.append('9) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(el,el,el,el,el,el,el,el,el,el,el,el))  # axis ratio (b/a) in each band
		if GRf['THETA'][i]>=0:
			Th=(GRf['THETA'][i]-90)
		else:
			Th=(GRf['THETA'][i]+90)
        
		Data.append('10) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(Th,Th,Th,Th,Th,Th,Th,Th,Th,Th,Th,Th)) # position angle (PA), same value in each band
	Data.append('Z) 0')   #  Skip this model in output image?  (yes=1, no=0)
	Data.append('#-------sky----------')
	Data.append('0) sky')
	Data.append(Median_sky(grupo,field)) # sky background       [ADU counts]
	Data.append('2) 0.000      0 ') # dsky/dx (sky gradient in x) 
	Data.append('3) 0.000      0 ') # dsky/dy (sky gradient in y)
	Data.append('Z) 0')   # Skip this model in output image?  (yes=1, no=0)
	# Guarda cada linea de data en un archivo
	fic = open("inputs/galfit_%s_%s.input"%(grupo,field), "w")
	for line in Data:
		print(line, file=fic)
	fic.close()
	return X,Y
    
Datos_S_field= S.group_by('Field')
Fields=Datos_S_field.groups.keys #numero de grupos

for f in Fields['Field']:
	print(size)	
	#por_campo(f)
	SF=S[S['Field']==f]
	Datos_SF= SF.group_by('Groups')
	GS=Datos_SF.groups.keys #grupos
	zp=Z[Z['Field']==f]
	for g in GS['Groups']:
		mask=Datos_SF.groups.keys['Groups'] == g#'cCGs-4007' #Gt['Groups'][i]
		#print(([i],G['Groups'][i]))
		GR=Datos_SF.groups[mask]
		print(GR)
		Da_f= GR.group_by('Field')
		X,Y=arh_galfit(GR,zp,f,size)
		X=np.array(X)
		Y=np.array(Y)
		mask_ge(f,g,X,Y)

