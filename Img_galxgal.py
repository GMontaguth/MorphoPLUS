import numpy as np
import matplotlib.pyplot as plt
from astropy.table import *
from astropy.io import fits
import splusdata
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from astropy.nddata import CCDData
conn = splusdata.connect('gpardo', 'gNGC5054')
Filters=np.array(['R','F378','F395','F410','F430','F515','F660','F861','G','I','Z','U'])

def img(Fi):
	for i in range(len(Fi)):
		for banda in Filters:
			hdulist = conn.get_cut(Fi['RA_1'][i], Fi['DEC_1'][i], 200, banda)
			hdu = hdulist[1].data
			hdr = hdulist[1].header #header
			im = fits.PrimaryHDU(hdu, header=hdr)
			im.writeto('Field_Img/Gal_%s_%s.fits'%(Fi['ID'][i],banda))
	return
              
def Median_sky(galaxy):
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
        image = CCDData.read('Field_Img/Gal_%s_%s.fits'%(galaxy,Filtros[i]),unit="adu") #Abre la imagen donde esta el grupo
        mean, median, std = sigma_clipped_stats(image.data, sigma=3.0)
        sk.append(mean)
    I='1) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(sk[0],sk[1],sk[2],sk[3],sk[4],sk[5],sk[6],sk[7],sk[8],sk[9],sk[10],sk[11])
    #print(I)
    return I
    
def arh_galfit(GRf,Z):
    galaxy=GRf['ID'] #nombre del grupo
    Table=GRf
    Field=GRf['Field']
    print(galaxy)
    z=Z[Z['Field']==Field] #print(z)
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
    B='B) Gal_%s.fits'%(galaxy) 
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
    F='F) Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits,Field_Img/mask/mask_%s.fits'%(galaxy,galaxy,galaxy,galaxy,galaxy,galaxy,galaxy,galaxy,galaxy,galaxy,galaxy,galaxy)
	# File with parameter constraints (ASCII file)
    G='G) none'
    # Image region to fit (xmin xmax ymin ymax)
    H='H) %f %f %f %f'%(31.0,171.0,31.0,171.0)
    # Size of the convolution box (x y)
    I='I) %i %i'%(20,190)
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
    #pp=Fi['Grupo_Gal'][i].split('_')
    #position=(int(pp[1]),int(pp[2]))
    for j in range(len(Filtros)):
        a.append('Field_Img/Gal_%s_%s.fits'%(galaxy,Filtros[j]))# 'CGs_FITS/R_%i_%s_%s_g_%s.fits'%(frame,Field,Filters[j],grupo))
        d.append('Field_Img/psf/psf_%s_%s.fits'%(Field,Filtros[j]))
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
    Data.append('#----------Galaxia -----------')
    Data.append('0) sersic') # Object type
    x1=(100)
    y1=(100)
    Data.append('1) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(x1,x1,x1,x1,x1,x1,x1,x1,x1,x1,x1,x1))  # position x
    Data.append('2) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(y1,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1,y1))  # position y
    Data.append('3) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 12'%(Table['u_auto'],Table['J0378_auto'],Table['J0395_auto'],Table['J0410_auto'],Table['J0430_auto'],Table['g_auto'],
		Table['J0515_auto'],Table['J0515_auto'],Table['J0660_auto'],Table['i_auto'],Table['J0861_auto'],Table['z_auto'])) # total magnitude in each band
    Data.append('4) 7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0 2') # R_e in each band
    Data.append('5) 3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0 2')  # Sersic exponent in each band
    el=1/GRf['ELONGATION']
    Data.append('9) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(el,el,el,el,el,el,el,el,el,el,el,el))  # axis ratio (b/a) in each band
    if GRf['THETA']>=0:
        Th=(GRf['THETA']-90)
    else:
        Th=(GRf['THETA1']+90)
        
    Data.append('10) %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f 1'%(Th,Th,Th,Th,Th,Th,Th,Th,Th,Th,Th,Th)) # position angle (PA), same value in each band
    Data.append('Z) 0')   #  Skip this model in output image?  (yes=1, no=0)
    Data.append('#-------sky----------')
    Data.append('0) sky')
    Data.append(Median_sky(galaxy)) # sky background       [ADU counts]
    Data.append('2) 0.000      0 ') # dsky/dx (sky gradient in x) 
    Data.append('3) 0.000      0 ') # dsky/dy (sky gradient in y)
    Data.append('Z) 0')   # Skip this model in output image?  (yes=1, no=0)
    # Guarda cada linea de data en un archivo
    fic = open("galfit_%s.input"%(galaxy), "w")
    for line in Data:
        print(line, file=fic)
    fic.close()
    return
    
Z= Table.read('Catalogos/ZPs.csv')
              


def galporgal(Fi):
	img(Fi)
	Data=[]
	for i in range(len(Fi)):
		arh_galfit(Fi[i],Z)      
		galaxy=Fi['ID'][i]
		Data.append("chmod 777 galfit_%s.input"%(galaxy))
		Data.append('./galfitm-1.4.4-linux-x86_64 galfit_%s.input'%(galaxy))
	fic = open('ejecutable_gal.sh','w')
	for line in Data:
		print(line, file=fic)
	fic.close()
	return
	

