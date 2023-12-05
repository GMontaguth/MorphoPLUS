### Este cogigo va en la carpeta de imagenes para crear las imagnes de los cluster members en cada filtro
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import *
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
import splusdata
conn = splusdata.connect(usuario,contraseña) ## from splus.cloud
from segmetation import main,main_gal
from astropy.table import Table,setdiff 
from astropy.io import ascii
from psf_new import make_psf
from table_generation import tables
from ejecutable import c,size,S
from Img_galxgal import galporgal

#from img_galxgal import galporgal


def recortar_mascara(field): #le ingreso el nombre de la región del stripe82 y la posicion central de cada grupo/cluster member
	Filtros=np.array(['R','F378','F395','F410','F430','F515','F660','F861','G','I','Z','U']) #filtros de splus
	#print('Filtros')
	#print(type(F[0]))
	for i in range(len(Filtros)):
		#print(field,Filtros[i])
		hdulist = conn.get_field(field, Filtros[i])
		hdu = hdulist[1].data
		hdr = hdulist[1].header
		fwhm, beta = (hdr["HIERARCH OAJ PRO FWHMMEAN"],hdr["HIERARCH OAJ PRO FWHMBETA"])
		make_psf(fwhm, beta, "Field_Img/psf/psf_"  + field + "_" + Filtros[i] + ".fits")
		for j in range(len(c)):
			for k in range(len(c)):
				position=(c[j],c[k])
				Tablef,Tabled=tables(S,field,position,size)
				if len(Tablef)>0:
					#print(len(Tablef))
					
					#print(j,k)
					#print(position)
					hdr['Cut'] ="%s %s / position center group and size"%(position,size)
					Recortada=Cutout2D(hdu, position, size)
					im_Re = fits.PrimaryHDU(Recortada.data, header=hdr) ##Convertir a fits la imagen how to convertrecortada
					im_Re.writeto('Field_Img/%i_%i_%s_%s.fits'%(position[0],position[1],field,Filtros[i]),overwrite=True)
			
				else:
					a=0
	main(field,c,size,S)
	
	return
    
def gal(field):
	ids=[]
	p=[]
	for j in range(len(c)):
		for k in range(len(c)):
			position=(c[j],c[k])
			Tablef,Tabled=tables(S,field,position,size)
			#print(len(Tablef),len(Tabled))
			if len(Tabled)>0:
				#print('hola')
				ids.append(Tabled['ID'][0])
				p.append(position)
			else:
				a=0
	#print(ids)
	return ids,p
				
S=Table.read('Catalogos/SPLUS_Table.csv')
Datos_S= S.group_by('Field')
Fields=Datos_S.groups.keys 
#Fields= Fields[0:1]

##################### PARA QUE GENERE LAS IMAGENES E INPUTS DE GALAXIAS INDIVIDUALES#################################
IDS_g=[]
p_g=[]
for f in Fields:
	recortar_mascara(f[0])
	I,P=gal(f[0])
	IDS_g.append(I)
	p_g.append(I)


IDS_g = np.concatenate(IDS_g)
p_g = np.concatenate(p_g)
if len(IDS_g)>0:
	g_S = Table([IDS_g,p_g], names=['ID', 'positon'])
	gal_S = join(S, g_S, keys='ID')
	ascii.write(gal_S, 'Catalogos/g_S.csv', format='csv', fast_writer=False,overwrite=True)
	galporgal(gal_S)
	main_gal(gal_S)

