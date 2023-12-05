from astropy.table import *
from astropy.io import ascii
from table_generation import tables
import numpy as np
import csv

S=Table.read('Catalogos/SPLUS_Table.csv')
#archivo_csv = "Catalogos/g_S.csv"

Datos_S= S.group_by('Field')
Fields=Datos_S.groups.keys 
Fields= Fields[0:1]
#Fields=Fields[0:1]
Data=['python Recortar.py']
Data2=[]
c=[550,1650,2750,3850,4950,6050,7150,8250,9350,10450]
size=1100
def ejecutable(Fields):
	for f in Fields:
		Data.append('chmod 777 dopsfex_mask_%s.sh'%(f[0]))
		Data.append('./dopsfex_mask_%s.sh'%(f[0]))
	#Data.append('cd inputs')
	archivo_gal=[]
	for f in Fields:	
		for j in range(len(c)):
			for k in range(len(c)):
				position=(c[j],c[k])
				mask= Datos_S.groups.keys['Field'] == f[0]
				N=Datos_S.groups[mask]
				Tablef,Tabled=tables(S,f[0],position,size)
				archivo_gal.append(len(Tabled))
				if len(Tablef)>0:
					Data.append('chmod 777 galfit_%i_%i_%s.input'%(position[0],position[1],f[0]))
					Data2.append('./galfitm-1.4.4-linux-x86_64 galfit_%i_%i_%s.input'%(position[0],position[1],f[0]))
				else:
					no=0
	if sum(archivo_gal)>0:
		Data.append('chmod 777 dopsfex_mask_gal.sh') #GAL INDIVIDUALES
		Data.append('./dopsfex_mask_gal.sh')
		Data.append('python mascara.py')
		for h in Data2:
			Data.append(h)
		Data.append('chmod 777 ejecutable_gal.sh')
		Data.append('./ejecutable_gal.sh')
	else:
		Data.append('python mascara.py')
		for h in Data2:
			Data.append(h)
	Data.append('python leer_header_output.py')
	fic = open('ejecutable.sh','w')
	for line in Data:
		print(line, file=fic)
	fic.close()
	return


ejecutable(Fields)
