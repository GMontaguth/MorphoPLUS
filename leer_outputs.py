#Este codigo lee los archivos de salidade de GalfitM para cada grupo,y los organiza en una tabla. Luegho esta tabla se une con la tabla de magnitudes de SPLUS
#Luego une cada tabla de cada grupo a una sola tabla, y la guarda como .csv

from astropy.table import *
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import os
import re

def tables(table,field,position_center,size):
	Datos_A=table.group_by('Field')
	mask= Datos_A.groups.keys['Field'] == field
	N=Datos_A.groups[mask]
	N2=N[(N['X']>position_center[0]-275) & (N['X']<=position_center[0]+size/2) & (N['Y']>position_center[1]-size/2) & (N['Y']<=position_center[1]+size/2)]
	N2['X']=N2['X']-position_center[0]+size/2
	N2['Y']=N2['Y']-position_center[1]+size/2
	return N2
 

def leer(position,field,N,GR):
    '''
    Este funcion lee los archivos de salida de galfit .band y los convierte en una tabla
    
    --------
    Input:
    grupo: Nombre del grupo
    field: El field donde se encuentra en grupo en splus
    N: El numero de galaxias en el grupo
    GR: Tabla que contienen toda la informacion de splus para las galaxias del grupo 
    
    -------
    Output:
    Tables2:Es una tabla que contiene los resultados del fit de galfitM y la informacion de splus para cada galaxia en el grupo 
    float(Ch.split('=')[1]): el chi cuadrado del ajuste con GalfitM para el grupo
    '''
    
    path='Galfitm_%i_%i_%s.galfit.01.band'%(position[0],position[1],field) #lee los archivos .band

    log_data=open(path,'r')
    D=[]
    for line in log_data: 
        D.append(line)
    nu,n378,n395,n410,n430,ng,n515,nr,n660,ni,n861,nz=[],[],[],[],[],[],[],[],[],[],[],[]
    Reu,Re378,Re395,Re410,Re430,Reg,Re515,Rer,Re660,Rei,Re861,Rez=[],[],[],[],[],[],[],[],[],[],[],[]
    Xu,X378,X395,X410,X430,Xg,X515,Xr,X660,Xi,X861,Xz=[],[],[],[],[],[],[],[],[],[],[],[]
    Yu,Y378,Y395,Y410,Y430,Yg,Y515,Yr,Y660,Yi,Y861,Yz=[],[],[],[],[],[],[],[],[],[],[],[]
    PAu,PA378,PA395,PA410,PA430,PAg,PA515,PAr,PA660,PAi,PA861,PAz=[],[],[],[],[],[],[],[],[],[],[],[]
    abu,ab378,ab395,ab410,ab430,abg,ab515,abr,ab660,abi,ab861,abz=[],[],[],[],[],[],[],[],[],[],[],[]
    Iu,I378,I395,I410,I430,Ig,I515,Ir,I660,Ii,I861,Iz=[],[],[],[],[],[],[],[],[],[],[],[]
    G=[]
    for i in range(N):
        #grupo_gal
        G.append('%s_%i_%i_%i'%(field,position[0],position[1],i))
        #X
        Xn=D[45+i*14].split(',')
        Xu.append(float(Xn[0].split(' ')[2])) #+np.mean(GR['X'])-250
        #Y
        Yn=D[46+i*14].split(',')
        Yu.append(float(Yn[0].split(' ')[2])) #+np.mean(GR['Y'])-250
        #Integrated magnitude
        In=D[47+i*14].split(',')
        #print(In)
        z=Z[Z['Field']==field]
        Iu.append(-2.5*np.log10(abs(float(In[0].split(' ')[2])/681))+z[0][1])
        I378.append(-2.5*np.log10(abs(float(In[1])))+z[0][2])
        I395.append(-2.5*np.log10(float(In[2]))+z[0][3])
        I410.append(-2.5*np.log10(float(In[3]))+z[0][4])
        I430.append(-2.5*np.log10(float(In[4]))+z[0][5])
        Ig.append(-2.5*np.log10(float(In[5]))+z[0][6])
        I515.append(-2.5*np.log10(float(In[6]))+z[0][7])
        Ir.append(-2.5*np.log10(float(In[7]))+z[0][8])
        I660.append(-2.5*np.log10(float(In[8]))+z[0][9])
        Ii.append(-2.5*np.log10(float(In[9]))+z[0][10])
        I861.append(-2.5*np.log10(float(In[10]))+z[0][11])
        Iz.append(-2.5*np.log10(float(In[11].split(' ')[0]))+z[0][12])
        #  R_e (effective radius)
        Ren=D[48+i*14].split(',')
        Reu.append(float(Ren[0].split(' ')[2]))
        Re378.append(float(Ren[1]))
        Re395.append(float(Ren[2]))
        Re410.append(float(Ren[3]))
        Re430.append(float(Ren[4]))
        Reg.append(float(Ren[5]))
        Re515.append(float(Ren[6]))
        Rer.append(float(Ren[7]))
        Re660.append(float(Ren[8]))
        Rei.append(float(Ren[9]))
        Re861.append(float(Ren[10]))
        Rez.append(float(Ren[11].split(' ')[0]))
        #Sersic index n
        nn=D[49+i*14].split(',')
        nu.append(float(nn[0].split(' ')[2]))
        n378.append(float(nn[1]))
        n395.append(float(nn[2]))
        n410.append(float(nn[3]))
        n430.append(float(nn[4]))
        ng.append(float(nn[5]))
        n515.append(float(nn[6]))
        nr.append(float(nn[7]))
        n660.append(float(nn[8]))
        ni.append(float(nn[9]))
        n861.append(float(nn[10]))
        nz.append(float(nn[11].split(' ')[0]))
        #Axis ratio (b/a)
        abn=D[53+i*14].split(',')
        abu.append(float(abn[0].split(' ')[2]))
        ab378.append(float(abn[1]))
        ab395.append(float(abn[2]))
        ab410.append(float(abn[3]))
        ab430.append(float(abn[4]))
        abg.append(float(abn[5]))
        ab515.append(float(abn[6]))
        abr.append(float(abn[7]))
        ab660.append(float(abn[8]))
        abi.append(float(abn[9]))
        ab861.append(float(abn[10]))
        abz.append(float(abn[11].split(' ')[0]))
        #Position angle (PA) [deg: Up=0, Left=90]
        PAn=D[54+i*14].split(',')
        #print(PAn[0].split(' '))
        PAu.append(float(PAn[0].split(' ')[1]))
        PA378.append(float(PAn[1]))
        PA395.append(float(PAn[2]))
        PA410.append(float(PAn[3]))
        PA430.append(float(PAn[4]))
        PAg.append(float(PAn[5]))
        PA515.append(float(PAn[6]))
        PAr.append(float(PAn[7]))
        PA660.append(float(PAn[8]))
        PAi.append(float(PAn[9]))
        PA861.append(float(PAn[10]))
        PAz.append(float(PAn[11].split(' ')[0]))
    Ch=D[3].split(',')[0]
    Tables= Table([G,Xu,Yu,Iu,I378,I395,I410,I430,Ig,I515,Ir,I660,Ii,I861,Iz,Reu,Re378,Re395,Re410,Re430,Reg,Re515,Rer,Re660,Rei,Re861,Rez,
                  nu,n378,n395,n410,n430,ng,n515,nr,n660,ni,n861,nz,abu,ab378,ab395,ab410,ab430,abg,ab515,abr,ab660,abi,ab861,abz,
                  PAu,PA378,PA395,PA410,PA430,PAg,PA515,PAr,PA660,PAi,PA861,PAz], names=('Grupo_Gal','X','Y','mag_u','mag_J0378','mag_J0395','mag_J0410','mag_J0430','mag_g','mag_J0515','mag_r','mag_J0660','mag_i','mag_J0861','mag_z',
                                                                                        'Re_u','Re_J0378','Re_J0395','Re_J0410','Re_J0430','Re_g','Re_J0515','Re_r','Re_J0660','Re_i','Re_J0861','Re_z',
                                                                                        'n_u','n_J0378','n_J0395','n_J0410','n_J0430','n_g','n_J0515','n_r','n_J0660','n_i','n_J0861','n_z',
                                                                                        'a/b_u','a/b_J0378','a/b_J0395','a/b_J0410','a/b_J0430','a/b_g','a/b_J0515','a/b_r','a/b_J0660','a/b_i','a/b_J0861','a/b_z',
                                                                                     'PA_u','PA_J0378','PA_J0395','PA_J0410','PA_J0430','PA_g','PA_J0515','PA_r','PA_J0660','PA_i','PA_J0861','PA_z'))

    Tables2=hstack([Tables, GR])
    return Tables2, float(Ch.split('=')[1])
    
Bands=np.array(['U','F378','F395','F410','F430','G','F515','R','F660','I','F861','Z'])
 
def grafico(position,field):
	f = fits.open('Galfitm_%i_%i_%s.fits'%(position[0],position[1],field))
	fig, axs = plt.subplots(3, 12,figsize=(20,8),  sharex=True, sharey=True)
	vmax=[0.03,0.03,0.03,.1,.1,0.7,0.5,0.7,0.5,0.7,0.7,0.7]
	for j  in range(12):
		axs[0, j].imshow(f[j].data,cmap='gray', vmin=-vmax[j],vmax=vmax[j])
		axs[0, j].set_title(f'Filter {Bands[j]}')
		axs[1, j].imshow(f[12+j].data,cmap='gray',vmin=-vmax[j],vmax=vmax[j])
		axs[2, j].imshow(f[24+j].data,cmap='gray',vmin=-vmax[j],vmax=vmax[j])
		axs[0,0].text( 60, 140, '5.5\"', fontsize=10,color='blue') 
		axs[0,0].plot([60, 80], [150, 150], 'b-', lw=3)
	plt.savefig('%i_%i_%s.svg'%(position[0],position[1],field),format='svg', dpi=1200)
	#plt.show()
	plt.close()
	return
 
S=Table.read('Catalogos/SPLUS_Table.csv')
Datos_S= S.group_by('Field')
Fields=Datos_S.groups.keys 
Z= Table.read('Catalogos/ZPs.csv')
#Fields=Fields[0:1]

Tabla=[]
ch=[]
GN=[]
F=[]
Tabla2=[]
ch2=[]
GN2=[]
F2=[]
#c=[600,1600,2600,3600,4600,5600,6600,7600,8600,9500]
c=[1375,1925,2475,3025,3575,4125,4675,5225,5775,6325,6875,7425,7975,8525,9075,9625,1175,10725]
size=550
Fields=["SPLUS-n14s32"]
PATH = "."
global output_pattern, template, table
 # # Extract all galaxy names
ff = glob.glob(os.path.join(PATH, "*.01.band"))
C1=[]
C2=[]
Fields2=[]
Name=[]
for f1 in ff:
	filename = os.path.basename(f1)	
	pieces = filename.split("_")
	C1.append(int(pieces[1]))
	C2.append(int(pieces[2]))
	pieces2=pieces[3].split(".")
	Fields2.append(pieces2[0])
	#print(filename,pieces)
	Name.append(filename)



for i in range(len(C1)):
	p=(C1[i],C2[i])
	f=Fields2[i]
	Tablef=tables(S,f,p,size)
	if len(Tablef)>9:
		print(p,len(Tablef),f)
	if len(Tablef)>0:
		T=leer(p,f,len(Tablef),Tablef)
		grafico(p,f)
		f float(T[1])<3:
			Tabla.append(T[0])
			ch.append(T[1])	
			F.append(f)
			GN.append('%i_%i'%(p[0],p[1]))
		else:
			Tabla2.append(T[0])
			ch2.append(T[1])	
			F2.append(f)
			GN2.append('%i_%i'%(p[0],p[1]))
	else:
		no=0
#	
Table_GalfitM=Tabla[0] #Selecciona los datos en la tabla solo los que tienen el grupo 
for i in range(1,len(Tabla)):
    Table_GalfitM=vstack([Table_GalfitM, Tabla[i]])
ascii.write(Table_GalfitM, 'GalfitM_output_A1644.csv', format='csv', fast_writer=False,overwrite=True) #guarda la tabla completa con los resultados de galfitm y splus
Chi_Table= Table([GN,F,ch],names=('Groups','Field','Chi^2/nu'))
ascii.write(Chi_Table, 'Chi2_A1644.csv', format='csv', fast_writer=False,overwrite=True) #guarda el chi cuadrado del ajuste para cada grupo


Table_GalfitM2=Tabla2[0] #Selecciona los datos en la tabla solo los que tienen el grupo 
for i in range(1,len(Tabla2)):
    Table_GalfitM2=vstack([Table_GalfitM2, Tabla2[i]])
ascii.write(Table_GalfitM2, 'GalfitM_output_malos_A1644.csv', format='csv', fast_writer=False,overwrite=True) #guarda la tabla completa con los resultados de galfitm y splus
Chi_Table2= Table([GN2,F2,ch2],names=('Groups','Field','Chi^2/nu'))
ascii.write(Chi_Table2, 'Chi2_malos_A1644.csv', format='csv', fast_writer=False,overwrite=True) #guarda el chi cuadrado del ajuste para cada grupo
