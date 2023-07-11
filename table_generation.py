import numpy as np
from astropy.io import ascii
from astropy.table import *


def tables(table,field,position_center,size):
	Datos_A=table.group_by('Field')
	mask= Datos_A.groups.keys['Field'] == field
	N=Datos_A.groups[mask]
	N2=N[(N['X']>position_center[0]-size/2) & (N['X']<=position_center[0]+size/2) & (N['Y']>position_center[1]-size/2) & (N['Y']<=position_center[1]+size/2)] 
	N2['X']=N2['X']-position_center[0]+size/2
	N2['Y']=N2['Y']-position_center[1]+size/2
	mascara = np.logical_or.reduce([ N2['X'] >= (size - 50), N2['X'] <= 50, N2['Y'] >= (size - 50), N2['Y'] <= 50])
	N3= N2[mascara]
	N4=N2[(N2['X']<(size-50)) & (N2['X']>50) & (N2['Y']<(size-50)) & (N2['Y']>50)] #las galaxias donde se hara el ajuste\
	#print(np.shape(N2),np.shape(N4),np.shape(N3))
	return N4,N3
