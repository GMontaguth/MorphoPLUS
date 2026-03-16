import numpy as np
import matplotlib.pyplot as plt
from astropy.table import *
from astropy.io import fits
import splusdata

# Credentials are loaded from config.py — fill in your own username and
# password there. That file is excluded from version control (.gitignore).
from config import SPLUS_USERNAME, SPLUS_PASSWORD
conn = splusdata.connect(SPLUS_USERNAME, SPLUS_PASSWORD)


Fi=Table.read('Field_L17_simple_random_sample_5ene_selec.csv')

for i in range(len(Fi)):
    print(i)
    Filters=np.array(['R','F378','F395','F410','F430','F515','F660','F861','G','I','Z','U'])
    for banda in Filters:
        hdulist = conn.get_cut(Fi['RA'][i], Fi['DEC'][i], 200, banda)
        hdu = hdulist[1].data
        hdr = hdulist[1].header #header
        im = fits.PrimaryHDU(hdu, header=hdr)
        #print(hdr,hdu)
        im.writeto('Field_Fits/Field_%s_%s.fits'%(Fi['ID'][i],banda))
