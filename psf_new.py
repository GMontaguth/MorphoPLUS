import numpy as np
from astropy.io import fits


radius = 10

def make_psf(fwhm, beta, outfile):
    import numpy as np
    from astropy.io import fits as pf
    fwhm = fwhm/0.5
    #beta = beta/0.5 
    alpha = fwhm / (2 * np.sqrt(np.power(2., 1/beta) - 1.))
    r = np.linspace(-radius, radius, 2 * radius + 1)
    print(r)
    X, Y = np.meshgrid(r, r)
    R = np.sqrt(X**2 + Y**2)
    I = (beta - 1.) / (np.pi * alpha**2) * \
        np.power(1. + np.power(R / alpha, 2), -beta)
    hdu = fits.PrimaryHDU(I)
    hdulist = pf.HDUList([hdu])
    hdulist.writeto(outfile)
    return

