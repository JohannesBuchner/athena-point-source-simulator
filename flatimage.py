import numpy
import astropy.io.fits as pyfits

npix = 9
image = numpy.ones((npix,npix))
image[0,0] = 10
image[-1,0] = 10
image[0,-1] = 10
image[-1,-1] = 10
hdu = pyfits.PrimaryHDU(data = image)
import datetime, time
now = datetime.datetime.fromtimestamp(time.time())
nowstr = now.isoformat()
nowstr = nowstr[:nowstr.rfind('.')]
hdu.header['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
hdu.header['DATE'] = nowstr
hdu.header['HDUCLASS'] = 'HEASARC/SIMPUT'
hdu.header['HDUCLAS1'] = 'IMAGE'
hdu.header['HDUVERS'] = '1.1.0'
hdu.header['CTYPE1'] = 'RA---TAN'
hdu.header['CTYPE2'] = 'DEC--TAN'
hdu.header['CRVAL1'] = 40.2
hdu.header['CRVAL2'] = 12.8
hdu.header['CRPIX1'] = 3
hdu.header['CRPIX2'] = 3
#hdu.header['CD1_1'] = 1
#hdu.header['CD1_2'] = 0
#hdu.header['CD2_1'] = 1
#hdu.header['CD2_2'] = 0
hdu.header['CUNIT1'] = 'deg'
hdu.header['CUNIT2'] = 'deg'
# wfi is 40x40arcmin
# we do 50% overlap
# 40*1.5/60 / npix ~= 0.1 degrees
hdu.header['CDELT1'] = 0.1
hdu.header['CDELT2'] = 0.1
hdu.writeto('flatimage.fits', clobber=True)

