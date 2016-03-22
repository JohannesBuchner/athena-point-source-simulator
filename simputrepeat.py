import astropy.io.fits as pyfits
import sys
import numpy

a = pyfits.open(sys.argv[1])

e = 'SRC_CAT'
# SRC_ID, SPECTRUM, IMAGE, TIMING need to be re-indexed
# for b-sources
center = a[e].data
b = 0.1759
tmax = 6 * 2 * numpy.pi
tmin = numpy.pi * 3 / 4
rmax = 0.2

alldata = [center]
for i, t in enumerate(numpy.linspace(0, tmax, 6*4*4+1)):
	if t < tmin: continue
	next = center.copy()
	next['SRC_ID'] += i
	next['RA']  += rmax * numpy.cos(t) * t/tmax
	next['DEC'] += rmax * numpy.sin(t) * t/tmax
	alldata.append(next)

# when updating data, some header information is dropped by pyfits
orig_header = dict(a[e].header)
a[e].data = numpy.hstack(tuple(alldata))
for k, v in orig_header.iteritems():
	if k not in a[e].header and not k.startswith('NAXIS'):
		a[e].header[k] = v

# save
a.writeto(sys.argv[2], clobber=True)

