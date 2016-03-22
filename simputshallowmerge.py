import astropy.io.fits as pyfits
import sys
import numpy

a = pyfits.open(sys.argv[1])
b = pyfits.open(sys.argv[2])

e = 'SRC_CAT'
# SRC_ID, SPECTRUM, IMAGE, TIMING need to be re-indexed
# for b-sources
b[e].data['SRC_ID'] += a[e].data['SRC_ID'].max() + 1
for k in 'SPECTRUM', 'IMAGE', 'TIMING':
	# modify the indexes in b to refer to its file
	b[e].data[k] = ['%s%s' % (sys.argv[2], v) if v.startswith('[') else v for v in b[e].data[k]]

# when updating data, some header information is dropped by pyfits
orig_header = dict(a[e].header)
a[e].data = numpy.hstack((a[e].data, b[e].data))
for k, v in orig_header.iteritems():
	if k not in a[e].header and not k.startswith('NAXIS'):
		a[e].header[k] = v

# save
a.writeto(sys.argv[3], clobber=True)

