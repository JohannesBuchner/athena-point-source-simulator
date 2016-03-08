import astropy.io.fits as pyfits
import sys, os
import numpy
import matplotlib.pyplot as plt

# open event file
f = pyfits.open(sys.argv[1])
#info = etree.parse(sys.argv[2])
evtheader = f[0].header
arf_filename = os.path.join(os.path.dirname(sys.argv[2]), evtheader['ANCRFILE'])
rmf_filename = os.path.join(os.path.dirname(sys.argv[2]), evtheader['RESPFILE'])
arf = pyfits.open(arf_filename)
rmf = pyfits.open(rmf_filename)
telescope = rmf['MATRIX'].header['TELESCOP']
instrument = rmf['MATRIX'].header['INSTRUME']
#print [e.name for e in arf]
#print [e.name for e in rmf]
#print rmf['EBOUNDS'].data.dtype

def create_pi(events, filename, **header):
	#bins = rmf['EBOUNDS'].data['E_MIN'].tolist() + [rmf['EBOUNDS'].data['E_MAX'][-1]]
	hist, bins = numpy.histogram(events, bins=numpy.arange(len(rmf['EBOUNDS'].data)+1))
	hdus = []
	hdu = pyfits.PrimaryHDU()
	import datetime, time
	now = datetime.datetime.fromtimestamp(time.time())
	nowstr = now.isoformat()
	nowstr = nowstr[:nowstr.rfind('.')]
	hdu.header['CREATOR'] = """Johannes Buchner <johannes.buchner.acad@gmx.com>"""
	hdu.header['DATE'] = nowstr
	hdus.append(hdu)

	counts = numpy.array(list(enumerate(hist)), dtype=[('CHANNEL', '>i2'), ('COUNTS', '>i4')])
	print counts[counts['COUNTS']>0]
	hdu = pyfits.BinTableHDU(data=counts)
	hdu.header['DATE'] = nowstr
	hdu.header['EXTNAME'] = 'SPECTRUM'
	hdu.header['CORRFILE'] = 'none'
	hdu.header['CORRSCAL'] = 1.0
	hdu.header['HDUCLASS'] = 'OGIP'
	hdu.header['HDUVERS'] = '1.2.0'
	hdu.header['LONGSTRN'] = 'OGIP 1.0'
	hdu.header['GROUPING'] = 0
	hdu.header['FILTER'] = 0
	hdu.header['POISSERR'] = True
	hdu.header['SYS_ERR'] = 0
	hdu.header['QUALITY'] = 0
	hdu.header['DETCHANS'] = len(rmf['EBOUNDS'].data)
	hdu.header['HDUCLAS1'] = 'SPECTRUM'
	for k, v in header.iteritems():
		hdu.header[k] = v
	hdu.header['CHANTYPE'] = 'PI'
	hdu.header['HDUCLAS2'] = 'TOTAL'
	hdu.header['HDUCLAS3'] = 'COUNT'
	hdus.append(hdu)
	hdus = pyfits.HDUList(hdus)
	hdus.writeto(filename, clobber=True)

# plot events
data = f['EVENTS'].data
#print data['PI']
x = data['RAWX'] 
y = data['RAWY']
xlo, xhi = x.min(), x.max()
ylo, yhi = y.min(), y.max()
hist, xbins, ybins = numpy.histogram2d(x, y, bins=(range(xlo, xhi+1), range(ylo, yhi+1)))
plt.imshow(hist.transpose(), extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal')
plt.savefig(sys.argv[1] + '.pdf', bbox_inches='tight')
xcenter, ycenter = (255.5, 255.5)
radius1 = 2
radius2 = 4
radius3 = 8
mask_src = (x - xcenter)**2 + (y - ycenter)**2 < radius1**2
mask_bkg2 = (x - xcenter)**2 + (y - ycenter)**2 > radius2**2
mask_bkg3 = (x - xcenter)**2 + (y - ycenter)**2 < radius3**2
mask_bkg = mask_bkg2 * mask_bkg3
header = dict(
	ANCRFILE=arf_filename,
	RESPFILE=rmf_filename,
	TELESCOP=telescope,
	INSTRUME=instrument,
	EXPOSURE=evtheader['TSTOP'] - evtheader['TSTART'],
	AREASCAL=radius1**2 * 1.,
	BACKSCAL=(radius3**2 - radius2**2) * 1. / radius1**2,
	)

bkg = create_pi(data['PI'][mask_bkg], 'obs_bkg.pi', 
	**header)
src = create_pi(data['PI'][mask_src], 'obs.pi', 
	BACKFILE='obs_bkg.pi', 
	**header)

sys.exit(0)

# filter events


bintable


