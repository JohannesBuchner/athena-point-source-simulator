import astropy.io.fits as pyfits
import sys, os
import numpy
import matplotlib.pyplot as plt
import scipy.stats, scipy.signal

import joblib
mem = joblib.Memory('.')

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
plt.imshow((hist > 0).transpose(), extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal')
plt.savefig(sys.argv[1] + '.pdf', bbox_inches='tight')
plt.close()

import cv2

def detect_sources(x, y):
	hist, xbins, ybins = numpy.histogram2d(x, y, bins=(range(xlo, xhi+1), range(ylo, yhi+1)))
	hist[hist > 3] = 3
	gray = (hist/hist.max()*255).astype('uint8')
	blurred = cv2.GaussianBlur(gray, (7, 7), 0)
	blurred = ((blurred > 10)*255).astype('uint8')
	kernel = numpy.ones((3,3),numpy.uint8)
	blurred = cv2.morphologyEx(blurred, cv2.MORPH_OPEN, kernel)
	plt.imshow(blurred.transpose(), extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal')
	plt.savefig(sys.argv[1] + '_detect_in.pdf', bbox_inches='tight')
	gray = blurred
	v = numpy.median(gray)
	sigma = 0.33
	lower = int(max(0, (1.0 - sigma) * v))
	upper = int(min(255, (1.0 + sigma) * v))
	edged = cv2.Canny(gray, lower, upper)
	plt.imshow(edged.transpose(), extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal')
	plt.savefig(sys.argv[1] + '_detect.pdf', bbox_inches='tight')
	return []

def detect_sources(x, y):
	hist, xbins, ybins = numpy.histogram2d(x, y, bins=(range(xlo, xhi+1), range(ylo, yhi+1)))
	hist[hist > 3] = 3
	gray = (hist/hist.max()*255).astype('uint8')
	blurred = cv2.GaussianBlur(gray, (7, 7), 0)
	blurred = ((blurred > 10)*255).astype('uint8')
	plt.imshow(blurred, extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal', origin='lower')
	plt.savefig(sys.argv[1] + '_detect.pdf', bbox_inches='tight')
	circles = cv2.HoughCircles(blurred, cv2.cv.CV_HOUGH_GRADIENT, dp=1.2, 
		param1=20, param2=5, 
		minDist=2, minRadius=0, maxRadius=5)

	if circles is not None:
		# convert the (x, y) coordinates and radius of the circles to integers
		circles = circles[0, :] #numpy.round(circles[0, :]).astype("int")
	 	
		# loop over the (x, y) coordinates and radius of the circles
		for (x, y, r) in circles:
			# draw the circle in the output image, then draw a rectangle
			# corresponding to the center of the circle
			#cv2.circle(output, (x, y), r, (0, 255, 0), 4)
			#cv2.rectangle(output, (x - 5, y - 5), (x + 5, y + 5), (0, 128, 255), -1)
			print 'detection', x, y, r
			plt.plot(x, y, 'x', color='r')
	
	#blurred = ((blurred > 10)*255).astype('uint8')
	
	#kernel = numpy.ones((3,3),numpy.uint8)
	#blurred = cv2.morphologyEx(blurred, cv2.MORPH_OPEN, kernel)
	plt.xlim(0, blurred.shape[1])
	plt.ylim(0, blurred.shape[0])
	plt.savefig(sys.argv[1] + '_detect.pdf', bbox_inches='tight')
	plt.close()
	return []

@mem.cache
def compute_sig2(hist, w=5, h=5):
	bkg = numpy.zeros(hist.shape)
	print 'sliding...'
	results = numpy.zeros(hist.shape, dtype=int)
	for i in range(h, hist.shape[0]-h):
		for j in range(w, hist.shape[1]-w):
			k = hist[i-h:i+h+1,j-w:j+w+1].sum()
			results[i,j] = k
	
	print 'testing...'
	for i in range(h, hist.shape[0]-h):
		for j in range(w, hist.shape[1]-w):
			k = results[i,j]
			if not k > 0: continue
			nearbys = numpy.array([results[i+di, j+dj] for di in [-2*h, -h, 0, h, 2*h] for dj in [-2*w, -w, 0, w, 2*w] if i+di > 0 and j+dj > 0 and i+di<hist.shape[0] and j+dj<hist.shape[1] and results[i+di, j+dj] >= 0 and not (di == 0 and dj == 0)])
			nearbys.sort()
			mid = scipy.stats.mstats.mquantiles(nearbys, 0.25)[0]
			nearbys2 = nearbys[nearbys <= 4*(mid+1)]
			if mid > 5:
				kmean = numpy.nanmedian(nearbys2[:-2])
			else:
				kmean = numpy.nanmean(nearbys2[:-2])
			if numpy.isnan(kmean) or kmean > 5:
				print nearbys
			bkg[i,j] = kmean
	return bkg, results


@mem.cache
def compute_bkg(hist, w=5, h=5):
	print 'sliding...'
	results = numpy.zeros(hist.shape)
	for i in range(h, hist.shape[0]-h):
		for j in range(w, hist.shape[1]-w):
			v = hist[i-h:i+h+1,j-w:j+w+1]
			if numpy.isnan(v).all():
				results[i,j] = numpy.nan
			else:
				results[i,j] = numpy.nanmean(v)
	
	print 'testing...'
	results2 = numpy.zeros(hist.shape) - 1
	for i in range(h, hist.shape[0]-h):
		for j in range(w, hist.shape[1]-w):
			nearbys = numpy.array([results[i+di, j+dj] for di in [-2*h, -h, 0, h, 2*h] for dj in [-2*w, -w, 0, w, 2*w] if i+di > 0 and j+dj > 0 and i+di<hist.shape[0] and j+dj<hist.shape[1] and results[i+di, j+dj] >= 0 and not (di == 0 and dj == 0)])
			nearbys.sort()
			mid = scipy.stats.mstats.mquantiles(nearbys, 0.25)[0]
			nearbys2 = nearbys[nearbys <= 4*(mid+1/(w*h))]
			kmean = numpy.nanmean(nearbys2[:-2])
			if numpy.isnan(kmean) or kmean > 5./(w*h):
				print nearbys
			results2[i,j] = kmean
	return results2, results

@mem.cache
def compute_sig(hist):
	w = h = 7
	kernel = numpy.zeros((3*w,3*h))
	kernel[w:-w,h:-h] = 1
	smallkernel = numpy.ones((w,h))
	assert smallkernel.sum() == kernel.sum(), (smallkernel.sum(), kernel.sum())
	print 'sliding...'
	blurred = scipy.signal.convolve(hist, smallkernel, mode='same')
	blurred = blurred.astype(int)
	plt.imshow(blurred, extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal', origin='lower')
	plt.savefig(sys.argv[1] + '_detect_in.pdf', bbox_inches='tight')
	plt.close()
	
	ul = numpy.zeros((3*w,3*h))
	ul[:w,:h] = 1
	um = numpy.zeros((3*w,3*h))
	um[:w,h:-h] = 1
	ur = numpy.zeros((3*w,3*h))
	ur[:w,-h:] = 1
	ll = numpy.zeros((3*w,3*h))
	ll[-w:,:h] = 1
	lm = numpy.zeros((3*w,3*h))
	lm[-w:,h:-h] = 1
	lr = numpy.zeros((3*w,3*h))
	lr[-w:,-h:] = 1
	ml = numpy.zeros((3*w,3*h))
	ml[w:-w,:h] = 1
	mr = numpy.zeros((3*w,3*h))
	mr[w:-w,-h:] = 1

	kernels = [ul, um, ur, ll, lm, lr, ml, mr]
	ul = numpy.zeros((5*w,5*h))
	ul[:w,:h] = 1
	um = numpy.zeros((5*w,5*h))
	um[:w,h:h*2] = 1
	ur = numpy.zeros((5*w,5*h))
	ur[:w,-h:] = 1
	ll = numpy.zeros((5*w,5*h))
	ll[-w:,:h] = 1
	lm = numpy.zeros((5*w,5*h))
	lm[-w:,h:h*2] = 1
	lr = numpy.zeros((5*w,5*h))
	lr[-w:,-h:] = 1
	ml = numpy.zeros((5*w,5*h))
	ml[w:w*2,:h] = 1
	mr = numpy.zeros((5*w,5*h))
	mr[w:w*2,-h:] = 1
	kernels += [ul, um, ur, ll, lm, lr, ml, mr]
	um = numpy.zeros((5*w,5*h))
	um[:w,2*h:h*3] = 1
	lm = numpy.zeros((5*w,5*h))
	lm[-w:,2*h:h*3] = 1
	ml = numpy.zeros((5*w,5*h))
	ml[2*w:w*3,:h] = 1
	mr = numpy.zeros((5*w,5*h))
	mr[2*w:w*3,-h:] = 1
	kernels += [um, lm, ml, mr]
	um = numpy.zeros((5*w,5*h))
	um[:w,3*h:h*4] = 1
	lm = numpy.zeros((5*w,5*h))
	lm[-w:,3*h:h*4] = 1
	ml = numpy.zeros((5*w,5*h))
	ml[3*w:w*4,:h] = 1
	mr = numpy.zeros((5*w,5*h))
	mr[3*w:w*4,-h:] = 1
	kernels += [um, lm, ml, mr]

	print 'nearbys...'
	nearbys = numpy.median([scipy.signal.convolve(hist, k, mode='same') for k in kernels], axis=0)
	nearbys[nearbys == 0] = 0.1
	print nearbys
	plt.imshow(nearbys, extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', 
		aspect='equal', origin='lower', vmin=0, vmax=nearbys.max())
	plt.savefig(sys.argv[1] + '_detect_bkg.pdf', bbox_inches='tight')
	plt.close()
	sig = scipy.stats.poisson.sf(blurred, nearbys)
	return sig

def detect_sources(x, y):
	"""
	Sliding box detection
	----------------------
	the counts in the box are compared to the background,
	which is from a 5x5 grid of boxes (without the middle, source box).
	The background rate is estimated using the median.
	Finally, the poisson survival rate is computed (1-cdf),
	and sources with >1e-6 accepted.
	"""
	hist, xbins, ybins = numpy.histogram2d(x, y, bins=(range(xlo, xhi+1), range(ylo, yhi+1)))
	plt.imshow(hist.transpose(), extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal')
	plt.savefig(sys.argv[1] + '_detect_in.pdf', bbox_inches='tight')
	plt.close()
	#sig = compute_sig(hist)
	#sig = compute_sig2(hist)
	w=7
	h=7
	bkg, blurred = compute_sig2(hist, w=w, h=h)
	plt.imshow(bkg, extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal', origin='lower')
	plt.savefig(sys.argv[1] + '_detect_bkg.pdf', bbox_inches='tight')
	plt.close()
	
	plt.imshow(blurred, extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal', origin='lower')
	plt.savefig(sys.argv[1] + '_detect_in2.pdf', bbox_inches='tight')
	plt.close()
	
	bkg[bkg < 0.1] = 0.1
	sig = scipy.stats.poisson.sf(blurred, bkg)
	print sig.max(), sig.min()
	sig[sig < 1e-300] = 1e-300
	sig = -numpy.log10(sig)
	gray = sig.copy()
	gray[gray < 3] = 3
	gray[gray > 6] = 6
	plt.figure()
	plt.imshow(gray, extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal', origin='lower')
	plt.savefig(sys.argv[1] + '_detect.pdf', bbox_inches='tight')
	plt.close()
	
	# go through pixels, ordered by significance
	hist_with_holes = hist.copy()
	hist_with_holes[sig > 3] = numpy.nan
	_, bkgmap = compute_bkg(hist_with_holes, h=11, w=11)
	hist_with_holes = hist.copy()
	print 'bkgmap zero: %.2f%%' % (bkgmap == 0).mean()
	print bkgmap.max(axis=0), bkgmap.max(axis=1)
	bkgmap[bkgmap == 0] = 0.1/(h*w)
	print bkgmap.min(), bkgmap.max(), bkgmap.max(axis=0), bkgmap.max(axis=1)
	plt.imshow(bkgmap, extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal', origin='lower', vmin=0)
	plt.savefig(sys.argv[1] + '_detect_bkg2.pdf', bbox_inches='tight')
	plt.close()
	Y, X = numpy.meshgrid(numpy.arange(hist.shape[1]), numpy.arange(hist.shape[0]))
	w, h = 4, 4
	xweights = X[:2*w+1,:2*h+1] - w
	yweights = Y[:2*w+1,:2*h+1] - h
	positions = []

	for s, x, y in sorted(zip(sig[sig > 3], X[sig > 3], Y[sig > 3]), reverse=True):
		#print 'source significance:', (x, y), s, sig[x,y]
		s = sig[x,y]
		if s < 3:
			continue
		#print 'considering source at', (x,y),' with significance', s
		for i in range(3):
			#    take all counts within radius 5,
			counts = hist_with_holes[x-w:x+w+1,y-h:y+h+1]
			if counts.sum() == 0:
				break
			if counts.shape != xweights.shape:
				# source at border
				break
			#    compute centroid
			xnext = int(round((xweights * counts).sum() / counts.sum()))
			ynext = int(round((yweights * counts).sum() / counts.sum()))
			if xnext == 0 and ynext == 0:
				break
			#print '   updating centroid', (xnext, ynext), (x+xnext, y+ynext)
			x += xnext
			y += ynext
		# compute significance
		nctsbkg = bkgmap[x,y]
		if nctsbkg == -1:
			print '   can not test here'
			continue
		ncts = counts.sum()
		s = scipy.stats.poisson.sf(ncts, nctsbkg * counts.size)
		print '   significance', s, ncts, nctsbkg * counts.size
		#    if significant:
		if s < 1e-5:
			print '   storing XXXXXXXXX'
			#      * store x,y, and extract spectrum
			#      * set number of counts to zero in that region
			hist_with_holes[x-w:x+w+1,y-h:y+h+1] = 0
			#      * set significance in that region to zero
			sig[x-h:x+h+1,y-h:y+h+1] = 0
			# 
			plt.plot(y, x, 'x', color='r', alpha=0.5)
			plt.text(y, x, ' %d' % ncts, color='r', alpha=0.2)
			positions.append([y, x])
	plt.imshow(hist, extent=[ylo, yhi, xlo, xhi], cmap='gray_r', interpolation='none', aspect='equal', origin='lower')
	plt.savefig(sys.argv[1] + '_detect2.pdf', bbox_inches='tight')
	plt.close()
	return positions

def extract(xcenter, ycenter, file_prefix, **extra_header):
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
		**extra_header)

	bkg = create_pi(data['PI'][mask_bkg], file_prefix+'_bkg.pi', 
		**header)
	src = create_pi(data['PI'][mask_src], file_prefix+'.pi', 
		BACKFILE='obs_bkg.pi', 
		**header)

mask = numpy.logical_and(data['PI'] > 0.5, data['PI'] < 10)
positions = detect_sources(x[mask], y[mask])
numpy.savetxt('detected.txt', positions)
for i, (x, y) in enumerate(positions):
	extract(x, y, file_prefix = 'det_%d' % i, DETX=x, DETY=y)




