# init SIXTE
# init XSPEC

XMLFILE := ${SIXTE}/share/sixte/instruments/athena/1190mm_wfi_wo_filter/depfet_b_1l_ff_large.xml

all: obs.pi

sphere0708.fits:
	wget http://www.mpe.mpg.de/~mbright/data/sphere0708.fits 

obs.simput: CTsphere.xcm constlightcurve.dat sphere0708.fits
	simputfile RA=40.2  Dec=12.8            XSPECFile="CTsphere.xcm"            LCFile=constlightcurve.dat            MJDREF=50800.0            Emin=2            Emax=10.0            srcFlux=1e-14            Simput="obs.simput"  clobber=yes

obs_repeated.simput: obs.simput
	python simputrepeat.py obs.simput obs_repeated.simput

flatimage.fits: 
	python flatimage.py

extbkg.simput: extbkg.xcm constlightcurve.dat flatimage.fits
	# srcFlux is total, with the image size being (0.1 deg)^2
	# I computed the source spectrum flux in 0.5-10keV and set it here
	simputfile RA=40.2  Dec=12.8            XSPECFile="extbkg.xcm" LCFile=constlightcurve.dat            MJDREF=50800.0            Emin=0.5            Emax=10            srcFlux=1e-13 	           Simput="extbkg.simput"  clobber=yes IMAGE=flatimage.fits 

obs+extbkg.simput: obs_repeated.simput extbkg.simput
	python simputshallowmerge.py obs_repeated.simput extbkg.simput obs+extbkg.simput

obs_events.fits: obs+extbkg.simput
	time runsixt EventList="obs_events.fits" PatternList="obs_pattern.fits"        \
		Mission="ATHENA"         Instrument="WFI"         Mode="1190mm_wfi_wo_filter"         XMLFile=${XMLFILE} \
		Simput="obs+extbkg.simput"         Exposure=10000.         RA=40.2 Dec=12.8         MJDREF=50814.0 clobber=yes

obs.pi obs_bkg.pi: obs_events.fits
	python createspectrum.py obs_events.fits ${SIXTE}/share/sixte/instruments/athena/1190mm_wfi_wo_filter/depfet_b_1l_ff_large.xml

obsplot.ps: obs.pi
	xspec < test.xspec

obsplot.pdf: obsplot.ps
	ps2pdf obsplot.ps

clean:
	rm -f obs.simput obs_events.fits 
	rm -f obs.pi obs_bkg.pi 
	rm -f obsplot.ps obsplot.pdf
