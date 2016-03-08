# init SIXTE
# init XSPEC

XMLFILE := ${SIXTE}/share/sixte/instruments/athena/1190mm_wfi_wo_filter/depfet_b_1l_ff_large.xml

all: obs.pi

obs.simput: CTsphere.xcm constlightcurve.dat
	simputfile RA=40.2  Dec=12.8            XSPECFile="CTsphere.xcm"            LCFile=constlightcurve.dat            MJDREF=50800.0            Emin=2            Emax=10.0            srcFlux=1e-11            Simput="obs.simput"  clobber=yes

obs_events.fits: obs.simput
	runsixt EventList="obs_events.fits" PatternList="obs_pattern.fits"        \
		Mission="ATHENA"         Instrument="WFI"         Mode="1190mm_wfi_wo_filter"         XMLFile=${XMLFILE} \
		Simput="obs.simput"         Exposure=10000.         RA=40.2 Dec=12.8         MJDREF=50814.0 clobber=yes

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
