=====================
Implementation notes
=====================

Source model
--------------

Murray Brightman's sphere model is used here to produce 1e24/cm² sources.

Background model
-----------------

The background spectrum has, when set up with the normalisations for /arcmin² instead of /cm²::

	flux 0.5 2
	Model Flux  0.000901 photons (1.8037e-12 ergs/cm^2/s) range (0.50000 - 2.0000 keV)
	XSPEC12>flux 0.5 10
	 Model Flux 1.6371e-06 photons (2.9029e-15 ergs/cm^2/s) range (0.50000 - 10.000 keV)

The per-square degree flux in the 2-10keV band is therefore (multiplying by 60x60) 1e-11 erg/deg²/s
The per-0.1 square degree flux is (multiplying by 6x6) 1e-13 erg/(0.1deg)²/s


