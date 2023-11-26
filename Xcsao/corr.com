
# correlation record

common/correlation/ templw1,templw2,objrms,ntmpl,nchop

real	templw1		# blue end of spectrum for best fit template
real	templw2		# red end of spectrum for best fit template
real	objrms		# rms of prepared object spectrum
short	ntmpl		# number of templates in record
short	nchop		# number of emission lines chopped from spectrum
data nchop/0/


#  chopped emission line structure (one per line chopped)

common/temchop/ cchop,rchop

short	cchop		# log lambda of line center
short	rchop		# half-width of pixels chopped


#  template correlation structure (one per template)

common/tempcor/ tcenter,theight,twidth,trmsa,trmss,tshft,tpw

real	tcenter		# center of correlation peak (ln lambda)
real	theight		# height of correlation peak
real	twidth		# width (fwhm) of correlation peak
real	trmsa		# antisymmetric rms
real	trmss		# symmetric rms
real	tshft		# template shift (km/sec)
real	tpw		# pixels per log wavelength
