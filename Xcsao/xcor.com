# File Xcsao/xcor.com
# September 17, 2002
# By Doug Mink

# Vectors for cross-correlation results in XCFIT

common /xcor/ ntmp, xcor, xvel, shspec, shtemp, wltemp

int ntmp	# Number of templates to cross-correlate from countemp()
pointer xcor	# cross-correlation returned from xcorfit
pointer xvel	# cross-correlation velocities from xcorfit
pointer shspec	# Spectrum pixels in log-lambda for xcorfit
pointer shtemp	# Template pixels in log-lambda for xcorfit
pointer wltemp	# Wavelength vector for template overlap
