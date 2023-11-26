# File Xcsao/xplt.com
# August 13, 2007
# By Doug Mink

# Vectors for plotting results in XCPLOT and XCORPLOT

common /xplt/ maxpix, scont, cspec, smspec, smcspec, maxpts4, xlev, fraclev

int maxpix		# Number of pixels allocated in following buffers
pointer scont		# Object spectrum continuum
pointer smspec          # Smoothed object spectrum
pointer cspec           # Continuum-subtracted object spectrum
pointer smcspec         # Smoothed continuum-subtracted object spectrum
int maxpts4		# Number of pixels allocated in the following buffers
pointer xlev		# X coordinate of peak fit pixel to be marked
pointer fraclev		# Fraction of peak for Y of marked peak fit pixel

# Aug 13 2007	New file
