# File rvsao/Util/subcon.x
# August 3, 1994
# By Doug Mink, Harvard-SMithsonian Center for Astrophysics
# After John Tonry and Guillermo Torres

# Copyright(c) 1994 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# Fit and subtract the continuum from a spectrum

include "rvsao.h"

procedure subcon (n,fdata,wltemp,curve,title,filename)

int	n		# Number of points in spectrum
real	fdata[ARB]	# Spectrum
real	wltemp[ARB]	# Wavelengths for spectrum
real	curve[ARB]	# Curve values subtracted (returned)
char	title[ARB]	# Title for plot
char	filename[ARB]	# Name of data file

real	lpoly(), scale[2]
double	coeff[10]
int	i, j, ncoeff
char	title_line[SZ_LINE]

include "rvsao.com"

begin

	do i = 1, n {
	    curve(I) = i
	    }

	scale[1] = curve[1]
	scale[2] = curve[n]
	ncoeff = contdim + 1

	call fitlpo (n,curve,fdata,scale,ncoeff,coeff)

	for (j = 1; j <= n; j = j + 1) {
	    curve[j] = lpoly (curve[j],scale,ncoeff,coeff)
	    fdata[j] = fdata[j] - curve[j]
	    }

	if (pltcon) {
	    call strcpy (title,title_line,SZ_LINE)
	    call strcat (filename, title_line, SZ_LINE)
	    call plotspec (n, fdata, title, wltemp, "Wavelength")
	    }

end
# Sep 29 1983	Written by John Tonry
# Jan    1989	Modified to run under IRAF by Guillermo Torres
# Jul 30 1990	Translated to SPP by Doug Mink
# Jul 31 1990	Modified to use fraction rather than percent
# Dec 16 1991	Use wavelength vector for region being cross-correlated

# Aug  3 1994	Change common and header from fquot to rvsao
