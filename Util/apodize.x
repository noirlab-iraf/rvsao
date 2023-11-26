# File rvsao/Util/apodize.x
# December 22, 1997
# By Doug Mink
# After J. Tonry (9/29/83) and G. Torres (Jan/1989)

# Copyright(c) 1997 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# APODIZE apodizes the ends of the spectrum with a cosine bell
# The cosine bell starts rising from zero at 1 and falls to zero at NPTS

include "rvsao.h"
include "emv.h"

procedure apodize (n,fdata,wltemp,fraction,title,filename)


int	n		# Total number of points in spectrum
real	fdata[ARB]	# Data array
real	wltemp[ARB]	# Wavelengths for data array
real	fraction	# Fraction to apodize at each end
char	title[ARB]
char	filename[ARB]

double	pi, factor, fact
int	i, nsquash, nsm
char	title_line[SZ_LINE]

include "rvsao.com"
include "emv.com"

begin
	pi = 3.1415926535897932d0

	nsquash = min (nint (float (n) * fraction), n / 2)
	if (nsquash <= 1) return
	fact = pi / double (nsquash - 1)

	do i = 0,  nsquash-1 {
	    factor = .5 * (1. - Cos (fact * double (i)))
	    fdata(1+i) = factor * fdata(1+i)
	    fdata(n-i) = factor * fdata(n-i)
	    }

	if (pltapo) {
	    call strcpy (title,title_line, SZ_LINE)
	    call strcat (filename, title_line, SZ_LINE)
	    nsm = nsmooth
	    call plotspec (n,fdata,title_line, wltemp,
		 "Wavelength in Angstroms",nsm)
	    }

end
# Jul 30 1990	SPP version by Doug Mink
# Dec 16 1991   Use wavelength vector for region being cross-correlated

# Apr 13 1994	Remove unused variables x1 and x2
# May  3 1994	Change arguments for PLOTSPEC call
# Jun 15 1994	Add smoothing argument to PLOTSPEC call
# Aug  3 1994	Change common and header from fquot to rvsao

# Sep 20 1995	Compute factor in double precision

# Dec 22 1997	Add units to wavelength label of graph
