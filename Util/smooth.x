# File rvsao/Util/smooth.x
# January 14, 1993
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1993 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

procedure smooth (fdata, npoints, nsmooth)
 
real	fdata[ARB]		# vector to be smoothed
int 	npoints			# number of points in vector
int 	nsmooth			# number of times to smooth it

real	val0, val1, val2
int	ismooth, ip, np

begin

	if (nsmooth <= 0) return

	do ismooth = 1, nsmooth {
	    val0 = fdata[1]
	    val1 = fdata[1]
	    np = npoints - 1
	    do ip = 1, np {
		val2 = fdata[ip+1]
		fdata[ip] = (val0 + val1 + val1 + val2) / 4.0
		val0 = val1
		val1 = val2
		}
	    val2 = fdata[npoints]
	    fdata[npoints] = (val0 + val1 + val1 + val2) / 4.0
	    }

end
# Dec 16 1991	Change variable name spectrum to fdata
