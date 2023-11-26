# File rvsao/Xcor/pknum.x
# July 25, 1990
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1991 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#--- Find number of peaks in cross-correlation function

procedure pknum (xcor, k1, k2, npeaks)

real	xcor[ARB]	# Cross-correlation vector
int	k1, k2		# Limits over which to count peaks
int	npeaks		# Number of peaks (returned)

int	k, k0, slope, last_slope

begin

	npeaks = 0
	k0 = k1 + 1
	if (xcor[k0] >= xcor[k1])
	    last_slope = 1
	else
	    last_slope = -1

	for (k = k0; k < k2; k = k+1) {
	    if (xcor[k+1] >= xcor[k])
		slope = 1
	    else
		slope = -1
	    if (slope != last_slope) {
		npeaks = npeaks + 1
		last_slope = slope
		}
	    }
end
