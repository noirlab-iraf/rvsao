# File rvsao/Xcor/pkwidth.x
# May 22, 1998
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1998 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#--- Find peak and full width at given power of cross-correlation function

include "rvsao.h"

procedure pkwidth (newpts,xcor,xvel,fracpeak,imax, height, width, il, ir)

int	newpts		# Number of points in cross-correlation vector
real	xcor[ARB]	# Cross-correlation vector
real	xvel[ARB]	# Cross-correlation vector velocities
double	fracpeak	# fraction of peak to use for fitting
int	imax		# Bin containing maximum (computed if zero)
double	height		# actual height of peak (returned)
double	width		# actual width of peak (returned)
int	il		# first bin for peak fit (returned)
int	ir		# last bin for peak fit (returned)

double	rr, rl, fracp, dindef
int	i,iwidth,pmin,pmax
real	velmin, velmax

include "rvsao.com"

begin
	dindef = INDEFD
	pmax = 0
	pmin = 0
	if (debug) {
	    call printf ("PKWIDTH: initial peak at %d / %d\n")
		call pargi (imax)
		call pargi (newpts)
	    call flush (STDOUT)
	    }

# Use the input peak position if it is set
	if (imax > 0) {
	    height = xcor[imax]
	    pmin = imax
	    pmax = imax
	    }

# Find the peak in the cross-correlation function if it's not already set.
	else {
	    height = 0.
	    imax = 0;
	    if (minvel == dindef)
		velmin = -100000.0
	    else
		velmin = minvel
	    if (maxvel == dindef)
		velmax = 200000.0
	    else
		velmax = maxvel
	    do i = 1, newpts {
		if ((minvel == dindef || xvel[i] >= velmin) &&
		    (maxvel == dindef || xvel[i] <= velmax)) {
		    if (pmin == 0) pmin = i
		    pmax = i
		    if (xcor[i] >= height) {
			if (i+1 <= newpts && i-1 > 0) {
			    if (xcor[i] >= xcor[i-1] && xcor[i] >= xcor[i+1]) {
				height = xcor[i]
				imax = i
#				if (debug) {
#				    call printf ("PKWIDTH: %d %8.2f %6.3f %6.3f %6.3f\n")
#					call pargi (i)
#					call pargr (xvel[i])
#					call pargr (xcor[i-1])
#					call pargr (xcor[i])
#					call pargr (xcor[i+1])
#				    }
				}
			    }
			}
   		    }
		}
	    }

	if (pmax == 0 && pmin == 0) {
	    rr = 0.d0
	    rl = 0.d0
	    height = 0.d0
	    }

# Set limits and width if number of points to fit is given
	else if (fracpeak > 1) {
	    iwidth = fracpeak * 0.5d0
	    il = imax - iwidth
	    ir = imax + iwidth
	    rl = il
	    rr = ir
	    }

# Set actual level for fraction of peak
	else {
	    if (fracpeak == 0)
		fracp = height * 0.5d0
	    else
		fracp = height * fracpeak

	    if (debug) {
		call printf ("PKWIDTH: fracpeak=%4.2f * %f, search from %d to %d\n")
		    call pargd (fracpeak)
		    call pargd (height)
		    call pargi (pmin)
		    call pargi (pmax)
		}

# Find left fractional power point
	    il = imax
	    while (xcor[il] > fracp && il > 0) {
		il = il - 1
		}
	    if (xcor[il+1] != xcor[il])
		rl = il + ((fracp-xcor[il]) / (xcor[il+1]-xcor[il]))
	    else
		rl = il

#  Find right fractional power point
	    ir = imax
	    while (xcor[ir] > fracp) {
		ir = ir + 1
		}
	    if (xcor[ir-1] != xcor[ir])
		rr = ir - ((fracp - xcor[ir]) / (xcor[ir-1] - xcor[ir]))
	    else
		rr = ir
	    }

	width = rr - rl

	if (debug) {
	    call printf ("PKWIDTH: %8.3f - %8.3f = %6.3f  c= %d, h= %6.4f\n")
	    call pargd (rl)
	    call pargd (rr)
	    call pargd (width)
	    call pargi (imax)
	    call pargd (height)
	    }

end

# Oct  3 1990	Add minpix and maxpix limits
# May 28 1991	Get minshift, maxshift, and pksrch from cl
# Sep 25 1991	Make height and width double precision
# Nov  4 1991	Make fracpeak double

# Feb 11 1992	When fracpeak > 1, treat it as number of points to fit
# Mar 27 1992	Use velocity instead of pixel limits; change argument order
# May 28 1992	Move limited peak search to cursor subroutine where it belongs
# Nov 30 1992	Avoid divide by zero by testing spectrum values for equality
# 		Pass debug as argument; drop inclusion of fquot common
# Dec  1 1992	Correct velocity limits if user puts them in backwards
# Aug  9 1993	Read minvel and maxvel instead of minshift and maxshift

# Mar 22 1995	Don't accept peak if it isn't max of surroundng points
# Jul 13 1995	Accept peak if it is equal to either of the adjacent points
# Aug 18 1995	Move MINVEL and MAXVEL acquisition to XCFIT

# May 12 1998	Return 0 width if no peak has been found
# May 22 1998	Fix bug so preset peak position works
