# File rvsao/Util/vcombine.x
# September 20, 2011
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1991-2011 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  Combine cross-correlation and emission line velocities using algorithm
#  derived from John Tonry's SAO Nova software by Bill Wyatt

include "rvsao.h"

procedure vcombine (velxc,vxerr,vr,velem,veerr,nlfit,vel,verr,debug)

double	velxc		# Radial velocity based on best template correlation
double	vxerr		# Error in cross-correlation radial velocity
double	vr		# R-value of best cross-correlation
double	velem		# Radial velocity based on emission lines
double	veerr		# Error in emission line radial velocity
int	nlfit		# Number of emission lines in velocity fit
double	vel		# Heliocentric radial velocity of object (returned)
double	verr		# Error in radial velocity (returned)
bool	debug		# True for printed information

double	disperr		# Dispersion error in km/sec
double fc,vdiff
double emerr,emerr2
double xcerr,xcerr2
double err,err2
double rlim, dindef,disperr2

begin
	dindef = INDEFD
	rlim = 15.d0
	disperr2 = 0.d0
#	disperr = 15.d0
#	if (disperr == 0.0 || disperr == dindef)
#	    disperr2 = 225.d0
#	else
#	    disperr2 = disperr * disperr

	if (velem == dindef && velxc == dindef) {
	    vel = dindef
	    verr = dindef
	    return
	    }
	else if (velem == dindef || veerr <= 0. || nlfit == 0) {
	    if (velxc != dindef) {
		vel = velxc
		verr = sqrt ((vxerr*vxerr) + disperr2)
		}
	    else {
		vel = 0.d0
		verr = 0.d0
		}
	    return
	    }
	else if (velxc == dindef || vxerr <= 0.) {
	    if (velem != dindef) {
		vel = velem
		verr = sqrt ((veerr*veerr) + disperr2)
		}
	    else {
		vel = 0.d0
		verr = 0.d0
		}
	    return
	    }
	else {
	    emerr = veerr
	    emerr2 = emerr * emerr
	    xcerr = vxerr
	    xcerr2 = xcerr * xcerr
	    vdiff = velxc - velem
	    fc = vdiff * vdiff / (xcerr2 + emerr2)
	    if (fc > 8) {
		if (vr > rlim) {
		    vel = velxc
		    err = vxerr
		    }
		else if (nlfit < 3 && vr > 4.) {
		    vel = velxc
		    err = vxerr
		    }
		else {
		    vel = velem
		    err = veerr
		    }
		}
	    else {
		if (xcerr2 > 0 && emerr2 > 0) {
		    err2 = 1.d0 / ((1.d0/emerr2) + (1.d0/xcerr2))
		    vel = ((velxc / xcerr2) + (velem / emerr2)) * err2
		    err = sqrt (err2)
		    }
		else if (xcerr2 > 0) {
		    err = xcerr
		    vel = velxc
		    }
		else if (emerr2 > 0) {
		    err = emerr
		    vel = velem
		    }
		}
	    }

	err2 = err * err
	err = sqrt (err2 + disperr2)
	verr = err

	if (debug) {
	    call printf ("VCOMBINE: vel = %.3f +- %.3f km/sec\n")
		call pargd (vel)
		call pargd (verr)
	    call printf ("VCOMBINE: xc vel = %.3f +- %.3f km/sec, R = %.2f\n")
		call pargd (velxc)
		call pargd (xcerr)
		call pargd (vr)
	    call printf ("VCOMBINE: em vel = %.3f +- %.3f km/sec, %d lines\n")
		call pargd (velem)
		call pargd (emerr)
		call pargi (nlfit)
	    }
end

# Dec  5 1991	Set vel = velxc and verr = vxerr if no emission velocity

# Aug 11 1992	Use nlfit from getim.com instead of nfit from emv common
# Nov 24 1992	Avoid divide by zero errors
# Nov 30 1992	Drop emv header file

# Jun 23 1994	Pass velocities as arguments, not in labelled common
# Aug  3 1994	Change header from fquot to rvsao
# Aug  3 1994	Remove quadrature addition of 30 km/sec error to emission error
# Aug 18 1994	Clean up code

# Jul 13 1995	Change R-value for cross-correlation-only from 10 to 15
# Jul 13 1995	Add debugging code
# Jul 19 1995	Deal with INDEF velocities

# Mar 27 1997	Back R-value for cross-correlation-only back to 10
# Mar 27 1997	Back R-value for cross-correlation-only back to 10
# Mar 27 1997	Only take emission alone if R < 4 and more than 2 lines

# Feb 12 1998	If only error is zero, return zero velocity

# Dec  9 2008	Make disperr dispersion error explicit to fix systematic error

# Sep 20 2011	Set disperr2 to zero because dispersion error is not
#		always known.  Instrument errors should still be added
#		to error
