# File rvsao/Emvel/emguess.x
# February 7, 2006
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Bill Wyatt and John Tonry

# Copyright(c) 1997-2006 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# Z Guessing routine - find largest emission line in certain
# regions, assume it's the most likely line for that region, and
# return the implied redshift.

# Returns a guess at 1+Z, or 0 if nothing found. 
# returns 1+Z = -10.0 if error 

include "emv.h"

procedure emguess (fdata,fcont,npts,blue,red,zvel,debug)

real	fdata[ARB]	# Spectrum with continuum removed
real	fcont[ARB]	# Continuum fit to spectrum
int	npts		# The size of the spectrum 
double	blue, red	# Wavelength limits
double	zvel		# Z velocity guess  (returned)
bool	debug		# print diagnostic statements if true

double	wavelength
double	hmax
double	blue_limit	# Spectrum blue limit in angstroms extended by WEXTRA
double	red_limit	# Spectrum red limit in angstroms extended by WEXTRA
int	i,ip,ip1, ip2
int	nlines,maxl, iguess
double	avheight
double	height
double	zval, c0
double	lcenter[30]
double	lheight[30]
double	wcs_p2w(), wcs_w2p()

include "emv.com"

begin
	c0 = 299792.5d0
	if (debug) {
	    call printf ("EMGUESS:  looking for a velocity among %d lines\n")
		call pargi (nsearch)
	    call printf ("EMGUESS:  wavelength range %7.3fa - %7.3fa %d pts \n")
		call pargd (blue)
		call pargd (red)
		call pargi (npts)
	    call printf ("EMGUESS:  line= %d points, min= %4.1f sigma\n")
		call pargi (npfit)
		call pargd (zsig)
	    }

	blue_limit = blue - (red - blue) * WEXTRA
	red_limit  = red  + (red - blue) * WEXTRA
	zvel = 0.0
	maxl = 30

#  Find the biggest line out of each wavelength range 
	hmax   = 0.0

	do iguess = nsearch, 1, -1 {

	    if (bluelim[iguess] > red)
		next
	    if (redlim[iguess] < blue)
		next

	# See if wavelength regions are within the spectrum
	    ip1 = idnint (wcs_w2p (bluelim[iguess]))
	    ip = idnint (wcs_w2p (blue_limit))
	    if (ip > ip1) ip1 = ip
	    if (ip1 < 1) ip1 = 1

	    ip2 = idnint (wcs_w2p (redlim[iguess]))
	    ip = idnint (wcs_w2p (red_limit))
	    if (ip < ip2) ip2 = ip
	    if (ip2 > npts) ip2 = npts

	    if (debug) {
		call printf ("EMGUESS:  looking for line %d = %s %7.2fA\n")
		call pargi (iguess)
		call pargstr (nmsearch[1,iguess])
		call pargd (restwave[iguess])
		call flush (STDOUT)
		}

	# Look for lines in this region
	    call emfind (fdata, fcont, npts, ip1, ip2, npfit, zsig, lcenter,
			 lheight, maxl, nlines, avheight, debug)

	# Select highest peak
	    wavelength = 0.0
	    height = 0.0
	    do i = 1, nlines {
		if (lheight[i] > height) {
		    height = lheight[i]
		    wavelength = wcs_p2w (lcenter[i])
		    }
		}	    

	# Compute z from wavelength
	    if (hmax < height) {
		hmax = height
		zval = (wavelength / restwave[iguess]) - 1.d0
		zvel = zval * c0
		if (debug) {
		    call printf ("EMGUESS: %s %8.3fa: h= %6.2f z= %6.4f = %9.2f km/sec\n")
			call pargstr (nmsearch[1,iguess])
			call pargd (wavelength)
			call pargd (height)
			call pargd (zval)
			call pargd (zvel)
		    }
		}
	    }

	return
end
# Dec  4 1991	Pass through continuum fit to entire spectrum

# Mar 26 1992	Use emission lines read from file for guess
# Mar 26 1992	Check spectrum limits extended by fractional amount

# Jun  2 1993	Move wavelength <-> pixel conversion to subroutines

# Mar 24 1994	Pass sigma limit in labelled common
# Apr 26 1994	Return velocity instead of Z
# May  2 1994	Fix velocity format

# Feb  4 1997	Declare and set C0; it is no longer in emv.com
# Sep 30 1997	Pass NPTS to emfind

# Feb  7 2006	Print z, not 1+z
