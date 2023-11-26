# File rvsao/Emvel/emsrch.x
# January 17, 2007
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Bill Wyatt and John Tonry

# Copyright(c) 1991-2007 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
#--- Emission Line Finder
#    Return a list of lines found in /emv/ common

include "emv.h"

# Order is important - If the line finder finds a duplicate line,
# it ignores the new one. Therefore, close lines should be in expected
# strength order.         
    

procedure emsrch (smdata,smcont,npts,blue,red,zval,debug)

real	smdata[ARB]	# Spectrum, smoothed in pixel space
real	smcont[ARB]	# Continuum fit to smoothed spectrum
int	npts		# Number of points in the spectrum
double	blue, red	# Wavelengths limits of the lines in the fit
double	zval		# 1+Z
bool	debug		# if true, print out lines found

int	iline, ifound
int	nlines		# Number of emission lines for which to search
double	center

double	blue_limit	# Spectrum blue limit in angstroms extended by WEXTRA
double	red_limit	# Spectrum red limit in angstroms extended by WEXTRA
double	zwave, zbcont, zrcont
double	pxobs0[MAXREF],wlrest0[MAXREF]
double	wlmax
double	wcs_p2w()
int	j, imax, irmax
int	iref[MAXREF]

include "emv.com"
 
begin

	nfound = 0
	nlines = nref
	blue_limit = blue - (red - blue) * WEXTRA
	red_limit  = red  + (red - blue) * WEXTRA
	if (debug) {
	    call printf ("EMSRCH: looking for emission lines at ")
	    call printf ("1+z= %6.4f from %7.3fa to %7.3fa\n")
		call pargd (zval)
		call pargd (blue_limit)
		call pargd (red_limit)
	    }

# Attempt to locate each candidate line
	do iline = 1, nlines {

	# Check that shifted continuum is within blue wavelength limit
	    zbcont = bcont[iline] * zval
	    if (zbcont < blue_limit) {
		if (debug) {
		    call printf ("EMSRCH: continuum %d: %7.3fa < %7.3fa\n")
			call pargi (iline)
			call pargd (zbcont)
			call pargd (blue_limit)
		    }
		next
		}

	# Check that shifted line is within blue wavelength limit
	    zwave  = wlref[iline] * zval
	    if (zwave < blue_limit) {
		if (debug) {
		    call printf ("EMSRCH: line %d: %7.3fa < %7.3fa\n")
			call pargi (iline)
			call pargd (zwave)
			call pargd (blue_limit)
		    }
		next
		}

	# Check that shifted line is within red wavelength limit
	    if (zwave > red_limit) {
		if (debug) {
		    call printf ("EMSRCH: line %d: %7.3fa > %7.3fa\n")
			call pargi (iline)
			call pargd (zwave)
			call pargd (red_limit)
		    }
		next
		}

	# Check that shifted continuum is within red wavelength limit
	    zrcont = rcont[iline] * zval
	    if (zrcont > red_limit) {
		if (debug) {
		    call printf ("EMSRCH: continuum %d: %7.3fa > %7.3fa\n")
			call pargi (iline)
			call pargd (zrcont)
			call pargd (red_limit)
		    }
		next
		}

	    if (debug) {
		call printf ("EMSRCH: line %d %s: %7.3fa   %7.3fa - %7.3fa\n")
		    call pargi (iline)
		    call pargstr (nmref[1,iline])
		    call pargd (zwave)
		    call pargd (zbcont)
		    call pargd (zrcont)
		}

	# Search for this line
	    call gxyline (smdata,smcont,npts,blue_limit,red_limit,
			  zwave,zbcont,zrcont,center,debug)

	    if (center > 0.d0) {

	# if previous line was the same, throw out this line
		if (nfound > 0 && abs (center - pxobs[nfound]) < 1.0) {
		    if (debug)
			call printf ("duplicates previous line -> ignored\n")
		    }

	# otherwise, add center wavelength and pixel to table
		else {
		    nfound = nfound + 1
		    iref[nfound] = iline
		    pxobs[nfound] = center
		    wlrest[nfound] = wlref[iline]
		    wlobs[nfound] = wcs_p2w (pxobs[nfound])
		    call strcpy (nmref[1,iline],nmobs[1,nfound],SZ_ELINE)
		    if (debug) {
			call printf ("EMSRCH: %d: pixel %7.2f = %7.2f %s\n")
			    call pargi (nfound)
			    call pargd (pxobs[nfound])
			    call pargd (wlref[iline])
			    call pargstr (nmref[1,iline])
			}
		    }
		}
	    if (nfound == MAXREF)
		break
	    }

#  Reorder lines by wavelength
	call amovd (wlrest,wlrest0,nfound)
	call amovd (pxobs,pxobs0,nfound)
	do ifound = nfound, 1, -1 {
	    wlmax = 0.d0
	    imax = 0
	    do j = 1, nfound {
		if (wlrest0[j] > wlmax) {
		    wlmax = wlrest0[j]
		    imax = j
		    }
		}
	    if (imax > 0) {
		pxobs[ifound] = pxobs0[imax]
		pxobs0[imax] = 0.d0
		wlrest[ifound] = wlrest0[imax]
		wlrest0[imax] = 0.d0
		wlobs[ifound] = wcs_p2w (pxobs[ifound])
		irmax = iref[imax]
		call strcpy (nmref[1,irmax],nmobs[1,ifound],SZ_ELINE)
		}
#	    if (debug) {
#		call printf ("EMSRCH: %d: pixel %7.2f = %7.2f <- %7.2f %s\n")
#		    call pargi (ifound)
#		    call pargd (pxobs[ifound])
#		    call pargd (wlobs[ifound])
#		    call pargd (wlrest[ifound])
#		    call pargstr (nmobs[1,ifound])
#		}
	    }
	if (debug) {
	    do ifound = 1, nfound {
		call printf ("EMSRCH: %d: pixel %7.2f = %7.2f <- %7.2f %s\n")
		    call pargi (ifound)
		    call pargd (pxobs[ifound])
		    call pargd (wlobs[ifound])
		    call pargd (wlrest[ifound])
		    call pargstr (nmobs[1,ifound])
		}
	    }

# Set pointer from reference list to list of found lines
	do iline = 1, nlines {
	    lfound[iline] = 0
	    do ifound = 1, nfound {
		if (wlrest[ifound] == wlref[iline])
		    lfound[iline] = ifound
		}
	    }

	if (debug) {
	    call printf ("EMSRCH: %d emission lines found\n")
		call pargi (nfound)
	    }
	return
end


# Attempt to locate a line at a specific wavelength, also given a
# region to estimate continuum and noise level from.
# RETURN is 0.0 for no line, pixel center otherwise.

procedure gxyline (fdata, fcont, npts, blue_limit, red_limit,
		    lineval, blue_cont, red_cont, center, debug)

real	fdata[ARB]	# Spectrum data in pixels
real	fcont[ARB]	# Continuum data in pixels
int	npts		# Number of pixels in spectrum
double	blue_limit	# Blue end limit, Angstroms
double	red_limit	# Red end limit
double	lineval		# estimated line center, Angstroms
double	blue_cont	# blue end of  continuum location, Angstroms
double	red_cont	# red end of continuum location, Angstroms
double	center		# Center of line in pixels (RETURNED)
bool	debug		# if true, print out lines found

double	height
#double	wavelength
double	vrms
int	blue_pix, red_pix
int	i, nlines, maxl
double	vave, val
double	wave1, wave2
double	lcenter[30]
double	lheight[30]
double	wcs_w2p()
#double	wcs_p2w()

int     ipl, ipr, ip1, ip2

include "emv.com"

begin
	maxl = 30

# Calculate a rms limit from the continuum region given
	blue_pix  = idnint (wcs_w2p (blue_cont))
	red_pix  = idnint (wcs_w2p (red_cont))
	if (red_pix < blue_pix) {
	    blue_pix  = idnint (wcs_w2p (red_cont))
	    red_pix  = idnint (wcs_w2p (blue_cont))
	    }

  # ERROR - continuum misidentified?
	if (red_pix < blue_pix+4) {
	    center = 0.d0
	    return
	    }

	vave = 0.d0
	do i = -16, 15  {
	    vave = vave + fdata[blue_pix + i]
	    }
	vrms = 0.d0
	do i = blue_pix, red_pix {
	    val  = fdata[i] - (vave / 32.d0)
	    vrms = vrms + (val * val)
	    vave = vave + fdata[i+16] - fdata[i-16]
	    }

	vrms = zsig * sqrt (vrms / (red_pix - blue_pix + 1))

#	if (debug) {
#	    call printf(" RMS of %d pts = %.2f\n")
#		call pargi (red_pix - blue_pix + 1)
#		call pargd (vrms/zsig)
#	    }

	center = wcs_w2p (lineval)

# Keep search within defined spectrum
	wave1 = lineval - wspan
	ip1 = idnint (wcs_w2p (wave1))
	wave2 = lineval + wspan
	ip2 = idnint (wcs_w2p (wave2))
	if (ip2 < ip1) {
	    ip1 = idnint (wcs_w2p (wave2))
	    ip2 = idnint (wcs_w2p (wave1))
	    }

	if (debug) {
	    call printf ("GXYLINE: %.4fA: %.4fA - %.4fA = %d - %d\n")
		call pargd (lineval)
		call pargd (wave1)
		call pargd (wave2)
		call pargi (ip1)
		call pargi (ip2)
	    }

	ipl = idnint (wcs_w2p (blue_limit))
	ipr = idnint (wcs_w2p (red_limit))
	if (ipl > ipr) {
	    ipr = idnint (wcs_w2p (blue_limit))
	    ipl = idnint (wcs_w2p (red_limit))
	    }
	if (ip1 < ipl) ip1 = ipl
	if (ip1 < 1) ip1 = 1
	if (ip2 > ipr) ip2 = ipr
	if (ip2 > npts) ip2 = npts

# Look for lines in this region
	call emfind (fdata,fcont,npts,ip1,ip2,npfit,zsig,lcenter,lheight,
		     maxl,nlines,height,debug)

# Select highest peak
#	wavelength = 0.0
	height = 0.0
	do i = 1, nlines {
	    if (lheight[i] > height) {
		height = lheight[i]
		center = lcenter[i]
#		wavelength = wcs_p2w (lcenter[i])
		}
	    }	    

#	if (debug) {
#	    call printf ("  Line near %.2f, C = %.1f H = %.1f W = %.1f\n")
#		call pargd (lineval)
#		call pargd (center)
#		call pargd (height)
#		call pargd (wavelength)
#	    }

	if (height > vrms) 
	    return
	else
	    center = 0.d0
	    return
end
# Dec  3 1991	Pass continuum fit vector through
# Dec  6 1991	Print observed wavelength in debug summary
# Dec 18 1991	Fix wavelength to pixel calculation

# Mar 23 1992	Pass continuum limits from file read by EMFIT
# Mar 27 1992	Get allowable pixel shift from parameter pspan
# Mar 30 1992	Set gxyline search window in velocity, not pixels
# Nov 24 1992	Dimension arrays MAXREF instead of 12

# Jun  1 1993	Drop width argument in gxyline
# Jun  2 1993	Move wavelength <-> pixel conversion to subroutines

# Mar 23 1994	Pass sigma limit in labelled common
# Jun  9 1994	Reject line unless height >, not >=, vrms
# Aug 22 1994	Consistently use sz_eline when copying lines

# Aug 10 2000	Fix bug in sorting found by Erwin de Blok, ATNF

# Mar 23 2001	Set line limits to work if spectrum is reversed

# Jan 17 2007	Add vector to line index of fit lines to emission line list
