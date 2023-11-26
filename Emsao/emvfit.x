# File rvsao/Emvel/emvfit.x
# February 15, 2006
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Bill Wyatt and John Tonry

# Copyright(c) 1991-2006 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

# Print out emission line velocities, then return the weighted velocity.

include "emv.h"

# Some constants for the acceptance criteria
define	WMIN  0.4	# line must be broader than this * mean comp. width
define	WMAX  1.7	# line must be narrower than this * mean comp. width
define	SMIN  2.0	# EQW/DEQW must be greater than this
define	PIXMAX 20	# Must no closer than this to ends of spectrum

procedure emvfit (blue, red, drms, meanwidth, verbose, zav, zrms, zchi, nz)

double	blue		# Blue wavelength limit of comparison
double	red		# Red wavelength limit of comparison
double	drms		# RMS of comparison line fit, Angstroms
double	meanwidth	# Mean comparison line width
bool	verbose		# if true, print lines out as well as return velocity
double	zav		# Weighted average 1+Z (returned)
double	zrms		# RMS of 1+Z (returned)
double	zchi		# Chi-square (returned)
int	nz		# Number of lines in zvel (returned)
double	cz,czerr	# "velocity" and error
double  verr		# Line velocity error
double  verr2		# Square of line velocity error
double  weight		# Sum of squared velocity errors
double	weight2
double	lwmin		# Minimum fraction of mean line width
double	lwmax		# Maximum fraction of mean line width
double	lsmin		# EQW/DEQW must be greater than this
double	ave_vel
double	disp
double	c0
double	zvel, zvel2
double	fiterr, fiterr2
double	wv_cen, wv_cen_plusone
double	dwdp
double	rms2		# Square of root mean square of residuals
bool	lineok[MAXREF]
int	ifound
double	wcs_p2w()
char	reason[SZ_LINE]
char	tstring[SZ_LINE]
real	clgetr()
int	clscan()

include "emv.com"

begin
	weight  = 0.d0
	ave_vel = 0.d0
	disp    = 0.d0
	nz = 0
	zav = 0.d0
	zrms = 0.d0
	c0 = 299792.5d0
	if (clscan ("lwmin") == EOF)
	    lwmin = WMIN
	else
	    lwmin = clgetr ("lwmin")
	if (clscan ("lwmax") == EOF)
	    lwmax = WMAX
	else
	    lwmax = clgetr ("lwmax")
	if (clscan ("lwmax") == EOF)
	    lwmax = WMAX
	else
	    lwmax = clgetr ("lwmax")
	if (clscan ("lsmin") == EOF)
	    lsmin = SMIN
	else
	    lsmin = clgetr ("lsmin")

	if (verbose) {
	    if (nfound == 1)
		call printf("\n EMVFIT: 1 line found\n")
	    else {
		call printf("\n EMVFIT: %d lines found\n")
		    call pargi (nfound)
		}
	    }

	if (nfound <= 0)
	   return

	if (verbose) {
	    call printf("\n Rest lam  Obs. lam     Pixel    vel      dvel     ")
	    call printf("EQW     DEQW  Chi**2/free\n")
	    }
	do ifound = 1, nfound {
	    verr2 = 0.d0
	    verr = 0.d0
	    linedrop[ifound] = 0

	# Acceptance criterion
	    if (emparams[4,1,ifound] > 0.0) {
		lineok[ifound] = TRUE
		reason[1] = EOS
		}
	    else {
		call sprintf (tstring,SZ_LINE," Too narrow: %6.3f <= 0")
		    call pargd (emparams[4,1,ifound])
		lineok[ifound] = FALSE
		linedrop[ifound] = 1
		call strcat (tstring, reason, SZ_LINE)
		}
	    if (lineok[ifound]) {
		wv_cen = wcs_p2w (emparams[4,1,ifound])
		wv_cen_plusone = wcs_p2w (emparams[4,1,ifound]+1.d0)
		wlobs[ifound] = wv_cen
		dwdp  = wv_cen_plusone - wv_cen
		if (wv_cen == wv_cen_plusone) {
		    wv_cen_plusone = wcs_p2w (emparams[4,1,ifound]-1.d0)
		    dwdp  = wv_cen - wv_cen_plusone
		    }

	# Tests for acceptance
	    # Width of Gaussian comp. to comparison
		if (meanwidth > 0.0) {
		    if ((emparams[6,1,ifound] * 2.d0*dwdp) < (lwmin * meanwidth)) {
			lineok[ifound] = FALSE
			linedrop[ifound] = 2
			call sprintf (tstring,SZ_LINE," Too narrow: %7.3f < %7.2f")
			    call pargd (emparams[6,1,ifound] * 2.d0*dwdp)
			    call pargd (lwmin * meanwidth)
			call strcat (tstring, reason, SZ_LINE)
			}
		    if ((emparams[6,1,ifound] * 2.d0*dwdp) > (lwmax * meanwidth)) {
			lineok[ifound] = FALSE
			linedrop[ifound] = 3
			call sprintf (tstring,SZ_LINE," Too wide: %7.3f < %7.2f")
			    call pargd (emparams[6,1,ifound] * 2.d0*dwdp)
			    call pargd (lwmax * meanwidth)
			call strcat (tstring, reason, SZ_LINE)
			}
		    }

		if (emparams[7,2,ifound] != 0.d0) {
		    if (lsmin > (emparams[7,1,ifound] / emparams[7,2,ifound])) {
			lineok[ifound] = FALSE
			linedrop[ifound] = 4
			call sprintf (tstring,SZ_LINE,"EQW/DEQW %6.3f / %6.3f < %6.3f")
			    call pargd (emparams[7,1,ifound])
			    call pargd (emparams[7,2,ifound])
			    call pargr (lsmin)
			call strcat (tstring, reason, SZ_LINE)
			}
		    }

	    # Too close to blue end of spectrum
		if (wv_cen - blue  < PIXMAX * dwdp) {
		    lineok[ifound] = FALSE
		    linedrop[ifound] = 5
		    call sprintf (tstring,SZ_LINE," Too blue: %7.2f-%7.2f < %7.2f")
			call pargd (wv_cen)
			call pargd (blue)
			call pargd (PIXMAX * dwdp)
		    call strcat (tstring, reason, SZ_LINE)
		    }

	    # Too close to red end of spectrum
		if (red - wv_cen < PIXMAX * dwdp) {
		    lineok[ifound] = FALSE
		    linedrop[ifound] = 6
		    call sprintf (tstring,SZ_LINE," Too red: %7.2f-%7.2f < %7.2f")
			call pargd (wv_cen)
			call pargd (red)
			call pargd (PIXMAX * dwdp)
		    call strcat (tstring, reason, SZ_LINE)
		    }

	    # Cursor mode overrides
		if (override[ifound] > 0) {
		    lineok[ifound] = TRUE
		    call strcpy (" Cursor keep",tstring,SZ_LINE)
		    call strcat (tstring, reason, SZ_LINE)
		    }
		if (override[ifound] < 0) {
		    lineok[ifound] = FALSE
		    call strcpy (" Cursor drop",tstring,SZ_LINE)
		    call strcat (tstring, reason, SZ_LINE)
		    }

	# Note this velocity does not include the HCV or BCV
		zvel  = (wlobs[ifound] / wlrest[ifound]) - 1.d0
		cz = zvel * c0
		emparams[9,1,ifound] = zvel

	# Formal error of center => error in velocity from Gaussian fit
		verr = 2 * emparams[4,2,ifound] * dwdp / wv_cen
		if (verr == 0.d0) {
		    lineok[ifound] = FALSE
		    linedrop[ifound] = 7
		    }
		verr2 = verr * verr

	# Add component of error due to fitting residuals
		fiterr = drms / wv_cen
		fiterr2 = fiterr * fiterr
		verr2 = verr2 + fiterr2
		emparams[9,2,ifound] = sqrt (verr2)
		czerr = sqrt (verr2) * c0

		if (lineok[ifound]) {
		    weight = weight + (1.0 / verr2)
		    ave_vel = ave_vel + (zvel / verr2)
		    nz = nz + 1
		    emparams[10,1,ifound] = sqrt (verr2)
		    }
		else
		    emparams[10,1,ifound] = 0.
		if (emparams[8,2,ifound] < 1.d0)
		    emparams[8,2,ifound] = -1.d0
		if (verbose) {
		    call printf("%9.1f %9.1f %9.1f %8.1f %7.1f %7.3f %7.3f %7.3f")
			call pargd (wlrest[ifound])
			call pargd (wv_cen)
			call pargd (emparams[4,1,ifound])
			call pargd (cz)
			call pargd (czerr)
			call pargd (emparams[7,1,ifound])
			call pargd (emparams[7,2,ifound])
			call pargd (emparams[8,1,ifound] / emparams[8,2,ifound])
			if (lineok[ifound]) {
			    call printf (" _")
			    if (override[ifound] > 0) {
				call printf ("%s")
				    call pargstr (reason)
				}
			    call printf ("\n")
			    }
			else {
			    call printf(" X %s\n")
				call pargstr (reason)
			    }
			call flush (STDOUT)
		    }
		}
	    else if (verbose) { 		# fit not done
		call printf("%10.1f%10s%10.1f\n")
		    call pargd (wlrest[ifound])
		    call pargstr (" ")
		    call pargd (pxobs[ifound])
		call flush (STDOUT)
		}
	    }

	if (nz > 0) {
	    ave_vel = ave_vel / weight
	    rms2 = 0.d0
	    weight2 = 1.d0 / weight
	    weight = sqrt (weight2)
	    do ifound = 1, nfound {
	    
		if (lineok[ifound]) {
		    wv_cen = wcs_p2w (emparams[4,1,ifound])
		    wv_cen_plusone = wcs_p2w (emparams[4,1,ifound]+1.d0)
		
		# Note this velocity does not include the HCV or BCV
		    zvel  = (wv_cen / wlrest[ifound]) - 1.d0 - ave_vel
		
		# Formal err of center => err in velocity from Gaussian fit
		    verr = 2 * emparams[4,2,ifound] * (wv_cen_plusone - wv_cen) / wv_cen
		    verr2 = verr * verr
		    verr2 = verr2 + fiterr2
		    if (zvel == 0.d0)
			zvel2 = 0.d0
		    else
			zvel2 = zvel * zvel
		    if (verr2 != weight2)
			disp = disp + (zvel2 / (verr2 - weight2))
		    rms2 = rms2 + zvel2
		    }
		}
	    disp = disp / double (nz)
	    rms2 = rms2 / double (nz)
	    fiterr2 = 0.d0
	    }
	if (nz == 1) {
#	    disp = emparams[8,1,ifound]
	    disp = 0.d0
	    rms2 = verr2

# Component of error due to fitting residuals
	    fiterr = (drms * 2.0) / (blue + red) # in terms of Z
	    fiterr2 = fiterr * fiterr
	    }

	zav  = ave_vel + 1.d0
#	zrms = sqrt (weight2 + fiterr2)
	zrms = sqrt (rms2 + fiterr2)
	zchi  = disp

	if (verbose) {
	    call printf("\n Average z = %.3f or %.1f +/- %.2f km/sec, chi**2 = %.2f\n")
		call pargd (ave_vel)
		call pargd (ave_vel * c0)
		call pargd (zrms * c0)
		call pargd (zchi)
	    }

	return
end
# Dec 18 1991	Fix wavelength calculation by subtracting 1 from pixel

# Mar 25 1992	Change GLEN to MAXREF in dimension of lineok and CEE to c0
# Dec  2 1992	Rename squared velocity error instead of overloading verr

# Jun  1 1993	Set lineok using if statements
# Jun  2 1993	Move wavelength <-> pixel conversion to subroutines

# May  2 1994	Clean up loop code
# May 18 1994	Drop line if fit center pixel is negative
# Jun  9 1994	Comment out all uses of fiterr as it is always 0
# Aug 15 1994	Save reasons for dropping lines
# Aug 15 1994	Fix bugs in computation of equivalent width and its error
# Aug 15 1994	Bring back fiterr to include dispersion error

# Jan 24 1996	Use RMS for error if >1 line found
# Feb  4 1997	Declare and set C0; it is no longer in emv.com
# Feb 19 1997	Add LWMIN, LWMAX, and LSMIN as parameters

# Feb 12 1998	If single line fit, pass its error as velocity error

# Feb 15 2006	Fix for Linux
