# File rvsao/Emvel/emlfit.x
# October 19, 2004
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Bill Wyatt and John Tonry

# Copyright(c) 2004 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# Fit Gaussians to emission lines, up to three at a time
 
include "emv.h"

procedure emlfit (fdata, fcont, sky, npix, debug)

real	fdata[ARB]	# Spectrum array with continuum removed
real	fcont[ARB]	# Continuum array 
real	sky[ARB]	# Sky spectrum (if non-NULL) 
int	npix		# number of pixels 
bool	debug		# print diagnostic messages if true

int	i, j, k, m
int	lp, rp
int	ierr		# MINI verbosity flag and error return 
double	wlfit		# Wavelength to fit on either side of peak(s)
int	iscombo()
int	pixleft, pixright, ngauss, nparams, npoly
double	gcenter[3]	# centers for multiple gaussians
double	gheight[3]	# heights for multiple gaussians
double	gwidth[3]	# widths for multiple gaussians
double	fitvec[20]	# working space for 3 Gaussians + 5 cont. params 
double	covar[20] 
double	cont, area, skylevel
double	wv_cen, wv_cen_plusone
double	px_cen, px_cen_plusone
double	dwdp, pi, s2pi
int	jj, ic, icenter, n, nit
int	ngfit		# Number of Gaussians to fit
bool	igauss[3]	# Flags for found lines in combination
real	hheight		# Half-height of gaussian, cfd, open()
double	height, sigma, relht, relsig
double	wlj
double	xi2, xi2dof	# Chi-squared and chi-squared per degree of freedom
double	wcs_p2w(), wcs_w2p()
double	meanback, meanpix, diffback, diffpix
double	eqw
char	gch

include "emv.com"

# Fit results are returned in the eparams array in emv.com
# emparams[1-3,1,line] are the continuum polynomial coefficients
#                      No midpoint is subtracted from the independent variable
# emparams[4-6,1,line] are the fit Gaussian center (in pixels),
#                      height (in counts), and half-width (in pixels)
# emparams[7,1,line] is the equivalent width
# emparams[8,1,line] is the Chi**2,
# emparams[8,2,line] is the number of degrees of freedom.
# emparams[9,1,line] is the velocity in km/sec, not computed in this routine
# emparams[1-9,2,line] are the errors in the above

begin
	ierr = 0
	pi = 3.14159265358979323846d0
	s2pi = sqrt (2.d0 * pi)
	npoly = nlcont

# for each line or set of lines, fit the Gaussians 
	i = 1
	m = 0
	if (debug) {
	    call printf ("EMLFIT: %d lines found\n")
		call pargi (nfound)
	    }
	if (nfound <= 0)
	    return
	while (i <= nfound) {

	# See if the line is part of a combination 
	    ngauss = iscombo (wlrest, i, nfound, ic, wlfit, ngfit, igauss)

	# Set number of terms to fit
	    npoly = 2
	    nparams = (ngfit * 3) + npoly

	    if (debug) {
		call printf ("EMLFIT: line %d %8.3fa:  ngauss= %d/%d  np= %d\n")
		    call pargi (i)
		    call pargd (wlrest[i])
		    call pargi (ngauss)
		    call pargi (ngfit)
		    call pargi (nparams)
		}

	# Set pixel limits to the fit 
	    pixleft = idnint (wcs_w2p (wlobs[i] - wlfit))
	    pixright = idnint (wcs_w2p (wlobs[i+ngauss-1] + wlfit))
	    if (pixright < pixleft) {
		pixright = idnint (wcs_w2p (wlobs[i] - wlfit))
		pixleft = idnint (wcs_w2p (wlobs[i+ngauss-1] + wlfit))
		}
	    if (pixleft < 1) pixleft = 1
	    if (pixright > npix) pixright = npix
	    if (debug) {
		call printf ("EMLFIT: fitting pixels %d to %d\n")
		    call pargi (pixleft)
		    call pargi (pixright)
		}

	# Set centers, heights, and widths of found gaussians to fit
	    n = i - 1
	    do j = 1, ngfit {
		if (igauss[j]) {
		    n = n + 1
		    gcenter[j] = pxobs[n]
		    gheight[j] = fdata[int (gcenter[j])]
		    hheight = 0.5 * gheight[j]

		# Find approx. FWHM of line 
		    lp = idnint (gcenter[j])
		    do k = lp - 1, pixleft {
			lp = lp - 1
			if (fdata[lp] < hheight) 
			    break
			}
		    rp = idnint (gcenter[j])
		    do k = rp + 1, pixright {
			rp = rp + 1
			if (fdata[rp] < hheight) 
			    break
			}
		    gwidth[j] = double (rp - lp - 1)
		    if (gwidth[j] < 1.d0)  {
			gwidth[j] = 1.d0
			gch = '*'
			}
		    else
			gch = ' '
		    if (debug) {
			call printf ("EMLFIT: fit %d %d/%d %8g %8g %8g%c\n")
			    call pargi (i)
			    call pargi (j)
			    call pargi (ngfit)
			    call pargd (gcenter[j])
			    call pargd (gheight[j])
			    call pargd (gwidth[j])
			    call pargc (gch)
			}
		    }
		}

	# Set centers, heights, and widths of merged gaussians to fit
	    do j = 1, ngfit {
		if (!igauss[j]) {
		    if (edrop[j,ic]) {
			gheight[j] = 0.d0
			gwidth[j] = 0.d0
			gcenter[j] = 0.d0
			}
		    else {
			do k = 1, ngfit {
			    if (igauss[k]) {
				gheight[j] = gheight[k] * eht[j,ic] / eht[k,ic]
				gwidth[j] = gwidth[k]
				wlj =  elines[j,ic] * wlobs[i] / wlrest[i]
				gcenter[j] = wcs_w2p (wlj)
				if (debug) {
				    call printf ("EMLFIT: fit %d %d/%d %8g %8g %8g using %d\n")
				    call pargi (i)
				    call pargi (j)
				    call pargi (ngfit)
				    call pargd (gcenter[j])
				    call pargd (gheight[j])
				    call pargd (gwidth[j])
				    call pargi (k)
				    }
				break
				}
			    }
			}
		    }
		}

	    # Drop third line in combination if not found 
	    if (ngfit > 2 && edrop[3,ic] && gcenter[3] == 0.d0)
		ngfit = 2

	    # Drop second line in combination if not found 
	    if (ngfit > 1 && edrop[2,ic] && gcenter[2] == 0.d0) {
		if (ngfit == 2)
		    ngfit = 1
		else {
		    ngfit = 2
		    gcenter[2] = gcenter[3]
		    gheight[2] = gheight[3]
		    gwidth[2] = gwidth[3]
		    igauss[2] = igauss[3]
		    }
		}

	    # Drop first line in combination if not found 
	    if (ngfit > 1 && edrop[1,ic] && gcenter[1] == 0.d0) {
		if (ngfit == 2) {
		    ngfit = 1
		    gcenter[1] = gcenter[2]
		    gheight[1] = gheight[2]
		    gwidth[1] = gwidth[2]
		    igauss[1] = igauss[2]
		    }
		else {
		    ngfit = 2
		    gcenter[1] = gcenter[2]
		    gheight[1] = gheight[2]
		    gwidth[1] = gwidth[2]
		    igauss[1] = igauss[2]
		    gcenter[2] = gcenter[3]
		    gheight[2] = gheight[3]
		    gwidth[2] = gwidth[3]
		    igauss[2] = igauss[3]
		    }
		}
	    nparams = (ngfit * 3) + npoly

	# Fit multiple gaussians
	    ierr = 0
	    do j = 1, 11 {
		fitvec[j] = 0.d0
		}
	    diffback = double (fdata[pixright] - fdata[pixleft])
	    diffpix = double (pixright - pixleft)
	    fitvec[2] = diffback / diffpix
	    meanback = 0.5d0 * double (fdata[pixleft] + fdata[pixright])
	    meanpix = 0.5d0 * double (pixright + pixleft)
	    fitvec[1] = meanback - (meanpix * fitvec[2])
	    fitvec[3] = gcenter[1]
	    fitvec[4] = gheight[1]
	    fitvec[5] = gwidth[1]
	    if (ngfit > 1) {
		fitvec[6] = gcenter[2]
		fitvec[7] = gheight[2]
		fitvec[8] = gwidth[2]
		}
	    if (ngfit > 2) {
		fitvec[9] = gcenter[3]
		fitvec[10] = gheight[3]
		fitvec[11] = gwidth[3]
		}
	    call emmin (fdata, npix, nparams, fitvec, covar,
			  pixleft, pixright, nit, xi2, xi2dof, ierr, debug)

	# Transfer the continuum poly and each gaussian 
	    do j = 1, ngfit {
		if (!igauss[j]) {
		    next
		    }
		m = m + 1
		do k = 1, 10 {
		    emparams[k,1,m] = 0.
		    emparams[k,2,m] = 0.
		    }

	# Error? 0=OK, 1=iter. exceeded, 2=Singular 
		if (ierr < 2) { 
		    if (ierr == 1 && debug) {
			call printf ("Line %d not converging after iteration %d\n")
			    call pargi (i)
			    call pargi (nit)
			}
		    jj = npoly + ((j - 1) * 3)

		# Continuum polynomial and formal errors
		    if (npoly > 0) {
			do k = 1, npoly {
			    emparams[k,1,m] = fitvec[k]
			    emparams[k,2,m] = covar[k]
			    }
			}

		# Gaussian and formal errors
		    do k = 1, 3 {  
			emparams[k+3,1,m] = fitvec[jj + k]
			emparams[k+3,2,m] = covar[jj + k]
			}

		# Find Angstroms to pixels relation at this line
		    px_cen = emparams[4,1,m]
		    px_cen_plusone = emparams[4,1,m] + 1.d0
		    wv_cen = wcs_p2w (px_cen)
		    wlobs[m] = wv_cen
		    wv_cen_plusone = wcs_p2w (px_cen_plusone)
		    dwdp  = wv_cen_plusone - wv_cen
		    if (dwdp < 0.d0) dwdp = -dwdp

		# Gaussian half-width in pixels
		    sigma = fitvec[jj+3]

		# Gaussian height or depth in counts per pixel
		    height = abs (fitvec[jj+2])

		# Area inside of line (pixels * counts/pixel)
		    area = s2pi * height * sigma

		# Mean continuum level under/above line (counts/pixel)
		    icenter = idnint (fitvec[jj+1])
		    cont = fcont[icenter] + fitvec[1] +
			   (fitvec[2] * (fitvec[jj+1] - pixleft))
		    if (cont <= 0)
			cont = 1.d0

		# Equivalent width is area of line divided by continuum
		    if (cont > mincont)
			eqw = area / cont

		# If continuum is less than minimum, set to 1
		    else
			eqw = area

		# Convert equivalent width from pixels to angstroms
		    eqw = eqw * dwdp
		    emparams[7,1,m] = eqw

		# Error in area if equivalent width cannot be computed
		    if (cont <= mincont) {
			if (sigma > 0.d0)
			    relsig = covar[jj+3] / sigma
			else
			    relsig = 1.d-4
			if (height > 0.d0)
			    relht = covar[jj+2] / height
			else
			    relht = 1.d-4
			emparams[7,2,m] = area * sqrt(relht*relht+relsig*relsig)
			}

		# Error in equivalent width
		    else {
			if (skyspec) {
			    skylevel = 0.0
			    do k = fitvec[jj+1] - 2, fitvec[jj+1] + 2 {
				skylevel = skylevel + sky[k]
				}
			    skylevel = skylevel / 5.0
			    }
			else
			    skylevel = 0.0
			area = area + 8.0 * (cont + skylevel) * sigma +
				(cont + 2.0 * skylevel) /
				(dwdp * double (pixright-pixleft+1) * (eqw**2))
			emparams[7,2,m] = sqrt (abs (area)) / cont
			}

		# Chi^2, degrees of freedom 
		    emparams[8,1,m] = xi2
		    emparams[8,2,m] = xi2 / xi2dof
		    }

		else if (debug) {
		    call printf ("Line %d not fit, singular matrix in iteration %d\n")
			call pargi (i)
			call pargi (nit)
		    }
		}
	    i = i + ngauss
	    }
	return
end


# input list MUST BE WAVELENGTH-SORTED! 

int procedure iscombo (restw, n, limit, icombo, wlfit, ngfit, igauss)

double	restw[ARB]	# Wavelength sorted list 
int	n		# Index of restw to test 
int	limit		# Size of list 
int	icombo		# Which combination is found
double	wlfit		# Wavelength to fit on either side of peak(s)
int	ngfit		# Number of Gaussians to fit
bool	igauss[3]	# Flags for found lines in combination

int i, j, k, m, ngauss

include "emv.com"

begin

	icombo = 0
	ngauss = 0
	m = n
	igauss[1] = FALSE
	igauss[2] = FALSE
	igauss[3] = FALSE
	do i = 1, ncombo {
	    if (ngauss != 0)
		break

	# Check for membership in i'th combination
	    do j = 1, numcom[i] {
		if (ngauss != 0)
		    break

	    # If member of combination, note that fact
		if (abs (restw[m] - elines[j,i]) < 0.1) {
		    ngauss = 1
		    igauss[j] = TRUE
		    wlfit = edwl[i]
		    ngfit = numcom[i]
		    icombo = i

		# See if another line in the combo is in the input list
		# by testing the next combo element against next restw line.
		# restw MUST BE SORTED! 
		    do k = j+1, ngfit {
			m = m + 1
			if (m > limit)
			    break
			else if (abs (restw[m] - elines[k,i]) < 0.1) {
			    ngauss = ngauss + 1
			    igauss[k] = TRUE
			    }
			else
			    break
			}
		    }
		}
	    }

	# If not a member of a combination, fit line as single gaussian
	if (ngauss == 0) {
	    ngauss = 1
	    ngfit = 1
	    wlfit = wgfit[n]
	    igauss[1] = TRUE
	    }
	return ngauss
end
# Dec  4 1991	Use fit continuum array
# Dec  6 1991	Fit all lines in combination whether identified or not

# Mar 24 1992	Use combination line information from file
# Mar 25 1992	Change GWIDTH to 10
# Mar 31 1992	Avoid dividing by cont if it is =0
# Apr  1 1992	Set fit width in wavelength, not pixels
# Nov  9 1992	Allow 0-3 polynomial continuum coefficients
# Nov 24 1992	Set fit half-width for single lines from emlines.dat file

# Jun  8 1993	Move wavelength <-> pixel conversion to subroutines

# Apr 22 1994	Drop length argument from GAUSS call
# Jun 13 1994	Add continuum to MINIFIT and GAUSS calls
# Aug 12 1994	Clean up code; set wlobs from Gaussian fit, not search

# Feb 19 1997	Drop unfound combination lines, if allowed
# Mar 14 1997	New version using Doug's least squares fitting subroutines
# Apr 17 1997	Compute equivalent width correctly for emission lines
# May  6 1997	Rework equivalent width once again
# May  9 1997	Add MINCONT minimum continuum for equivalent width else area
# May  9 1997	Correct computation of area error if no equivalent width
# Jul  8 1997	Return if no lines found

# Mar 23 2001	Fix limits to work with reversed spectrum

# Oct 19 2004	Fix equivalent width computation so units match
