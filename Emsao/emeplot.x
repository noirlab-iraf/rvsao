# File rvsao/Emsao/emeplot.x
# November 18, 2008
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1991-2008 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  EMEPLOT uses IRAF plot routines to plot a spectrum
#  and label its emission and absorption lines

include <gset.h>
include <gio.h>
include <pkg/gtools.h>
include <ttyset.h>
include <fset.h>
include "rvsao.h"
include "emv.h"

procedure emeplot (gfd, specfile, specim, mspec, pix, wpix, npix)
 
pointer	gfd		# Pointer to graphics control structure
char	specfile[ARB]	# Name of spectrum file
pointer	specim		# IRAF image descriptor structure
int	mspec		# Spectrum number in multispec, else 0
real	pix[ARB]	# Data vector to plot
real	wpix[ARB]	# Wavelengths for data vector
int	npix		# Number of points in pix

real	xmin, xmax	# Beginning and ending x coordinates
real	ymin, ymax	# Beginning and ending y coordinates
bool	alldata		# True if entire spectrum is being plotted
bool	dispvert
real	x, y, y1, y2, dy, dy4
double	wx
int	i,j, iy, ny, npix1, jm1, jp1, ifound
char	text[SZ_LINE+1]
#char	xlab[SZ_LINE+1]
char	labform[32], labname[32], wlname[16]
double	bcz,velplot, dindef

include	"emv.com"
include	"rvsao.com"

common/emp/ shift
double	shift
 
begin
	dindef = INDEFD
	call gclear (gfd)
	dispvert = TRUE

# If npix is passed as negative, view is zoomed in on portion of data
	if (npix < 0) {
	    alldata = FALSE
	    npix = -npix
	    }
	else
	    alldata = TRUE
	npix1 = npix - 1

# Set format for line labels
	if (dispem && !dispabs)
	    dispvert = TRUE
	else
	    dispvert = FALSE
	if (dispvert)
	    call strcpy ("u=180;p=r;q=h;f=i;h=c;s=0.6", labform, 32)
	else
	    call strcpy ("q=h;f=i;h=c;s=0.8", labform, 32)

# Set velocity to which to shift line labels
	switch (vplot) {
	    case VEMISS:
		velplot = spevel
	    case VCORREL:
		velplot = spxvel
	    case VCOMB:
		velplot = spvel
	    case VGUESS:
                velplot = gvel
	    default:
		vplot = VCOMB
		velplot = spvel
	    }
	if (velplot == dindef)
	    velplot = 0.d0

# If not zoomed, plot information about velocities
	if (alldata) {
#	    call emiplot (gfd, specfile, specim, mspec)
	    text[1] = EOS
	    }

	# Print specfile
	else {

	    switch (vplot) {
		case VEMISS:
		    call sprintf (text,SZ_LINE,
			"%s  em vel= %.2f (%.2f) km/sec")
			call pargstr (specid)
			call pargd (velplot)
			call pargd (speerr)
		case VCORREL:
		    call sprintf (text,SZ_LINE,
			"%s  corr vel= %.2f (%.2f) km/sec R= %.2f")
			call pargstr (specid)
			call pargd (velplot)
			call pargd (spxerr)
			call pargd (spxr)
		case VCOMB:
		    call sprintf (text,SZ_LINE,
			"%s  comb vel= %.2f (%.2f) km/sec R= %.2f")
			call pargstr (specid)
			call pargd (velplot)
			call pargd (sperr)
			call pargd (spxr)
		case VGUESS:
		    if (velplot < 1.d0)
			velplot = velplot * c0
		    call sprintf (text,SZ_LINE,"%s vel= %8.2f km/sec")
			call pargstr (specid)
			call pargd (velplot)
		default:
		}
	    }
	bcz = 1.d0 + (spechcv / c0)
	shift = (1.d0 + (velplot / c0)) / bcz

#  Set axis limits
	xmin = wpix[1]
	xmax = wpix[npix]
	if (xmin > xmax) {
	    xmin = wpix[npix]
	    xmax = wpix[1]
	    }

# Set minimum and maximum Y values - limits = +- 10%
	ymin = pix[1]
	ymax = pix[1]
	do i = 2, npix {
	    if (pix[i] < ymin)
		ymin = pix[i]
	    if (pix[i] > ymax)
		ymax = pix[i]
	    }
	if (ymin > 0.) {
	    ymin = 0.
	    dy = ymax / 100.
	    if (dispvert)
		ymax = ymax + (dy * 30.)
	    else
		ymax = ymax + (dy * 10.)
	    }
	else if (ymin == ymax) {
	    dy = .01
	    ymin = ymin - (dy * 100.)
	    ymax = ymax + (dy * 100.)
	    }
	else {
	    ymin = ymin - dy
	    dy = (ymax - ymin) / 100.
	    if (dispvert)
		ymax = ymax + (dy * 30.)
	    else
		ymax = ymax + (dy * 10.)
	    }
	dy4 = dy * 0.25

	call gseti (gfd, G_WCS, 1)
	call gswind (gfd, xmin, xmax, ymin, ymax)

# Set viewport to make space for results at right if they are requested
	if (alldata)
	    call gsview (gfd, 0.08, 0.66, 0.10, 0.90)
	else
	    call gsview (gfd, 0.0, 1.0, 0.0, 1.0)

	call gseti (gfd,G_TXQUALITY,GT_HIGH)

#  Plot and label axes
	if (mspec < 0)
	    call glabax (gfd,"","Wavelength in angstroms","")
	else
	    call glabax (gfd,text,"Wavelength in angstroms","")

#  Plot spectrum
	call gpline (gfd, wpix, pix, npix)
	call gflush (gfd)

#  Label emission lines
	if ((dispem || plotem) && nref > 0) {
	    do i = 1, nref {
		ifound = lfound[i]
		if (ifound > 0)
		    wx = wlrest[ifound] * shift
		else
		    wx = wlref[i] * shift
		x = wx
		if (x > xmin && x < xmax) {
		    do j = 1, npix1 {
			if ((wx >= wpix[j] && wx < wpix[j+1]) ||
			    (wx <= wpix[j] && wx > wpix[j+1])) {
			    jm1 = j - 1
			    if (jm1 < 1) jm1 = 1
			    jp1 = j + 1
			    if (jp1 > npix) jp1 = npix
			    y1 = pix[j]
			    if (y1 < pix[jm1])
				y1 = pix[jm1]
			    if (y1 < pix[jp1])
				y1 = pix[jp1]
			    y1 = y1 + (2.5 * dy)
#			    break
			    }
			}
		    if (dispvert) {
			y2 = y1 + (10.0 * dy)
			if (y2 > ymax - (30.0 * dy)) y2 = ymax - (25.0 * dy)
			}
		    else {
			y2 = y1 + (30.0 * dy)
			if (y2 > ymax - (5.0 * dy)) y2 = ymax - (5.0 * dy)
			}
		    if (y2 < y1) y2 = y1 + dy
		    if (y1 < ymax) {
			if (ifound > 0 &&
			    (wtobs[ifound] > 0 || override[ifound] >= 0)) {
			    call gline (gfd,x,y1,x,y2)
			    }
			else {
			    ny = 1. + ((y2 - y1) / dy)
			    y = y1
			    do iy = 1, ny {
				call gline (gfd, x, y-dy4, x, y+dy4)
				y = y + dy
				}
			    }
			call strcpy (nmref[1,i], labname, 32)
			if (dispvert) {
			    call strcat (" ", labname, 32)
			    call sprintf (wlname, 16,"%8.3f")
				call pargd (wlref[i])
			    call strcat (wlname, labname, 32)
			    }
			call gtext (gfd, x,y2+dy,labname,labform)
			}
		    }
		}
	    }
	call gflush (gfd)

#  Label absorption lines
	if (dispabs && nabs > 0) {
	    do i = 1, nabs {
		wx = wlabs[i] * shift
		x = wx
		if (x > xmin && x < xmax) {
		    do j = 1, npix1 {
			if ((wx >= wpix[j] && wx < wpix[j+1]) ||
			    (wx <= wpix[j] && wx > wpix[j+1])) {
			    jm1 = j - 1
			    if (jm1 < 1) jm1 = 1
			    jp1 = j + 1
			    if (jp1 > npix) jp1 = npix
			    y1 = pix[j]
			    if (y1 > pix[jm1])
				y1 = pix[jm1]
			    if (y1 > pix[jp1])
				y1 = pix[jp1]
			    y1 = y1 - (2.5 * dy)
			    }
			}
		    y2 = y1 - (30.0 * dy)
		    if (y2 < ymin + (5.0 * dy)) y2 = ymin + (5.0 * dy)
		    if (y2 > y1) y2 = y1 - dy
		    ny = 1. + ((y1 - y2) / dy)
		    y = y1
		    do iy = 1, ny {
			call gline (gfd,x,y+dy4,x,y-dy4)
			y = y - dy
			}
		    call strcpy (nmabs[1,i], labname, 32)
		    if (dispvert) {
			call strcat (" ", labname, 32)
			call sprintf (wlname, 16,"%8.3f")
			    call pargd (wlabs[i])
			call strcat (wlname, labname, 32)
			}
		    call gtext (gfd, x,y2-3*dy,labname,labform)
		    }
		}
	    }
	call gflush (gfd)

# If not zoomed, plot information about velocities
	if (alldata) {
	    call emiplot (gfd, specfile, specim, mspec)
	    text[1] = EOS
	    }
end

# Sep 17 1991	Add log-lambda plotting
# Dec  6 1991	Use shifted rest wavelength, not observed wavelength for line y
# Dec 18 1991	Use separate variable for plotted velocity, not VEL
# Dec 19 1991	Position line markers better

# Mar 23 1992	Read absorption line information from a file

# Dec  2 1993	Print multispec spectrum number

# Apr 13 1994	Drop image header pointer as argument; fix call to gset
# May 18 1994	Set space between markers and spectrum better
# Jun  9 1994	If dy is 0, set it to .001
# Jun 23 1994	Get velocities from fquot instead of getim labelled common
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug 10 1994	Print version and date on plot
# Aug 17 1994	Print summary page if spectrum is not zoomed
# May 15 1995	Change all sz_fname to sz_pathname, which is longer
# Jun 20 1995	Do not set specid; use setting from getspec
# Aug 25 1995	If plotting at corr. velocity and DISPEM, plot all em. lines
# Oct  2 1995	Use GVEL instead of reading parameter CZ_GUESS

# Apr  7 1997	Fix bugs in calling programs
# Jul  8 1997	If VELPLOT is INDEF, set it to 0
# Sep 25 1997	Draw dashed lines instead of points

# Jul 31 1998	Make printing of heading optional

# Aug 19 1999	Plot emission lines for combination as well as correlation vel

# Mar 28 2001	Make plotting work if spectrum is reversed

# May 12 2004	Drop dw and dlw from common/emp/

# Jan 17 2007	Always mark emission lines, keeping track of which were fit
# Jan 31 2007	Drop unused variable nelines
