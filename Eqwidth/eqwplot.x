# File rvsao/Eqwidth/eqwplot.x
# May 14, 2009
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 2005-2009 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  EWQPLOT uses IRAF plot routines to plot a spectrum
#  and label regions where equivalent widths are being measured

include <gset.h>
include <gio.h>
include <pkg/gtools.h>
include <ttyset.h>
include <fset.h>
include "rvsao.h"
include	"eqw.h"

procedure eqwplot (specfile, specim, pix, wpix, npix, wlab, bands, nbands)
 
pointer	gfd		# Pointer to graphics control structure
char	specfile[ARB]	# Name of spectrum file
pointer	specim		# IRAF image descriptor structure
real	pix[ARB]	# Data vector to plot
real	wpix[ARB]	# Wavelengths for data vector
int	npix		# Number of points in pix
char	wlab[ARB]	# Label for wavelength axis
pointer	bands		# Band data structure
int	nbands		# Number of bands

real	xmin, xmax	# Beginning and ending x coordinates
real	ymin, ymax	# Beginning and ending y coordinates
real	x, y, y1, y2, dy, dy4
double	w1, w2
int	i,j, k, iy, ny, npix1, jm1, jp1
double	dindef
pointer	band1

include	"rvsao.com"

common/emp/ shift
double	shift
 
begin
	dindef = INDEFD
	call gclear (gfd)
	npix1 = npix - 1

	# Print specfile

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
	    ymax = ymax + (dy * 10.0)
	    }
	else if (ymin == ymax) {
	    dy = .01
	    ymin = ymin - (dy * 100.)
	    ymax = ymax + (dy * 100.)
	    }
	else {
	    ymin = ymin - dy
	    dy = (ymax - ymin) / 100.
	    ymax = ymax + (dy * 10.0)
	    }
	dy4 = dy * 0.25

	call gseti (gfd, G_WCS, 1)
	call gswind (gfd, xmin, xmax, ymin, ymax)

# Set viewport
	call gsview (gfd, 0.0, 1.0, 0.0, 1.0)
	call gseti (gfd,G_TXQUALITY,GT_HIGH)

#  Plot and label axes
	call glabax (gfd,"",wlab,"")

#  Plot spectrum
	call gpline (gfd, wpix, pix, npix)
	call gflush (gfd)

#  Label equivalent width regions
	if (nbands > 0)
	    do i = 1, nbands {
		band1 = BAND(bands,i,BAND1)
		w1 = BAND_W1(band1)
		w2 = BAND_W2(band1)
		x = w1
		do k = 1, 2 {
		if (x > xmin && x < xmax) {
		    do j = 1, npix1 {
			if ((w1 >= wpix[j] && w1 < wpix[j+1]) ||
			    (w1 <= wpix[j] && w1 > wpix[j+1])) {
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
		    y2 = y1 + (30.0 * dy)
		    if (y2 > ymax - (5.0 * dy)) y2 = ymax - (5.0 * dy)
		    if (y2 < y1) y2 = y1 + dy
		    if (y1 < ymax)
			call gline (gfd,x,y1,x,y2)
		    else {
			ny = 1. + ((y2 - y1) / dy)
			y = y1
			do iy = 1, ny {
			    call gline (gfd, x, y-dy4, x, y+dy4)
			    y = y + dy
			    }
			}
		    call gtext (gfd, x,y2+dy,Memc[BAND_ID(band1)],
				"q=h;f=i;h=c;s=0.8")
		    }
		}
	    }
	call gflush (gfd)
	call gclose (gfd)

end

# Sep 21 2005	New program based on Emsao/emeplot.x

# May 14 2009	Read wavelength limits, not center/width from BAND structure
