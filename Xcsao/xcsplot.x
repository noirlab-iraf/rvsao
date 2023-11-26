# File rvsao/Xcsao/xcsplot.x
# September 20, 2000
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Stephen Levine and Jon Morse (major changes)

# Copyright(c) 2000 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

include	<gset.h>
include	<gio.h>
include	<smw.h>
include	"rvsao.h"
 
#  Plot spectrum on XCSAO summary page
 
procedure xcsplot (gp, npoints, xpts, ypts, minwav, maxwav, dispmode)
 
pointer gp		# Graphics structure pointer
int	npoints		# Number of points in spectrum to plot
real	xpts[ARB]	# Array of wavelengths to plot
real	ypts[ARB]	# Array of intensities to plot
double	minwav		# Minimum wavelength for cross-correlation
double	maxwav		# Maximum wavelength for cross-correlation
int	dispmode	# Display mode, negate for zero-minimum plot

int	npw		# Number of points in spectrum to plot
int	ipw		# First wavelength point to plot
int	i
real	wx, wy, wy1, wy2, dy, dy4, dwl
pointer	w
char	xlab[SZ_LINE+1]
int	strlen()

include "rvsao.com"
include "results.com"

begin

	call gseti (gp, G_WCS, 1)

# Set font and size of tick labels
	call gseti (gp,G_TXQUALITY,GT_HIGH)
	call gsetr (gp,G_TICKLABELSIZE,0.65)

	ipw = 1
	npw = npoints

# Find points within wavelength limits (ignored for now)
	if (dispmode == 0) {
	    if (minwav != maxwav) {
		do i = 1, npoints {
		    if (xpts[i] < minwav)
			ipw = i + 1
		    if (xpts[i] > maxwav) {
			npw = i - 1
			break
			}
		    }
		}
	    npw = npw - ipw + 1
	    }
 
# Set scale for plot of data array
	call gswind (gp, INDEF, INDEF, 0., INDEF)

# Set the window x (wavelength) limits.
	w = GP_WCSPTR (gp, GP_WCS(gp))
	if (dispmode == 0) {
	    WCS_WX1(w) = minwav
	    WCS_WX2(w) = maxwav
	    }
	else {
	    dwl = (xpts[npoints] - xpts[1]) / real (npoints - 1)
	    WCS_WX1(w) = xpts[1] - dwl
	    WCS_WX2(w) = xpts[npoints] + dwl
	    }

	GP_WCSSTATE(gp) = MODIFIED
	call gpl_reset()

# Get label for dispersion axis
	if (correlate == COR_PIX)
	    call strcpy ("Pixels",xlab,SZ_LINE)
	else {
	    call strcpy (LABEL(specsh),xlab,SZ_LINE)
	    if (strlen (UNITS(specsh)) > 0) {
		call strcat (" in ",xlab,SZ_LINE)
		call strcat (UNITS(specsh),xlab,SZ_LINE)
		}
	    }
	if (npw > 0) {
	    if (dispmode < 0)
		call rvscale (gp, ypts[ipw], npw, -2)
	    else
		call rvscale (gp, ypts[ipw], npw, 2)
	    }

# Set viewport to half page         
	call gsview (gp, 0.06, 0.66, 0.60, 0.90)

# Plot axes
       	call glabax (gp, "", xlab, "")
	call gflush (gp)

# Plot spectrum
	if (npw > 0) {
	    call gpline (gp, xpts, ypts, npoints)
	    call gflush (gp)
	    }

# Plot wavelength limits (when plotting full spectrum)
	if (dispmode != 0 && minwav != maxwav) {
	    w = GP_WCSPTR(gp,GP_WCS(gp))
	    wy1 = WCS_WY1(w)
	    wy2 = WCS_WY2(w)
	    dy = (wy2 - wy1) / 50.0
	    dy4 = dy * 0.25
	    wy = wy1
	    wx = minwav
	    do i = 1, 50 {
		call gline (gp, wx, wy-dy4, wx, wy+dy4)
		wy = wy + dy
		}
	    wy = wy1
	    wx = maxwav
	    do i = 1, 50 {
		call gline (gp, wx, wy-dy4, wx, wy+dy4)
		wy = wy + dy
		}
	    }
	call gflush (gp)
end
# May 10 1995	Flush graphcs buffer
# Jun 19 1995	If no points inside wavelength limits, do not try to plot

# Jan 10 1996	Add option to plot spectrum from zero if DISPMODE is -1

# Mar 14 1997	Set dispersion axis label from spectrum header
# Mar 14 1997	Plot all points if no wavelength limits
# Apr  8 1997	Always plot full spectrum; mark limits for displayed template
# Sep 25 1997	Plot dashed lines for limits instead of points

# Jan 12 1998	Add option to plot part of spectrum within wavelength limits

# Sep 20 2000	Label x-axis as pixels if correlating in pixel space
