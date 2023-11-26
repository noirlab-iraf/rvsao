# File rvsao/Xcsao/xcorplot.x
# July 22, 2010
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Stephen Levine and Jon Morse

# Copyright(c) 1995-2010 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

include	<gset.h>
include	<gio.h>
include	"rvsao.h"
 
#  Plot cross-correlation peak on XCSAO summary page

procedure xcorplot (gp, npoints, xpts, ypts, it)
 
pointer gp		# Graphics structure pointer
int 	npoints		# Number of points in spectrum
real	xpts[ARB]	# Array of x-coordinates to plot
real	ypts[ARB]	# Array of y-coordinates to plot
int	it		# Template in R-value order for which to plot cross-correlation

real	xcr, dxcr, msize, y01[2]
real	rindef
double	fracpeak,ph,pw
int	il, ir, pkindx, i, npts4, itemp
char	gtitle[SZ_LINE]	# Graph title
#int	npfit

include "rvsao.com"
include "results.com"
include "xplt.com"

begin
	itemp = itr[it]
	rindef = INDEFR
	call sprintf (gtitle,SZ_LINE,"%d Corr. Template: %s")
	    call pargi (it)
	    call pargstr (tempname[1,itemp])
	call gseti (gp, G_WCS, 2)

# Set font and size of tick labels
	call gseti (gp,G_TXQUALITY,GT_HIGH)
	call gsetr (gp,G_TICKLABELSIZE,0.65)
 
# Set scale for plot of correlation
	if (xcr0 == rindef)
	    xcr = zvel[itemp]
	else
	    xcr = xcr0
	if (xcrdif == rindef) {
	    dxcr = 20. * tvw[itemp]
	    if (dxcr < 1.d0) dxcr = 300.d0
	    }
	else
	    dxcr = xcrdif
	if (dxcr == 0) return
	call gswind (gp, xcr-dxcr, xcr+dxcr, INDEF, INDEF)

	if (npoints > 0)
	    call rvscale (gp, ypts, npoints, 2)
	else {
	    y01[1] = 0.
	    y01[2] = 1.
	    call rvscale (gp, y01, 2, 2)
	    }

# Set viewport to half page         
	call gsview (gp, 0.06, 0.66, 0.15, 0.45)

# Plot axes
	if (correlate == COR_PIX)
	    call glabax (gp, gtitle,"Pixel Shift", "")
	else if (correlate == COR_WAV)
	    call glabax (gp, gtitle,"Wavelength Shift [Angstroms]", "")
	else
	    call glabax (gp, gtitle,"Velocity [km/s]", "")
	call gflush (gp)

# Plot cross-correlation
	if (npoints > 0)
	    call gpline (gp, xpts, ypts, npoints)
	call gflush (gp)

# Mark points used in peak fit if PKFRAC is negated
	if (npoints > 0 && pkfrac < 0.d0) {
	    fracpeak = -pkfrac
	    call pkwidth (npts,ypts,xpts,fracpeak,pkindx,ph,pw,il,ir)
#	    npfit = ir - il + 1
#	    msize = 0.005
#	    call gpmark (gp, xpts[il], ypts[il], npfit, GM_CIRCLE, msize, msize)

	    npts4 = npts / 4
	    if (npts4 > maxpts4) {
		if (fraclev != NULL) {
		    call mfree (fraclev, TY_REAL)
		    fraclev = NULL
		    }
		if (xlev != NULL) {
		    call mfree (xlev, TY_REAL)
		    xlev = NULL
		    }
		maxpts4 = npts4
		call malloc (fraclev, npts4, TY_REAL)
		call malloc (xlev, npts4, TY_REAL)
		}
	    for (i = 0; i < npts4; i = i + 1) {
		Memr[fraclev+i] = fracpeak * ph
		Memr[xlev+i] = xpts[(4*i)+2]
		}
	    msize = 0.001
	    call gpmark (gp, Memr[xlev], Memr[fraclev], npts4, GM_CIRCLE, msize, msize)
	    call gflush (gp)
	    }

end

# May 10 1995	Flush graphics buffer
# Aug 22 1995	Mark points used in peak fit if PKFRAC is negated
# Aug 22 1995	Mark minimum fit level if PKFRAC is negated
# Sep 20 1995	Don't plot data if there are no correlation points to plot
# Oct  5 1995	Change plotting characters for peak

# Apr 14 1997	Get R-value-sorted template number from XCPLOT instead of title
# May  2 1997	Always test against rindef, not INDEFR

# Sep 19 2000	Add labels for correlations in pixel and wavelength space

# Aug 13 2007	Allocate buffers only when needed; include xplt.com

# Jul 22 2010	Free individual pointers only if already allocated
