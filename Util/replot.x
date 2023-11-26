# File rvsao/Util/replot.x
# September 18, 2006
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss, Johns Hopkins University

# Copyright(c) 2006 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  REPLOT -- Replot the current array (Adapted from IRAF onedspec.splot)
 
procedure replot (gfd, gt, npts, pix, xvec)
 
pointer	gfd, gt		# Pointers for graphics control
int	npts		# Number of points in pix
real	pix[ARB]	# The data array to plot
real	xvec[ARB]	# Beginning and ending x coordinates

real	x1, x2		# Beginning and ending x coordinates
 
begin
	call gclear (gfd)
	x1 = xvec[1]
	x2 = xvec[npts]
	if (x1 > x2) {
	    x1 = xvec[npts]
	    x2 = xvec[1]
	    }
	call gswind (gfd, x1, x2, INDEF, INDEF)
	call gascale (gfd, pix, npts, 2)
	call gt_swind (gfd, gt)
	call gt_labax (gfd, gt)
	if (npts > 1)
	    call gpline (gfd, xvec, pix, npts)
	return
end

#  REPLOTY -- Replot the current array (Adapted from IRAF onedspec.splot)
 
procedure reploty (gfd, gt, npts, pix, xvec, y1, y2)
 
pointer	gfd, gt		# Pointers for graphics control
int	npts		# Number of points in pix
real	pix[ARB]	# The data array to plot
real	xvec[ARB]	# X coordinates of data array
real	y1, y2		# Minimum and maximum y coordinates

real	x1, x2		# Beginning and ending x coordinates
real	y1x, y2x, rindef, diff
int	i
 
begin
	rindef = INDEFR
	call gclear (gfd)
	x1 = xvec[1]
	x2 = xvec[npts]
	if (x1 > x2) {
	    x1 = xvec[npts]
	    x2 = xvec[1]
	    }
	if (y1 == rindef) {
	    y1x = pix[1]
	    do i = 1, npts {
		if (pix[i] < y1x)
		    y1x = pix[i]
		}
	    }
	else
	    y1x = y1
	if (y2 == rindef) {
	    y2x = pix[1]
	    do i = 1, npts {
		if (pix[i] > y2x)
		    y2x = pix[i]
		}
	    diff = (y2x - y1x) * 0.05
	    y1x = y1x - diff
	    y2x = y2x + diff
	    }
	else
	    y2x = y2

	call gswind (gfd, x1, x2, y1x, y2x)
#	call gascale (gfd, pix, npts, 2)
	call gt_swind (gfd, gt)
	call gt_labax (gfd, gt)
	if (npts > 1)
	    call gpline (gfd, xvec, pix, npts)
	return
end

# June	1987	Gerard Kriss
# Oct	1991	Doug Mink	Use x vector instead of end points
# Dec 19 1991	Call GPLINE instead of GAMOVE and GADRAW

# Feb 14 1992	Removed unused variable i

# May  4 1994	Switch order of number of points and spectrum arguments

# Mar 26 2001	Always plot spectra with wavelength increasing left to right

# Mar 23 2005	reploty adapted from replot with added y limits

# Sep 18 2006	Added reploty to replot.x
