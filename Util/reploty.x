# File rvsao/Util/reploty.x
# March 23, 2005
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss, Johns Hopkins University

# Copyright(c) 2005 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#  REPLOTY -- Replot the current array (Adapted from IRAF onedspec.splot)
 
procedure reploty (gfd, gt, npts, pix, xvec, y1, y2)
 
pointer	gfd, gt		# Pointers for graphics control
int	npts		# Number of points in pix
real	pix[ARB]	# The data array to plot
real	xvec[ARB]	# X coordinates of data array
real	y1, y2		# Minimum and maximum y coordinates

real	x1, x2		# Beginning and ending x coordinates
 
begin
	call gclear (gfd)
	x1 = xvec[1]
	x2 = xvec[npts]
	if (x1 > x2) {
	    x1 = xvec[npts]
	    x2 = xvec[1]
	    }
	call gswind (gfd, x1, x2, y1, y2)
#	call gascale (gfd, pix, npts, 2)
	call gt_swind (gfd, gt)
	call gt_labax (gfd, gt)
	if (npts > 1)
	    call gpline (gfd, xvec, pix, npts)
	return
end

# Mar 23 2005	Adapted from replot by Doug Mink
