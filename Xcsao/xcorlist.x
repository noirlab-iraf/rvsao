#  List cross-correlation peak on XCSAO summary page

include "rvsao.h"

procedure xcorlist (gp, npoints, xpts, ypts, it)
 
pointer gp		# Graphics structure pointer
int 	npoints		# Number of points in spectrum
real	xpts[ARB]	# Array of x-coordinates to plot
real	ypts[ARB]	# Array of y-coordinates to plot
int	it		# Template in R-value order for which to plot cross-correlation

real	xcr, dxcr, msize, y01[2]
real	rindef
pointer	fraclev, xlev, sp
double	fracpeak,ph,pw
int	il, ir, pkindx, i, npts4, itemp
char	gtitle[SZ_LINE]	# Graph title
#int	npfit

include "rvsao.com"
include "results.com"

begin
	itemp = itr[it]
	rindef = INDEFR
	call printf ("XCORLIST: %d Corr. Template: %s\n")
	    call pargi (it)
	    call pargstr (tempname[1,itemp])

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
	call printf ("XCORLIST: velocity range %f to %f\n")
	    call pargr (xcr-dxcr)
	    call pargr (xcr+dxcr)
	call flush (STDOUT)

	call printf ("XCORLIST: vel[1] = %.3f, vel[%d] =  %.3f\n")
	    call pargi (ixi)
	    call pargr (xpts[1])
	    call pargi (npoints)
	    call pargr (xpts[npoints])
	call flush (STDOUT)
	call printf ("XCORLIST: xcor = %.3f - %.3f\n")
	    call pargr (xcor[1])
	    call pargr (xcor[npts])
	call flush (STDOUT)

#	if (npoints > 0)
#	    call rvscale (gp, ypts, npoints, 2)
#	else {
#	    y01[1] = 0.
#	    y01[2] = 1.
#	    call rvscale (gp, y01, 2, 2)
#	    }

	return
end
