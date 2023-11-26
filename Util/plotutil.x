# File rvsao/Util/plotutil.x
# September 7, 1999
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

include	<syserr.h>
include <gset.h>
include	<gio.h>

# Open output device and return graphics descriptor pointer

pointer procedure openplot (device)

char	device[ARB]

pointer gp
pointer gopen()

begin             

# Open plotter
	gp = gopen (device, NEW_FILE, STDGRAPH)
 	call gclear (gp)

# Add x- and y-axis minor ticks.
	call gseti (gp, G_YNMINOR, YES)
        call gseti (gp, G_XNMINOR, YES)

	return (gp)

end


#  Close graphics output device

procedure closeplot (gp)

pointer gp

begin

# Close up graph window
	call gflush (gp)
	call gclose (gp)

end
 
# Copyright(c) 1986 Association of Universities for Research in Astronomy Inc.
# GASCALE modified to add 10% above and below graph by Doug Mink, CfA

# RVSCALE -- Scale the world coordinates of either the X or Y axis to fit the
# data vector.  This is done by setting the WCS limits to the minimum and
# maximum pixel values of the data vector.  The original WCS limits are
# overwritten.

procedure rvscale (gp, v, npts, axis)

pointer	gp			# graphics descriptor
real	v[ARB]			# data vector
int	npts			# length of data vector
int	axis			# axis to be scaled (1=X, 2=Y)

int	start, i
real	minval, maxval, pixval, diff
pointer	w

begin
	# Find first definite valued pixel.  If entire data vector is
	# indefinite we cannot perform our function and must abort.

	for (start=1;  start <= npts;  start=start+1)
	    if (!IS_INDEF (v[start]))
		break
	if (start > npts)
	    call syserr (SYS_GINDEF)

	minval = v[start]
	maxval = minval

	# Compute min and max values of data vector.
	do i = start+1, npts {
	    pixval = v[i]
	    if (!IS_INDEF(pixval))
		if (pixval < minval)
		    minval = pixval
		else if (pixval > maxval)
		    maxval = pixval
	}

	w = GP_WCSPTR (gp, GP_WCS(gp))
	diff = (maxval - minval) * 0.1

	# Set the window limits.
	switch (axis) {
	case 1:
	    WCS_WX1(w) = minval
	    WCS_WX2(w) = maxval
	case 2:
	    WCS_WY1(w) = minval - diff
	    WCS_WY2(w) = maxval + diff
	case -2:
	    WCS_WY1(w) = 0.0
	    WCS_WY2(w) = maxval + diff
	default:
	    call syserr (SYS_GSCALE)
	}

	GP_WCSSTATE(gp) = MODIFIED
	call gpl_reset()
end

# Aug 15 1994	New program

# Jan 10 1996	Add possibility of 0 lower limit for spectrum plot

# Sep  7 1999	Put all includes at start of file
