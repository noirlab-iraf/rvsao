# File rvsao/Xcor/pkfitp.x
# August 16, 1996
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss

# Copyright(c) 1995 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#--- Fit peak of cross-correlation

define MAXPWIDTH 100

procedure pkfitp (npt,xind,xcor,xvel,fracpeak,pkindx,npfit,center,width,
		  height,debug)

int	npt		# Number of points in cross-correlation vector
real	xind[ARB]	# Cross-correlation vector indices
real	xcor[ARB]	# Cross-correlation vector
real	xvel[ARB]	# velocities corresponding to cross-correlation
double	fracpeak	# fraction of peak to use for fitting
int	pkindx		# Bin containing maximum
int	npfit		# number of points used in peak fit (returned)
double	center		# center of peak (returned)
double	width		# width of peak (returned)
double	height		# height of peak (returned)
bool	debug		# True for diagnostic printout

int	il, ir, ilu, iru
double	cu, cl
double	wu, wl
double	hu, hl
double	pwidth
double	pwu, pwl, pweight
double	co[3]
int	ncoef

begin

# Find the peak in the cross-correlation function.
# The peak found here can be altered in plotspec with the 'p' keystroke.
 
	call pkwidth (npt,xcor,xvel,fracpeak,pkindx,height,pwidth,il,ir)
	if (debug) {
	    call printf ("PKFITP: imax = %d, il= %d, ir = %d\n")
		call pargi (pkindx)
		call pargi (il)
		call pargi (ir)
	    call printf ("PKFITP: height = %f, width = %f\n")
		call pargd (height)
		call pargd (width)
	    }

	if (il == ir) {
	    npfit = 0;
	    if (pkindx > 0) {
		center = double (xind[pkindx])
		height = double (xcor[pkindx])
		}
	    else {
		center = 0.d0
		height = 0.d0
		}
	    width = 0.d0
	    if (debug)
		call printf ("PKFITP: No peak found, returning zeroes\n");
	    return;
	    }

	if (fracpeak == 0) fracpeak = 0.5d0

# Fit the peak down to the specified fraction
	if (fracpeak < 1.d0) {

	# Compute fit for endpoints just above half of peak
	    ilu = il + 1
	    iru = ir - 1
	    npfit = iru - ilu + 1
	    ncoef = 3
	    call earlin (xind[ilu],xcor[ilu],npfit,co,-ncoef)
	    call parabolapeak (co, cu, hu, wu)
	    pwu = double (iru - ilu)
	    if (debug) {
		call printf ("PKFITP U %d-%d %d pts c: %8.3f p: %6.4f w: %6.3f\n")
		    call pargi  (ilu)
		    call pargi  (iru)
		    call pargi  (npfit)
		    call pargd (cu)
		    call pargd (hu)
		    call pargd (wu)
		}

	# Compute fit for endpoints just below half of peak
	    npfit = ir - il + 1
	    ncoef = 3
	    call earlin (xind[il],xcor[il],npfit,co,-ncoef)
	    call parabolapeak (co,cl, hl, wl)
	    pwl = double (ir - il)
	    if (debug) {
		call printf ("PKFITP L %d-%d %d pts c: %8.3f p: %6.4f w: %6.3f\n")
		    call pargi (il)
		    call pargi (ir)
		    call pargi (npfit)
		    call pargd (cl)
		    call pargd (hl)
		    call pargd (wl)
		call printf ("PKFITP pwu: %6.3f pwl: %6.3f pw: %6.3f\n")
		    call pargd (pwu)
		    call pargd (pwl)
		    call pargd (pwidth)
		}

#	Compute weighted average of peak parameters
	    pweight = (pwidth - pwl) / 2.0
	    center = cl + ((cl - cu) * pweight)
	    height = hl + ((hl - hu) * pweight)
	    width = wl + ((wl - wu) * pweight)
	    }

# Fit a fixed-width parabola to the peak
# Try 11 pts for normal peaks, 5 pts for very narrow peaks (late type stars)

	else if (fracpeak > 1.d0) {
	    npfit = fracpeak
	    ir = il + npfit - 1
	    ncoef = 3
	    call earlin (xind[il],xcor[il],npfit,co,-ncoef)
	    call parabolapeak (co,center, height, width)
	    }

# If fit finds a peak which is outside of the fit area, assume it is wrong
	if (center > xind[ir] || center < xind[il] || height <= 0.) {
	    if (pkindx > 0) {
		center = double (xind[pkindx])
		height = double (xcor[pkindx])
		}
	    else {
		center = 0.d0
		height = 0.d0
		}
	    width = 0.d0
	    }
end


#--- Routine to compute the center x0, height h, width w of a parabola

procedure parabolapeak (a,x0,h,w)

double	a[3], h, w, x0

begin

	if (a[3] == 0) {
	    x0 = 0
	    h = 1
	    w = 1
	    return
	    }
	x0 = -a[2] / (2.d0 * a[3])
	h = a[1] - a[2]*a[2] / (4.d0 * a[3])
	w = -h / (2.d0 * a[3])
	if (w >= 0)
	    w = dsqrt (w)
	else
	    w = -dsqrt (-w)
end
# Sep 25 1991	Make all computations double instead of real, if possible

# Mar 27 1992	Pass velocity vector to PKWIDTH (and change argument order)
# Nov 30 1992	Pass debug as argument; drop inclusion of fquot common

# Apr 13 1994	Remove unused variable nhalf

# Jul 3 1995	If fit out of bounds, return initial peak and zero width
# Aug 21 1995	Drop debug from PKWIDTH call

# Aug 16 1996	Zero returned values if no peak is found
