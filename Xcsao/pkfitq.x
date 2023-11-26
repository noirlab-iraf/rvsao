# File rvsao/Xcor/pkfitq.x
# August 21, 1995
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss

# Copyright(c) 1994 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#--- Fit quartic to peak of cross-correlation

procedure pkfitq (npt,xind,xcor,xvel,fracpeak,imax,npfit,center,width,height,debug)

int	npt		# Number of points in cross-correlation vector
real	xind[ARB]	# Cross-correlation vector indices
real	xcor[ARB]	# Cross-correlation vector
real	xvel[ARB]	# velocities corresponding to cross-correlation
double	fracpeak	# Fraction of peak height above which fit is done
int	imax		# Bin containing maximum
int	npfit		# number of points used in peak fit (returned)
double	center		# center of peak (returned)
double	width		# width of peak (returned)
double	height		# height of peak (returned)
bool	debug		# True for diagnosting listing

int	il, ir, ilu, iru, ilh, irh
double	clh, crh
double	cu, cl
double	wu, wl
double	hu, hl
double	pwidth
double	pwu, pwl, pweight
double	co[5]
int	ncoef

begin

# Find the peak in the cross-correlation function.
# The peak found here can be altered in plotspec with the 'p' keystroke.
 
	call pkwidth (npt,xcor,xvel,0.5d0,imax,height,pwidth,ilh,irh)
	call pkwidth (npt,xcor,xvel,fracpeak,imax,height,pwidth,il,ir)

	if (fracpeak == 0) fracpeak = 0.5

# Fit the peak down to the specified fraction

	if (fracpeak < 1.) {

	# Compute fit for endpoints just above half of peak
	    ilu = il + 1
	    iru = ir - 1
	    npfit = iru - ilu + 1
	    ncoef = 5
	    call earlin (xind[ilu],xcor[ilu],npfit,co,-ncoef)
	    clh = xind[ilh]
	    crh = xind[irh]
	    cu = xind[imax]
	    call quarticpeak (co,clh,crh,cu,hu,wu)
	    pwu = dble (iru - ilu)
	    if (debug) {
		call printf ("PKFITP U %d pts c: %f p: %f w: %f\n")
		    call pargi (npfit)
		    call pargd (cu)
		    call pargd (hu)
		    call pargd (wu)
		}

	# Compute fit for endpoints just below half of peak
	    npfit = ir - il + 1
	    ncoef = 5
	    call earlin (xind[il],xcor[il],npfit,co,-ncoef)
	    clh = xind[ilh]
	    crh = xind[irh]
	    cl = xind[imax]
	    call quarticpeak (co,clh,crh,cl,hl,wl)
	    pwl = dble (ir - il)
	    if (debug) {
		call printf ("PKFITP L %d pts c: %f p: %f w: %f\n")
		    call pargi  (npfit)
		    call pargd (cl)
		    call pargd (hl)
		    call pargd (wl)
		call printf ("PKFITP pwu: %f pwl: %f pw: %f\n")
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

# Fit the peak to a fixed number of points
# Try 11 pts for normal peaks, 5 pts for very narrow peaks (late type stars)

	else if (fracpeak > 1) {
	    npfit = fracpeak
	    ncoef = 3
	    call earlin (xind[il],xcor[il],npfit,co,-ncoef)
	    clh = xind[ilh]
	    crh = xind[irh]
	    center = xind[imax]
	    call quarticpeak (co,clh,crh,center,height,width)
	    }
end


#--- Routine to compute the center and half max points of a quartic
#    x0 and xl, xr are provided with initial guesses

procedure quarticpeak (a,xl,xr,x0,h,w)

double	a[5], x0, xl, xr, h, w

double	acc, p, z, x, x1, step, s, xd, xc, diff
int	i, iroot

begin

	p(z) = a[1] + z*(a[2] + z*(a[3] + z*(a[4] + z*a[5])))
#	d(z) = a[2] + z*(2*a[3] + z*(3*a[4] + z*4*a[5]))
#	c(z) = 2*a[3] + z*(6*a[4] + z*12*a[5])
	acc = 1.d-3
	h = p(x0)
#	call printf ("PKFITQ: initial center= %f height= %f\n")
#	    call pargd (x0)
#	    call pargd (h)

#  Get the center with Newton

	diff = acc + 1.
	for (i = 1; i <= 50 && abs(diff) > acc; i = i + 1) {
	    xd = a[2] + x0*(2*a[3] + x0*(3*a[4] + x0*4*a[5]))
	    xc = 2*a[3] + x0*(6*a[4] + x0*12*a[5])
	    diff = xd / xc
	    x = x0 - diff
	    if (abs (diff) > acc)
		x0 = x
	    }

	x0 = .5 * (x + x0)
	h = p(x0)
#	call printf ("PKFITQ: center= %f height= %f\n")
#	    call pargd (x0)
#	    call pargd (h)
#	call printf ("PKFITQ: coeff = %f %f %f %f %f\n")
#	    call pargd (a[1])
#	    call pargd (a[2])
#	    call pargd (a[3])
#	    call pargd (a[4])
#	    call pargd (a[5])

#  Now test to see whether the polynomial ever reaches half max
#  and if it does, compute where it is (fudge a0 and look for roots)

	a[1] = a[1] - (0.5d0 * h)
	for (iroot = 1; iroot <= 2; iroot = iroot + 1) {
	    if (iroot == 1)
		x1 = xl
	    else
		x1 = xr
	    step = .1
	    s = 0.1
	    while (s < 3. && p(x0+s*(x1-x0)) >= 0) {
		s = s + step
		}

# No sign change and we've probably got no root

	    if (s >= 3) {
		call printf ("QUARTICPEAK: poly never reaches half maximum\n")
		x = x1
		}

#  Get the zero with Newton's method

	    else {
		i = 1
		x = (s - step/2) * (x1 - x0) + x0
		while (i <= 50 && abs (diff) >= acc) {
		    xd = a[2] + x*(2*a[3] + x*(3*a[4] + x*4*a[5]))
		    diff = p(x) / xd
		    i = i + 1
		    x = x - diff
		    }
		x = .5 * (x1 + x)
		}

	    if (iroot == 1)
		xl = x
	    else
		xr = x
	    }

	a[1] = a[1] + 0.5d0*h
	w = xr - xl
	return
end
# Sep 25 1991	Make height, width, and peakfrac double for PKWIDTH call

# Feb 11 1992	Clarify comments on fixed number of point fit
# Mar 27 1992	Pass velocity vector to PKWIDTH (and change argument order)
# Nov 30 1992	Pass debug as argument; drop inclusion of fquot common

# Apr 13 1994	Drop unused variable nhalf

# Aug 21 1995	Drom DEBUG from PKWIDTH call
