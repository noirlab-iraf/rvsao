# File rvsao/Xcor/pkfitc.x
# August 21, 1995
# By Doug Mink,Harvard-Smithsonian Center for Astrophysics
# After Guillermo Torres

# Copyright(c) 1992 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#**************************************************************
# Fit the following function to the observations [xind,xcor]: *
#                                                             *
#                  A cos(Bz)                                  *
#         f(x) = ------------    where  z = x - C             *
#                 1 + D z**2                                  *
#                                                             *
# Array a[4] receives initial values for the fit and returns  *
# the adjusted parameters A, B, C, D.                         *
#                                                             *
# G. Torres (Jan/1989)                                        *
#**************************************************************

procedure pkfitc (npt,xind,xcor,xvel,peakfrac,pkindx,npfit,center,width,height,debug)

int	npt		# Number of data points to fit
real	xind[ARB]	# index vector
real	xcor[ARB]	# cross-correlation vector
real	xvel[ARB]	# velocity vector
double	peakfrac	# Fraction of peak height above which fit is done
int	pkindx		# Bin containing maximum
int	npfit		# number of points used in peak fit (returned)
double	center		# center of peak (returned)
double	width		# width of peak (returned)
double	height		# height of peak (returned)
bool	debug		# True for diagnostic listing

double	a[4],cov[4],denom,ff,dff,pi, dv, psi0, dpsi,g
int	il, ir, ierr, itry

begin
	pi = 3.1415926535897932d0

#  Find peak and limits for fit

	call pkwidth (npt,xcor,xvel,peakfrac,pkindx,height,width,il,ir)

	npfit = ir - il + 1

#  First approximations for B and D

	ierr = 1
	itry = 0
	while (ierr > 0 && itry < 2) {
	    dv = xind[ir] - xind[il]
	    a[1] = height
	    a[2] = pi / dv
	    a[3] = xind[pkindx]
	    a[4] = 16.d0 * (dsqrt(2.d0) - 1.d0) / (dv * dv)

#  This division by 2 is for better convergence (at least in some cases...)

	    a[4] = a[4] / 2.

#  Try different values for parameter D, according to itry: start
#  with D as computed above, and if it fails then try with twice
#  the starting value, and finally with half of the starting value

	    g =  -1.25 * (itry * itry) + 2.25 * itry + 1.
	    a[4] = g * a[4]

#	    if (debug) {
#		call printf ("PKFITC: calling MINI; try= %d\n")
#		    call pargi (itry)
#		}

#  Do fit
	    call xcmin (npfit,xind[il],xcor[il],4,a,debug,cov,ierr)

	    if (ierr > 0) {
		call printf ("PKFITC Error:  singular matrix\n")
		itry = itry + 1
		if (itry > 2) {
		    call printf ("PKFITC Error:  max. iteration exceeded\n")
		    if (ierr > 1) return
		    }
		}
	    }

#  Calculate FWHM of peak
#  (Newton - Raphson method to solve non - linear equation)
#  Maybe this can be improved

	psi0 = dv / 4.
	dpsi = 1.
	while (dpsi > 1.e-4) {
	    denom = 1. + a[4] * psi0 * psi0
	    ff = (cos (a[2] * psi0) / denom) - .5
	    dff =  -a[2] * sin (a[2] * psi0) / denom
	    dff = dff - (2. * a[4] * psi0 * cos (a[2] * psi0) / (denom * denom))
	    dpsi = ff / dff
	    psi0 = abs (psi0 - dpsi)
#	    call printf ("PKFITC: psi= %f dpsi= %f\n")
#		call pargd (psi0)
#		call pargd (dpsi)
	    }
	width = 2. * psi0

	height = a[1]
	center = a[3]

end


#--- Function to be minimized

procedure mfunk (n, x, y, a, chi)

int	n
real	x[ARB]
real	y[ARB]
double	a[4]
double	chi

double psi,f
int	i

begin

	chi = 0.
	for (i = 1; i <=n; i = i + 1) {
	    psi = x[i] - a[3]
	    f = a[1] * cos(a[2] * psi) / (1. + a[4] * psi**2)
	    chi = chi + (y[i] - f) ** 2
	    }
	chi = chi / (n - 4)

end
# Sep 25 1991	Make height, width, and peakfrac double for PKWIDTH call
# Oct 23 1991	Add debug argument to xcmin call which now calls SPP program

# Mar 27 1992	Pass velocity vector to PKWIDTH (and change argument order)
# Nov 30 1992	Pass debug as argument; drop inclusion of fquot common

# Aug 21 1995	Drop DEBUG from PKWIDTH call
