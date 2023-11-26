# File rvsao/Emvel/emfind.x
# September 30, 1997
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Bill Wyatt's C and Forth code

# Copyright(c) 1997 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.

#--- Find all lines within section of spectrum which meet criteria

procedure emfind (fdata, fcont, npts, ipix1, ipix2, npfit, lsigma, center
		  height, maxlines, nlines, avheight, debug)

real	fdata[ARB]		# Spectrum (continuum subtracted)
real	fcont[ARB]		# Continuum fit to spectrum
int	npts			# Number of pixels in spectrum
int	ipix1, ipix2		# First and last pixels to search
int	npfit			# Number of points to fit
double	lsigma			# Number of sigma for line acceptance
double	center[ARB]		# Centers for lines
double	height[ARB]		# Heights for lines
int	maxlines		# Maximum number of lines that can be returned
int	nlines			# Number of lines returned
double	avheight		# Average line height (returned)
bool	debug 

int     i
int	nl			# number of lines found
double	avh
bool	isline()		# determine whether at center of line
bool	linepar()		# compute line parameters
double	lcenter,lheight,lwidth	# single line center, height, width
int	ipmin, ipmax		# First and last fittable line centers

begin
	if (debug) {
	    call printf("EMFIND: searching pixels %d to %d, max= %d \n")
		call pargi (ipix1)
		call pargi (ipix2)
		call pargi (maxlines)
	    call flush (STDOUT)
	    }
	nl = 0
	avh = 0.d0

	ipmin = ipix1 + npfit
	ipmax = ipix2 - npfit
	do i = ipmin, ipmax {
	    if (isline (i,fdata,fcont,npts,npfit,lsigma,debug)) {
		if (linepar (i,fdata,npts,npfit,lcenter,lheight,lwidth,debug)) {
		    if (nl < maxlines) {
			nl = nl + 1
			center[nl] = lcenter
			height[nl] = lheight
			avh = avh + lheight
			}
		    if (debug) {
			if (nl == 1)
			    call printf("    Center    Height   Width\n")
			call printf(" %9.2f %9.2f %7.2f\n")
			    call pargd (lcenter)
			    call pargd (lheight)
			    call pargd (lwidth)
			}
		    }
		}
	    }
	nlines = nl
	if (nl > 0)
	    avheight = avh / double (nl)
	else
	    avheight = 0.d0
	return
end


#  Given a pixel index, returns 1 if:
#    1) the pixel value is the max of the npts+1 neighbors on each side.
#    2) the sqrt(average counts in the 3 peak pixels) * sigma are <= height.

bool procedure isline (xval, fdata, fcont, npts, npfit, lsigma, debug)

int	xval		# index of center pixel
real	fdata[ARB]	# Spectrum with continuum removed
real	fcont[ARB]	# Continuum fit to spectrum
int	npts		# number of pixels in spectrum
int	npfit		# number of pixels to check for peak (/2)
double	lsigma		# number of sigma for line acceptance
double	sigma		# sigma at line
bool	debug

int	i, i1, i2
real	ymax		# Pixel value at line being checked
real	ycont		# Continuum for line +-1
real	dp		# Number of pixels for sigma
double	sqrt()

begin

	ymax = fdata[xval]

#	if (debug) {
#	    call printf ("ISLINE: checking line at pixel %d = %f\n")
#		call pargi (xval)
#		call pargr (ymax)
#	    }

#  If peak is to the left of this pixel, return
	do i = 1, npfit {
	    if (fdata[xval-i] >= ymax) {
		return FALSE
		}
	    }

#  If peak is to the right of this pixel, return
	do i = 1, npfit {
	    if (fdata[xval+i] >= ymax) {
		return FALSE
		}
	    }

#  Compute sigma for region = sqrt (total counts)
        i1 = xval - npfit
        if (i1 < 1)
            i1 = 1
        i2 = xval + npfit
        if (i2 > npts)
	    i2 = npts
	ymax = 0.0
	ycont = 0.0
	dp = 0
        do i = i1, i2 {
	    ymax = ymax + fdata[i]
	    ycont = ycont + fdata[i]
	    dp = dp + 1.0
	    }
	if (dp > 0.d0) {
	    ymax = ymax / dp
	    ycont = ycont / dp
            sigma = sqrt (abs (ymax + ycont))
	    }
	else {
#	    if (debug) {
#		call printf ("ISLINE: No points at %d +- %d (1-%d)\n")
#		    call pargi (xval)
#		    call pargi (npfit)
#		    call pargi (npts)
#		}
	    return FALSE
	    }

#	if (debug) {
#	    call printf ("ISLINE %d: line value= %8.3f continuum= %8.3f nsigma= %8.3f\n")
#		call pargi (xval)
#		call pargr (ymax)
#		call pargr (ycont)
#		call pargd (sigma*lsigma)
#	    }
	
        if (ymax >= sigma * lsigma)
	    return TRUE
	else if (ymax > abs (ycont * 1.e10))
	    return TRUE
	else
	    return FALSE
end


#  from the FORTH program by J. Tonry: 
# 
#  If n = npfit (so 2n + 1 points are fit)
#    f(x) = h + a (x-x')**2 
#         = h [ 1 - 2 { (x-x')/w }**2 ]
# 
#  also,  a, b, c are:
#    a = -45/{n(n+1)(2n-1)(2n+1)(2n+3)} * Sum(i=-n, n) {n(n+1)/3 - i**2}yi
#    b = -3/{n(n+1)(2n+1)} * Sum(i=-n,n) {i * yi}
#    c = 1/(2n+1) * Sum(i=-n,n) yi
# 
#  Simplifying, if:
#    u = 3 / { n (n+1) }
#    v =  2n + 1 
#    w = (2n-1)(2n+3) = 4n**2 + 4n - 3
# 
#  Then:
#    a = -15 u / { v (2n-1)(2n+3)} * Sum(i=-n,n) { yi/u - yi * i**2 }
#    a = -u/v * Sum(-n, n){ i * yi }
#    a = Sum(-n, n){ yi } / v
# 
#  Then:
#    x' = center = x0 + b/2a
#    h  = height = c - a [ n(n+1)/3 + (b/2a)**2 ]
#    w  = width  = sqrt(-2h/a)

# returns true for success, false for error

bool procedure linepar (xzero, fdata, npts, npfit, center, height, width, debug)

int	xzero		# Initial guess of line (should be peak!)
real	fdata[ARB]	# Spectrum data to fit
int	npts		# Number of pixels in spectrum
int	npfit		# +/- number of pixels to fit
double	center		# line center in pixels (returned)
double	height		# line height (returned)
double	width		# line width in pixels (returned)
bool	debug		# true for diagnostic printing

int	i, i1, i2
double	a, b, c
double	u, v, w
double	sumy, sumyi, sumyii
double	y, dn, di

begin

	if (npfit <= 0 || xzero < npfit) {
#	    if (debug) {
#		call printf ("LINEPAR: npfit= %d xzero= %d\n")
#		    call pargi (npfit)
#		    call pargi (xzero)
#		}
	    return FALSE
	    }

	sumy   = 0.d0
	sumyi  = 0.d0
	sumyii = 0.d0
	i1 = xzero - npfit
	if (i1 < 1)
	    i1 = 1
	i2 = xzero + npfit
	if (i2 > npts)
	    i2 = npts
	dn = 0.d0
	do i = i1, i2 {
	    y = fdata[i]
	    di = double (i - xzero)
	    sumy   = sumy + y
	    sumyi  = sumyi + (y * di)
	    sumyii = sumyii + (y * di * di)
	    dn = dn + 1
	    }

	if (dn > 0.d0) {
	    u = 3.d0 / (dn * dn + dn)
	    v = 2.d0 * dn + 1.d0
	    w = ((4.d0 * dn) + 4.d0) * dn - 3.d0

	    a = (-15.d0 * u) / (v * w) * (sumy / u - sumyii)
	    b = (-u * sumyi) / v
	    c = sumy / v

	    w = b / ( 2.d0 * a )
	    center    = xzero + w
	    height    = c - a * (1.d0 / u + (w * w))
	    width     = sqrt (abs (-2.d0 * height / a))
	    }
	else {
#	    if (debug) {
#		call printf ("LINEPAR: No points at %d +- %d (1-%d)\n")
#		    call pargi (xzero)
#		    call pargi (npfit)
#		    call pargi (npts)
#		}
	    return FALSE
	    }

#	if (debug) {
#	    call printf ("LINEPAR: line found at pixel %d -> %8.3f\n")
#		call pargi (xzero)
#		call pargd (center)
#	    }

	return TRUE
end

# Dec  4 1991	Use interactive continuum fit to whole spectrum instead
#		of fitting local continuum for each line

# Mar 27 1992	Print diagnostic as first executable code

# Apr 13 1994	Drop avwidth as it is not used
# Apr 13 1994	Use npts instead of constant 3 in isline
# Apr 22 1994	Drop avw as it is not used
# Apr 25 1994	Limit search to fittable pixels only

# Sep 30 1997	Pass total number of points in spectrum as argument
# Sep 30 1997	Compute sigma from 2*NPFIT+1 points, not 3 in ISLINE
# Sep 30 1997	Check spectrum limits when doing fit in LINEPAR
