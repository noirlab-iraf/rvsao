# File Util/rebin.x
# July 2, 2008
# Adapted by Doug Mink to read MWCS data
# Adapted by Stephen Levine , University of Wisconsin (madraf::levine)
# from AL_REBIN in T_REBIN (in onedspec)

# REBIN rebins a spectrum to a specified range in wavelength
# REBINL rebins a spectrum to a specified range in log wavelength

include <mach.h>
include <ctype.h>
include <error.h>
include <imhdr.h>
include <fset.h>
include <math/iminterp.h>
include <math/curfit.h>
include <smw.h>

#rebinning defs	
define RB_NEAREST	1	# nearest neighbor
define RB_LINEAR	2	# linear
define RB_POLY3		3	# 3rd order polynomial
define RB_POLY5		4	# 5th order polynomial
define RB_SPLINE3	5	# cubic spline
define RB_SINC		6	# sinc
define RB_LSINC		7	# look-up table sinc
define RB_DRIZZLE	8	# drizzle
define RB_SUMS		9
define RB_FUNCTIONS	"|nearest|linear|poly3|poly5|spline3|sinc|lsinc|drizzle|sums|"

procedure rebinl (imin, sh, indlog, imout, outpix, outw0, outdw, pixshift)

real	imin[ARB]	# Input spectrum
pointer	sh		# Spectrum header structure
double	indlog		# Log wavelength shift for input spectrum
real	imout[ARB]	# Output spectrum
int	outpix		# Number of pixels in output spectrum
double	outw0		# Output starting log wavelength
double	outdw		# Output log wavelength increment
double	pixshift	# Pixel shift in input image

pointer	invert, sp
char	interp_mode[SZ_LINE]
int	inpix		# Number of pixels in input spectrum
double	lw, lw0		# Log wavelength
double	px		# Pixel number (1=center of first)
int	mode, user_mode
int	i
int	clgwrd()
double	wcs_l2p()
#int	fd, open()
#double	inw0		# Input starting wavelength
#double	inw1		# Input final wavelength

begin

# Get rebinning method
	user_mode = clgwrd ("interp_mode", interp_mode, SZ_LINE, RB_FUNCTIONS)
	call cfit_mode (user_mode,mode)

# Make room for inverted solution
	call smark (sp)
	call salloc (invert, outpix+1, TY_DOUBLE)

# Compute pixel position as a function of log lambda.
	call wcs_set (sh)
	call wcs_pixshift (pixshift)

# Interpolate input spectrum to output log lambda spectrum
	inpix = SN(sh)

# For SUMS interpolation, compute upper and lower limits of each pixel
	if (mode == RB_SUMS) {
	    do i = 1, outpix+1 {
		if (indlog != 0.d0) {
		    lw0 = outw0 + (outdw * (double (i) - 1.5))
		    lw = lw0 - indlog
		    }
		else {
		    lw = outw0 + (outdw * (double (i) - 1.5))
		    }
		px = wcs_l2p (lw)
		Memd[invert+i-1] = px
		}
	    call resum (imin, imout, Memd[invert], inpix, outpix)
	    }

# For other interpolations, compute pixel centers
	else {
	    do i = 1, outpix {
		if (indlog != 0.d0) {
		    lw0 = outw0 + (outdw * double (i - 1))
		    lw = lw0 - indlog
		    }
		else {
		    lw = outw0 + (outdw * double (i - 1))
		    }
		px = wcs_l2p (lw)
		Memd[invert+i-1] = px
		}
	    call reinterp (imin, imout, Memd[invert], inpix, outpix, mode)
	    }

	call sfree (sp)

end


procedure rebin (imin, sh, indlog, imout, outpix, outw0, outdw, pixshift, ispix)

real	imin[ARB]	# Input spectrum
pointer	sh		# Spectrum header structure
double	indlog		# Log wavelength shift for input spectrum
real	imout[ARB]	# Output spectrum
int	outpix		# Number of pixels in output spectrum
double	outw0		# Output starting wavelength
double	outdw		# Output wavelength increment
double	pixshift	# Pixel shift in input image
bool	ispix		# TRUE if rebinning in pixel space, else FALSE

pointer	invert, sp
char	interp_mode[SZ_LINE]
int	inpix		# Number of pixels in input spectrum
double	lw		# Log wavelength
double	wl,wl0		# Wavelength in Angstroms
double	px		# Pixel number (1=center of first)
int	mode, user_mode
int	i
int	clgwrd()
double	wcs_w2p()

begin

# Get rebinning method
	user_mode = clgwrd ("interp_mode", interp_mode, SZ_LINE, RB_FUNCTIONS)
	call cfit_mode (user_mode,mode)

# Make room for inverted solution
	call smark (sp)
	call salloc (invert, outpix+1, TY_DOUBLE)

# Compute wavelength at each pixel position
	call wcs_set (sh)
	call wcs_pixshift (pixshift)
	wl0 = 0.d0

# Interpolate input spectrum to output spectrum
	inpix = SN(sh)

# For SUMS interpolation, compute upper and lower limits of each pixel
	if (mode == RB_SUMS) {
	    do i = 1, outpix+1 {
		if (indlog != 0.d0) {
		    wl0 = outw0 + (outdw * (double (i) - 1.5))
		    lw = dlog10 (wl0) - indlog
		    wl = 10.d0 ** lw
		    }
		else
		    wl = outw0 + (outdw * (double (i) - 1.5))
		if (ispix)
		    px = wl
		else
		    px = wcs_w2p (wl)
		Memd[invert+i-1] = px
		}
	    call resum (imin, imout, Memd[invert], inpix, outpix)
	    }

# For other interpolations, compute pixel centers
	else {
	    do i = 1, outpix {
		if (indlog != 0.d0) {
		    wl0 = outw0 + (outdw * double (i - 1))
		    lw = dlog10 (wl0) - indlog
		    wl = 10.d0 ** lw
		    }
		else
		    wl = outw0 + (outdw * double (i - 1))
		if (ispix)
		    px = wl
		else
		    px = wcs_w2p (wl)
		Memd[invert+i-1] = px
#		call printf ("%d: %.4fA -> %.4fA = %.4f: %.4f\n")
#		    call pargi (i)
#		    call pargd (wl0)
#		    call pargd (wl)
#		    call pargd (lw)
#		    call pargd (px)
		}
	    call reinterp (imin, imout, Memd[invert], inpix, outpix, mode)
	    }

	call sfree (sp)

end


# CFIT_MODE -- Transform users input mode to CURFIT package mode.

procedure cfit_mode (user_mode, interp_mode)

int	user_mode, interp_mode

begin
	switch (user_mode) {
	    case RB_NEAREST:
		interp_mode = II_NEAREST

	    case RB_LINEAR:
		interp_mode = II_LINEAR

	    case RB_SPLINE3:
		interp_mode = II_SPLINE3

	    case RB_POLY3:
		interp_mode = II_POLY3

	    case RB_POLY5:
		interp_mode = II_POLY5

	    case RB_SINC:
		interp_mode = II_SINC

	    case RB_DRIZZLE:
		interp_mode = II_DRIZZLE

	    case RB_SUMS:
		interp_mode = RB_SUMS
	    }
end


# RESUM -- Rebinning using a partial pixel summation technique to 
#          preserve total flux.

procedure resum (pixin, pixout, invert, npin, npout)

real	pixin[ARB]	# Input spectrum
real	pixout[ARB]	# Output spectrum (returned)
double	invert[ARB]	# Output to input pixel limit map
int 	npin		# Number of pixels in input spectrum
int 	npout		# Number of pixels in output spectrum

int 	i
double	x1, x2, xlo, xhi, xmin, xmax

real	pixel_parts()

begin
	# Initialize
	xmin = 0.5d0
	xmax = double (npin) + 0.5d0
	x2 = invert [1]
	xlo = 0.0
	xhi = 0.0

	do i = 1, npout {
	    x1 = x2
	    x2 = invert[i + 1]
	    if (x2 < x1) {
		xhi = x1
		xlo = x2
		}
	    else {
		xlo = x1
		xhi = x2
		}

	    # Integrate over the partial pixels from xlo->xhi
	    pixout[i] = pixel_parts (pixin, xlo, xhi)
	    }
end


# PIXEL_PARTS -- Integrate over partial pixels to obtain total flux
# 	         over specified region.

real procedure pixel_parts (y, x1, x2)

real	y[ARB]		# Spectrum vector
double	x1		# Starting fractional pixel
double	x2		# Ending fractional pixel

int 	i, i1, i2, npix
double	sum, frac1, frac2, xi1, xi2

begin
	# Remember that pixel centers occur at intgral values
	# so a pixel extends from i-0.5 to i+0.5

#	call printf ("%.5f - %.5f: \n")
#	    call pargd (x1)
#	    call pargd (x2)
	frac1 = 0.d0
	frac2 = 0.d0

	if (x1 >= x2 || x2 <= x1) {
#	    call printf ("0.0\n")
	    return (0.d0)
	    }

	# Compute first and last input pixel numbers
	i1 = int (x1 + 0.5d0)
	i2 = int (x2 + 0.5d0)

	# Flux from part of single pixel
	if (i1 == i2 && i1 > 0) {
	    frac1 = x2 - x1
	    sum = frac1 * y[i1]
	    }

	else {

	# Flux from first pixel
	    if (i1 > 0) {
		xi1 = double (i1) + 0.5d0
		frac1 = xi1 - x1
		sum = frac1 * y[i1]
		}

	# Flux from last pixel
	    if (i2 > 0) {
		xi2 = double (i2) - 0.5d0
		frac2 = x2 - xi2
		sum = sum + frac2 * y[i2]
		}
	    }

	# Include entire pixels between first and last
	npix = 0
	if (i2 > i1+1) {
	    do i = i1+1, i2-1 {
		sum = sum + y[i]
		npix = npix + 1
		}
	    }
#	call printf ("%fx%d %fx%d %d = %f\n")
#	    call pargd (frac1)
#	    call pargi (i1)
#	    call pargd (frac2)
#	    call pargi (i2)
#	    call pargi (npix)
#	    call pargd (sum)

	return (sum)
end

# REINTERP -- Rebin the vector by interpolation.

procedure reinterp (pixin, pixout, invert, npin, npout, mode)

real	pixin[ARB]	# Input spectrum
real	pixout[ARB]	# Output spectrum (returned)
double	invert[ARB]	# Output to input pixel map
int 	npin		# Number of pixels in input spectrum
int 	npout		# Number of pixels in output spectrum
int 	mode		# Interpolation mode

int	j
real	xpos

real	arieval()

begin
	do j = 1, npout {
	    xpos = invert[j]
	    if (xpos < 1.0 || xpos > real (npin))
		pixout[j] = 0.0
	    else
		pixout[j] = arieval (xpos, pixin, npin, mode)
	    }
end

# Jun  4 1993	Transform wavelength to pixel in subroutine
# Jun 15 1993	Clarify subroutines; add nearest value option
# Jun 16 1993	Change computations to double
# Jun 30 1993	Add spectrum header structure argument
# Jun 30 1993	Perform log-wavelength shift (for templates)
# Jul  1 1993	Keep subroutine arguments same type for Decstations
# Oct  6 1995	Change SHDR_* calls to SPHD_*

# Aug  7 1996	Use WCS_L2P instead of SHDR_PL

# May 16 1997	Add REBIN to rebin to linear wavelength; change file name
# Jul 17 1997	Fix REBIN linear wavelength bug: log10, not log
# Oct  6 1997	Use smw.h from IRAF, not local version

# Feb  9 2001	Update interpolation modes to match IRAF 2.11.3 iminterp
# Feb 12 2001   Make SUMS interpolation option work with negative dispersion
# Feb 14 2001	Rewrite rebin(), resum() and pixel_parts() to use pixel edges

# Jun  7 2005	Rewrite pixel_parts to work correctly on first and last pixels

# Mar 10 2008	Add ispix argument to rebin() for rebinning in pixel space
# Jul  2 2008	Fix rebinl() to compute pixel edges for SUMS interpolation
