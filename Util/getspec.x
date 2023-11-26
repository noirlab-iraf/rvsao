# File rvsao/Util/getspec.x
# July 9, 2010
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss, Johns Hopkins University

# Copyright(c) 1989-2010 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
#  GETSPEC opens the image specified by specfile and returns pointers
#  to the data and the image descriptor.  Data is shared in "rvsao.com"
 
include	"rvsao.h"
include	<smw.h>

procedure getspec (specfile, mspec, mband, spectrum, specim, spwrite)
 
char	specfile[ARB]	# Data file name
int	mspec		# Number of spectrum to read from multispec file
int	mband		# Multispec band
pointer	spectrum	# Spectrum data (returned)
pointer	specim		# Image header structure (returned)
bool	spwrite		# TRUE if writing to header and/or image

int	world		# 1=wavelength 2=pixel
 
int	npix, ap, nzpix
char	tstr[SZ_PATHNAME]
char	slash
 
real	compbcv()
double	wcs_p2l(), wcs_p2w()
double	bcvel,hcvel,dindef,w1,w2
double	pixshift
int	mode, i, ipix1, ipix2, ip1, ip2
int	imaccess(), strldx()
 
include	"rvsao.com"
 
begin
	dindef = INDEFD
	slash = '/'

# Assume access to first extension of multi-extension FITS file

	if (debug) {
	    call printf ("GETSPEC: %s %d\n")
		call pargstr (specfile)
		call pargi (mspec)
	    }

	if (correlate == COR_PIX)
	    world = 2
	else
	    world = 1

# Check spectrum file for writeability before opening it
	mode = READ_ONLY
	if (spwrite) {
	    mode = READ_WRITE
	    if (imaccess (specfile, READ_WRITE) == NO) {
		call eprintf ("GETSPEC: cannot write to %s; not saving results\n")
		    call pargstr (specfile)
        	mode = READ_ONLY
		}
	    }

# Read spectrum
	call getimage (specfile,mext,mspec,mband,spectrum,specim,specsh,npix,
		       specname,mode,world,debug)
	if (specim == ERR)
	    return

# Set up spectrum line and aperture in specid
	i = strldx (slash, specfile)
	if (i > 0)
	    i = i + 1
	else
	    i = 1
	call strcpy (specfile[i],specid,SZ_PATHNAME)
	if (mspec > 0) {
	    ap = AP(specsh)
	    if (ap != mspec) {
		call sprintf (tstr,SZ_PATHNAME,"[%d ap%d]")
		    call pargi (mspec)
		    call pargi (ap)
		call strcat (tstr, specid, SZ_PATHNAME)
		}
	    else {
		call sprintf (tstr,SZ_PATHNAME,"[%d]")
		    call pargi (mspec)
		call strcat (tstr, specid, SZ_PATHNAME)
		}
	    }

# Set number of spectrum pixels in labelled common
	specpix = npix
	if (debug) {
	    call printf ("GETSPEC: %s  %d pixels read\n")
		call pargstr (specid)
		call pargi (specpix)
	    }

# Set up pixel <-> wavelength transformation
	call wcs_set (specsh)

# Adjust pixel <-> wavelength transformation
        call imgdpar (specim, "PIXSHIFT", pixshift)
	call wcs_pixshift (pixshift)

# Set wavelength and log-wavelength limits of spectrum
	specdc = DC(specsh)
        call imgipar (specim,"DC-FLAG ",specdc)
	w1 = wcs_p2w (double (NP1(specsh)))
	if (w1 > 0.d0)
	    spstrt = wcs_p2l (double (NP1(specsh)))
	else {
	    call eprintf ("GETSPEC: Illegal starting wavelength %.2f\n")
		call pargd (w1)
	    call close_image (specim, specsh)
	    specim = ERR
	    return
	    }
	w2 = wcs_p2w (double (NP2(specsh)))
	if (w2 > 0.d0)
	    spfnsh = wcs_p2l (double (NP2(specsh)))
	else {
	    call eprintf ("GETSPEC: Illegal ending wavelength %.2f\n")
		call pargd (w2)
	    call close_image (specim, specsh)
	    specim = ERR
	    return
	    }
	if (debug) call printf ("GETSPEC: about to read velocity info\n")

# Set bad pixels at start of spectrum to zero and change wavelength limit
	ipix1 = 1
	do i = NP1(specsh), NP1(specsh)+300 {
	    if (Memr[spectrum+i-1] < MIN_PIXEL_VALUE) {
		Memr[spectrum+i-1] = 0.0
		ipix1 = i + 1
		spstrt = wcs_p2l (double (ipix1))
		if (debug) {
		    call printf ("GETSPEC:  Pixel %d = %.2f = %.7g set to zero\n")
			call pargi (i)
			call pargd (wcs_p2w (double(i)))
			call pargr (Memr[spectrum+i-1])
		    call printf ("GETSPEC:  Start of spectrum is now %.2f\n")
			call pargd (wcs_p2w (double(ipix1)))
		    }
		}
	    else
		break
	    }

# Set bad pixels at end of spectrum to zero and change wavelength limit
	ipix2 = NP2(specsh)
	do i = NP2(specsh), NP2(specsh)-300, -1 {
	    if (Memr[spectrum+i-1] < MIN_PIXEL_VALUE) {
		Memr[spectrum+i-1] = 0.0
		ipix2 = i - 1
		spfnsh = wcs_p2l (double (ipix2))
		if (debug) {
		    call printf ("GETSPEC:  Pixel %d = %.2f = %.7g set to zero\n")
			call pargi (i)
			call pargd (wcs_p2w (double(i)))
			call pargr (Memr[spectrum+i-1])
		    call printf ("GETSPEC:  End of spectrum is now %.2f\n")
			call pargd (wcs_p2w (double(ipix2)))
		    }
		}
	    else
		break
	    }

# Interpolate across bad pixels anywhere else in the spectrum
	ip1 = 0
	ip2 = 0
	nzpix = 0
	do i = ipix1, ipix2 {
	    if (Memr[spectrum+i-1] < MIN_PIXEL_VALUE) {
		if (ip1 == 0)
		    ip1 = i
		ip2 = i
		if (debug) {
		    call printf ("GETSPEC:  Fill bad pixel %d = %.2f = %.7g\n")
			call pargi (i)
			call pargd (wcs_p2w (double(i)))
			call pargr (Memr[spectrum+i-1])
		    }
		}
	    else if (ip1 > 0) {
		call fillpix (npix, Memr[spectrum], ip1, ip2)
		if (debug) {
		    call printf ("GETSPEC:  Fill bad pixels from %d to %d ")
			call pargi (ip1)
			call pargi (ip2)
		    call printf ("= %.2f to %.2f\n")
			call pargd (wcs_p2w (double(ip1)))
			call pargd (wcs_p2w (double(ip2)))
		    }
		ip1 = 0
		ip2 = 0
		}

# Add to count of non-zero pixels in spectrum
	    if (Memr[spectrum+i-1] != 0.0) {
		nzpix = nzpix + 1
		}
	    }

# If there are no non-zero pixels in the spectrum give up on it
	if (nzpix == 0) {
	    call eprintf ("GETSPEC: Entire spectrum is zeroes\n")
	    call close_image (specim, specsh)
	    specsh = NULL
	    specim = ERR
	    return
	    }

# Get radial velocity and error
	call vrhead (mspec, specim, bcvel, hcvel)
	if (debug) {
	    call printf ("GETSPEC: hcvel = %.3f, bcvel = %.3f, svcor = %d\n")
		call pargd (hcvel)
		call pargd (bcvel)
		call pargi (svcor)
	    }

# Calculate heliocentric velocity, or use one from file
	spechcv = dindef
	switch (svcor) {
	    case NONE:
		spechcv = 0.d0
	    case HCV:
		spechcv = compbcv (specim,0,debug)
	    case BCV:
		spechcv = compbcv (specim,1,debug)
		specvb = TRUE
	    case FHCV:
		if (hcvel != dindef) {
		    spechcv = hcvel
		    specvb = FALSE
		    }
		else if (bcvel != dindef) {
		    spechcv = bcvel
		    specvb = TRUE
		    }
		else {
		    spechcv = compbcv (specim,0,debug)
		    specvb = FALSE
		    }
	    case FBCV:
		if (bcvel != dindef) {
		    spechcv = bcvel
		    specvb = TRUE
		    }
		else if (hcvel != dindef) {
		    spechcv = hcvel
		    specvb = FALSE
		    }
		else {
		    spechcv = compbcv (specim,1,debug)
		    specvb = TRUE
		    }
	    }
	if (spechcv == dindef)
	    spechcv = 0.d0

# Debugging information about this spectrum
	if (debug) {
	    call printf ("GETSPEC: vel: %.3f hcv: %.3f bcv? %b\n")
		call pargd (spvel)
		call pargd (spechcv)
		call pargb (specvb)
	    call flush (STDOUT)

	    call printf ("GETSPEC: lambda %9.4f - %9.4f by %9.4f")
		call pargd (w1)
		call pargd (w2)
		call pargr (WP(specsh))
	    if (DC(specsh) != DCLOG) {
		call printf (" = %d points\n")
		    call pargi (specpix)
		}
	    else
		call printf ("\n")
	    call flush (STDOUT)

	    dlogw = (spfnsh - spstrt) / double (specpix - 1)
	    call printf ("GETSPEC: log lambda %9.7f - %9.7f by %9.7f")
		call pargd (spstrt)
		call pargd (spfnsh)
		call pargd (dlogw)
	    if (DC(specsh) == DCLOG) {
		call printf (" = %d points\n")
		    call pargi (specpix)
		}
	    else
		call printf ("\n")
	    call flush (STDOUT)
	    }

end
 
# June	1987	Gerard Kriss
# Sept	1989	Doug Mink--add parameter for 2-D files
# July	1990	Doug Mink--add HCV
# Sept	1990	Doug Mink--add velocity parameters
# Oct 	1990	Doug Mink--add hcv or bcv from file
# June	1991	Doug Mink--open image READ_WRITE when needed

# Feb 14 1992	Call new getimage which can handle multispec files
# Apr 20 1992	Pass mspec as argument

# May 24 1993	Add MWCS pixel<->wavelength transformation
# May 25 1993	Get radial velocity information directly from the image header
# Jun  2 1993	Use default velocity from mwcs
# Jun  3 1993	Use onedspec spectrum header structure
# Jun 16 1993	Fix use of onedspec spectrum header
# Jul  7 1993	Add spectrum header to getimage arguments
# Nov 23 1993	Read cross-correlation velocity error BEFORE R-value
# Dec  1 1993	Fix bug setting multispec keyword names

# Mar  8 1994	Fix bug listing wavelength increment
# Apr 12 1994	Return MWCS header pointer
# May 23 1994	Use header-save flag from labelled common
# Jun 15 1994	Remove unused variables
# Jun 22 1994	Save all previous velocity results
# Jun 23 1994	Keep MWCS pointer in SHDR pointer
# Aug  3 1994	Change common and header from fquot to rvsao
# Dec  7 1994	Read velocity information using XCRHEAD and EMRHEAD
# Dec  8 1994	Initialize mwcs here instead of in main programs

# Mar 28 1995	Return name of file if error occurs
# May 10 1995	Print velocities to meters/sec
# May 26 1995	Handle barycentric correction for multispec spectra
# Jun 19 1995	Set ID string with file, line, and aperture (if not =line)
# Jul  7 1995	Fix bug switching between BCV and HCV header keywords
# Jul 13 1995	Add DEBUG to COMPBCV calls
# Aug  4 1995	If BCV and HCV aren't in header, compute BCV or HCV
# Sep 21 1995	Open file read_write unless it cannot be written

# Aug  7 1996	Use smw.h

# Apr  9 1997	EPRINTF error messages; do not print message if GETIMAGE error
# Apr 18 1997	Fix handling of velocity and bcv flag
# May  1 1997	Assign INDEFD, not INDEF to dindef
# Aug 27 1997	Add argument for multispec band
# Oct  6 1997	Flush stdout after each diagnostic line
# Dec 19 1997	Use BCV and HCV keywords even if mspec > 1

# Feb 13 1998	Fix use of BCV and HCV in multispec files
# Apr 22 1998	Drop use of getim common; pass everything in rvsao.com
# May  5 1998	Use SPECID to trick compiler into working
# Jun 11 1998	Use pixel limits from SHDR

# Mar 11 1999	Add argument spwrite to indicate whether file will be written
# Sep 15 1999	Do not read old emission line and correlation results

# Sep 13 2000	Add option to work in pixel as well as wavelength space

# Aug  7 2001	Drop directory from specid

# Feb  5 2002	For fiber spectra with same WCS, drop end points with bad values
# Feb  6 2002	Interpolate across any other pixels with bad values

# Jul 10 2003	Read 0th extension of multi-extension FITS files (for now)
# Aug  1 2003	Add argument for image extension to getimage()

# Aug 24 2004	Double integer argument to wcs_p2w()

# May 23 2005	Add PIXSHIFT to pixel when converting to and from wavelength

# Jul 08 2010	Use mext from rvsao.com instead of setting it to zero
# Jul 09 2010	Return error if spectrum has no non-zero pixels
