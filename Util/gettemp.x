# File rvsao/Util/gettemp.x
# July 08, 2010
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss, Johns Hopkins University

# Copyright(c) 1992-2010 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# GETTEMP opens one template spectrum, extracts wavelength and velocity
# information from the header, and returns the template spectrum.
# Data is shared in common in "rvsao.com"
 
include <imhdr.h>
include	"rvsao.h"
include	<smw.h>
define  LEN_USER_AREA   2880
 
procedure gettemp (tempfile,mspec,mband,tempdir,tempim,tempsh,itemp)
 
char	tempfile[ARB]	# Template file name
int	mspec		# Spectrum aperture to read from multispec or 2-D file
int	mband		# Spectrum band to read from multispec file
char	tempdir[ARB]	# Template directory
pointer	tempim		# Template image header structure (returned)
pointer	tempsh		# Template spectrum header structure (returned)
int	itemp		# Template number for saved information

int	world		# 1=wavelength 2=pixel
int	mextt		# Image extension number

char	temppath[SZ_PATHNAME]
int	mode, nline, npix
double	bcvel,hcvel,dindef,w1,w2, pixshift
char	tempnm[SZ_PATHNAME]
real	compbcv()
double	wcs_p2l(), wcs_p2w()
pointer	tempspec	# Template spectrum

include	"rvsao.com"

begin
	mode = READ_ONLY
	mextt = 0

	call strcpy (tempdir,temppath,SZ_PATHNAME)
	call strcat (tempfile,temppath,SZ_PATHNAME)
	if (debug) {
	    call printf ("GETTEMP: %d %s %d %d\n")
		call pargi (itemp)
		call pargstr (temppath)
		call pargi (mspec)
		call pargi (mband)
	    }
	if (correlate == COR_PIX)
	    world = 2
	else
	    world = 1

# Get template spectrum and image header
	call getimage (temppath,mextt,mspec,mband,tempspec,tempim,tempsh,npix,
		       tempnm, mode, world, debug)

	if (tempim == ERR) {
	    call printf ("GETTEMP: Cannot read %s\n")
		call pargstr (tempname[1,itemp])
	    return
	    }
	call strcpy (tempnm, tempname[1,itemp], SZ_PATHNAME)
	if (debug) {
	    call printf ("GETTEMP: %s  %d pixels read\n")
		call pargstr (tempname[1,itemp])
		call pargi (npix)
	    call printf ("GETTEMP: %f %f %f\n")
		call pargr (Memr[tempspec])
		call pargr (Memr[tempspec+1])
		call pargr (Memr[tempspec+2])
	    }

	temppix[itemp] = npix

# Set filter flag
	call imgipar (tempim, "FI-FLAG ", tempfilt[itemp])

#  Get radial velocity and error
	tempvel[itemp] = 0.d0
	tempshift[itemp] = 0.d0
	dindef = INDEFD
	bcvel = dindef
	hcvel = dindef
	call vthead (nline, tempim, tempvel[itemp], bcvel, hcvel, debug)
	call imgdpar (tempim, "TSHIFT  ", tempshift[itemp])

# Calculate template heliocentric velocity, or use one from file
	temphcv[itemp] = 0.d0
	switch (tvcor) {
	    case HCV:
		temphcv[itemp] = compbcv (tempim,0,debug)
	    case BCV:
		temphcv[itemp] = compbcv (tempim,1,debug)
	    case FHCV:
		if (hcvel != dindef)
		    temphcv[itemp] = hcvel
		else if (bcvel != dindef)
		    temphcv[itemp] = bcvel
		else
                    temphcv[itemp] = compbcv (tempim,0,debug)
	    case FBCV:
		if (bcvel != dindef)
		    temphcv[itemp] = bcvel
		else if (hcvel != dindef)
		    temphcv[itemp] = hcvel
		else
                    temphcv[itemp] = compbcv (tempim,1,debug)
	    case NONE:
		temphcv[itemp] = 0.d0
	    }

	if (debug) {
	    call printf ("GETTEMP: vel: %.3f tsh: %.3f hcv: %.3f\n")
		call pargd (tempvel[itemp])
		call pargd (tempshift[itemp])
		call pargd (temphcv[itemp])
	    }


#  Set wavelength and log-wavelength limits for template spectrum
	call wcs_set (tempsh)
	pixshift = 0.0
	call imgdpar (tempim, "PIXSHIFT  ", pixshift)
	call wcs_pixshift (pixshift)
	tempdc[itemp] = DC(tempsh)
	w1 = wcs_p2w (1.d0)
	if (w1 > 0.d0)
	    testrt[itemp] = wcs_p2l (1.d0)
	else {
	    call printf ("GETTEMP: Illegal starting wavelength %.3f in %s\n")
		call pargd (w1)
		call pargstr (tempfile)
	    call close_image (tempim, tempsh)
	    tempim = ERR
	    return
	    }

	w2 = wcs_p2w (double (npix))
	if (w2 > 0.d0)
	    tefnsh[itemp] = wcs_p2l (double (npix))
	else {
	    call printf ("GETTEMP: Illegal ending wavelength %.3f in %s\n")
		call pargd (w2)
		call pargstr (tempfile)
	    call close_image (tempim, tempsh)
	    tempim = ERR
	    return
	    }
	if (debug) {
	    call printf ("GETTEMP: lambda %9.4f - %9.4f by %9.4f")
		call pargd (w1)
		call pargd (w2)
		call pargr (WP(tempsh))
	    if (DC(tempsh) != DCLOG) {
		call printf (" = %d points\n")
		    call pargi (temppix[itemp])
		}
	    else
		call printf ("\n")
	    dlogw = (tefnsh[itemp] - testrt[itemp]) / double (temppix[itemp] - 1)
	    call printf ("GETTEMP: log lambda %9.7f - %9.7f by %9.7f")
		call pargd (testrt[itemp])
		call pargd (tefnsh[itemp])
		call pargd (dlogw)
	    if (DC(tempsh) == DCLOG) {
		call printf (" = %d points\n")
		    call pargi (temppix[itemp])
		}
	    else
		call printf ("\n")
	    }
	return
end
# May 22 1992	Do not change mspec
# Dec  2 1992	Print log wavelength limits in debug mode

# May 24 1993	Set up MWCS pixel<->wavelength transformations
# May 25 1993	Get all velocity information directly from image header
# May 25 1993	Read spectrum using getimage
# Jun  4 1993	Use onedspec spectrum header structure
# Jun 15 1993	Fix use of onedspec spectrum header
# Jun 30 1993	Return spectrum header as argument
# Jul  7 1993	Add spectrum header to getimage arguments
# Aug 20 1993	Add more diagnostic output

# Mar  9 1994	Fix bug listing wavelength incremen
# Apr 12 1994	Return MWCS header pointer
# Apr 13 1994	Fix parg calls to match types
# Apr 20 1994	Fix SSCAN call
# May  3 1994	Use no temporary storage
# Jun 15 1994	Set TEMPFILT from spectrum image header
# Jun 23 1994	Keep MWCS pointer in SHDR structure
# Aug  3 1994	Change common and header from fquot to rvsao
# Dec  7 1994	Make TEMPFILT integer instead of boolean

# Feb 27 1995	Fix null directory
# Mar 13 1995	Remove unused STRLEN function declaration
# May 10 1995	Print velocities to meters/sec
# May 15 1995	Change all sz_fname to sz_pathname, which is longer
# Jul 13 1995	Add DEBUG to COMPBCV calls

# Jun 20 1996	Print error message if template cannot be read
# Jun 24 1996	Pass string argument to GETIMAGE, not array
# Aug  7 1996	Use smw.h
# Aug 16 1996	Make error messages more informative

# Jan 28 1997	Fix bug which overwrote template titles
# Aug 27 1997	Add band argument for multispec spectra

# Feb  9 1998	Always return BCV/HCV of 0.0 if not set
# Feb 13 1998	Call VTHEAD to extract velocity from template spectrum header
# Apr 22 1998	Drop use of getim.com; extracted needed parameters from  header locally

# Sep 13 2000	Add option to work in pixel as well as wavelength space

# Aug  1 2003	Add image extension argument to getimage()

# May 25 2005	Set pixel shift from header keyword PIXSHIFT

# Mar 13 2009	Stop adding .imh to file name debugging listing

# Jul 08 2010	Use mextt for FITS extension instead of mext, which is now in rvsao.com
