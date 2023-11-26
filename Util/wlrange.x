# File rvsao/Util/wlrange.x
# June 15, 2009
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1999-2009 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# WLRANGE finds the starting and stopping wavelengths for the region of overlap
# between a list of spectra
 
include	<imhdr.h>
include <smw.h>
include	<imio.h>
include	<fset.h>
include "rvsao.h"

procedure wlrange (spectra, wl1, wl2, npix)

char	spectra[ARB]	# List of spectrum files
double	wl1		# Maximum first pixel wavelength of spectra in list
double	wl2		# Minimum final pixel wavelength of spectra in list
int	npix		# Maximum number of pixels of spectra in list

pointer	speclist
char	specfile[SZ_PATHNAME]	# Object spectrum file name
char	specpath[SZ_PATHNAME]	# Object spectrum path name
char	specdir[SZ_PATHNAME]	# Directory for object spectra
char	svel_corr[SZ_LINE]	# Type of velocity correction for spectrum
				# (none | file | heliocentric | barycentric)
int	mspec		# Object spectrum to read from multispec file
int	mband		# Object band to read from multispec file

pointer spectrum	# Object spectrum
pointer specim		# Object image header structure
int	ldir
double	tz,zshift

int	nmspec0		# Number of object multispec spectra
int	mspec_range[3,MAX_RANGES]
int	ip,jp,lfile	# Limits for multispec aperture decoding
char	lbracket[3]	# "[({"
char	rbracket[3]	# "])}"
double	wli1, wli2	# Wavelengths of first and last pixels in spectrum
double	wlr1, wlr2	# Rest wavelengths of first and last pixels in spectrum
double	zcomp
double	dindef
double	wcs_p2w()
 
bool	clgetb()
int	clgeti()
double	clgetd()
int	strdic(), stridx(), stridxs()
int	decode_ranges(),get_next_number()
int	imtgetim(), imaccess(), strlen()
int	ispec
bool	rdebug
pointer	imtopen()

define	newspec_ 10
define	newap_	 20
define	endspec_	 90

include	"rvsao.com"
include	"sum.com"
 
begin
	dindef = INDEFD
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	rdebug = debug
	ispec = 0
	npix = 0

# Get task parameters.

# Multispec spectrum numbers (use only first if multiple files)
	if (decode_ranges (specnums, mspec_range, MAX_RANGES, nmspec0) == ERR)
	    call error (1, "WLRANGE: Illegal multispec list")
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	ldir = strlen (specdir)
	if (specdir[1] != EOS && specdir[ldir] != '/') {
	    specdir[ldir+1] = '/'
	    specdir[ldir+2] = EOS
	    }
	mband = clgeti ("specband")

# Print processing information
	rdebug  = clgetb ("debug")
	debug = FALSE

# Type of heliocentric velocity correction to be used
	call clgstr ("svel_corr",svel_corr,SZ_LINE)
	svcor = strdic (svel_corr,svel_corr,SZ_LINE, HC_VTYPES)

# Template object velocity
	velocity = clgetd ("velcomp")
        zcomp = clgetd ("zcomp")
	if (zcomp != dindef)
	    velocity = c0 * zcomp

# Initialize wavelength range
	wl1 = 0.d0
	wl2 = 0.d0

# Get next object spectrum file name from the list
	speclist = imtopen (spectra)
newspec_
#	if (debug) {
#	    call printf ("WLRANGE: About to get next image from list\n")
#	    call flush (STDOUT)
#	    }
	if (imtgetim (speclist, specfile, SZ_PATHNAME) == EOF) {
	    if (debug) {
		call printf ("WLRANGE: No more images in input list\n")
		call flush (STDOUT)
		}
	    call imtclose (speclist)
	    if (debug) {
		call printf ("\n")
		call flush (STDOUT)
		}
	    return
	    }
#	if (debug) {
#	    call printf ("WLRANGE: Next image is %s\n")
#		call pargstr (specfile)
#	    call flush (STDOUT)
#	    }

# Check for specified apertures in multispec spectrum file
	ip = stridxs (lbracket,specfile)
	if (ip > 0) {
	    lfile = strlen (specfile)
	    specfile[ip] = EOS
	    ip = ip + 1
	    jp = 1
	    while (stridx (specfile[ip],rbracket) == 0 && ip <= lfile) {
		specnums[jp] = specfile[ip]
		specfile[ip] = EOS
		ip = ip + 1
		jp = jp + 1
		}
	    specnums[jp] = EOS
	    if (decode_ranges (specnums,mspec_range,MAX_RANGES,nmspec) == ERR)
		call error (1, "WLRANGE: Illegal multispec list")
	    }
	else
	    nmspec = nmspec0
	if (debug) {
	    call printf ("WLRANGE: next file is %s [%s] = %d aps\n")
		call pargstr (specfile)
		call pargstr (specnums)
		call pargi (nmspec)
	    call flush (STDOUT)
	    }

# Check for readability of object spectrum
	call strcpy (specdir,specpath,SZ_PATHNAME)
	call strcat (specfile,specpath,SZ_PATHNAME)
	if (imaccess (specpath, READ_ONLY) == NO) {
	    call eprintf ("WLRANGE: cannot read spectrum file %s \n")
		call pargstr (specpath)
	    go to newspec_
	    }

# Get next multispec number from list
	mspec = -1
newap_
	if (nmspec <= 0)
	    go to newspec_
	if (get_next_number (mspec_range, mspec) == EOF)
	    go to newspec_

	ispec = ispec + 1

# Load spectrum
	ldir = strlen (specdir)
	if (ldir > 0) {
	    if (specdir[ldir] != '/') {
		specdir[ldir+1] = '/'
		specdir[ldir+2] = EOS
		}
	    call strcpy (specdir,specpath,SZ_PATHNAME)
	    call strcat (specfile,specpath,SZ_PATHNAME)
	    }
	else
	    call strcpy (specfile,specpath,SZ_PATHNAME)
	specsh = NULL
	call getspec (specpath, mspec, mband, spectrum, specim, FALSE)
	if (specim == ERR) {
	    call eprintf ("OVERLAP: Error reading spectrum %s\n")
		call pargstr (specpath)
	    go to endspec_
	    }

# Add known velocity components as 1+z's
	if (spvel == dindef)
	    spvel = 0.d0
	if (spechcv == dindef)
	    spechcv = 0.d0
	tz = (1.d0 + (spvel / c0)) /
	     (1.d0 + (spechcv / c0))

# Set z offset for rebinning
	if (velocity == dindef) {
	    zshift = 1.d0
	    }
	else
	    zshift = (1.d0 + (velocity / c0)) / tz
#	if (rdebug) {
#	    call printf ("WLRANGE: spvel= %.2f, spechcv= %.2f, velocity= %.2f")
#		call pargd (spvel)
#		call pargd (spechcv)
#		call pargd (velocity)
#	    call printf (", zshift is %.5f\n")
#		call pargd (zshift)
#	    call flush (STDOUT)
#	    }

# Get wavelengths at first and last pixels
	wli1 = wcs_p2w (double (NP1(specsh)))
	wli2 = wcs_p2w (double (NP2(specsh)))
	if (wli1 < wli2) {
	    wlr1 = wli1 * zshift
	    wlr2 = wli2 * zshift
	    }
	else {
	    wlr1 = wli2 * zshift
	    wlr2 = wli1 * zshift
	    wli1 = wlr1 / zshift
	    wli2 = wlr2 / zshift
	    }

# Check for maximum minimum wavelength
	if (ispec == 1 || wlr1 > wl1)
	    wl1 = wlr1

# Check for minimum maximum wavelength
	if (ispec == 1 || wlr2 < wl2)
	    wl2 = wlr2

	if (SN(specsh) > npix)
	    npix = SN(specsh)

	if (rdebug) {
	    call printf ("WLRANGE: %3d:  im: %.2fA - %.2fA, rest: %.2fA - %.2fA\r")
		call pargi (ispec)
		call pargd (wli1)
		call pargd (wli2)
		call pargd (wlr1)
		call pargd (wlr2)
#	    call printf ("WLRANGE: %3d: lim %.2fA - %.2fA %s\n")
#		call pargi (ispec)
#		call pargd (wl1)
#		call pargd (wl2)
#		call pargstr (specfile)
	    call flush (STDOUT)
	    }

# Move on to next aperture or next image
endspec_
	call close_image (specim, specsh)
	nmspec = nmspec - 1
	if (nmspec > 0)
	    go to newap_
	go to newspec_
end

# Apr 21 1997	New subroutine
# May  2 1997	Always test against dindef, not INDEFD
# May  5 1997	Pass file template for spectrum list
# May 20 1997	Do not shift if output velocity is INDEF
# Jun 18 1997	Use local double ZCOMP instead of real Z0 in common
# Jul 28 1997	Fix bug by initializing the spectrum counter ISPEC

# Mar 20 1998	Fix error messages
# Apr  7 1998	Fix bug by adding mband to GETSPEC call
# Jun 12 1998	Use only portion of spectrum with WCS

# Mar 11 1999	Add write argument to getspec()

# Jun  4 2004	Return maximum number of pixels as well as wavelength range

# Jun 10 2009	Fix range setting for reversed spectra
# Jun 15 2009	Drop most debugging printouts
