# File rvsao/Util/t_wl2pix.x
# October 7, 2009
# By Doug Mink, SAO Telescope Data Center

include	<imhdr.h>
include	<smw.h>

define	LEN_USER_AREA	100000
define	SZ_FORMAT	8

# T_WL2PIX - Compute the pixel for a specific wavelength in a 
# one dimensional wavelength calibrated image

procedure t_wl2pix()

char	imname[SZ_PATHNAME]	# image name
char	objname[SZ_PATHNAME]	# object name
int	mext			# FITS image extension
int	mspec			# Number of spectrum to read from multispec file
int	mband			# Number of band to read from multispec file
pointer	im			# Image header structure
pointer	sh			# Spectrum header structure
pointer	pix			# pointer to pixel values
int	npts			# Number of pixels in spectrum
double	wl			# Wavelength at pixel in angstroms
double	pixel			# Pixel at wavelength
int	world			# wavelength or pixel
int	mode
double	p1, p2, w1, w2, dw
char	waveform[SZ_FORMAT]
char	pixform[SZ_FORMAT]
bool	debug, verbose, db

double	clgetd()
int	clgeti()
bool	clgetb()
int	strlen()
double	wcs_w2p(), wcs_p2w()

begin


	call clgstr ("spectrum", imname, SZ_PATHNAME)
	if (strlen (imname) < 1) {
	    call printf ("WL2PIX returns the pixel at a given wavelength in a spectrum\n")
	    return
	    }
	mext = clgeti ("specext")
	mspec = clgeti ("specnum")
	mband = clgeti ("specband")
	world = 1

	wl = clgetd ("wavelength")

	verbose = FALSE
	verbose = clgetb ("verbose")
	debug = FALSE
	debug = clgetb ("debug")
	db = FALSE

	call strcpy ("%9.3f", waveform, SZ_FORMAT)
	call clgstr ("waveform", waveform, SZ_FORMAT)

	call strcpy ("%8.3f", pixform, SZ_FORMAT)
	call clgstr ("pixform", pixform, SZ_FORMAT)

	mode = READ_ONLY

	# Open the image
	call getimage (imname, mext, mspec, mband, pix, im, sh, npts, objname,
		       mode, world, db)
	if (im == ERR) {
	    if (debug)
		call printf ("PIX2WL: Cannot read %s\n")
		    call pargstr (imname)
	    return
	    }

	call wcs_set (sh)
	pixel = wcs_w2p (wl)
	call clputd ("pixel", pixel)

	if (verbose) {
	    call printf (pixform)
		call pargd (pixel)
	    call printf ("\n")
	    }
	if (debug) {
	    p1 = pixel - 0.5d0
	    w1 = wcs_p2w (p1)
	    p2 = pixel + 0.5d0
	    w2 = wcs_p2w (p2)
	    dw = w2 - w1
	    call printf ("%s: ")
		call pargstr (imname)
	    call printf (waveform)
		call pargd (wl)
	    call printf (" -> ")
	    call printf (pixform)
		call pargd (pixel)
	    call printf (" (%.4f/pix)\n")
		call pargd (dw)
	    }

	return
end

# Aug 12 2008	New task based on listspec

# Oct  7 2009	Change image file name length from SZ_FNAME to SZ_PATHNAME
