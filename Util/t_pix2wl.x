# File rvsao/Util/t_pix2wl.x
# October 7, 2009
# By Doug Mink, SAO Telescope Data Center

include	<imhdr.h>
include	<smw.h>

define	LEN_USER_AREA	100000
define	SZ_FORMAT	8

# T_PIX2WL - Compute the wavelength for a specific pixel in a 
# one dimensional wavelength calibrated image

procedure t_pix2wl()

char	imname[SZ_PATHNAME]	# image name
char	objname[SZ_PATHNAME]	# object name
int	mext			# FITS image extension
int	mspec			# Number of spectrum to read from multispec file
int	mband			# Number of band to read from multispec file
pointer	im			# Image header structure
pointer	sh			# Spectrum header structure
pointer	pix			# pointer to pixel values
int	npts			# Number of pixels in spectrum
double	pixel			# Pixel for wavelength in angstroms
double	wl			# Wavelength at pixel in angstroms
int	world			# wavelength or pixel
char	waveform[SZ_FORMAT]
char	pixform[SZ_FORMAT]
double	p1, p2, w1, w2, dw, dindef
bool	debug, verbose, db
int	mode

double	clgetd()
int	clgeti()
bool	clgetb()
int	strlen()
double	wcs_p2w()


begin
	dindef = INDEFD

	call clgstr ("spectrum", imname, SZ_PATHNAME)
	if (strlen (imname) < 1) {
	    call printf ("PIX2WL returns the wavelength at a given pixel in a spectrum\n")
	    return
	    }
	mext = clgeti ("specext")
	mspec = clgeti ("specnum")
	mband = clgeti ("specband")
	world = 1

	pixel = clgetd ("pixel")

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

# If pixel greater than 0, print wavelength at that pixel
	if (pixel > 0.0) {
	    wl = wcs_p2w (pixel)
	    call clputd ("wavelength", wl)
	    call clputd ("wave2", dindef)
	    if (verbose) {
		call printf (waveform)
		    call pargd (wl)
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
		call printf (pixform)
		call pargd (pixel)
		    call printf (" -> ")
		call printf (waveform)
		    call pargd (wl)
		call printf (" (%.4f/pix)\n")
		    call pargd (dw)
		}
	    }

# If pixel less than or equal to 0, print wavelength range of spectrum
	else {
	    p1 = 0.5d0
	    w1 = wcs_p2w (p1)
	    call clputd ("wavelength", w1)
	    p2 = double (npts) + 0.5d0
	    w2 = wcs_p2w (p2)
	    call clputd ("wave2", w2)
	    if (verbose) {
		call printf (waveform)
		    call pargd (w1)
		call printf (" - ")
		call printf (waveform)
		    call pargd (w2)
		call printf ("\n")
		}
	    if (debug) {
		dw = (w2 - w1) / (p2 - p1 + 1.d0)
		call printf ("%s: ")
		    call pargstr (imname)
		call printf (waveform)
		    call pargd (w1)
		call printf ("-")
		call printf (waveform)
		    call pargd (w2)
		call printf (" (%d pixels, %.4f/pix)\n")
		    call pargi (npts)
		    call pargd (dw)
		}
	    }

	return
end

# Aug 12 2008	New task based on listspec

# Oct  7 2009	Change image file length from SZ_FNAME to SZ_PATHNAME
