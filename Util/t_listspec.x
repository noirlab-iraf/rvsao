# File procdata/Src/t_listspec.x
# June 10, 2011
# By Doug Mink, SAO Telescope Data Center

include	<imhdr.h>
include	<smw.h>
include "rvsao.h"

define	LEN_USER_AREA	100000
define	SZ_FORMAT	8

# T_LISTSPEC - Compute and print the pixel value and corresponding wavelength
# for all the pixels in a list of one dimensional wavelength calibrated
# images.

procedure t_listspec()

char	images[SZ_PATHNAME]	# Image file list
char	specfile[SZ_PATHNAME]	# image file name (possible with extensions)
char	specpath[SZ_PATHNAME]	# image file pathname
char	specdir[SZ_PATHNAME]	# image file directory
char	image[SZ_PATHNAME]		# image file name
char	imname[SZ_PATHNAME]	# image name
char	outname[SZ_PATHNAME]	# Output file name
int	mspec			# Number of spectrum to read from multispec file
int	mband			# Number of band to read from multispec file
int	i, icol, ipix, ixpix
pointer	im			# Image header structure
pointer	sh			# Spectrum header structure
pointer	pix			# pointer to pixel values
pointer	wav			# pointer to pixel wavelengths
int	npts			# Number of pixels in spectrum
real	wavelength		# Wavelength at pixel in angstroms
double	wl			# Wavelength at pixel in angstroms
real	value			# Image value at pixel
real	minwav			# Minimum wavelength in angstroms
real	maxwav			# Maximum wavelength in angstroms
int	minpix, maxpix, npix
int	mspec_range[3,MAX_RANGES]
int	mext			# FITS image extension
int	nmspec, nmspec0
char	lbracket[3]		# "[({"
char	rbracket[3]		# "])}"
int	speclist		# List of spectrum files
char	specnums[SZ_LINE]
char	str[SZ_LINE]
int	world			# wavelength or pixel

real	fpix, fmin, fmax, fnorm, fmean
int	outfd
double	dpix, wl1, wl2
 
bool	debug, verbose, heading, logwav, printlim
int	mode, nc, ns, ncol,ip, jp, lfile
double	fsum
real	renorm
bool	outfile
char	cdot, cslash
char	fluxform[SZ_FORMAT]
char	waveform[SZ_FORMAT]
char	pixform[SZ_FORMAT]
char	numform[SZ_FORMAT]
char	velform[SZ_FORMAT]
char	columns[8]
double	c0

real	clgetr()
int	clpopnu()
int	clgeti(), clscan()
bool	clgetb()
int	open(), imtgetim(), imaccess()
int	stridxs(), stridx()
int	strldx(), strncmp(), strlen()
double	wcs_p2w(), wcs_p2l()
int	decode_ranges(),get_next_number()

define  newspec_ 10
define  newap_   20
define	endxc_   30

begin
	cdot = '.'
	cslash = '/'
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	c0 = 299792.5d0

	# Open input list of images
	speclist = clpopnu ("input")
	call clgstr ("input", images, SZ_PATHNAME)
	if (strlen (images) < 1) {
	    call printf ("LISTSPEC lists a spectrum in an arbitrary ASCII format\n")
	    call printf (" n=ap p=pixel w=wavelength f=flux v=velocity d=delta wavelength\n")
	    call printf (" can be listed in any order with specified formats.\n")
	    return
	    }
	call clgstr ("specdir", specdir, SZ_PATHNAME)

	mext = clgeti ("specext")

# Multispec spectrum numbers (use only first if multiple files)
	call clgstr ("specnum",specnums,SZ_LINE)
	if (decode_ranges (specnums, mspec_range, MAX_RANGES, nmspec0) == ERR){
            call sprintf (str, SZ_LINE, "LISTSPEC: Illegal multispec list <%s>")
                call pargstr (specnums)
            call error (1, str)
            }

	mspec = clgeti ("specnum")
	mband = clgeti ("specband")
	world = 1

	call clgstr ("columns",columns,7)
	ncol = strlen (columns)
	printlim = clgetb ("printlim")

	minwav = clgetr ("lambda1")
	maxwav = clgetr ("lambda2")

	minpix = clgeti ("pix1")
	maxpix = clgeti ("pix2")

	verbose = clgetb ("verbose")
	

	debug = FALSE
	if (clscan ("debug") != EOF)
	    debug = clgetb ("debug")

	heading = FALSE
	if (clscan ("heading") != EOF)
	    heading = clgetb ("heading")

	renorm = 0.0
	if (clscan ("renormalize") != EOF)
	    renorm = clgetr ("renormalize")

	outfile = FALSE
	if (clscan ("outfile") != EOF)
	    outfile = clgetb ("outfile")

	call strcpy ("%g", fluxform, SZ_FORMAT)
	if (clscan ("fluxform") != EOF)
	    call clgstr ("fluxform", fluxform, SZ_FORMAT)

	logwav = FALSE
	if (clscan ("logwav") != EOF)
	    logwav = clgetb ("logwav")

	call strcpy ("%9.3f", waveform, SZ_FORMAT)
	if (clscan ("waveform") != EOF)
	    call clgstr ("waveform", waveform, SZ_FORMAT)

	call strcpy ("%3d", numform, SZ_FORMAT)
	if (clscan ("numform") != EOF)
	    call clgstr ("numform", numform, SZ_FORMAT)

	call strcpy ("%9.3f", velform, SZ_FORMAT)
	if (clscan ("velform") != EOF)
	    call clgstr ("velform", velform, SZ_FORMAT)

	call strcpy ("%3d", pixform, SZ_FORMAT)
	if (clscan ("pixform") != EOF)
	    call clgstr ("pixform", pixform, SZ_FORMAT)

	mode = READ_ONLY

	if (debug) {
	    call printf ("LISTSPEC: parameters all read.\n")
	    call flush (STDOUT)
	    }

# Loop over images

# Get next object spectrum file name from the list
newspec_
        if (imtgetim (speclist, specfile, SZ_PATHNAME) == EOF)
           go to endxc_

	if (debug) {
	    call printf ("LISTSPEC: Spectrum %s to be read\n")
		call pargstr (specfile)
	    call flush (STDOUT)
	    }

# Check for specified apertures in multispec spectrum file
	ip = stridxs (lbracket,specfile)
	if (ip > 0) {
	    lfile = strlen (specfile)
	    specfile[ip] = EOS
	    jp = 0
	    ip = ip + 1
	    while (stridx (specfile[ip],rbracket) == 0 && ip <= lfile) {
		jp = jp + 1
		specnums[jp] = specfile[ip]
		specfile[ip] = EOS
		ip = ip + 1
		}
	    if (jp > 0)
		specnums[jp+1] = EOS
	    else
		call strcpy ("0",specnums,SZ_LINE)
	    if (decode_ranges (specnums,mspec_range,MAX_RANGES,nmspec) == ERR){
		call sprintf (str, SZ_LINE, "T_XCSAO: Illegal multispec list <%s>")
		    call pargstr (specnums)
		call error (1, str)
		}
	    }
	else
	    nmspec = nmspec0
	if (debug) {
	    call printf ("LISTSPEC: next file is %s [%s] = %d aps\n")
		call pargstr (specfile)
		call pargstr (specnums)
		call pargi (nmspec)
	    }
# Check for readability of object spectrum
	call strcpy (specdir,specpath,SZ_PATHNAME)
	call strcat (specfile,specpath,SZ_PATHNAME)
	if (imaccess (specpath, READ_ONLY) == NO) {
	    call eprintf ("XCSAO: cannot read spectrum file %s \n")
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

	# Try to open the image
	call getimage (specpath, mext, mspec, mband, pix, im, sh, npts, imname,
		       mode, world, debug)
	if (im == ERR)
	    go to newspec_

# Open output file
	if (outfile) {
	    nc = strldx (cdot, image)
	    ns = strldx (cslash, image)
	    if (nc > 0) {
		if (strncmp (image[nc+1],"imh",3) == 0) {
		    call strcpy (image[ns+1], outname, nc-ns-1)
		    outname[nc] = EOS
		    }
		else if (strncmp (image[nc+1],"fit",3) == 0) {
		    call strcpy (image[ns+1], outname, nc-ns-1)
		    outname[nc] = EOS
		    }
		else
		    call strcpy (image[ns+1], outname, 32)
		}
	    if (nmspec0 > 1) {
		call sprintf (str, SZ_LINE, "_%d")
		    call pargi (mspec)
		call strcat (str, outname, SZ_PATHNAME)
		}
	    call strcat (".wav", outname, SZ_PATHNAME)
	    outfd = open (outname, APPEND, TEXT_FILE)
	    }
	else
	    outfd = STDOUT

	if (heading) {
	    if (mspec > 0 && mband > 0) {
		call fprintf (outfd, "# %s[%d] %s band %d: %d pixels\n")
		    call pargstr (image)
		    call pargi (mspec)
		    call pargstr (imname)
		    call pargi (mband)
		    call pargi (npts)
		}
	    else if (mspec > 0) {
		call fprintf (outfd, "# %s[%d] %s : %d pixels\n")
		    call pargstr (image)
		    call pargi (mspec)
		    call pargstr (imname)
		    call pargi (npts)
		}
	   else if (mband > 0) {
		call fprintf (outfd, "# %s %s band %d: %d pixels\n")
		    call pargstr (image)
		    call pargi (mspec)
		    call pargstr (imname)
		    call pargi (npts)
		}
	    else {
		call fprintf (outfd, "# %s %s : %d pixels\n")
		    call pargstr (image)
		    call pargstr (imname)
		    call pargi (npts)
		}
	    }

# Get image wavelengths to list
	wav = SX(sh)

# Print image name
	if (verbose) {
	    if (outfile) {
		call printf ("Spectrum %s %d - %d -> %s\n")
		    call pargstr (imname)
		    call pargi (NP1(sh))
		    call pargi (NP2(sh))
		    call pargstr (outname)
		}
	    else {
		call printf ("Spectrum %s %d - %d\n")
		    call pargstr (imname)
		    call pargi (NP1(sh))
		    call pargi (NP2(sh))
		}
	    call flush (STDOUT)
	    }

# Loop over pixels
	if (minpix == INDEFI || minpix < NP1(sh))
	    minpix = NP1(sh)
	if (maxpix == INDEFI || maxpix > NP2(sh))
	    maxpix = NP2(sh)

	fsum = 0.d0
	npix = 0
	do i = minpix, maxpix {
	    fpix = Memr[pix + i - 1]
	    if (i == minpix) {
		fmax = fpix
		fmin = fpix
		}
	    else if (fpix > fmax)
		fmax = fpix
	    else if (fpix < fmin)
		fmin = fpix
	    fsum = fsum + double (fpix)
	    npix = npix + 1
	    }

	fmean = fsum / double (npix)
	if (debug) {
	    call printf ("LISTSPEC: pixels %d-%d: %.8g, %.8g-%.8g\n")
		call pargi (minpix)
		call pargi (maxpix)
		call pargr (fmean)
		call pargr (fmin)
		call pargr (fmax)
	    }

# Compute normalization factors so maximum value is renorm
	if (renorm > 0.0) {
	    fnorm = renorm / (fmax - fmin)
	    if (debug) {
		call printf ("LISTSPEC: fnorm = %.2f\n")
		    call pargr (fnorm)
		}
	    }

# Compute normalization factors so mean value is renorm
	else if (renorm < 0.0) {
	    fnorm = -renorm / fmean
	    fmin = 0.0
	    if (debug) {
		call printf ("LISTSPEC: fnorm = %.2f\n")
		    call pargr (fnorm)
		}
	    }
	else {
	    fnorm = 1.0
	    fmin = 0.0
	    }

	if (debug) {
	    call printf ("Output formats are w=%s, f=%s\n")
		call pargstr (waveform)
		call pargstr (fluxform)
	    call flush (STDOUT)
	    }

	call wcs_set (sh)
	if (printlim)
	    ixpix = maxpix - minpix
	else
	    ixpix = 1
	do ipix = minpix, maxpix, ixpix {
	    dpix = double (ipix)
	    wavelength = Memr[wav+i-1]
	    if (logwav) {
		wl = wcs_p2l (dpix)
		wl1 = wcs_p2l (dpix-0.5d0)
		wl2 = wcs_p2l (dpix+0.5d0)
		}
	    else {
		wl = wcs_p2w (dpix)
		wl1 = wcs_p2w (dpix-0.5d0)
		wl2 = wcs_p2w (dpix+0.5d0)
		}
	    value = Memr[pix+ipix-1]
	    if ((minwav == INDEFR || wavelength >= minwav) &&
		(maxwav == INDEFR || wavelength <= maxwav)) {

		if (renorm != 0.0)
		    value = (value - fmin) * fnorm

	    # Write values to the standard output
		do icol = 1, ncol {
		    if (columns[icol] == 'n') {
			call fprintf (outfd, numform)
			    call pargi (mspec)
			}
		    else if (columns[icol] == 'p') {
			call fprintf (outfd, pixform)
			    call pargi (ipix)
			}
		    else if (columns[icol] == 'w') {
			call fprintf (outfd, waveform)
			    call pargd (wl)
			}
		    else if (columns[icol] == 'f') {
			call fprintf (outfd, fluxform)
			    call pargr (value)
			}
		    else if (columns[icol] == 'd') {
			call fprintf (outfd, waveform)
			if (wl1 > wl2)
			    call pargd (wl1 - wl2)
			else
			    call pargd (wl2 - wl1)
			}
		    else if (columns[icol] == 'v') {
			call fprintf (outfd, velform)
			if (wl1 > wl2)
			    call pargd (((wl1 - wl2) / wl) * c0)
			else
			    call pargd (((wl2 - wl1) / wl) * c0)
			}
		    if (icol < ncol)
			call fprintf (outfd, " ")
		    }
		call fprintf (outfd, "\n")
		}
	    }

# Flush standard output before opening next image
	call flush (outfd)

# Close image
	call close_image (im, sh)

# Close output file
	if (outfile)
	    call close (outfd)

# Move on to next aperture or next image
	nmspec = nmspec - 1
	if (nmspec > 0)
	    go to newap_
	go to newspec_

# Close input list
endxc_
	call imtclose (speclist)
	return
end

# Mar 25 2008	New task based on listwav
# Apr 23 2008	Fix delta wavelength bug

# Apr 15 2009	Add ranges of apertures and write them to separate files
# Apr 15 2009	Add specdir parameter
# Oct  7 2009	Make images length SZ_PATHNAME instead of SZ_FNAME

# Jun 10 2011	Add velform format for velocity per pixel output
