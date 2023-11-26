#*** File rvsao/Util/getimage.x
#*** August 11, 2010
#*** By Doug Mink

#  GETIMAGE -- Opens IRAF images for gettemp and getspec and returns the
#  image file descriptors and extracted header information.
#  Adapted from the onedspec package in IRAF.
 
include	<imhdr.h>
include	<smw.h>
 
define	LEN_USER_AREA	100000
 
procedure getimage (image,mext,mspec,band,pix,im,sh,npix,name,mode,world,debug)
 
char	image[ARB]	# Spectrum image file name
int	mext		# FITS image extension
int	mspec		# Number of spectrum to read from multispec or 2-D file
int	band		# Multispec band
pointer	pix		# Spectrum [returned]
pointer	im		# Image header structure [returned]
pointer	sh		# Spectrum header structure [returned]
int	npix		# Number of pixels in spectrum [returned]
char	name[ARB]	# Name of object [returned]
int	mode		# Mode (READ_ONLY, READ_WRITE)
int	world		# wavelength or pixel
bool	debug
 
double	ra, dec, w0, wpc
pointer mw, smw_openim()
int	nline, ap, ndim
int	nband		# Band number for 3-D images
int	ip, nc, ncat, nfib, i
int	npix2, npix3
int	ltitle
int	ibeam
char	temp[32]
# char	units[32]
char	apstring[80]
char	impath[SZ_PATHNAME]
char	imname[SZ_PATHNAME]
char	tname[SZ_PATHNAME]
char	rastr[32]
char	decstr[32]

int	imaccf(), imaccess()
int	ctoi(), ctowrd()
int	strlen(), strldxs(), strncmp()
pointer	immap()
errchk	immap()
 
begin
	im = NULL
	pix = NULL
	sh = NULL
	call strcpy (image, impath, SZ_PATHNAME)
#	if (mext < 1) {
#	    call strcpy (image, imname, SZ_PATHNAME)
#		call pargstr (image)
#	    }
#	else {
	    call sprintf (imname, SZ_PATHNAME, "%s[%d]")
		call pargstr (image)
		call pargi (mext)
#	    }
	if (debug) {
	    call printf ("GETIMAGE: image file %s about to be read\n")
		call pargstr (imname)
	    call flush (STDOUT)
	    }

	if (imaccess (impath, mode) == NO) {
	    if (imaccess (impath, READ_ONLY) == NO) {
		if (debug) {
		    call printf ("GETIMAGE: image file %s cannot be read\n")
			call pargstr (imname)
		    call flush (STDOUT)
		    }
		im = ERR
		return
		}
	    else {
		call eprintf ("GETIMAGE: Cannot write to image %s\n")
		    call pargstr (image)
		}
	    }

#  Map the image
	iferr (im = immap (imname, mode, LEN_USER_AREA)) {
	    if (im != NULL) call imunmap (im)
	    if (debug) {
		call printf ("GETIMAGE: image file %s cannot be opened\n")
		    call pargstr (imname)
		call flush (STDOUT)
		}
	    if (mext == 0) {
		iferr (im = immap (image, mode, LEN_USER_AREA)) {
		    if (mode == READ_ONLY)
			call eprintf ("GETIMAGE: Cannot read image %s\n")
		    else
			call eprintf ("GETIMAGE: Cannot write to image %s\n")
			call pargstr (image)
		    im = ERR
		    return
		    }
		}
	    else {
		if (mode == READ_ONLY)
		    call eprintf ("GETIMAGE: Cannot read image %s\n")
		else
		    call eprintf ("GETIMAGE: Cannot write to image %s\n")
		    call pargstr (image)
		im = ERR
		return
		}
	    }
	ndim = IM_NDIM(im)
	if (debug) {
	    call printf ("GETIMAGE: %d x %d x %d %d-D image: spectrum %d\n")
		call pargi (IM_LEN(im,1))
		call pargi (IM_LEN(im,2))
		call pargi (IM_LEN(im,3))
		call pargi (ndim)
		call pargi (mspec)
	    call flush (STDOUT)
	    }

	if (ndim > 1) {
	    npix2 = IM_LEN(im,2)
	    npix3 = IM_LEN(im,3)
	    if (ndim > 2 && npix3 > 1) {
		if (band > npix3) {
		    call eprintf ("GETIMAGE: Band %d > %d in image %s\n")
			call pargi (band)
			call pargi (npix3)
			call pargstr (image)
		    call imunmap (im)
		    im = ERR
		    return
		    }
		else if (mspec > npix2) {
		    call eprintf ("GETIMAGE: Ap %d > %d(3) in image %s\n")
			call pargi (mspec)
			call pargi (npix2)
			call pargstr (image)
		    call imunmap (im)
		    im = ERR
		    return
		    }
		else {
		    nline = mspec
		    nband = band
		    }
		}
	    else if (ndim == 2 && mspec > npix2) {
		call eprintf ("GETIMAGE: Ap %d > %d(2) in image %s\n")
		    call pargi (mspec)
		    call pargi (npix2)
		    call pargstr (image)
		    call imunmap (im)
		    im = ERR
		    return
		    }
	    else {
		nband = 1
		if (mspec > 0)
		    nline = mspec
		else
		    nline = 1
		}
	    }
	else {
	    nline = 1
	    nband = 1

# Set up FITS World Coordinate System keywords if not present in 1-D spectrum
	    if (imaccf (im, "W0") == YES && imaccf (im, "CRVAL1") == NO) {
		call imgdpar (im, "W0",w0)
		call imaddi (im,"CRPIX1",1)
		call imaddd (im,"CRVAL1",w0)
		}
	    if (imaccf (im, "WPC") == YES && imaccf (im, "CDELT1") == NO) {
        	call imgdpar (im, "WPC",wpc)
		call imaddd (im,"CDELT1",wpc)
		}
	    }

# If stacked spectra and CD2_2 is not set, set it to 1
	if (ndim > 1) {
	    if (imaccf (im, "CD2_2") == NO) {
		call imaddd (im, "CD2_2", 1.d0)
		}
	    }

# If working in pixel space, set WCS accordingly
#	if (world == 2) {
#	    call strcpy ("angstroms", units, 32)
#	    call imgspar (im, "CRUNIT1", units, 32)
#	    if (strncmp (units, "pixel", 5) != 0) {
#		call imaddi (im, "DC-FLAG", 0)
#		call imaddi (im, "CRPIX1", 1)
#		call imaddd (im, "CRVAL1", 1.d0)
#		call imaddd (im, "CDELT1", 1.d0)
#		call imaddd (im, "CD1_1", 1.d0)
#		if (ndim > 1) {
#		    call imaddi (im, "CRPIX2", 1)
#		    call imaddd (im, "CRVAL2", 1.d0)
#		    call imaddd (im, "CDELT2", 1.d0)
#		    call imaddd (im, "CD2_2", 1.d0)
#		    }
#		}
#	    }

        ap = INDEFI

#  Set dispersion axis for 1-D spectra
	if (ndim < 3 && imaccf (im, "DISPAXIS") == NO) {
	    if (ndim < 2) {
		call imaddi (im,"DISPAXIS",1)
		if (debug) {
		    call eprintf ("GETIMAGE: Dispersion axis set to 1 in image %s\n")
			call pargstr (image)
		    }
		}
	    else {
		if (IM_LEN(im,1) == 1) {
		    call imaddi (im,"DISPAXIS",2)
		    if (debug) {
			call eprintf ("GETIMAGE: Dispersion axis set to 2 in image %s\n")
			    call pargstr (image)
			}
		    }
		else if (IM_LEN(im,2) == 1) {
		    call imaddi (im,"DISPAXIS",1)
		    if (debug) {
			call eprintf ("GETIMAGE: Dispersion axis set to 1 in image %s\n")
			    call pargstr (image)
			}
		    }
		else {
		    call imaddi (im,"DISPAXIS",1)
		    if (debug) {
			call eprintf ("GETIMAGE: No dispersion axis in %d-dim image %s, 1 assumed\n")
			    call pargi (ndim)
			    call pargstr (image)
			}
		    }
		}
	    }

#  Open spectrum world coordinate system
	if (debug) {
	    call printf ("GETIMAGE: Ready to open MWCS\n")
	    call flush (STDOUT)
	    }
	mw = smw_openim (im)
	call smw_daxis (mw, im, 0, INDEFI, INDEFI)
	call smw_saxes (mw, NULL, im)

	if (mw == ERR) {
	    call eprintf ("GETIMAGE: MWCS error in image %s\n")
		call pargstr (image)
	    call imunmap (im)
	    im = ERR
	    return
	    }

#  Set position variables to 0:00 if missing
	if (debug) {
	    call printf ("GETIMAGE: Checking positions in header\n")
	    call flush (STDOUT)
	    }
	temp[0] = EOS
	call imgspar (im, "RA", temp, 32)
	if (strlen (temp) < 1) {
	    call imastr (im, "RA","00:00:00.000")
	    if (debug) {
		call imgspar (im, "RA", temp, 32)
		call printf ("GETIMAGE: RA set to %s\n")
		    call pargstr (temp)
		call flush (STDOUT)
		}
	    }
	temp[0] = EOS
	call imgspar (im, "DEC", temp, 32)
	if (strlen (temp) < 1) {
	    call imastr (im, "DEC","00:00:00.00")
	    if (debug) {
		call imgspar (im, "DEC", temp, 32)
		call printf ("GETIMAGE: DEC set to %s\n")
		    call pargstr (temp)
		call flush (STDOUT)
		}
	    }
	temp[0] = EOS
	call imgspar (im, "HA", temp, 32)
	if (strlen (temp) < 1) {
	    call imastr (im, "HA","00:00:00.000")
	    if (debug) {
		call imgspar (im, "HA", temp, 32)
		call printf ("GETIMAGE: HA set to %s\n")
		    call pargstr (temp)
		call flush (STDOUT)
		}
	    }
	temp[0] = EOS
	call imgspar (im, "ST", temp, 32)
	if (strlen (temp) < 1) {
	    call imastr (im, "ST","00:00:00.000")
	    if (debug) {
		call imgspar (im, "ST", temp, 32)
		call printf ("GETIMAGE: ST set to %s\n")
		    call pargstr (temp)
		call flush (STDOUT)
		}
	    }

#  Open spectrum header
	if (debug) {
	    call printf ("GETIMAGE: Ready to open SH\n")
	    call flush (STDOUT)
	    }
	call shdr_open (im, mw, nline, nband, ap, SHDATA, sh)
	if (sh == ERR) {
	    call eprintf ("GETIMAGE: Spectrum header error in image %s\n")
		call pargstr (image)
	    call imunmap (im)
	    im = ERR
	    return
	    }

	if (world < 2 && (strncmp (LABEL(sh),"Pix",3) == 0 ||
	    strncmp (LABEL(sh),"pix",3) == 0)) {
	    if (debug) {
	    call eprintf ("GETIMAGE: Non-spectral dispersion %s in image %s\n")
		call pargstr (LABEL(sh))
		call pargstr (image)
		}
#	    call imunmap (im)
#	    im = ERR
#	    return
	    }

	npix = IM_LEN (im, 1)

#  Get image line to plot
	if (debug) call printf ("GETIMAGE: Ready to get image line\n")
	pix = SY(sh)

	w0 = W0(sh)
	wpc = WP(sh)
	if (w0 <= 1.d0) {
	    call imgipar (im,"DC-FLAG ",DC(sh))
	    call imgdpar (im, "W0",w0)
            call imgdpar (im, "WPC",wpc)
	    if (DC(sh) == 1) {
		W0(sh) = 10.d0 ** w0
		WP(sh) = 10.d0 ** wpc
		W1(sh) = 10.d0 ** (w0 + (wpc * (npix - 1)))
		}
	    else {
		W0(sh) = w0
		WP(sh) = wpc
		W1(sh) = w0 + (wpc * (npix - 1))
		}
	    CTWL(sh) = NULL
	    }

	if (debug) {
	    call printf ("GETIMAGE: line = %d, aperture = %d %d\n")
		call pargi (nline)
		call pargi (AP(sh))
		call pargi (ap)
	    }

#  Title is object name (and spectrum number if multispec)
	call strcpy (TITLE(sh), name, SZ_PATHNAME)
	ltitle = strlen (name)
	if (ltitle < 1) {
	    call strcpy (IM_TITLE(im), name, SZ_PATHNAME)
	    ltitle = strlen (name)
	    }
	if (ltitle < 1) {
	    ip = strldxs ("$/",image)
	    if (ip > 0)
		call strcpy (image[ip+1],name,SZ_PATHNAME)
	    else
		call strcpy (image, name, SZ_PATHNAME)
	    }
	if (mspec > 0) {
	    call sprintf (temp, 32, "_%d")
		call pargi (mspec)
	    call strcat (temp, name, SZ_PATHNAME)
	    }

# If APIDi keyword is present, read position information from it
	call sprintf (temp, 32, "APID%d")
	    call pargi (mspec)
	if (imaccf (im, temp) == YES) {
	    call imgspar (im, temp, apstring, 80)
	    do i = 1, 80 {
		if (apstring[i] != ' ')
		    tname[i] = apstring[i]
		else {
		    tname[i] = EOS
		    ip = i + 1
		    break
		    }
		}
	    nc = ctowrd (apstring, ip, rastr, 31)
	    if (nc > 0)
		call imastr (im, "RA",rastr)
	    nc = ctowrd (apstring, ip, decstr, 31)
	    if (nc > 0)
		call imastr (im, "DEC",decstr)
	    nc = ctoi (apstring, ip, ncat)
	    ncat = -10
	    if (nc > 0)
		call imaddi (im, "CATNUM", ncat)
	    nc = ctoi (apstring, ip, nfib)
	    if (nc > 0)
		call imaddi (im, "FIBER", nfib)
	    ibeam = -10
	    if (ncat == 0) {
		ibeam = 2
		call strcpy (tname, name, SZ_PATHNAME)
		}
	    else if (ncat < 0) {
		ibeam = 0
		call strcpy (tname, name, SZ_PATHNAME)
		}
	    else {
		ibeam = 1
		if (debug) {
		    call printf ("GETIMAGE: field name is \"%s\", object name is \"%s\"\n")
			call pargstr (IM_TITLE(im))
			call pargstr (tname)
		    }
		if (strncmp (tname, "target", 6) == 0) {
		    call sprintf (name, SZ_PATHNAME, "%s_%d")
			call pargstr (IM_TITLE(im))
			call pargi (ncat)
		    }
		else
		    call strcpy (tname, name, SZ_PATHNAME)
		}
	    if (debug) {
		    call printf ("GETIMAGE: ap %d: fiber=%d beam= %d ra=%h dec=%h name=%s\n")
		    call pargi (mspec)
		    call pargi (nfib)
		    call pargi (ibeam)
		    call pargd (ra)
		    call pargd (dec)
		    call pargstr (name)
		}
	    if (ncat > -10)
		call imaddi (im, "BEAM", ibeam)
	    }

#  Use telescope pointing direction RA and DEC if RA keyword not there
	if (imaccf (im, "RA") == NO) {
	    if (imaccf (im, "POSTN-RA") == YES) {
		call imgdpar (im, "POSTN-RA",ra)
		ra = ra / 15.d0
		call imaddd (im, "RA",ra)
		call imgdpar (im, "POSTN-DEC",dec)
		call imaddd (im, "DEC",dec)
		}
	    }

end
# Feb 14 1992	Rewrite to read any multispec file, not just echelle
# 		Read velocities from strings in multispec files
# Feb 18 1992	Read RA and Dec from telescope direction in multifiber files
# Apr 22 1992	Get mspec as input argument
# May 22 1992	Do not change mspec
# Aug 11 1992	Set number of emission lines fit correctly

# Jan 22 1993	Check APNUM* only if mspec > 0
# May 25 1993	Move all velocity information out of subroutine
# May 28 1993	Set up MWCS pixel<->wavelength transformations
# Jun  4 1993	Use shdr subroutine from onedspec package to open mwcs
# Jun 16 1993	Fix shdr subroutine use
# Jul  7 1993	Add spectrum header to getimage arguments
# Jul 12 1993	Turn off debugging
# Aug 11 1993	Reset mspec if it exceeds limit found in shdr_open
# Aug 24 1993	Fix error checking
# Dec  2 1993	Increase user area length from 2880 to 100000
# Dec  3 1993	Add spectrum number to spectrum name

# Mar 23 1994	Use filename as object name if no object name is present
# Apr  8 1994	Add FITS WCS info if not present in 1-D spectrum
# Apr 12 1994	Handle 3-D multispec files correctly
# Apr 12 1994	Return mwcs structure
# Apr 13 1994	Drop unused variable naxes; fix shdr close argument list
# May  3 1994	Do not free pix 
# May 23 1994	Improve image open error message
# Jun 15 1994	Add flag for xc filtering; drop log lambda flag from common
# Jun 23 1994	Keep mwcs pointer in shdr

# Jan 19 1995	Return error if spectrum header error occurs
# Mar 29 1995	Return error if dispaxis not set in multispec image
# May 15 1995	Change all sz_fname to sz_pathname, which is longer
# Jun  6 1995	Set dispersion axis for 1 or 2 dimensions
# Jun 14 1995	Assume dispersion axis is 1 if none is set
# Jun 19 1995	Print aperture number if different that multispec sequence
# Jun 19 1995	Get spectrum name from spectrum header OR image title
# Jun 19 1995	Add aperture number to name in SHDR_OPEN, not here
# Oct  6 1995	Call SPHD_OPEN instead of SHDR_OPEN

# Jul  8 1996	Pass debugging flag in to this subroutine
# Aug  7 1996	Use smw.h; add smw_openim
# Aug 15 1996	Add smw_daxis to set summing to one

# Apr  9 1997	Return error if aperture not 0 or within image dimension
# Apr  9 1997	Return error if dispersion units are pixels
# Aug 27 1997	Add argument for multispec band

# Apr 22 1998	Drop getim.com

# Aug  2 2000	Add spectrum number to name of spectrum
# Sep 26 2000	Add world argument for wavelength, pixel

# Sep  7 2001	Fix units comparison
# Sep 19 2001	Test for pixel units correctly

# Aug  1 2003	Add argument for image extension number
# Oct 20 2003	Pass spectrum even if no dispersion 

# Jan 16 2004	Print nonspectral dispersion message only in debug mode

# Dec  8 2005	Do not add aperture/band to image name if 0

# Jun  1 2006	Read position, name, and beam from APID if present

# Apr  6 2007 If no name for spectrum, make file_ap, not file [ap]
# Jun 13 2007	If APIDi doesn't contain complete info, do not overwrite header
# Jun 14 2007	Read RA and Dec as strings from APIDi

# Mar 10 2008	Ignore world argument; resetting header confuses shdr_open()
# Jul  2 2008	Check for file existence using imaccess(), not access()

# Apr 10 2009	If no pointing direction is present, set it to 0:00 0:00
# Oct  7 2009	Make all file names length SZ_PATHNAME instead of SZ_FNAME

# Jul 16 2010	If stacked spectra and CD2_2 is not set, set it to 1
# Aug 11 2010	Always add extension number to name
