# File Util/t_velset.x
# Modified by Doug Mink from stsdas.playpen.newredshift
# October 1, 1997

include	<imhdr.h>
include	<error.h>

# VELSET -- Changes the redshift of spectra.
#
# The input spectra are given by an image template list. The output is either 
# a matching list of spectra or a directory. The number of input spectra may 
# be either one or match the number of output spectra. Image sections are 
# ignored, since the user wants a exact copy of the input, however with the
# wavelength scale modified.
#
#							Ivo Busko  7/17/89

procedure t_velset()

char	imtlist1[SZ_LINE]		# Input spectra list
char	imtlist2[SZ_LINE]		# Output spect. list/directory
double	newvel				# New redshift
bool	verbose				# Print operations?
bool	velz				# New redshift in Z (yes) or km/sec

char	image1[SZ_PATHNAME]		# Input image name
char	image2[SZ_PATHNAME]		# Output image name
char	dirname1[SZ_PATHNAME]		# Directory name
char	dirname2[SZ_PATHNAME]		# Directory name

int	list1, list2, root_len
double	c0

int	imtopen(), imtgetim(), imtlen()
int	fnldir(), isdirectory()
bool	clgetb()
double	clgetd()

begin
	c0 = 299792.5d0

	# Get task parameters
	call clgstr ("input", imtlist1, SZ_LINE)
	velz	= clgetb ("velz")
	newvel  = clgetd ("newvel")
	if (velz)
	    newvel = c0 * (newvel - 1.d0)
	call clgstr ("output", imtlist2, SZ_LINE)
	verbose = clgetb ("verbose")

	# If the output string is a directory, generate names for
	# the new images accordingly.

	if (isdirectory (imtlist2, dirname2, SZ_PATHNAME) > 0) {
	    list1 = imtopen (imtlist1)
	    while (imtgetim (list1, image1, SZ_PATHNAME) != EOF) {

		# Strip an eventual image section first.  Place the input 
		# image name, without a directory or image section, in string 
		# dirname1.
		call imgimage (image1, image2, SZ_PATHNAME)
		root_len = fnldir (image2, dirname1, SZ_PATHNAME)
		call strcpy (image2[root_len + 1], dirname1, SZ_PATHNAME)

		# Assemble output image name. Strip again image section from
		# input image name.
		call strcpy (dirname2, image2, SZ_PATHNAME)
		call strcat (dirname1, image2, SZ_PATHNAME)
		call imgimage (image1, image1, SZ_PATHNAME)
		call printf ("VELSET: %s -> %s\n")
		    call pargstr (image1)
		    call pargstr (image2)
		call flush (STDOUT)

		# Do it.
		call vels_image (image1, image2, velz, newvel, verbose)
	    }
	    call imtclose (list1)

	} else {

	    # Expand the input and output image lists.
	    list1 = imtopen (imtlist1)
	    list2 = imtopen (imtlist2)
	    if (imtlen (list1) != imtlen (list2)) {
	        call imtclose (list1)
	        call imtclose (list2)
	        call error (0, "Number of input and output images not the same")
	    }

	    # Do each set of input/output images. First strip any sections.
	    while ((imtgetim (list1, image1, SZ_PATHNAME) != EOF) &&
		(imtgetim (list2, image2, SZ_PATHNAME) != EOF)) {
		call imgimage (image1, image1, SZ_PATHNAME)
		call imgimage (image2, image2, SZ_PATHNAME)
		call printf ("VELSET: %s -> %s\n")
		    call pargstr (image1)
		    call pargstr (image2)
		call flush (STDOUT)
		call vels_image (image1, image2, velz, newvel, verbose)
	    }

	    call imtclose (list1)
	    call imtclose (list2)
	}
end

# VELS_IMAGE -- Copy a spectrum, changing the wavelength scale accordingly
# to the value of newz relative to the VELOCITY parameter. 
#
# If the input and output image names are equal, just 
# issue a warning message. Keywords that describe the wavelength axis are 
# looked for in the header, in the following order: first, ONEDSPEC keywords 
# W0 and WPC; second, CD keywords CRVALn and CDn_1, and at last FITS keywords 
# CRVALn and CDELTn. Here n is the value of parameter AXIS, unless the ONEDSPEC 
# header keyword DISPAXIS is found, in which case it takes precedence over the 
# parameter. Logarithmic wavelength scale is treated correctly either if the 
# LOG task parameter is set to yes or if the ONEDSPEC header keyword DC-FLAG is
# found in the header with value 1. HISTORY records are appended to the 
# header.

procedure vels_image (image1, image2, velz, newvel, verbose)

char	image1[ARB]	# Input spectrum
char	image2[ARB]	# Output spectrum
bool	velz		# Velocity in Z (yes) or km/sec (no)
double	newvel		# New redshift in km/sec
bool	verbose		# Print the operation

pointer	im, sp, str
int	dispaxis, dcflag
double	corre, crval2, cdelt2, crval1, cdelt1, c0, oldvel, oldz, newz
double	hcv
double	imgetd()

pointer immap()
int	imgeti()
int	imaccf()
bool	streq()

begin
	if (streq (image1, image2))
	    call error (0, "Same input and output images.")

	# Check keywords on input image header
	im = immap (image1, READ_ONLY, 0)

	# First check which one is the dispersion axis.
	if (imaccf (im, "DISPAXIS") == YES)
	    dispaxis = imgeti (im, "DISPAXIS")
	else {
	    call printf ("VELSIMAGE: No DISPAXIS for %s\n")
		call pargstr (image1)
            call imunmap (im)
	    return
	    }

	# Next check image dimensionality
	if (dispaxis < 1 || dispaxis > IM_NDIM(im)) {
	    call eprintf ("VELSIMAGE: non-existent axis %d for %s.\n")
		call pargi (dispaxis)
	        call pargstr (image1)
	    call imunmap (im)
	    return
	    }

	# Check if wavelength scale is linear or logarithmic.
	# If no indication is found on the image, assume linear.
	if (imaccf (im, "DC-FLAG") == YES)
	    dcflag = imgeti (im, "DC-FLAG")
	else
	    dcflag = 0
	if (dcflag == -1) {
	    call eprintf ("%s : not calibrated to wavelength.\n")
	        call pargstr (image1)
	    call imunmap (im)
	    return
	    }
	else if (dcflag == 0) {
	    call eprintf ("%s : not log wavelength.\n")
	        call pargstr (image1)
	    call imunmap (im)
	    return
	    }

	# Convert velocity from km/sec to Z = 1 + v/C
	c0 = 299792.5d0
	if (imaccf (im, "VELOCITY") == YES) {
	    oldvel = imgetd (im, "VELOCITY")
	    if (imaccf (im,"BCV") == YES)
		hcv = imgetd (im, "BCV")
	    else if (imaccf (im,"HCV") == YES)
		hcv = imgetd (im, "HCV")
	    else
		hcv = 0.d0
	    oldz = (oldvel - hcv) / c0
	    if (verbose) {
		call printf ("VELSET: Old velocity is %f + %f = %f/n")
		    call pargd (oldvel)
		    call pargd (hcv)
		    call pargd (oldz)
		}
	    }

	call imunmap (im)

	# Create and open the new image.
#	call imcopy (image1, image2)
	call vels_imcopy (image1, image2)
	im = immap (image2, READ_WRITE, 28800)

	if (imaccf (im,"VELOCITY") != YES) {
	    call imaddd (im, "VELOCITY", newvel)
	    call imunmap (im)
	    return
	    }

	# Update keywords with new values corrected by the redshift
	# difference between oldz and newz.
	
	newz = newvel / c0
	corre = (1.d0 + newz) / (1.d0 + oldz)

	call vels_get (im, dispaxis, crval1, cdelt1)

	if (dcflag == 1) {
	    crval2 = crval1 + log10 (corre)
	    cdelt2 = cdelt1
	    }
	else {
	    crval2 = crval1 * corre
	    cdelt2 = cdelt1 * corre
	    }

	call vels_put (im, dispaxis, crval2, cdelt2)

	# If verbose print the operation.
	if (verbose) {
	    call eprintf ("%s %32t->%40t %s\n")
	        call pargstr (image1)
	        call pargstr (image2)
	    call eprintf (" origin = %.8g %40t  origin = %.8g\n")
	        call pargd (crval1)
	        call pargd (crval2)
	    call eprintf (" step   = %.8g %40t  step   = %.8g\n")
	        call pargd (cdelt1)
	        call pargd (cdelt2)
	    }

	# Update velocity in image header
	call imputd (im, "VELOCITY", newvel)
	if (imaccf (im,"BCV") == YES)
	    call imputd (im, "BCV", 0.d0)
	if (imaccf (im,"HCV") == YES)
	    call imputd (im, "HCV", 0.d0)

	# Update image HISTORY.
	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)
	if (velz) {
	    call sprintf (Memc[str], SZ_LINE, " VELSET: oldz = %.8f, newz = %.8f ")
		call pargd (oldz)
		call pargd (newz)
	    }
	else {
	    call sprintf (Memc[str], SZ_LINE, " VELSET: old vel = %.3f, new vel = %.3f ")
		call pargd (oldvel)
		call pargd (newvel)
	    }
	call vels_timelog (Memc[str], SZ_LINE)
	call imputh (im, "HISTORY", Memc[str])
	call sprintf (Memc[str], SZ_LINE, "input file to VELSET: %s ")
	    call pargstr (image1)
	call imputh (im, "HISTORY", Memc[str])
	call sfree (sp)

	call imunmap (im)
end

# VELS_GET -- Get header parameters which describe the wavelength axis.

procedure vels_get (im, axis, crval, cdelt)

pointer	im			#i: image pointer
int	axis			#i: dispersion axis
double	crval			#o: reference wavelength
double	cdelt			#o: wavelength step

char	keyword[SZ_FNAME]
double	imgetd()
int	imaccf()

begin
	# FITS or IRAF keyword
	# Wavelength scale zero is common for ST and FITS formats.
	call sprintf (keyword, SZ_FNAME, "CRVAL%1d")
	    call pargi (axis)
	if (imaccf (im, keyword) == YES)
	    crval = imgetd (im, keyword)

	# STSDAS keyword
	call sprintf (keyword, SZ_FNAME, "CD%1d_1")
	    call pargi (axis)
	if (imaccf (im, keyword) == YES)
	    cdelt = imgetd (im, keyword)

	# FITS keyword
	call sprintf (keyword, SZ_FNAME, "CDELT%1d")
	    call pargi (axis)
	if (imaccf (im, keyword) == YES)
	    cdelt = imgetd (im, keyword)

	# Old IRAF keywords
	if (imaccf (im, "W0") == YES) {
	    crval = imgetd (im, "W0")
	    cdelt = imgetd (im, "WPC")
	    }
end

# VELS_PUT -- Put header parameters which describe the wavelength axis.

procedure vels_put (im, axis, crval, cdelt)

pointer	im			#i: image pointer
int	axis			#i: dispersion axis
double	crval			#i: reference wavelength
double	cdelt			#i: wavelength step

int	imaccf()
char	keyword[SZ_FNAME]

begin
	# FITS or IRAF keyword
	# Wavelength scale zero is common for ST and FITS formats.
	call sprintf (keyword, SZ_FNAME, "CRVAL%1d")
	    call pargi (axis)
	if (imaccf (im, keyword) == YES)
	    call imputd (im, keyword, crval)

	# STSDAS keyword
	call sprintf (keyword, SZ_FNAME, "CD%1d_1")
	    call pargi (axis)
	if (imaccf (im, keyword) == YES)
	    call imputd (im, keyword, cdelt)

	# FITS keyword
	call sprintf (keyword, SZ_FNAME, "CDELT%1d")
	    call pargi (axis)
	if (imaccf (im, keyword) == YES)
	    call imputd (im, keyword, cdelt)

	# Old IRAF keywords
	if (imaccf (im, "W0") == YES) {
	    call imputd (im, "W0", crval)
	    call imputd (im, "WPC", cdelt)
	    }
end

# VELS_TIMELOG -- Prepend a time stamp to the given string.
#
# For the purpose of a history logging prepend a short time stamp to the
# given string.  Note that the input string is modified.

procedure vels_timelog (str, max_char)

char	str[ARB]		# String to be time stamped
int	max_char		# Maximum characters in string

pointer	sp, time, temp

begin
	call smark (sp)
	call salloc (time, SZ_LINE, TY_CHAR)
	call salloc (temp, max_char, TY_CHAR)

	call logtime (Memc[time], SZ_LINE)
	call sprintf (Memc[temp], max_char, "%s %s")
	    call pargstr (Memc[time])
	    call pargstr (str)
	call strcpy (Memc[temp], str, max_char)
	call sfree (sp)
end


# VELS_IMCOPY -- Copy an image.  Use sequential routines to permit copying
# images of any dimension.  Perform pixel i/o in the datatype of the image,
# to avoid unnecessary type conversion.
#
# This routine is basicaly task images.imcopy whith verbose option and
# output image section handling removed.

procedure vels_imcopy (image1, image2)

char	image1[ARB]			# Input image
char	image2[ARB]			# Output image

int	npix, junk
pointer	buf1, buf2, im1, im2
pointer	sp, imtemp
long	v1[IM_MAXDIM], v2[IM_MAXDIM]

int	imgnls(), imgnll(), imgnlr(), imgnld(), imgnlx()
int	impnls(), impnll(), impnlr(), impnld(), impnlx()
pointer	immap()

begin
	call smark (sp)
	call salloc (imtemp, SZ_PATHNAME, TY_CHAR)

	# Map the input image.
	im1 = immap (image1, READ_ONLY, 0)

	# Get a temporary output image name and map it as a copy of the 
	# input image.
	# Copy the input image to the temporary output image and unmap
	# the images.  Release the temporary image name.

	call xt_mkimtemp (image1, image2, Memc[imtemp], SZ_PATHNAME)
	im2 = immap (image2, NEW_COPY, im1)

	# Setup start vector for sequential reads and writes.

	call amovkl (long(1), v1, IM_MAXDIM)
	call amovkl (long(1), v2, IM_MAXDIM)

	# Copy the image.

	npix = IM_LEN(im1, 1)
	switch (IM_PIXTYPE(im1)) {
	case TY_SHORT:
	    while (imgnls (im1, buf1, v1) != EOF) {
		junk = impnls (im2, buf2, v2)
		call amovs (Mems[buf1], Mems[buf2], npix)
	    }
	case TY_USHORT, TY_INT, TY_LONG:
	    while (imgnll (im1, buf1, v1) != EOF) {
		junk = impnll (im2, buf2, v2)
		call amovl (Meml[buf1], Meml[buf2], npix)
	    }
	case TY_REAL:
	    while (imgnlr (im1, buf1, v1) != EOF) {
		junk = impnlr (im2, buf2, v2)
		call amovr (Memr[buf1], Memr[buf2], npix)
	    }
	case TY_DOUBLE:
	    while (imgnld (im1, buf1, v1) != EOF) {
		junk = impnld (im2, buf2, v2)
		call amovd (Memd[buf1], Memd[buf2], npix)
	    }
	case TY_COMPLEX:
	    while (imgnlx (im1, buf1, v1) != EOF) {
	        junk = impnlx (im2, buf2, v2)
		call amovx (Memx[buf1], Memx[buf2], npix)
	    }
	default:
	    call error (1, "unknown pixel datatype")
	}

	# Unmap the images.

	call imunmap (im2)
	call imunmap (im1)
	call xt_delimtemp (image2, Memc[imtemp])
	call sfree (sp)
end


# VELS_CONVERT -- Copy an image, translating in wavelength.  Use sequential
# routines to permit copying images of any dimension.  Perform pixel i/o in
# the datatype of the image, to avoid unnecessary type conversion.

# This routine is basicaly task images.imcopy with verbose option and
# output image section handling removed.

procedure vels_convert (image1, image2, inwl1, indwl, outwl1, outdwl)

char	image1[ARB]	# Input image
char	image2[ARB]	# Output image
double	inwl1		# Input starting wavelength
double	indwl		# Input delta wavelength
double	outwl1		# Output starting wavelength
double	outdwl		# Output delta wavelength

int	npix, junk
pointer	buf1, buf2, im1, im2
pointer	sp, imtemp, tbuf
long	v1[IM_MAXDIM], v2[IM_MAXDIM]

int	imgnls(), imgnll(), imgnlr(), imgnld(), imgnlx()
int	impnls(), impnll(), impnlr(), impnld(), impnlx()
pointer	immap()

begin
	call smark (sp)
	call salloc (imtemp, SZ_PATHNAME, TY_CHAR)

	# Map the input image.
	im1 = immap (image1, READ_ONLY, 0)

	# Get a temporary output image name and map it as a copy of the 
	# input image.
	# Copy the input image to the temporary output image and unmap
	# the images.  Release the temporary image name.

	call xt_mkimtemp (image1, image2, Memc[imtemp], SZ_PATHNAME)
	im2 = immap (image2, NEW_COPY, im1)

	# Setup start vector for sequential reads and writes.

	call amovkl (long(1), v1, IM_MAXDIM)
	call amovkl (long(1), v2, IM_MAXDIM)

	# Copy the image.

	npix = IM_LEN(im1, 1)
	switch (IM_PIXTYPE(im1)) {
	case TY_SHORT:
	    call salloc (tbuf, TY_SHORT, npix)
	    while (imgnls (im1, buf1, v1) != EOF) {
		junk = impnls (im2, buf2, v2)
		call amovs (Mems[buf1], Mems[buf2], npix)
	    }
	case TY_USHORT, TY_INT, TY_LONG:
	    call salloc (tbuf, TY_LONG, npix)
	    while (imgnll (im1, buf1, v1) != EOF) {
		junk = impnll (im2, buf2, v2)
		call amovl (Meml[buf1], Meml[buf2], npix)
	    }
	case TY_REAL:
	    call salloc (tbuf, TY_REAL, npix)
	    while (imgnlr (im1, buf1, v1) != EOF) {
		junk = impnlr (im2, buf2, v2)
		call amovr (Memr[buf1], Memr[buf2], npix)
	    }
	case TY_DOUBLE:
	    call salloc (tbuf, TY_DOUBLE, npix)
	    while (imgnld (im1, buf1, v1) != EOF) {
		junk = impnld (im2, buf2, v2)
		call amovd (Memd[buf1], Memd[buf2], npix)
	    }
	case TY_COMPLEX:
	    while (imgnlx (im1, buf1, v1) != EOF) {
	        junk = impnlx (im2, buf2, v2)
		call amovx (Memx[buf1], Memx[buf2], npix)
	    }
	default:
	    call error (1, "unknown pixel datatype")
	}

	# Unmap the images.

	call imunmap (im2)
	call imunmap (im1)
	call xt_delimtemp (image2, Memc[imtemp])
	call sfree (sp)
end
# Apr 15 1994	Modify STSDAS playpen.newredshift
# Apr 15 1994	Change multiple types of WCS parameters at once
# Apr 15 1994	Change parameters to use only header velocity and axis
# Apr 15 1994	Modify STSDAS playpen.newredshift
# Apr 15 1994	Change multiple types of WCS parameters at once
# Apr 15 1994	Change parameters to use only header velocity and axis
# Apr 20 1994	Fix IMACCF calls
# Apr 22 1994	Fix sprintf str

# Mar 28 1995	Free all stack pointers

# OCt  1 1997	Always print 4-digit date in time stamp
