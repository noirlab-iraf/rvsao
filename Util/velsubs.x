# File Util/velsubs.x
# October 1, 1997
# By Doug Mink

include <imhdr.h>
include <error.h>

# VELS_GET -- Get header parameters which describe the wavelength axis.

procedure vels_get (im, axis, crval, cdelt)

pointer	im			# Image structure pointer
int	axis			# Dispersion axis
double	crval			# Reference wavelength (returned)
double	cdelt			# Wavelength step (returned)

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

pointer	im			# Image structure pointer
int	axis			# Dispersion axis
double	crval			# Reference wavelength
double	cdelt			# Wavelength step

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
# Jul 13 1994	New file

# Mar 28 1995	Free all stack pointers

# Oct  1 1997	Print 4-digit year in time stamp
