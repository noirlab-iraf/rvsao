# File Util/t_velset.x
# Modified by Doug Mink from stsdas.playpen.newredshift
# July 13, 1994

include	<imhdr.h>
include	<error.h>
include <time.h>

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
double	newvel				# Velocity or shift
bool	verbose				# Print operations?
bool	velz				# New redshift in Z (yes) or km/sec
bool	shift				# Velocity shift (yes) or velocity (no)

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
	shift = clgetb ("shift")
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
		call printf ("VELSHIFT: %s -> %s\n")
		    call pargstr (image1)
		    call pargstr (image2)
		call flush (STDOUT)

		# Do it.
		call shiftimage (image1, image2, velz, shift, newvel, verbose)
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
		call printf ("VELSHIFT: %s -> %s\n")
		    call pargstr (image1)
		    call pargstr (image2)
		call flush (STDOUT)
		call shiftimage (image1, image2, velz, shift, newvel, verbose)
	    }

	    call imtclose (list1)
	    call imtclose (list2)
	}
end

# SHIFTIMAGE -- Copy a spectrum, changing the wavelength scale accordingly
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

procedure shiftimage (image1, image2, velz, shift, velocity, verbose)

char	image1[ARB]	# Input spectrum
char	image2[ARB]	# Output spectrum
bool	velz		# Velocity in Z (yes) or km/sec (no)
bool	shift		# Velocity shift (yes) or velocity (no)
double	velocity	# Velocity or shift in km/sec
bool	verbose		# Print the operation

pointer	im, sp, str
int	dispaxis, dcflag
double	corre, crval2, cdelt2, crval1, cdelt1, c0, oldvel, oldz, newz
double	hcv,newvel
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
		call printf ("VELSHIFT: Old velocity is %f + %f = %f/n")
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
	
	if (shift)
	    newvel = oldvel + velocity
	else
	    newvel = velocity
	newz = (newvel - hcv) / c0
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
	    call sprintf (Memc[str], SZ_LINE, " VELSHIFT: oldz = %.8f, newz = %.8f ")
		call pargd (oldz)
		call pargd (newz)
	    }
	else {
	    call sprintf (Memc[str], SZ_LINE, " VELSHIFT: old vel = %.3f, new vel = %.3f ")
		call pargd (oldvel)
		call pargd (newvel)
	    }
	call vels_timelog (Memc[str], SZ_LINE)
	call imputh (im, "HISTORY", Memc[str])
	call sprintf (Memc[str], SZ_LINE, "input file to VELSHIFT: %s ")
	    call pargstr (image1)
	call imputh (im, "HISTORY", Memc[str])
	call sfree (sp)

	call imunmap (im)
end
