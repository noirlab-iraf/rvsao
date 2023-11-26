# File rvsao/Util/putspec.x
# May 2, 2014
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics

include <error.h>
include <imhdr.h>
include <imset.h>
include <fio.h>
include <mach.h>
include <fset.h>


# PUTSPEC --  Write spectrum pixels to IRAF image file

procedure putspec (im, ispec, specbuf, spectype, debug)

pointer	im		# IRAF image descriptor
int	ispec		# Number of spectrum in 2D spectrum to write
pointer	specbuf		# Spectrum data buffer
int	spectype	# Spectrum data type
bool	debug		# True for diagnostic output

int	npix
long	v[IM_MAXDIM], nlines
pointer	outline, specline
real	linemax, linemin
long	clktime()

errchk	malloc
errchk	spec_change_pix, spec_put_image_line, spec_pix_lim

begin
	if (IM_NDIM(im) == 0) {
	    call printf ("WPIXTEMP: No pixel file created\n")
	    return
	}

	npix = IM_LEN(im, 1)
	nlines = 1

	IM_MAX(im) = -MAX_REAL
	IM_MIN(im) = MAX_REAL

	call amovkl (long(1), v, IM_MAXDIM)
	if (ispec > 1)
	    v[2] = long (ispec)

	if (debug) {
	    call printf ("PUTSPEC:  %d pixels, type %d\n")
	    call pargi (npix)
	    call pargi (IM_PIXTYPE(im))
	    }

# Set output image line buffer
	specline = specbuf
	call spec_put_image_line (im, outline, v, IM_PIXTYPE(im))

# Change pixels to appropriate type
call spec_change_pix (specline, outline, npix, spectype, IM_PIXTYPE(im))

# Calculate image maximum and minimum
	call spec_pix_lim (outline, npix, IM_PIXTYPE(im), linemin, linemax)
	IM_MAX(im) = max (IM_MAX(im), linemax)
	IM_MIN(im) = min (IM_MIN(im), linemin)

	if (debug) {
	    call printf ("PUTSPEC:  minimum is %f, maximum is %f\n")
	    call pargr (IM_MIN(im))
	    call pargr (IM_MAX(im))
	    }
	IM_CTIME(im) = clktime (long(0))
	return
end


# SPEC_PUT_IMAGE_LINE -- Get a buffer pointer to output a line to an IRAF file

procedure spec_put_image_line (im, buf, v, data_type)

pointer	im			# IRAF image descriptor
pointer	buf			# Pointer to output image line
long	v[ARB]			# imio pointer
int	data_type		# output pixel type

int	impnll(), impnlr(), impnld(), impnlx(), impnls()
errchk	impnll, impnlr, impnld, impnlx, impnls

begin
	switch (data_type) {
	case TY_SHORT, TY_USHORT:
	    if (impnls (im, buf, v) == EOF)
		call error (3, "SPEC_PUT_IMAGE_LINE: Error writing IRAF file")
	case TY_INT, TY_LONG:
	    if (impnll (im, buf, v) == EOF)
		call error (3, "SPEC_PUT_IMAGE_LINE: Error writing IRAF file")
	case TY_REAL:
	    if (impnlr (im, buf, v) == EOF)
		call error (3, "SPEC_PUT_IMAGE_LINE: Error writing IRAF file")
	case TY_DOUBLE:
	    if (impnld (im, buf, v) == EOF)
		call error (3, "SPEC_PUT_IMAGE_LINE: Error writing IRAF file")
	case TY_COMPLEX:
	    if (impnlx (im, buf, v) == EOF)
		call error (3, "SPEC_PUT_IMAGE_LINE: Error writing IRAF file")
	default:
	    call error (10, "SPEC_PUT_IMAGE_LINE: Unsupported IRAF image type")
	}
	return
end


# SPEC_CHANGE_PIX -- Change a line of numbers to the IRAF image type

procedure spec_change_pix (inbuf, outbuf, npix, in_type, out_type)

pointer inbuf			# array of archive data
pointer	outbuf			# pointer to IRAF image line
int	npix			# number of pixels
int	in_type			# input pixel type
int	out_type		# output pixel type

begin
	switch (out_type) {
	    case TY_SHORT, TY_USHORT:
		switch (in_type) {
	    	    case TY_SHORT, TY_USHORT:
			Call achtss (Mems[inbuf], Mems[outbuf], npix)
	    	    case TY_INT, TY_LONG:
			Call achtls (Memi[inbuf], Mems[outbuf], npix)
		    case TY_REAL:
			Call achtrs (Memr[inbuf], Mems[outbuf], npix)
		    case TY_DOUBLE:
			Call achtds (Memd[inbuf], Mems[outbuf], npix)
		    default:
			call error (10, "SPEC_CHANGE_LINE: Illegal archive type")
		    }
	    case TY_INT, TY_LONG:
		switch (in_type) {
	    	    case TY_SHORT, TY_USHORT:
			Call achtsl (Mems[inbuf], Meml[outbuf], npix)
	    	    case TY_INT, TY_LONG:
			Call achtll (Memi[inbuf], Meml[outbuf], npix)
		    case TY_REAL:
			Call achtrl (Memr[inbuf], Meml[outbuf], npix)
		    case TY_DOUBLE:
			Call achtdl (Memd[inbuf], Meml[outbuf], npix)
		    default:
			call error (10, "SPEC_CHANGE_LINE: Illegal archive type")
		    }
	    case TY_REAL:
		switch (in_type) {
	    	    case TY_SHORT, TY_USHORT:
			Call achtsr (Mems[inbuf], Memr[outbuf], npix)
	    	    case TY_INT, TY_LONG:
			Call achtlr (Memi[inbuf], Memr[outbuf], npix)
		    case TY_REAL:
			Call achtrr (Memr[inbuf], Memr[outbuf], npix)
		    case TY_DOUBLE:
			Call achtdr (Memd[inbuf], Memr[outbuf], npix)
		    default:
			call error (10, "SPEC_CHANGE_LINE: Illegal archive type")
		    }
	    case TY_DOUBLE:
		switch (in_type) {
	    	    case TY_SHORT, TY_USHORT:
			Call achtsd (Mems[inbuf], Meml[outbuf], npix)
	    	    case TY_INT, TY_LONG:
			Call achtld (Memi[inbuf], Meml[outbuf], npix)
		    case TY_REAL:
			Call achtrd (Memr[inbuf], Meml[outbuf], npix)
		    case TY_DOUBLE:
			Call achtdd (Memd[inbuf], Meml[outbuf], npix)
		    default:
			call error (10, "SPEC_CHANGE_LINE: Illegal archive type")
		    }
	    case TY_COMPLEX:
		switch (in_type) {
	    	    case TY_SHORT, TY_USHORT:
			Call achtsx (Mems[inbuf], Memx[outbuf], npix)
	    	    case TY_INT, TY_LONG:
			Call achtlx (Memi[inbuf], Memx[outbuf], npix)
		    case TY_REAL:
			Call achtrx (Memr[inbuf], Memx[outbuf], npix)
		    case TY_DOUBLE:
			Call achtdx (Memd[inbuf], Memx[outbuf], npix)
		    default:
			call error (10, "SPEC_CHANGE_LINE: Illegal archive type")
		    }
	    default:
		call error (10, "SPEC_CHANGE_LINE: Illegal IRAF image type")
	    }
	return
end


# SPEC_PIX_LIMITS -- Determine the maximum and minimum values in a line

procedure spec_pix_lim (buf, npix, pixtype, linemin, linemax)

pointer	buf			# pointer to IRAF image line
int	npix			# number of pixels
int	pixtype			# output data type
real	linemax, linemin	# min and max pixel values

short	smax, smin
long	lmax, lmin
real	rmax, rmin
double	dmax, dmin
complex	xmax, xmin

begin
	switch (pixtype) {
	case TY_SHORT, TY_USHORT:
	    call alims (Mems[buf], npix, smin, smax)
	    linemax = smax
	    linemin = smin
	case TY_INT, TY_LONG:
	    call aliml (Meml[buf], npix, lmin, lmax)
	    linemax = lmax
	    linemin = lmin
	case TY_REAL:
	    call alimr (Memr[buf], npix, rmin, rmax)
	    linemax = rmax
	    linemin = rmin
	case TY_DOUBLE:
	    call alimd (Memd[buf], npix, dmin, dmax)
	    linemax = dmax
	    linemin = dmin
	case TY_COMPLEX:
	    call alimx (Memx[buf], npix, xmin, xmax)
	    linemax = xmax
	    linemin = xmin
	default:
	    call error (30, "SPEC_PIX_LIMITS: Unknown IRAF type")
	}
	return
end

# May  2 2014	New program baed on writetemp.x
