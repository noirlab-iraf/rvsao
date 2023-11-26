# File rvsao/Util/shdr.x
# October 6, 1995
# IRAF onedspec sphd subroutines
# Modified by Doug Mink, Harvard-SMithsonian Center for Astrophysics

include	<error.h>
include <imhdr.h>
include	<imset.h>
include	<mwset.h>
include <math/iminterp.h>
include	"shdr.h"
include	"mwcs.h"


# SPHD_OPEN    -- Open spectrum header structure.
# SPHD_CLOSE   -- Close and free spectrum header structure.
# SPHD_2D      -- Set/get physical dispersion axis and number of lines to sum.
# SPHD_COPY    -- Make a copy of an SHDR structure.
# SPHD_SYSTEM  -- Set or change the system.
# SPHD_LW      -- Logical to world coordinate transformation
# SPHD_WL      -- World to logical coordinate transformation
# SPHD_REBIN   -- Rebin spectrum to dispersion of reference spectrum
# SPHD_LINEAR  -- Rebin spectrum to linear dispersion
# SPHD_EXTRACT -- Extract a specific wavelength region
# SPHD_GI      -- Load an integer value from the header
# SPHD_GR      -- Load a real value from the header
# SPHD_GWATTRS -- Get spectrum attribute parameters
# SPHD_SWATTRS -- Set spectrum attribute parameters


# SPHD_OPEN -- Open spectrum header structure.
# This routine sets header information, WCS transformations, and extracts the
# spectrum from MULTISPEC and TWODSPEC format images.  The
# spectrum from a 2D/3D format is specified by a logical line and band
# number.  Optionally a MULTISPEC spectrum may be selected by it's aperture
# number.  The physical dispersion axis and summing parameter in TWODSPEC
# images are obtained by a call to SHDR_2D.  The access modes are header only
# or header and data.  Special checks are made to avoid repeated setting of
# the header and WCS information common to all spectra in a 2D format
# provided the previously set structure is input.  Note that the logical to
# world and world to logical transformations require that the MWCS pointer
# not be closed.

define  MW_NLOGDIM      Memi[$1+12]             # dimension of logical system

procedure sphd_open (im, index1, index2, ap, mode, sh)

pointer	im			# IMIO pointer
int	index1			# Image index desired
int	index2			# Image index desired
int	ap			# Aperture number desired
int	mode			# Access mode
pointer	sh			# SHDR pointer

int	format, daxisl, daxisp, np, nsum
int	i, j, k, l, aaxis, pndim, np1, np2, axno[3], axval[3]
real	apmin, apmax, amax
double	r[3], w[3], mw_c1trand()
real	asumr()
char	temp[16]
double	dval, aplow, aphigh, z, c0, sphd_lw()
bool	newim, streq(),debug
int	mw_stati(), strncmp()
pointer	sp, key, str, coeff, ct, mw_sctran(), imgs3r(), un_open()
pointer	mw, mw_openim()		# MWCS pointer
errchk	mw_sctran, imgstr, imgeti, imgetr, sphd_2d, un_open, sphd_gwattrs

define	data_	90

begin
	call smark (sp)
	call salloc (key, SZ_FNAME, TY_CHAR)
	call salloc (str, SZ_LINE, TY_CHAR)
	coeff = NULL
	c0 = 299792.5d0
	debug = FALSE

	# Allocate basic structure or check if the same spectrum is requested
	if (sh == NULL) {
	    call calloc (sh, LEN_SHDR, TY_STRUCT)
	    newim = true
	    }
	else {
	    call imstats (im, IM_IMAGENAME, Memc[str], SZ_LINE)
	    newim = !streq (Memc[str], SPECTRUM(sh))
	    if (!newim) {
		if (INDEX1(sh)==index1 && max(1,INDEX2(sh))==index2) {
		    if (IS_INDEFI(ap) || AP(sh)==ap) {
			if (CTLW(sh) != NULL && CTWL(sh) != NULL &&
			    ((mode==SHDATA && SY(sh)!=NULL) ||
			    (mode==SHHDR && SY(sh)==NULL))) {
			    call sfree (sp)
			    return
			    }
			else {
			    np1 = NP1(sh)
			    np2 = NP2(sh)
			    np = np2 - np1 + 1
			    goto data_
			    }
			}
		    }
		}
	    }

	# Set parameters common to an entire image
	if (newim) {
	    call imstats (im, IM_IMAGENAME, SPECTRUM(sh), LEN_SHDRS)
	    call strcpy (IM_TITLE(im), TITLE(sh), LEN_SHDRS)
	    IM(sh) = im

	#  Open MWCS descriptor
	    mw = mw_openim (im)
	    MW(sh) = mw

	    # Get standard parameters
	    call sphd_gi (im, "OFLAG", OBJECT, TYPE(sh))
	    call sphd_gr (im, "EXPOSURE", 1.0, IT(sh))
	    call sphd_gr (im, "ITIME", IT(sh), IT(sh))
	    call sphd_gr (im, "EXPTIME", IT(sh), IT(sh))
	    call sphd_gr (im, "RA", 0.0, RA(sh))
	    call sphd_gr (im, "DEC", 0.0, DEC(sh))
	    call sphd_gr (im, "UT", 0.0, UT(sh))
	    call sphd_gr (im, "ST", 0.0, ST(sh))
	    call sphd_gr (im, "HA", 0.0, HA(sh))
	    call sphd_gr (im, "AIRMASS", 0.0, AM(sh))
	    call sphd_gi (im, "DC-FLAG", DCNO, DC(sh))
	    call sphd_gi (im, "EX-FLAG", ECNO, EC(sh))
	    call sphd_gi (im, "CA-FLAG", FCNO, FC(sh))
	    call sphd_gd (im, "VELOCITY", 0.d0, VEL(sh))

	    # Flag bad airmass value; i.e. 0
	    if (!IS_INDEF (AM(sh)) && AM(sh) < 1.)
		AM(sh) = INDEF

	    # Determine the format and dispersion axis
	    call mw_gwattrs (mw, 0, "system", Memc[key], SZ_FNAME)
	    call mw_seti (mw, MW_USEAXMAP, NO)
	    pndim = mw_stati (mw, MW_NDIM)

#	    call printf ("SHDR_OPEN: %s ndim = %d, useaxmap = %d, ldim = %d\n")
#		call pargstr (Memc[key])
#		call pargi (pndim)
#		call pargi (mw_stati (mw, MW_USEAXMAP))
#		call pargi (MW_NLOGDIM(mw))

	    call mw_gaxmap (mw, axno, axval, pndim)

	    IF (STREq (Memc[key], "multispec"))
		format = MULTISPEC
	    else
		format = TWODSPEC

	    if (debug) {
		if (format == MULTISPEC)
		    call printf ("SHDR_OPEN: Multispec format\n")
		else
		    call printf ("SHDR_OPEN: Twodspec format\n")
		}

	    switch (format) {
	    case MULTISPEC:
		daxisp = 1
		daxisl = axno[daxisp]
		nsum = 1

		if (daxisl == 0) {
		    if (axval[daxisp] == 0)
			daxisl = daxisp
		    else
			call error (1, "No dispersion axis")
		    }

		CTLW1(sh) = mw_sctran (MW(sh), "logical", "multispec", 3)
		CTWL1(sh) = mw_sctran (MW(sh), "multispec", "logical", 3)
	    case TWODSPEC:
		nsum = 1
		call sphd_2d (im, daxisp, nsum)
		daxisl = max (1, axno[daxisp])
		if (IM_LEN(im,daxisl) == 1)
		    daxisl = mod (daxisl, 2) + 1

		i = daxisp
		do daxisp = 1, pndim
		    if (axno[daxisp] == daxisl)
			break
		if (i != daxisp) {
		    call eprintf (
		      "WARNING: Dispersion axis %d not found. Using axis %d.\n")
		    call pargi (i)
		    call pargi (daxisp)
		    }
		if (debug) {
		    call printf ("SHDR_OPEN: dispersion axis %d\n")
			call pargi (daxisp)
		    }

		CTLW1(sh) = mw_sctran (MW(sh), "logical", "world", daxisp)
		if (debug)
		    call printf ("SHDR_OPEN: logical -> world transform done\n")
		CTWL1(sh) = mw_sctran (MW(sh), "world", "logical", daxisp)
		if (debug)
		    call printf ("SHDR_OPEN: world -> logical transform done\n")

		# Check that the dispersion type makes sense.
		if (DC(sh) == DCLOG) {
		    w[1] = mw_c1trand (CTLW1(sh), 1.d0)
		    w[2] = mw_c1trand (CTLW1(sh), double (IM_LEN[im,daxisl]))
		    if (abs(w[1]) > 20. || abs(w[2]) > 20.)
			DC(sh) = DCLINEAR
		    if (debug) {
			call printf ("SHDR_OPEN: log wavelength from %.6f to %.6f\n")
			    call pargd (w[1])
			    call pargd (w[2])
			}
		    }
	        }

	    # Convert physical dispersion axis to logical dispersion axis
	    daxisl = axno[daxisp]
	    if (daxisl == 0) {
		if (axval[daxisp] == 0)
		    daxisl = daxisp
		else
		    call error (1, "No dispersion axis")
	        }
	    aaxis = 3 - daxisl

	    # Set labels
	    iferr (call mw_gwattrs (mw, daxisp, "label", LABEL(sh), LEN_SHDRS))
		call strcpy ("", LABEL(sh), LEN_SHDRS)
	    if (streq (LABEL(sh), "multispe"))
		call strcpy ("", LABEL(sh), LEN_SHDRS)
	    iferr (call mw_gwattrs (mw, daxisp, "units", UNITS(sh), LEN_SHDRS))
		call strcpy ("", UNITS(sh), LEN_SHDRS)
	    if (strncmp (LABEL(sh),"Pixel",5) == 0) {
		call sfree (sp)
		call mfree (sh, TY_STRUCT)
		sh = ERR
		return
		}

	    # Set units
	    UN(sh) = un_open (UNITS(sh))
	    MWUN(sh) = un_open (UNITS(sh))

	    FORMAT(sh) = format
	    NSUM(sh) = nsum
	    AAXIS(sh) = aaxis
	    DAXISP(sh) = daxisp
	    DAXIS(sh) = daxisl
	    NDIM(sh) = IM_NDIM(im)
	    PNDIM(sh) = pndim
	    if (NDIM(sh) < 3)
		IM_LEN(im,3) = 1
	    if (NDIM(sh) < 2)
		IM_LEN(im,2) = 1
	    }
	else {
	    format = FORMAT(sh)
	    aaxis = AAXIS(sh)
	    daxisp = DAXISP(sh)
	    daxisl = DAXIS(sh)
	    }

	# Set WCS parameters for spectrum type
	INDEX1(sh) = max (1, min (IM_LEN(im,aaxis), index1))
	if (index1 > IM_LEN(im,aaxis))
	    index1 = IM_LEN(im,aaxis)
	INDEX2(sh) = max (1, min (IM_LEN(im,3), index2))

	if (debug) {
	    if (format == MULTISPEC)
		call printf ("SHDR_OPEN: Multispec format\n")
	    else
		call printf ("SHDR_OPEN: Twodspec format\n")
	    }

	switch (format) {
	case MULTISPEC:
	    # If an aperture is specified first try and find it.
	    # If not specified or not found then use the index.

	    np = IM_LEN(im,1)
	    np1 = 1
	    ct = mw_sctran (mw, "logical", "physical", 2)
	    AP(sh) = 0
	    if (!IS_INDEFI(ap)) {
		do i = 1, IM_LEN(im,2) {
		    j = mw_c1trand (ct, double(i))
		    call sphd_gwattrs (mw, j, AP(sh), BEAM(sh), DC(sh), dval,
			dval, np2, z, aplow, aphigh, coeff) 
		    VEL(sh) = (z - 1.d0) * c0
		    APLOW(sh) = aplow
		    APHIGH(sh) = aphigh
		    if (AP(sh) == ap) {
			INDEX1(sh) = i
			break
			}
		    }
		}
	    if (AP(sh) != ap) {
		i = INDEX1(sh)
		j = mw_c1trand (ct, double(i))
		call sphd_gwattrs (mw, j, AP(sh), BEAM(sh), DC(sh), dval,
		    dval, np2, z, aplow, aphigh, coeff) 
		APLOW(sh) = aplow
		APHIGH(sh) = aphigh
		VEL(sh) = (z - 1.d0) * c0
		}

	    PINDEX1(sh) = j
	    call sprintf (Memc[key], SZ_LINE, "APID%d")
		call pargi (j)
	    iferr (call imgstr (im, Memc[key], TITLE(sh), LEN_SHDRS)) {
		call strcpy (IM_TITLE(im), TITLE(sh), LEN_SHDRS)
		if (AP(sh) > 0 && AP(sh) != INDEX1(sh)) {
		    call sprintf (temp,16,"[%d ap%d]")
			call pargi (INDEX1(sh))
			call pargi (AP(sh))
		    }
		else {
		    call sprintf (temp,16,"[%d]")
			call pargi (INDEX1(sh))
		    }
		call strcat (temp,TITLE(sh),LEN_SHDRS)
		}

	    call mw_ctfree (ct)
	case TWODSPEC:
	    np = IM_LEN(im,daxisl)

	    ct = mw_sctran (mw, "logical", "physical", 3B)
	    r[daxisp] = 1.d0
	    r[aaxis] = INDEX1(sh)
	    call mw_ctrand (ct, r, w, 2)
	    i = w[daxisp]
	    r[daxisp] = np
	    call mw_ctrand (ct, r, w, 2)
	    j = w[daxisp]
	    call mw_ctfree (ct)

	    np1 = min (i, j)
	    np2 = max (i, j)
	    #AP(sh) = w[aaxis]
	    #BEAM(sh) = w[aaxis]
	    AP(sh) = INDEX1(sh)
	    BEAM(sh) = INDEX1(sh)
	    apmin = AP(sh) - NSUM(sh) / 2
	    APLOW(sh) = max (1., apmin)
	    apmax = APLOW(sh) + NSUM(sh) - 1
	    amax = IM_LEN(im,aaxis)
	    APHIGH(sh) = min (amax, apmax)
	    apmin = APHIGH(sh) - NSUM(sh) + 1
	    APLOW(sh) = max (1., apmin)
	    NSUM(sh) = nint (APHIGH(sh)) - nint (APLOW(sh)) + 1
	    PINDEX1(sh) = w[aaxis]
	    call strcpy (IM_TITLE(im), TITLE(sh), LEN_SHDRS)
	    }
	
	# Set NP1 and NP2 in logical coordinates.
	ct = mw_sctran (mw, "physical", "logical", daxisp)
	i = max (1, min (int (mw_c1trand (ct, double (np1))), np))
	j = max (1, min (int (mw_c1trand (ct, double (np2))), np))
	call mw_ctfree (ct)
	np1 = min (i, j)
	np2 = max (i, j)
	np = np2 - np1 + 1

	NP1(sh) = np1
	NP2(sh) = np2
	SN(sh) = np

	if (debug) {
	    call printf ("SHDR:  about to read data %d - %d\n")
	    call pargi (np1)
	    call pargi (np2)
	    }

data_	CTLW(sh) = CTLW1(sh)
	CTWL(sh) = CTWL1(sh)

	# Set linear approximation.
	W0(sh) = sphd_lw (sh, double(np1))
	W1(sh) = sphd_lw (sh, double(np2))
	WP(sh) = (W1(sh) - W0(sh)) / (np2 - np1)
	SN(sh) = np2 - np1 + 1

	if (mode == SHDATA) {
	    # Set WCS array
	    if (SX(sh) != NULL)
		call mfree (SX(sh), TY_REAL)
	    call malloc (SX(sh), np, TY_REAL)
	    do i = np1, np2
		Memr[SX(sh)+i-np1] = real (sphd_lw (sh, double(i)))

	    # Set spectrum array
	    if (SY(sh) == NULL)
		call mfree (SY(sh), TY_REAL)
	    call malloc (SY(sh), np, TY_REAL)

	    i = max (1, INDEX1(sh))
	    j = max (1, INDEX2(sh))
	    switch (FORMAT(sh)) {
	    case MULTISPEC:
		call amovr (Memr[imgs3r(im,np1,np2,i,i,j,j)], Memr[SY(sh)], np)
	    case TWODSPEC:
		apmin = AP(sh) - NSUM(sh) / 2
	        APLOW(sh) = max (1., apmin)
		apmax = APLOW(sh) + NSUM(sh) - 1
		amax = IM_LEN(im,aaxis)
	        APHIGH(sh) = min (amax, apmax)
		apmin = APHIGH(sh) - NSUM(sh) + 1
		APLOW(sh) = max (1., apmin)
	        NSUM(sh) = nint (APHIGH(sh)) - nint (APLOW(sh)) + 1
		k = nint (APLOW(sh))
		l = nint (APHIGH(sh))
		nsum = l - k + 1
		if (daxisl == 1) {
		    do i = k, l {
			if (i == k)
		    	    call amovr (Memr[imgs3r(im,np1,np2,i,i,j,j)],
				Memr[SY(sh)], np)
			else
			     call aaddr (Memr[imgs3r(im,np1,np2,i,i,j,j)],
				 Memr[SY(sh)], Memr[SY(sh)], np)
			}
		    }
		else if (daxisl == 2) {
		    do i = np1, np2
			Memr[SY(sh)+i-np1] =
			    asumr (Memr[imgs3r(im,k,l,i,i,j,j)], nsum)
		    }
		}
	    }
	else {
	    call mfree (SX(sh), TY_REAL)
	    call mfree (SY(sh), TY_REAL)
	    }

	#if (PNDIM(sh) < 2) {
	#    INDEX1(sh) = 0
	#    PINDEX1(sh) = 0
	#    }
	#if (IM_NDIM(im) < 3) {
	#    INDEX2(sh) = 0
	#    PINDEX2(sh) = 0
	#    }

	call mfree (coeff, TY_CHAR)
	call sfree (sp)
end


# SHDR_CLOSE -- Close and free spectrum header structure.

procedure sphd_close (sh)

pointer	sh			# SHDR structure

begin
	if (sh != NULL) {
	    if (SX(sh) != NULL)
		call mfree (SX(sh), TY_REAL)
	    if (SY(sh) != NULL)
		call mfree (SY(sh), TY_REAL)
	    if (MW(sh) != NULL)
		call mw_close (MW(sh))
	    if (UN(sh) != NULL)
		call un_close (UN(sh))
	    if (MWUN(sh) != NULL)
		call un_close (MWUN(sh))
	    call mfree (sh, TY_STRUCT)
	    }
	return
end


# SHDR_2D -- Set/get physical dispersion axis and number of lines to sum.
# If the IMIO pointer is NULL then the values are set otherwise
# the values are returned.  If the default values are zero (the initial
# values) and they are not in the image header then they are queried
# from the CL.

procedure sphd_2d (im, daxisp, nsum)

pointer	im			# IMIO pointer (get/set flag)
int	daxisp			# Physical dispersion axis
int	nsum			# Number of lines to sum

int	da, ns, imgeti()
#int	clgeti()
data	da/0/, ns/0/
#errchk	clgeti

begin
	if (im == NULL) {
	    if (!IS_INDEFI (daxisp))
		da = daxisp
	    if (!IS_INDEFI (nsum))
		ns = nsum
	    return
	    }

	daxisp = da
	if (daxisp == 0) {
	    iferr (daxisp = imgeti (im, "DISPAXIS"))
		daxisp = 1
#		daxisp = clgeti ("dispaxis")
	    }
#	nsum = ns
#	if (nsum == 0)
#	    nsum = clgeti ("nsum")
end


# SHDR_LW -- Logical to world coordinate transformation
# The transformation pointer is generally NULL only after SHDR_LINEAR

double procedure sphd_lw (sh, l)

pointer	sh			# SHDR pointer
double	l			# Logical coordinate
double	w			# World coordinate

double	l1, l2, w1, mw_c1trand()

begin
	if (CTLW(sh) != NULL) {
	    switch (FORMAT(sh)) {
	    case MULTISPEC:
		call mw_c2trand (CTLW(sh), l, double (INDEX1(sh)), w, w1)
	    case TWODSPEC:
		w = mw_c1trand (CTLW(sh), l)
		if (DC(sh) == DCLOG)
		    w = 10. ** max (-20D0, min (20D0, w))
	    }
	    }
	else {
	    switch (DC(sh)) {
	    case DCLINEAR:
		w = W0(sh) + (l - 1) * WP(sh)
	    case DCLOG:
		w = W0(sh) * 10. ** (log10(W1(sh)/W0(sh)) * (l-1) / (SN(sh)-1))
	    case DCFUNC:
		w = W0(sh)
		call mw_c2trand (CTWL1(sh), w, double (INDEX1(sh)), l1, w1)
		w = W1(sh)
		call mw_c2trand (CTWL1(sh), w, double (INDEX1(sh)), l2, w1)
		if (SN(sh) > 1)
		    l1 = (l2 - l1) / (SN(sh) - 1) * (l - 1) + l1
		else
		    l1 = l - 1 + l1
		call mw_c2trand (CTLW1(sh), l1, double (INDEX1(sh)), w, w1)
	    }
	    }

	iferr (call un_ctrand (MWUN(sh), UN(sh), w, w, 1))
	    ;
	return (w)
end


# SHDR_WL -- World to logical coordinate transformation
# The transformation pointer is generally NULL only after SHDR_LINEAR

double procedure sphd_wl (sh, w)

pointer	sh			# SHDR pointer
double	w			# World coordinate
double	l			# Logical coordinate

double	w1, l1, l2, mw_c1trand()
int	fd, open()
pointer	ct

begin
	iferr (call un_ctrand (UN(sh), MWUN(sh), w, w1, 1))
	    w1 = w
	fd = open ("sphdwl.log", APPEND, TEXT_FILE)
	ct = CT_D(CTWL(sh))
	if (CTWL(sh) != NULL) {
	    switch (FORMAT(sh)) {
	    case MULTISPEC:
		call mw_c2trand (CTWL(sh), w1, double (INDEX1(sh)), l, l1)
	    case TWODSPEC:
		if (DC(sh) == DCLOG)
		    w1 = log10 (w1)
		l = mw_c1trand (CTWL(sh), w1)
	    }
	call fprintf (fd,"%11.5f: format %d, index %d, type %d: %11.5f %11.5f\n")
	    call pargd (w1)
	    call pargi (FORMAT(sh))
	    call pargi (INDEX1(sh))
	    call pargi (CT_TYPE(ct))
	    call pargd (l)
	    call pargd (l1)
	    }
	else {
	    switch (DC(sh)) {
	    case DCLINEAR:
		l = (w1 - W0(sh)) / WP(sh) + 1
	    case DCLOG:
		l = log10(w1/W0(sh)) / log10(W1(sh)/W0(sh)) * (SN(sh)-1) + 1
	    case DCFUNC:
		call mw_c2trand (CTWL1(sh), w1, double (INDEX1(sh)), l, l1)

		w1 = W0(sh)
		call mw_c2trand (CTWL1(sh), w1, double (INDEX1(sh)), l1, w1)
		w1 = W1(sh)
		call mw_c2trand (CTWL1(sh), w1, double (INDEX1(sh)), l2, w1)
		if (l1 != l2)
		    l = (SN(sh) - 1) / (l2 - l1) * (l - l1) + 1
		else
		    l = l - l1 + 1
	    }
	call fprintf (fd,"%11.5f -> %11.5f: DC %d, index %d: %11.5f\n")
	    call pargd (w)
	    call pargd (w1)
	    call pargi (DC(sh))
	    call pargi (INDEX1(sh))
	    call pargd (l)
	    }

	call close (fd)
	return (l)
end


# SHDR_GI -- Load an integer value from the header

procedure sphd_gi (im, field, default, ival)

pointer	im
char	field[ARB]
int	default
int	ival

int	dummy, imaccf(), imgeti()

begin
	ival = default
	if (imaccf (im, field) == YES) {
	    iferr (dummy = imgeti (im, field))
		call erract (EA_WARN)
	    else
		ival = dummy
	    }
end


# SHDR_GR -- Load a real value from the header

procedure sphd_gr (im, field, default, rval)

pointer	im
char	field[ARB]
real	default
real	rval

int	imaccf()
real	dummy, imgetr()

begin
	rval = default
	if (imaccf (im, field) == YES) {
	    iferr (dummy = imgetr (im, field))
		call erract (EA_WARN)
	    else
		rval = dummy
	    }
end


# SHDR_GD -- Load a double value from the header

procedure sphd_gd (im, field, default, dval)

pointer	im
char	field[ARB]
double	default
double	dval

int	imaccf()
double	dummy, imgetd()

begin
	dval = default
	if (imaccf (im, field) == YES) {
	    iferr (dummy = imgetd (im, field))
		call erract (EA_WARN)
	    else
		dval = dummy
	    }
end




# SHDR_GWATTRS -- Get spectrum attribute parameters

procedure sphd_gwattrs (mw, line, ap, beam, dtype, w1, dw, nw, z, aplow, aphigh,
	coeff)

pointer	mw				# MWCS pointer
int	line				# Physical line number
int	ap				# Aperture number
int	beam				# Beam number
int	dtype				# Dispersion type
double	w1				# Starting coordinate
double	dw				# Coordinate interval
int	nw				# Number of valid pixels
double	z				# Redshift factor
double	aplow, aphigh			# Aperture limits
pointer	coeff				# Nonlinear coeff string (input/output)

int	i, j, sz_coeff, strlen(), ctoi(), ctod()
pointer	sp, key
errchk	mw_gwattrs

data	sz_coeff /SZ_LINE/

begin
	call smark (sp)
	call salloc (key, SZ_FNAME, TY_CHAR)

	if (coeff != NULL)
	    call mfree (coeff, TY_CHAR)
	call malloc (coeff, sz_coeff, TY_CHAR)

	call sprintf (Memc[key], SZ_FNAME, "spec%d")
	    call pargi (line)

	call mw_gwattrs (mw, 2, Memc[key], Memc[coeff], sz_coeff)
	while (strlen (Memc[coeff]) == sz_coeff) {
	    sz_coeff = 2 * sz_coeff
	    call realloc (coeff, sz_coeff, TY_CHAR)
	    call mw_gwattrs (mw, 2, Memc[key], Memc[coeff], sz_coeff)
	}

	i = 1
	j = ctoi (Memc[coeff], i, ap)
	j = ctoi (Memc[coeff], i, beam)
	j = ctoi (Memc[coeff], i, dtype)
	j = ctod (Memc[coeff], i, w1)
	j = ctod (Memc[coeff], i, dw)
	j = ctoi (Memc[coeff], i, nw)
	j = ctod (Memc[coeff], i, z)
	j = ctod (Memc[coeff], i, aplow)
	j = ctod (Memc[coeff], i, aphigh)
	if (Memc[coeff+i-1] != EOS)
	    call strcpy (Memc[coeff+i], Memc[coeff], sz_coeff)
	else
	    Memc[coeff] = EOS

	if (j == 0)
	    call error (1, "Syntax error in spectrum attribute parameter")

	call sfree (sp)
end


# SHDR_SWATTRS -- Set spectrum attribute parameters

procedure sphd_swattrs (mw, line, ap, beam, dtype, w1, dw, nw, z, aplow, aphigh,
	coeff)

pointer	mw				# MWCS pointer
int	line				# Physical line number
int	ap				# Aperture number
int	beam				# Beam number
int	dtype				# Dispersion type
double	w1				# Starting coordinate
double	dw				# Coordinate interval
int	nw				# Number of valid pixels
double	z				# Redshift factor
double	aplow, aphigh			# Aperture limits
char	coeff[ARB]			# Nonlinear coeff string

int	sz_val, strlen()
pointer	sp, key, val

begin
	sz_val = strlen (coeff) + SZ_LINE

	call smark (sp)
	call salloc (key, SZ_FNAME, TY_CHAR)
	call salloc (val, sz_val, TY_CHAR)

	# We can't use SPRINTF for the whole string because it can only
	# handle a limited length and trucates long coefficient strings.
	# Use STRCAT instead.

	call sprintf (Memc[key], SZ_FNAME, "spec%d")
	    call pargi (line)
	call sprintf (Memc[val], sz_val, "%d %d %d %g %g %d %g %.2f %.2f")
	    call pargi (ap)
	    call pargi (beam)
	    call pargi (dtype)
	    call pargd (w1)
	    call pargd (dw)
	    call pargi (nw)
	    call pargd (z)
	    call pargd (aplow)
	    call pargd (aphigh)
	if (coeff[1] != EOS) {
	    call strcat (" ", Memc[val], sz_val)
	    call strcat (coeff, Memc[val], sz_val)
	}
	call mw_swattrs (mw, 2, Memc[key], Memc[val])

	call sfree (sp)
end
# Jul  1 1993	Match arguments to min and max for Decstations
# Jul 12 1993	Turn off debugging
# Aug 11 1993	If selected aperture exceeds max, set it to max
# Aug 20 1993	Make more variables double

# Apr 25 1994	Malloc instead of realloc SX and SY, when possible
# Apr 25 1994	Make SX single rather than double
# Jun 23 1994	Keep MWCS within SH structure

# Jan 10 1995	Set default positions and times to 0, not INDEF, exposure to 1
# Jan 19 1995	Drop out if multispec label is Pixel
# Jan 23 1995	Always set useaxmap to no (this may not be a good idea)
# Mar 28 1995	Close only opened structures
# Jun 19 1995	Add aperture and line to title if not separately named
# Oct  6 1995	Change SHDR_* subroutines to SPHD_*
