# File rvsao/Makespec/addspec.x if more than one file
# August 22, 2013
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics

# Copyright(c) 1995-2013 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# ADDSPEC adjusts one spectrum and adds it to a composite spectrum
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include <smw.h>
include "contin.h"
include "rvsao.h"

#  Emission line chopping flags
define  XC_CTYPES       "|no|yes|specfile|"
define  NOCHOP          1       # no correction
define  CHOP            2       # Always remove emission lines
define  SPECFILE        3       # Remove emission lines if spectrum CHOPEM=T


#rebinning defs
define RB_NEAREST	1	# nearest neighbor
define RB_LINEAR	2	# linear
define RB_POLY3		3	# 3rd order polynomial
define RB_POLY5		4	# 5th order polynomial
define RB_SPLINE3	5	# cubic spline
define RB_SINC		6	# sinc
define RB_LSINC		7	# look-up table sinc
define RB_DRIZZLE	8	# drizzle
define RB_SUMS		9
define RB_FUNCTIONS     "|nearest|linear|poly3|poly5|spline3|sinc|lsinc|drizzle|sums|"

define	LN2		0.69314718

procedure addspec (specfile,specdir,specap,mband,nspecap, fproc,
		   compim,comppath,compspec,ispec,nspec,perspec)

char	specfile[ARB]	# Object spectrum file name
char	specdir[ARB]	# Directory for object spectra
int	specap		# Object aperture to read from multispec file
int	mband	# Object band to read from multispec file
int	nspecap		# Number of input apertures to read in total
bool	fproc		# true if processing each input file and writing out
pointer compim		# Composite image header structure
char	comppath[ARB]	# Composite spectrum file path
pointer compspec	# Composite spectrum
int	ispec		# Number of spectra combined in this composite so far
int	nspec		# Number of spectra stacked in output file
bool	perspec		# True if limits set per spectrum

real	continuum	# Continuum level for output composite spectrum
pointer spectrum	# Object spectrum
pointer specim		# Object image header structure
pointer	work		# Temporary working storage for continuum removal
int	i, j
char	specpath[SZ_PATHNAME]	# Object spectrum path name
char	instrument[SZ_LINE]
double	tz		# Corrections to correlation Z
char	wtitle[20]	# Wavelength axis caption for plots of spectrum
char	title[SZ_LINE]	# Title for plots of spectrum
char	chopstr[SZ_LINE]	# Chop out emission lines
				# (yes | no | specfile)
char	keyword[16]	# Keyword for velocity shift
char	ctype1[16]	# Keyword value for x axis

double	zshift		# Z for composite spectrum
bool	normin,normout	# renormalization flag
char	tstring[SZ_LINE]
char	fstring[16]
char	temp[16]
char	compname[SZ_IMTITLE]	# Title for output spectrum
char	outpath[SZ_PATHNAME]	# Name of output file

pointer	shspec		# Rebinned spectrum pixels
pointer	wlcomp		# Wavelength vector for spectrum overlap
pointer	wlspec		# Wavelength vector for complete input spectrum
double 	logmin		# Log-10 of minimum wavelength for cross-correlation
double 	logmax		# Log-10 of maximum wavelength for cross-correlation
double	exp1		# Exposure time in seconds
double	dnormin, dnormout
double	dindef
double	pixshift
real	rindef
double	ra0,ra,dec0,dec
double	zcomp
double	vshift, vshift0
double	spshft		# Spectrum shift in log wavelength
bool	complog		# True if rebinning in log wavelength
int	lcomp		# Length of composite pathname
bool	pixfill, ispix
real	emrej1,abrej1
int	nsp, npts0
char	dstr[SZ_LINE]
int	npwcs
int	iord
int	nsig
int	ssmooth		# Number of times to smooth input spectrum
pointer	wlwcs
pointer	specwcs

double	dlog10()
double	wcs_p2w(), wcs_getshift()
bool	clgetb()
real	clgetr()
double	clgetd()
int	clscan()
int	clgeti()
#int	imaccess()
int	strlen(), ldir, imaccf(), strdic(), strcmp(), strncmp()
pointer sp

define	endspec_	 90
include	"rvsao.com"
include	"contin.com"
include	"sum.com"
 
begin
	dindef = INDEFD
	rindef = INDEFR
	call smark (sp)
	c0 = 299792.5
	specsh = NULL
	call sprintf (wtitle,20,"Wavelength")
	emrej1 = clgetr ("em_reject")
	abrej1 = clgetr ("abs_reject")
	iord = 0
	nsig = -1

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
	npts0 = npts
	call getspec (specpath, specap, mband, spectrum, specim, FALSE)
	if (specim == ERR) {
	    call eprintf ("ADDSPEC: Error reading spectrum %s\n")
	 	call pargstr (specpath)
	    go to endspec_
	    }
	npts = npts0

# Allocate rebinned spectra and wavelength vectors
	call salloc (shspec, npts, TY_REAL)
	call salloc (wlcomp, npts, TY_REAL)
	call salloc (work, npts, TY_REAL)
	call salloc (wlspec,specpix,TY_REAL)

# Set flag if input spectrum is in pixels instead of wavelength
	call imgspar (specim, "CTYPE1", ctype1, 16)
	if (strncmp (ctype1, "PIXEL", 5) == 0) {
	    ispix = TRUE
	    do j = 1, specpix {
		Memr[wlspec+j-1] = real (j)
		}
	    }
	else {
	    ispix = FALSE
	    do j = 1, specpix {
		Memr[wlspec+j-1] = real (wcs_p2w (double (j)))
		}
	    }

# Redshift velocity for composite spectrum
	velocity = clgetd ("velcomp")
	zcomp = clgetd ("zcomp")
	if (zcomp != dindef)
	    velocity = c0 * zcomp
	if (debug) {
	    call printf ("ADDSPEC: Z = %.5f, CZ = %.2f\n")
		call pargd (zcomp)
		call pargd (velocity)
	    call flush (STDOUT)
	    }

# If plot enabled, show the object spectrum.
	npwcs = NP2(specsh) - NP1(specsh) + 1
	specwcs = spectrum + NP1(specsh) - 1
	wlwcs = wlspec + NP1(specsh) - 1
	if (debug) {
	    call printf ("ADDSPEC: %d points file %s from %f-%f\n")
		call pargi (npwcs)
		call pargstr (specfile)
		call pargr (ymin)
		call pargr (ymax)
	    call flush (STDOUT)
	    }
	if (pltspec)
#	    call plotsum (npwcs,Memr[specwcs],specfile,
#			  Memr[wlwcs],wtitle,nsmooth, ymin, ymax)
	    call plotsum (npwcs,Memr[spectrum],specfile,
			  Memr[wlspec],wtitle,nsmooth, ymin, ymax)

# Eliminate bad lines from spectrum
	if (clscan("fixbad") != EOF)
	    call clgstr ("fixbad", fstring, 16)
	else
	    call strcpy ("no", fstring, 16)

# Fill in bad lines before removing continuum
	if (strcmp (fstring, "no") != 0 ) {
	    call filllist (npwcs, Memr[specwcs], Memr[wlwcs], iord, pixfill, debug)
	    if (pltspec) {
		call sprintf (title,SZ_LINE,"%s filling bad lines")
		    call pargstr (specfile)
		call plotsum (npwcs,Memr[specwcs],title,
			      Memr[wlwcs],wtitle,nsmooth,ymin,ymax)
		}
	    }

# Spectrum line chopping and continuum removal
	emrej[1] = emrej1
	abrej[1] = abrej1
	schop = clgetb ("reject")
	call clgstr ("cont_remove",chopstr,SZ_LINE)
        conproc = strdic (chopstr, chopstr, SZ_LINE, CONT_TYPE)
	if (debug) {
	    call printf ("ADDSPEC: Input spectra continuum ")
	    call printf ("(%d) ")
		call pargi (conproc)
	    if (conproc == DIVCONT)
		call printf ("divided\n")
	    else if (conproc == SUBCONT)
		call printf ("subtracted\n")
	    else if (conproc == ZEROCONT)
		call printf ("divided - 1\n")
	    else if (conproc == CONTFIT)
		call printf ("returned\n")
	    else
		call printf ("not removed\n")
	    }

# Renormalize spectrum if requested
	dnormin = clgetd ("normin")
	if (dnormin <= 0.0)
	    normin = FALSE
	else
	    normin = TRUE
	if (normin) {
	    call renormalize (Memr[spectrum], specpix, Memr[specwcs], npwcs,
			      dnormin, spmean, nsp, nsig, debug)
#	    if (nsp < 1)
#		go to endspec_

	    if (pltspec) {
		call sprintf (title,SZ_LINE,"%s renormalized")
		    call pargstr (specfile)
		call plotsum (npwcs,Memr[specwcs],title,
			      Memr[wlwcs],wtitle,nsmooth,ymin,ymax)
		}
	    }
	tz = 0.d0

# Add known velocity components as 1+z's
	if (spvel != dindef)
	    tz = (1.d0 + (spvel / c0))
	else
	    tz = 1.d0
	if (spechcv != 0.d0 && spechcv != dindef) {
	    tz = tz / (1.d0 + (spechcv / c0))
	    if (velhc == dindef)
		velhc = spechcv
	    }

# Read previous velocity shift from header for rebinning
	vshift = dindef
	vshift0 = dindef
	call sprintf (keyword, 16, "VSKY%d")
	    call pargi (ispec)
	call imgdpar (specim, keyword, vshift)

# Set z offset for rebinning

# If no desired velocity, eliminate barycentric correction if present
	if (velocity == dindef) {
	    if (spechcv != 0.0 && spechcv != dindef && velhc != dindef)
		zshift = tz * (1.d0 + (velhc / c0))

# If no desired velocity and no bcv remove vshift, if present
	    else if (vshift != dindef)
		zshift = tz * (1.d0 - (vshift / c0))

# If no desired velocity and no bcv and no vshift, do not shift at all
	    else
		zshift = 0.d0
	    }

# If there is a desired velocity, shift to it, taking out bcv, too
	else {
	    zshift = (1.d0 + (velocity / c0)) / tz
#	    if (spechcv != 0.0 && spechcv != dindef && velhc != dindef)
#		zshift = zshift * (1.d0 + (velhc / c0))
	    vshift0 = (zshift - 1.d0) * c0
	    }
	if (debug) {
	    if (spvel == dindef)
		call printf ("ADDSPEC: spvel= INDEF, ")
	    else {
		call printf ("ADDSPEC: spvel= %.2f, ")
		    call pargd (spvel)
		}
	    if (spechcv == dindef)
		call printf ("spechcv= INDEF, ")
	    else {
		call printf ("spechcv= %.2f, ")
		    call pargd (spechcv)
		}
	    if (velocity == dindef)
		call printf ("velocity= INDEF, ")
	    else {
		call printf ("velocity= %.2f, ")
		    call pargd (velocity)
		}
	    call printf ("zshift is %.5f\n")
		call pargd (zshift)
	    call flush (STDOUT)
	    }
	if (zshift == 0.d0)
	    spshft = 0.d0
	else
	    spshft = dlog10 (zshift)

# Set wavelength limits if writing one output file per input file
	if (perspec) {
	    if (spshft == 0.0) {
		minwav0 = 10.d0 ** spstrt
		maxwav0 = 10.d0 ** spfnsh
		}
	    else {
		minwav0 = 10.d0 ** (spstrt + spshft)
		maxwav0 = 10.d0 ** (spfnsh + spshft)
		}
	    }

# Rebin in log wavelength
	complog = clgetb ("complog")
	if (complog) {

# Figure lambda and log-lambda limits and increment
	    logmin = dlog10 (minwav0)
	    logmax = dlog10 (maxwav0)

# Fill in globally available values of wavelength
	    logw0 = logmin
	    dlogw = (logmax - logmin) / double (npts - 1)

	    if (debug) {
		call printf("ADDSPEC: spectrum %.7f to %.7f = %d %d, %d npts ->%.3f\n")
		    call pargd (spstrt)
		    call pargd (spfnsh)
		    call pargi (SN(specsh))
		    call pargi (DC(specsh))
		    call pargi (npwcs)
		    call pargd (10.d0 ** spshft)
		call flush (STDOUT)
		call printf("ADDSPEC: binned to %.7f to %.7f by %.9f, %d pts\n")
		    call pargd (logmin)
		    call pargd (logmax)
		    call pargd (dlogw)
		    call pargi (npts)
		call flush (STDOUT)
		}

#  Rebin spectrum to match desired region
	    pixshift = wcs_getshift()
	    call rebinl (Memr[spectrum],specsh, spshft, Memr[shspec], npts,
		     logw0, dlogw, pixshift)

# Set up wavelength vector for intermediate graphs
	    do j = 1, npts {
		Memr[wlcomp+j-1] = 10 ** (logw0 + (j-1) * dlogw)
		}
	    call sprintf (title,SZ_LINE,"%s rebinned in log wavelength")
		call pargstr (specfile)
	    }

# Rebin in linear wavelength
	else {
	    if (ispix) {
		dw = 1.d0
		minwav0 = 1.d0
		}
	    else
		dw = (maxwav0 - minwav0) / double (npts - 1)
	    if (debug) {
		call printf("ADDSPEC: spectrum %.3fA to %.3fA = %d %d, %d npts ->%.3f\n")
		    call pargd (10.d0 ** spstrt)
		    call pargd (10.d0 ** spfnsh)
		    call pargi (SN(specsh))
		    call pargi (DC(specsh))
		    call pargi (npwcs)
		    call pargd (10.d0 ** spshft)
		call flush (STDOUT)
		call printf("ADDSPEC: binned to %.3fA to %.3fA by %.3fA, %d pts\n")
		    call pargd (minwav0)
		    call pargd (maxwav0)
		    call pargd (dw)
		    call pargi (npts)
		call flush (STDOUT)
		}

#  Rebin spectrum to match desired region
	    pixshift = wcs_getshift()
	    call rebin (Memr[spectrum],specsh, spshft, Memr[shspec], npts,
		     minwav0, dw, pixshift, ispix)

# Set up wavelength vector for intermediate graphs
	    do j = 1, npts {
		Memr[wlcomp+j-1] = minwav0 + (double (j-1) * dw)
		}
	    call sprintf (title,SZ_LINE,"%s rebinned in wavelength")
		call pargstr (specfile)
	    }

# Renormalize rebinned spectrum if requested
	dnormout = clgetd ("normout")
	if (dnormout <= 0.0)
	    normout = FALSE
	else
	    normout = TRUE
	if (normout) {
	    call renormalize (Memr[shspec], npts, Memr[shspec], npts,
			      dnormout, spmean, nsp, nsig, debug)
#	    if (nsp < 1)
#		go to endspec_

	    if (pltspec) {
		call sprintf (title,SZ_LINE,"%s renormalized")
		    call pargstr (specfile)
		call plotsum (npwcs,Memr[shspec],title,
			      Memr[wlspec],wtitle,nsmooth,ymin,ymax)
		}
	    }

# Plot rebinned spectrum, if requested
	if (pltspec)
	    call plotsum (npts,Memr[shspec],title,Memr[wlcomp],wtitle,nsmooth,ymin,ymax)

# Number of times to smooth (1-2-1) final data plot
	ssmooth = clgeti ("spec_smooth")

# Smooth and plot spectrum, if requested
	if (ssmooth > 0) {
	    call smooth (Memr[shspec], npts, ssmooth)
	    call sprintf (title,SZ_LINE,"%s smoothed %d times")
		    call pargstr (specfile)
		    call pargi (ssmooth)
	    call plotsum (npts,Memr[shspec],title,Memr[wlcomp],wtitle,0,ymin,ymax)
	    }

# Remove continuum from rebinned spectrum if requested
	emrej[1] = emrej1
	abrej[1] = abrej1
	schop = clgetb ("reject")
	call clgstr ("contout",chopstr,SZ_LINE)
        conproc = strdic (chopstr, chopstr, SZ_LINE, CONT_TYPE)
	if (debug) {
	    call printf ("ADDSPEC: Output spectrum continuum ")
	    call printf ("(%d) ")
		call pargi (conproc)
	    if (conproc == DIVCONT)
		call printf ("divided\n")
	    else if (conproc == SUBCONT)
		call printf ("subtracted\n")
	    else if (conproc == ZEROCONT)
		call printf ("divided - 1\n")
	    else if (conproc == CONTFIT)
		call printf ("returned\n")
	    else
		call printf ("not removed\n")
	    }
	if (conproc != NOCONT)
	    call icsubcony (npts,Memr[shspec],Memr[wlcomp],compname,1,
			   nsmooth,Memr[work], ymin, ymax)

	call logtime (dstr,SZ_LINE)

# Set up path for output file
	lcomp = strlen (comppath)
	if (lcomp == 0) { 
	    call strcpy (specfile, outpath, SZ_PATHNAME)
	    call mkfits (outpath)
	    }
	else if (comppath[lcomp] == '/') {
	    call strcpy (comppath, outpath, SZ_PATHNAME)
	    call strcat (specfile, outpath, SZ_PATHNAME)
	    call mkfits (outpath)
	    }
	else {
	    call strcpy (comppath, outpath, SZ_PATHNAME)
	    }
	if (fproc && strlen (compname) < 1)
	    call strcpy (specname, compname, SZ_IMTITLE)

	if (debug) {
	    call printf ("ADDSPEC: Output to %s -> %s\n")
		call pargstr (comppath)
		call pargstr (outpath)
	    }

# If output file is not yet open, open it
	if (compim == NULL) {
	    if (copyhead)
		call tmp_open (compim, outpath, compspec, npts, nspec, specim)
	    else
		call tmp_open (compim, outpath, compspec, npts, nspec, NULL)
	    if (compim == ERR) {
		call eprintf ("ADDSPEC: Error writing composite spectrum %s\n")
	 	    call pargstr (outpath)
		go to endspec_
		}
	    call compinit (specim,compim,compname,ra0,dec0,normin,normout,
			   nspecap, nspec, dstr)
	    }

	if (ispec == 1 || fproc) {

	# Subtract continuum level from input spectrum
	    if (conproc != NOCONT)
		call icsubcony (npts,Memr[shspec],Memr[wlcomp],specname,1,
			       nsmooth,Memr[work], ymin, ymax)

	    if (strcmp (fstring, "zero") == 0) {
		call zerolist (npts, Memr[shspec], Memr[wlcomp], debug)
		if (pltspec) {
		    call sprintf (title,SZ_LINE,"%s zeroing bad lines")
			call pargstr (specfile)
		    call plotsum (npts,Memr[shspec],title,
				   Memr[wlcomp],wtitle,nsmooth,ymin,ymax)
		    }
		}

	# Move input spectrum to composite spectrum
	    call amovr (Memr[shspec],Memr[compspec],npts)

	# Add continuum level to composite spectrum
	    continuum = clgetr ("cont_add")
	    if (continuum != 0.0)
		call aaddkr (Memr[compspec], continuum, Memr[compspec], npts)

	# Enter in log files
	    do i = 1, nlogfd {
		call fprintf (logfd[i],"ADDSPEC %s making %s\n")
		    call pargstr (dstr)
		    call pargstr (outpath)
		call fprintf (logfd[i],"ADDSPEC %d: %s")
		    call pargi (ispec)
		    call pargstr (specpath)
		if (specap > 0) {
		    if (mband > 0) {
			call fprintf (logfd[i],"[%d,%d]")
			    call pargi (specap)
			    call pargi (mband)
			}
		    else {
			call fprintf (logfd[i],"[%d]")
			    call pargi (specap)
			}
		    }
		if (normin) {
		    call fprintf (logfd[i]," / %g7")
			call pargd (spmean)
		    }
		call fprintf (logfd[i]," added")
		if (velocity != dindef) {
		    if (spvel != dindef && spechcv != 0.d0) {
			call fprintf (logfd[i],", shifted %.2f - %.2f")
			    call pargd (spvel)
			    call pargd (spechcv)
			}
		    else 
			call fprintf (logfd[i]," shifted")
		    call fprintf (logfd[i]," to %.2f")
			call pargd (velocity)
		    }
		else if (spvel != dindef && spechcv != 0.d0) {
		    call fprintf (logfd[i],", shifted - %.2f")
			call pargd (spechcv)
		    }
		call fprintf (logfd[i],"\n")
		}
	    }

	else {

	# Compare positions; if different from last, delete position from header
	    exp1 = 1.0
	    if (nspecap != nspec) {
		if (imaccf (specim, "RA") == YES) {
		    call imgstr (specim, "RA", tstring,SZ_LINE)
		    call sscan (tstring)
			call gargd (ra)
		    }
		else
		    ra = ra0
		if (imaccf (specim, "DEC") == YES) {
		    call imgstr (specim, "DEC", tstring,SZ_LINE)
		    call sscan (tstring)
			call gargd (dec)
		    }
		else
		    dec = dec0
		if (ra != ra0 || dec != dec0) {
		    if (imaccf (compim, "RA") == YES)
			call imdelf (compim, "RA")
		    if (imaccf (compim, "DEC") == YES)
			call imdelf (compim, "DEC")
		    if (imaccf (compim, "EPOCH") == YES)
			call imdelf (compim, "EPOCH")
		    if (imaccf (compim, "EQUINOX") == YES)
			call imdelf (compim, "EQUINOX")
		    }

	# Read exposure time from input spectrum header
		exp1 = 0.d0
		call imgdpar (specim, "EXPTIME",exp1)
		if (exp1 .le. 0.d0)
		    call imgdpar (specim, "EXPOSURE", exp1)
		if (exp1 .le. 0.d0)
		    call imgdpar (specim, "ITIME", exp1)
		if (exp1 < .001d0)
		    exp1 = 1.d0
		}

	# If non-zero exposure time, add this spectrum
	    if (exp1 > 0.d0) {

	    # Remove continuum
		if (conproc != NOCONT)
		    call icsubcony (npts,Memr[shspec],Memr[wlcomp],specname,1,
				   nsmooth,Memr[work], ymin, ymax)

		if (strcmp (fstring, "zero") == 0) {
		    call zerolist (npts, Memr[shspec], Memr[wlcomp], debug)
		    if (pltspec) {
			call sprintf (title,SZ_LINE,"%s zeroing bad lines")
			    call pargstr (specfile)
			call plotsum (npts,Memr[shspec],title,
				      Memr[wlcomp],wtitle,nsmooth,ymin,ymax)
			}
		    }

		if (nspec > 1) {
		    call amovr (Memr[shspec],Memr[compspec],npts)

	    # Add continuum level to composite spectrum
		    continuum = clgetr ("cont_add")
		    if (continuum != 0.0)
			call aaddkr (Memr[compspec], continuum, Memr[compspec], npts)
		    }
		else {
		    call aaddr (Memr[compspec],Memr[shspec],Memr[compspec],npts)
		    texp = texp + exp1
		    call imaddd (compim, "EXPTIME", texp)
		    }

	# Enter in log files
		do i = 1, nlogfd {
		    call fprintf (logfd[i],"ADDSPEC %d: %s")
			call pargi (ispec)
			call pargstr (specpath)
		    if (specap > 0) {
			if (mband > 0) {
			    call fprintf (logfd[i],"[%d,%d]")
				call pargi (specap)
				call pargi (mband)
			    }
			else {
			    call fprintf (logfd[i],"[%d]")
				call pargi (specap)
			    }
			}
		    if (normin) {
			call fprintf (logfd[i]," / %g10")
			    call pargd (spmean)
			}
		    call fprintf (logfd[i]," added")
		    if (velocity != dindef) {
			call fprintf (logfd[i],", shifted")
			if (spechcv != 0.d0 && spechcv != dindef){
			    call fprintf (logfd[i]," %.2f - %.2f")
				call pargd (spvel)
				call pargd (spechcv)
			    }
		call fprintf (logfd[i]," to %.2f")
			    call pargd (velocity)
			}
		    else if (spechcv != 0.d0 && spechcv != dindef){
			call fprintf (logfd[i],", shifted - %.2f")
			    call pargd (spechcv)
			}

		    if (ispec >= nspec)
			call fprintf (logfd[i],"\n")
		    else
			call fprintf (logfd[i],"\r")
		    }
		}
	    else {
		do i = 1, nlogfd {
		    call fprintf (logfd[i],"ADDSPEC: %s skipped EXP is %.2f\n")
			call pargstr (specpath)
			call pargd (exp1)
		    }
		}
	    }

	if (ispec == 1 && normin)
	    call imaddd (compim, "SPNORM", spmean)

# Add APID for this spectrum
	if (nspec > 1)
	    call addapnum (specim, compim, ispec)

# Note addition of this spectrum to the composite spectrum
	if (imaccf (specim, "INSTRUME") == YES)
	    call imgstr (specim, "INSTRUME", instrument, SZ_LINE)
	else
	    call strcpy ("none", instrument, SZ_LINE);
	call sprintf (tstring,SZ_LINE,"Add %s %s")
	    call pargstr (instrument)
	    call pargstr (specpath)
	if (specap > 0) {
	    if (mband > 0) {
		call sprintf (temp,SZ_LINE,"[%d,%d]")
		    call pargi (specap)
		    call pargi (mband)
		}
	    else {
		call sprintf (temp,SZ_LINE,"[%d]")
		    call pargi (specap)
		}
	    call strcat (temp, tstring, SZ_LINE)
	    }
	if (spvel != dindef || (spechcv != dindef && spechcv != 0.d0)) {
	    if (velocity == dindef) {
		call strcat (" at indef CZ",tstring, SZ_LINE)
		if (spechcv != dindef && spechcv != 0.d0) {
		    call sprintf (temp,SZ_LINE," %.2f")
			call pargd (spechcv)
		    call strcat (temp, tstring, SZ_LINE)
		    }
		}
	    else {
		call strcat (" at CZ",tstring, SZ_LINE)
		if (spvel != dindef) {
		    call sprintf (temp,SZ_LINE," %.2f")
			call pargd (spvel)
		    call strcat (temp, tstring, SZ_LINE)
		    }
		if (spechcv != dindef && spechcv != 0.d0) {
		    call sprintf (temp,SZ_LINE," %.2f")
			call pargd (spechcv)
		    call strcat (temp, tstring, SZ_LINE)
		    }
		}
	    }

	if (save_names) {
	    if (vshift0 != dindef) {
		call sprintf (keyword, 16, "VSKY%d")
		    call pargi (ispec)
		call imaddd (compim, keyword, vshift0)
		}
	    else
		call imputh (compim, "HISTORY", tstring)
	    }

# Check for writability
#	if (imaccess (outpath, WRITE_ONLY) == NO) {
#	    call eprintf ("ADDSPEC: cannot write to %s; not saving results\n")
#		call pargstr (outpath)
#	    IM_UPDATE(compim) = NO
#	    }
#	else
	    IM_UPDATE(compim) = YES

# If plot enabled, show the current composite spectrum
	if (plttemp) {
	    if (nspec > 1) {
		call sprintf (title,SZ_LINE,"%s spectrum %d of %d")
		    call pargstr (compname)
		    call pargi (ispec)
		    call pargi (nspec)
		    }
	    else if (ispec > 1) {
		call sprintf (title,SZ_LINE,"%s using %d spectra")
		    call pargstr (compname)
		    call pargi (ispec)
		}
	    else if (strlen (compname) > 0) {
		call sprintf (title,SZ_LINE,"%s using %s")
		    call pargstr (compname)
		    call pargstr (specfile)
		}
	    else {
		call sprintf (title,SZ_LINE,"%s using %s")
		    call pargstr (outpath)
		    call pargstr (specfile)
		}
	    call plotsum (npts,Memr[compspec],title,
			   Memr[wlcomp],wtitle,tsmooth,ymin,ymax)
	    }

# Free rebinned spectra, wavelength, and correlation vectors
endspec_
	call sfree (sp)

# Close the object spectrum image
	IM_UPDATE(specim) = NO
	call close_image (specim, specsh)
	return

end


# COMPINIT initializes header of composite spectrum

procedure compinit (specim, compim, compname, ra0, dec0, normin, normout,
		    nspecap, nspec, dstr)

pointer	specim		# Object image header structure
pointer	compim		# Composite image header structure
char	compname[ARB]	# Title for output spectrum
double	ra0		# Right ascension of first object
double	dec0		# Declination of first object
bool	normin,normout	# renormalization flag
int	nspecap		# Number of spectra read from input file
int	nspec		# Number of spectra written to output file
char	dstr[ARB]	# Date string for file time stamp

double	exp
double	alt
double	dindef, dtemp
char	tstring[SZ_LINE]
char	interp_mode[SZ_LINE]
int	bin_mode
char	bintype[16]
char	ctype1[16]	# Keyword value for x axis
int	strlen()
bool	complog
bool	ispix
bool	clgetb()
int     clgwrd(), strncmp()

include "rvsao.com"
include "contin.com"
include "sum.com"

begin
	dindef = INDEFD

	call clgstr ("compname",compname,SZ_IMTITLE)
	if (strlen (compname) > 0)
	    call strcpy (compname,IM_TITLE(compim),SZ_IMTITLE)
	else {
	    call strcpy (IM_TITLE(specim),IM_TITLE(compim),SZ_IMTITLE)
	    call strcpy (IM_TITLE(compim),compname,SZ_IMTITLE)
	    }

# Set flag if input spectrum is in pixels instead of wavelength
	call imgspar (specim, "CTYPE1", ctype1, 16)
	if (strncmp (ctype1, "PIXEL", 5) > 0)
	    ispix = FALSE
	else
	    ispix = TRUE

# Get pointing direction
	ra0 = dindef
	call imgdpar (specim, "RA",ra0)
	dec0 = dindef
	call imgdpar (specim, "DEC",dec0)

# Delete IRAF WCS keywords if not same number of input and output apertures
	if (nspecap != nspec || nspecap == 1) {
	    call imgxpar (compim,"WCSDIM")
	    call imgxpar (compim,"WAT0")
	    call imgxpar (compim,"WAT1")
	    call imgxpar (compim,"WAT2")

	    call strcpy ("0:00:00", tstring, SZ_LINE)
	    call imcpstr (specim, compim, "RA", tstring)

	    call strcpy ("0:00:00", tstring, SZ_LINE)
	    call imcpstr (specim, compim, "DEC", tstring)

	    dtemp = 2000.0
	    call imcpd (specim, compim, "EPOCH", dtemp)
	    call imcpd (specim, compim, "EQUINOX", dtemp)
	    call imcpstr (specim, compim, "INSTRUME", tstring)
	    call imcpstr (specim, compim, "SITENAME", tstring)
	    call imcpstr (specim, compim, "SITELONG", tstring)
	    call imcpstr (specim, compim, "SITELAT", tstring)
	    alt = 0.d0
	    call imcpd (specim, compim, "SITEELEV", alt)
	    }

	exp = 0.d0
	call imgdpar (specim, "EXPTIME",exp)
	if (exp <= 0.d0)
	    call imgdpar (specim, "EXPOSURE", exp)
	if (exp <= 0.d0)
	    call imgdpar (specim, "ITIME", exp)
	if (exp < .001d0)
	    exp = 1.d0
	call imaddd (compim, "EXPTIME", exp)
	call imaddi (compim, "DISPAXIS", 1)
	texp = 0.0

	if (conproc == DIVCONT)
	    call sprintf (tstring, SZ_LINE, "Cont div, %s order %d, %d iter  %d -%d")
	else if (conproc == ZEROCONT)
	    call sprintf (tstring, SZ_LINE, "Cont div-1, %s order %d, %d iter  %d -%d")
	else if (conproc == SUBCONT)
	    call sprintf (tstring, SZ_LINE, "Cont sub, %s order %d, %d iter -%d %d")
	if (conproc != NOCONT) {
	    if (function == CN_SPLINE3)
		call pargstr ("spline3")
	    else if (function == CN_LEGENDRE)
		call pargstr ("legendre")
	    else if (function == CN_CHEBYSHEV)
		call pargstr ("chebyshev")
	    else
		call pargstr ("spline1")
	    call pargi (order)
	    call pargi (niterate)
	    call pargr (lowrej[1])
	    call pargr (hirej[1])
	    if (normin || normout)
		call strcat (" normalized", tstring, SZ_LINE)
	    call imastr (compim, "HISTORY", tstring);
	    }
	complog = clgetb ("complog")
	if (complog) {
	    call imaddi (compim, "DC-FLAG", 1)
	    call imastr (compim, "CTYPE1","LINEAR")
	    call imaddr (compim, "CRPIX1", 1.)
	    call imaddd (compim, "CRVAL1", logw0)
	    call imaddd (compim, "CDELT1", dlogw)
	    call imaddd (compim, "CD1_1", dlogw)
	    }
	else {
	    call imaddi (compim, "DC-FLAG", 0)
	    if (ispix)
		call imastr (compim, "CTYPE1","PIXEL")
	    else
		call imastr (compim, "CTYPE1","LINEAR")
	    call imaddr (compim, "CRPIX1", 1.)
	    call imaddd (compim, "CRVAL1", minwav0)
	    call imaddd (compim, "CDELT1", dw)
	    call imaddd (compim, "CD1_1", dw)
	    call imastr (compim, "CTYPE2","LINEAR")
	    call imaddi (compim, "CRPIX2", 1)
	    call imaddi (compim, "CRVAL2", 1)
	    call imaddi (compim, "CDELT2", 1)
	    call imaddi (compim, "CD2_2", 1)
	    }

	call imaddi (compim, "NP1", 1)
	call imaddi (compim, "NP2", npts)
	call imaddd (compim, "BCV", 0.d0)
	if (velocity != dindef)
	    call imaddd (compim, "VELOCITY", velocity)

# Get rebinning method
	bin_mode = clgwrd ("interp_mode",interp_mode,SZ_LINE,RB_FUNCTIONS)
	switch (bin_mode) {
	    case RB_NEAREST:
		call strcpy ("nearest", bintype, 16)
	    case RB_LINEAR:
		call strcpy ("linear", bintype, 16)
	    case RB_POLY3:
		call strcpy ("poly3", bintype, 16)
	    case RB_POLY5:
		call strcpy ("poly5", bintype, 16)
	    case RB_SPLINE3:
		call strcpy ("spline3", bintype, 16)
	    case RB_SINC:
		call strcpy ("sinc", bintype, 16)
	    case RB_LSINC:
		call strcpy ("lsinc", bintype, 16)
	    case RB_DRIZZLE:
		call strcpy ("drizzle", bintype, 16)
	    case RB_SUMS:
		call strcpy ("sums", bintype, 16)
	    }
	call sprintf (tstring,SZ_LINE,"wavelength interpolated %s")
		call pargstr (bintype)
	call imastr (compim, "INTERP", tstring)

	call sprintf (tstring,SZ_LINE,"rvsao.sumspec %s run %s")
	    call pargstr (VERSION)
	    call pargstr (dstr)
	call imputh (compim, "HISTORY", tstring)

	return
end


procedure imcpstr (imin, imout, keyword, tstring)

pointer imin		# Input image header structure
pointer imout		# Output image header structure
char	keyword[ARB]	# Keyword to be transferred
char	tstring[SZ_LINE] # value of keyword (returned)

int	imaccf()
begin
	if (imaccf (imin, keyword) == YES) {
	    call imgstr (imin, keyword, tstring,SZ_LINE)
	    call imastr (imout, keyword, tstring)
	    }
	return
end


procedure imcpsd (imin, imout, keyword, dtemp)

pointer imin		# Input image header structure
pointer imout		# Output image header structure
char	keyword[ARB]	# Keyword to be transferred
double	dtemp		# value of keyword (returned)

char	tstring[SZ_LINE] # value of keyword (returned)
int	imaccf()
begin
	if (imaccf (imin, keyword) == YES) {
	    call imgstr (imin, keyword, tstring,SZ_LINE)
	    call imastr (imout, keyword, tstring)
	    call sscan (tstring)
		call gargd (dtemp)
	    }
	return
end


procedure imcpd (imin, imout, keyword, dtemp)

pointer imin		# Input image header structure
pointer imout		# Output image header structure
char	keyword[ARB]	# Keyword to be transferred
double	dtemp		# value of keyword (returned)

int	imaccf()
begin
	if (imaccf (imin, keyword) == YES) {
	    call imgdpar (imin, keyword, dtemp)
	    call imaddd (imout, keyword, dtemp)
	    }
	return
end


procedure renormalize (specwcs,npwcs,spectrum,npix,dnorm,spmean,nsp,nsig,debug)

real	specwcs[ARB]
int	npwcs
real	spectrum[ARB]
int	npix
double	dnorm		
double	spmean
int	nsp		# Number of pixels used in computation
int	nsig		# Number of sigma beyond which to reject pixels
			# If < 0, scale by median instead
bool	debug

real	spsig, sphigh, splow, rsig, snorm
int	i, np

int	awvgr()
real	amedr()

begin

	snorm = real (dnorm)

# Divide by median value of spectrum if nsig < 0
	if (nsig < 0) {
	    spmean = amedr (specwcs, npwcs)
		    nsp = 0
	    do i = 1, npix {
		if (spectrum[i]!= 0.0) {
		    spectrum[i] = snorm * spectrum[i] / spmean
		    nsp = nsp + 1
		    }
		}
	    if (debug) {
		call printf ("RENORMALIZE: %.2f * spectrum / median = %.2f\n")
		    call pargd (snorm)
		    call pargd (spmean)
		}
	    }

# Else divide by mean value of spectrum rejecting outside of nsig sigma
	else {
	    call aavgr (specwcs, npwcs, spmean, spsig)

# If nsig > 0, iterate twice with limits
	    if (nsig > 0) {
		rsig = real (nsig)
		sphigh = spmean + rsig * spsig
		splow = spmean - rsig * spsig
		np = awvgr (specwcs, npwcs, spmean, spsig, splow, sphigh)
		sphigh = spmean + rsig * spsig
		splow = spmean - rsig * spsig
		nsp = awvgr (specwcs, npwcs, spmean, spsig, splow, sphigh)
		if (nsp < 1) {
		    spmean = 1.d0
		    if (debug)
			call printf ("RENORMALIZE: Spectrum all zeroes, 0 mean used\n")
		    }
		}
	    nsp = 0
	    do i = 1, npix {
		if (spectrum[i] != 0.0) {
		    spectrum[i] = snorm * spectrum[i] / spmean
		    nsp = nsp + 1
		    }
		}
	    if (nsp < 1 && debug) {
		call printf ("RENORMALIZE: Spectrum all zeroes, 0 mean used\n")
		}
	    if (debug) {
		call printf ("RENORMALIZE: %.2f * spectrum / mean = %.2f\n")
		    call pargd (snorm)
		    call pargd (spmean)
		}
	    }
	return
end


# Change filename to explicit FITS file

procedure mkfits (filename)

char	filename[ARB]	# File name to fix

int	iext, strsearch()

begin
	iext = strsearch (filename, ".imh")
	if (iext > 0) {
	    filename[iext-4] = EOS
	    call strcat (".fits", filename, SZ_PATHNAME)
	    return
	    }
	iext = strsearch (filename, ".fit")
	if (iext == 0)
	    call strcat (".fits", filename, SZ_PATHNAME)

	return
end

procedure addapnum (imin, imout, i)

pointer	imin	# Input image structure
pointer	imout	# Output image structure
int	i	# aperture

int ap, beam, dtype, nw
real aplow[2], aphigh[2]
double w1, dw, z
pointer coeff, sp, key, str1, str2
pointer smw, smw_openim()
int imaccf()

begin
	if (imaccf (imin, "WAT2_002") == NO)
	    return
	call smark (sp)
	call salloc (key, SZ_FNAME, TY_CHAR)
	call salloc (str1, SZ_LINE, TY_CHAR)
	call salloc (str2, SZ_LINE, TY_CHAR)

	smw = smw_openim (imin)
	call smw_gwattrs (smw, i, 1, ap, beam, dtype, w1, dw, nw, z,
                    aplow, aphigh, coeff)
	call sprintf (Memc[str1], SZ_LINE, "%d %d")
	    call pargi (ap)
	    call pargi (beam)
	if (!IS_INDEF(aplow[1]) || !IS_INDEF(aphigh[1])) {
	    call sprintf (Memc[str2], SZ_LINE, " %.2f %.2f")
		call pargr (aplow[1])
		call pargr (aphigh[1])
	    call strcat (Memc[str2], Memc[str1], SZ_LINE)
	    if (!IS_INDEF(aplow[2]) || !IS_INDEF(aphigh[2])) {
		call sprintf (Memc[str2], SZ_LINE, " %.2f %.2f")
		    call pargr (aplow[2])
		    call pargr (aphigh[2])
		call strcat (Memc[str2], Memc[str1], SZ_LINE)
		}
	    }
	if (i < 1000)
	    call sprintf (Memc[key], SZ_FNAME, "APNUM%d")
	else
	    call sprintf (Memc[key], SZ_FNAME, "AP%d")
	    call pargi (i)
	call imastr (imout, Memc[key], Memc[str1])
	call mfree (coeff, TY_CHAR)
	call sfree (sp)

	return
end


# Jul 21 1995	New program
# Oct  6 1995	Do not add spectrum to composite unless exposure time is > 0
# Oct 11 1995	Add continuum removal

# Jan 10 1997	Change parameter names to show that em and abs lines are cut
# Jan 14 1997	Change tempvel to veltemp
# Jan 14 1997	Adjust continuum processing to current RVSAO methods
# Feb  3 1997	Drop unused parameter t_chop
# Feb  4 1997	Drop declaration of C0; it is now in rvsao.com
# Feb  4 1997	Fix calls to plotspec
# Feb  5 1997	Fix HISTORY statement; fix EPOCH transfer
# Mar 14 1997	Drop SPECSH from GETSPEC call; more it to common
# Apr 18 1997	Fix unrebinned graph of input spectra
# Apr 21 1997	Move wavelength limit setting to T_SUMTEMP
# Apr 21 1997	Add barycentric velocity correction to HISTORY
# Apr 25 1997	Add composite title from parameter tempobj
# Apr 25 1997	Add composite continuum from parameter continuum
	# Apr 29 1997	Change name from ADDTEMP to ADDSPEC and variable names, too
# May  1 1997	Include SMW.H instead of SHDR.H
# May  2 1997	Fix title in composite graph
# May  2 1997	Always test against dindef, not INDEFD
# May  5 1997	Drop S_ from REJECT parameters and change SCHOP to REJECT
# May  5 1997	Change S_CONT parameter to CONT_REMOVE
# May  5 1997	Change CONTINUUM parameter to CONT_ADD
# May  6 1997	Drop position from header if different objects added
# May  6 1997	Add NSMOOTH argument to ICSUBCON
# May 16 1997	Add option to rebin to linear wavelength
# May 16 1997	If VELOCITY is INDEF, do not shift or write to header
# May 21 1997	Add renormalization after continuum removal and rebinning
# May 21 1997	Fix option to rebin to linear wavelength; improve logging
# Jun  5 1997	Add FIXBAD and BADLINES to remove specific wavelength regions
# Jun 11 1997	Add end of line to log listing
# Jun 17 1997	Clean up handling of undefined velocities
# Jul 21 1997	Reset NPTS so it isn't change by the value in the input header
# Jul 22 1997	Make logging of input filenames in header optional
# Jul 22 1997	Move normalization from before to after continuum removal
# Aug 27 1997	Add multspec band as argument
# Dec 17 1997   Use EQUINOX if it is present as well as EPOCH
#
# Apr 17 1998	Change ctype from pixel to linear
# Apr 22 1998	Drop getim.com; extract needed parameters from header locally
# Jun 12 1998	Use only portions of spectra with WCS
# Dec 16 1998	Write normalization factor for first spectrum
# Dec 16 1998	Write indefinite velocity to history if velocity not being used
# Dec 18 1998	Add option to rebin 2D spectrum to 2D spectrum

# Mar 11 1999	Add write argument to getspec()
# Mar 18 1999	When renormalizing summed string, multiply by 1000
# Mar 19 1999	Fix logic for multiple spectrum summing
# May 11 1999	Suppress line feeds if adding multiple spectra to a file
# May 11 1999	Check for output file WRITE_ONLY, not READ_WRITE
# May 11 1999	Fix spectrum stacking
# Jun 10 1999	Add header HISTORY line for continuum removal
# Jun 10 1999	Keep average exposure if stacking spectra rather than adding
# Jun 10 1999	Add option of zeroing bad lines after continuum removal
# Jun 16 1999	Add normalization to HISTORY
# Jul 27 1999	Add option to copy entire input header to output file
# Jul 28 1999	Change normalization parameters to numbers
# Aug 18 1999	Add divide and subtract one continuum removal option
# Sep  1 1999	Set CD1_1 as well as CDELT1

# Jul  6 2000	Drop spmean declaration; it is in rvsao.com
# Jul  6 2000	Move renormalization to subroutine
# Jul 21 2000	Add barycentric velocity correction if requested for indef vel.
# Oct 24 2000	Use none if no instrument keyword in header

# Feb 13 2001	Add binning type to title
# Mar 19 2001	Delete IRAF WCSDIM and WAT keywords from output file header
# Apr  3 2001	Fix renormalization argumnet mismatch
# Apr 27 2001	Add option to rebin individual files to output directory

# Mar 29 2002	Add per spectrum wavelength limit option
# May 30 2002	Add option to read and write VSKYn in header

# Jun  3 2003	Add option to smooth spectrum, not just for graphing

# May  7 2004	Add options to return fit continuum in place of spectrum
# May 25 2004	Only add exposure times if single output spectrum
# Aug  3 2004	Remove continuum from input spectra using cont_remove, not contout
# Aug 25 2004	Change all uses of spmean to double precision
# Aug 25 2004	make sure strcpy always has 3 arguments

# Mar 23 2005	Change calls from plotspec() to plotsum() with scaling
# May 25 2005	Set pixel shift before rebinning
# Aug 30 2005	Add copyhead in sum.com to copy header from first spectrum
# Sep  7 2005	Renormalize by square root of sum of squares of intensity

# Mar  9 2006	Always set BCV in output file to zero; it's always corrected
# Apr 25 2006	Normalize all-zero spectra, too (divide by 1.0)
# Sep 18 2006	Call icsubcony() instead of icsubcon() to set graph flux/count limits
# Dec 18 2006	Add number of input apertures as argument
# Dec 18 2006	Do not change header if multispec to multispec copy

# Feb 14 2007	Add APIDi to multispec output file from MWCS structure
# Mar 27 2007	Add exposure times only if spectra are being added

# Jan 11 2008	Add fproc to argument list so it can be set correctly
# Jan 16 2008	Remove barycentric velocity correction when summing to V
# Mar  5 2008	Add pixfill argument to filllist() call
# May  8 2008	Add order=0 argument to filllist() call
# May  9 2008	Use IRAF subroutines to compute renormalization mean or median
# May  9 2008	Fix bug so constant continuum is correctly added
# May 12 2008	Fix bug so plotsum() not plotspec() is called for composite
# May 30 2008	Drop declarations of unused variables in renormalize
# Sep 16 2008	Keep position if number of input apertures = output apertures 

# Jan  7 2009	Do not add back velh except for velocity=INDEF
# Jan  9 2009	Delete IRAF WCS keywords if nspecap=1
# Apr 24 2009	Fix bug so that data with PIXEL WCS is dealt with correctly
# Jun 10 2009	Fix pixel flag setting

# Aug 22 2013	Change specband to mband to avoid conflict with rvsao.com
