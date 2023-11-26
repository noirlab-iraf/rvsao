# File rvsao/Subpix/pixshift.x
# June 9, 2008
# By Doug Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss, Johns Hopkins University and others

# Copyright(c) 2008 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# PIXSHIFT computes the fractional pixel shift and velocity for one spectrum
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include <smw.h>
include "rvsao.h"
include	"xcv.h"
include	"contin.h"

procedure pixshift (specfile,specdir,mspec,mband)

char	specfile[ARB]	# Object spectrum file name
char	specdir[ARB]	# Directory for object spectra
int	mspec		# Object aperture to read from multispec file
int	mband		# Object band to read from multispec file

pointer spectrum	# Object spectrum
pointer specim		# Object image header structure

bool	arcflag		# archive record writing flag
bool	newfit		# true to recompute pixel shift
bool	ispix		# TRUE if returning result in pixel space, else FALSE

int	i, j, itm0
char	specpath[SZ_PATHNAME]	# Object spectrum path name
char	tempfiles[SZ_PATHNAME]	# List of template spectra
char	tempfile[SZ_PATHNAME]	# Template spectrum file name
char	temppath[SZ_PATHNAME]	# Template spectrum path name
char	tempdir[SZ_PATHNAME]	# Directory for template spectra
char	tempnums0[SZ_LINE]	# List of multispec spectra to read
char	tempnums[SZ_LINE]	# List of multispec spectra to read
char	txfile[SZ_PATHNAME,MAXTEMPS]	# Correlation file names
int	rmode			# Report format (1=normal,2=one-line)
double	tz			# Corrections to correlation Z
char	xlab[SZ_LINE]		# Title for wavelength plots of spectrum

pointer	specwcs			# Pointer to first data point with WCS
pointer	wlwcs			# Pointer to first wavelength with WCS
int	npspec			# Number of spectrum points with WCS
pointer	tempwcs			# Pointer to first template data point with WCS
pointer	twlwcs			# Pointer to first template wavelength with WCS
int	nptemp			# Number of template points with WCS

int	ncvel[MAXTEMPS]	# Number of points in each correlation

double	zshift		# guessed Z for template offset
double	czguess		# initial velocity or Z for template offset
double	tmstrt,tmfnsh	# template limits after zguess applied
int	izpass		# Count for template shifts
int	nzpass		# Number of template shifts
bool	zinit		# TRUE if initial value for z is set
int	tspec		# Template spectrum to read from multispec file
int	ntpts		# ntemp * npts2
bool	divcon		# TRUE if template header says to divide by continuum

pointer	templist	# List of template files

double	tcz		# Velocity returned from xcorfit
double	tczerr		# Velocity error returned from xcorfit
double	tczr		# R returned from xcorfit
double	rmax		# Maximum R value for spectrum and best template
double 	minwav0,maxwav0	# Starting,ending wavelengths from parameter file
double 	minwav1,maxwav1	# Overlap wavelength limits from spectrum and template
double 	minwav,maxwav	# Wavelength limits for pixel shift computation
double 	logmin,logmax	# Log-10 of wavelength limits for pixel shift computation
bool	filexcor	# Flag to save correlation to a file
double	spwl1,spwl2	# Spectrum first and last wavelengths in angstroms
double	tmwl1,tmwl2	# Template first and last wavelengths in angstroms
int	ntx		# Number of templates actually used
int	tband		# Multispec band for templates
double	specd, spmin, spmax, spfloor
pointer	speci
double 	sptot
double	spshft, tmshft	# Spectrum shift in log wavelength
double	wcs_getshift(), wcs_p2w()
double	dlog10()
int	itemp, nsp, ncor, npts2, iord
double	tmpvel
double	pixshift, tpixshift

pointer tempim		# Template image header structure
pointer	tempsh		# Template spectrum header structure

int	ntspec,ntspec0	# Number of template multispec spectra
int	tspec_range[3,MAX_RANGES]
bool	echelle		# If true, template multispec numbers track object
int	ip,jp,lfile	# Limits for multispec aperture decoding
char	lbracket[3]	# left bracket, parenthesis, brace
char	rbracket[3]	# right bracket, parenthesis, brace
char	corstring[SZ_LINE] #
double	dindef
double	factor		# additional factor for renormalized data
char	cdot, cslash, cdollar, cdash, colon
char	tempsec[SZ_LINE]
bool	wenable		# If yes, open object spectrum for updating
bool	pixfill
int	ilb, iw
int	np1, np2
int	lo0, toplo0, tpnrn0, nrun0
 
bool	clgetb()
int	clgeti(), clscan()
real	clgetr(), max(), min()
double	clgetd()
int	strdic(), stridx(), stridxs()
int	decode_ranges(),get_next_number()
int	imtgetim(), imaccess(), strlen(), ldir
pointer	imtopen()

include	"rvsao.com"
include "results.com"
include "contin.com"
include "xcor.com"

define	refit_	 30
define	newtemp_ 40
define	newtap_	 50
define	nextemp_ 60
define	endtemp_ 65
define	endspec_ 70
define	endcorr_ 80
define	endxc_	 90
 
begin
	c0 = 299792.5
	cdot = '.'
	cslash = '/'
	cdash = '-'
	colon = ':'
	cdollar = '$'
	dindef = INDEFD
	specsh = NULL
	tempsh = NULL
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	wenable = savevel
	pixshift = 0.0
	tpixshift = 0.0

	correlate = COR_VEL
	call clgstr ("correlate", corstring, SZ_LINE)
	correlate = strdic (corstring, corstring, SZ_LINE, COR_TYPES)

	if (correlate == COR_NO) {
	    if (savevel) {
		wenable = TRUE
		savevel = FALSE
		}
	    else
		wenable = FALSE
	    }
	newresults = FALSE
	newfit = FALSE

# Echelle switch (template spectrum tracks object spectrum number if true)
	echelle = clgetb ("echelle")

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

# Check for readability of object spectrum
	if (imaccess (specpath, READ_ONLY) == NO) {
	    call eprintf ("SHIFTPIX: cannot read spectrum path %s \n")
		call pargstr (specpath)
	    return
	    }

	call getspec (specpath, mspec, mband, spectrum, specim, wenable)
	if (specim == ERR) {
	#   call eprintf ("SHIFTPIX: Error reading spectrum %s\n")
	# 	call pargstr (specpath)
	    go to endxc_
	    }
	pixshift = wcs_getshift()
	call emrhead (mspec, specim)
	if (correlate == COR_NO)
	    call xcrhead (mspec, specim)

	call strcpy (LABEL(specsh), xlab, SZ_LINE)
	if (strlen (xlab) <= 0)
	    call strcpy ("Wavelength", xlab, SZ_LINE)
	if (strlen (UNITS(specsh)) > 0) {
	    call strcat (" in ",xlab,SZ_LINE)
	    call strcat (UNITS(specsh),xlab,SZ_LINE)
	    }
	else
	    call strcat (" in Angstroms",xlab,SZ_LINE)

	np1 = NP1(specsh)
	np2 = NP2(specsh)
	if (debug) {
	    call printf ("SHIFTPIX: %d-pixel spectrum from %.3fA(%d) - %.3fA(%d)\n")
		call pargi (specpix)
		call pargd (wcs_p2w (double (np1)))
		call pargi (np1)
		call pargr (Memr[SX(specsh)+NP2(specsh)-1])
		call pargd (wcs_p2w (double (np2)))
	    call flush (STDOUT)
	    }

# Turn off spectrum continuum subtraction, if not needed
	tcont = clgetb ("subcont")
	if (tcont) {
	    contfit = TRUE
	    conproc = SUBCONT
	    }
	else {
	    conproc = NOCONT
	    contfit = FALSE
	    }
#	if (debug) {
#	    if (contfit)
#		call printf ("SHIFTPIX: Spectrum continuum will be subtracted\n")
#	    else
#		call printf ("SHIFTPIX: Spectrum continuum already subtracted\n")
#	    }

# Eliminate zeroes at left end of spectrum
#	while (Memr[spectrum+np1-1] == 0.0) {
#	    np1 = np1 + 1
#	    }
#	if (debug) {
#	    call printf ("SHIFTPIX: %d pixels were zero at blue end of spectrum\n")
#		call pargd (np1 - NP1(specsh))
#	    call flush (STDOUT)
#	    }

# Eliminate zeroes at right end of spectrum
#	while (Memr[spectrum+np2-1] == 0.0)
#	    np2 = np2 - 1
#	if (debug) {
#	    call printf ("SHIFTPIX: %d pixels were zero at red end of spectrum\n")
#		call pargd (NP2(specsh) - np2)
#	    call flush (STDOUT)
#	    }

# Set pixel limits of spectrum WCS
	specwcs = spectrum + np1 - NP1(specsh)
#	wlwcs = SX(specsh) + np1 - 1
	npspec = np2 - np1 + 1
	call malloc (wlwcs, npspec, TY_REAL)
	do i = np1, np2 {
	    Memr[wlwcs+i-np1] = real (wcs_p2w (double(i)))
	    }

# If plot enabled, show the object spectrum.
	if (pltspec)
	    call plotspec (npspec, Memr[specwcs], specid,
			   Memr[wlwcs], xlab, nsmooth)

# Eliminate bad lines from spectrum
	pixfill = FALSE
	if (clscan("fixbad") != EOF) {
	    if (clgetb ("fixbad")) {
		if (echelle)
		    iord = mspec
		else
		    iord = 0
		call filllist (npspec,Memr[specwcs],Memr[wlwcs],iord,pixfill,debug)
		if (pltspec)
		    call plotspec (npspec,Memr[specwcs],
				   specid,Memr[wlwcs],xlab,nsmooth)
		}
	    }

# Compute minimum and maximum values in spectrum
	spmax = Memr[specwcs]
	spmin = Memr[specwcs]
	do i = 2, npspec {
	    specd = Memr[specwcs+i-1]
	    if (specd > spmax)
		spmax = specd
	    if (specd < spmin)
		spmin = specd
	    }
	if (debug) {
	    call printf ("SHIFTPIX: %8g < counts < %8g from %.3fA(%d) - %.3fA(%d)\n")
		call pargd (spmin)
		call pargd (spmax)
		call pargr (Memr[wlwcs+np1-1])
		call pargi (np1)
		call pargr (Memr[wlwcs+np2-1])
		call pargi (np2)
	    call flush (STDOUT)
	    }
	if (debug) {
	    call printf ("SHIFTPIX: maximum count in spectrum is %f\n")
		call pargd (spmax)
	    call printf ("         minimum count in spectrum is %f\n")
		call pargd (spmin)
	    call flush (STDOUT)
	    }

# If RENORM is not set, set it if maximum counts are less than 1
	renorm = clgetb ("renormalize")
	if (!renorm) {
	    if (spmax < 1.0 && spmin > -1.0)
		renorm = TRUE
	    else
		renorm = FALSE
	    }

# Renormalize spectrum if requested
	if (renorm) {
	    sptot = 0.d0
	    nsp = 0
	    factor = 1000.0
	    do i = 1, npspec {
		speci = specwcs+i-1
		specd = Memr[speci]
		if (specd != 0.) {
		    sptot = sptot + specd
		    nsp = nsp + 1
		    }
		}
	    if (nsp > 0) {
		spmean = sptot / double (nsp)
		if (debug) {
		    call printf ("SHIFTPIX: mean count for %d pixels is %f\n")
			call pargi (nsp)
			call pargd (spmean)
		    call flush (STDOUT)
		    }
		if (spmean == 0.d0)
		    spmean = 1.d0
		if (spmin < 0.0) {
		    spfloor = factor * spmin / spmean
		    }
		else
		    spfloor = 0.0
		if (debug) {
		    call printf ("SHIFTPIX: minimum added count is %f\n")
			call pargd (spfloor)
		    call flush (STDOUT)
		    }
		do i = 1, specpix {
		    speci = spectrum+i-1
		    specd = Memr[speci]
		    if (specd != 0.) {
			Memr[speci] = factor * (specd / spmean)
			}
		    specd = Memr[speci] + spfloor
		    }
		}
	    else {
		call eprintf ("SHIFTPIX: Spectrum is all zeroes\n")
		go to endxc_
		}

#	If plot enabled, show the renormalized object spectrum.
	    if (pltspec)
		call plotspec (npspec, Memr[specwcs], specid,
			       Memr[wlwcs], xlab, nsmooth)
	    }
 
# Optional wavelength limits for pixel shift computation
	minwav0 = clgetd ("st_lambda")
	if (minwav0 == 0.d0)
	    minwav0 = dindef
	maxwav0 = clgetd ("end_lambda")
	if (maxwav0 == 0.d0)
	    maxwav0 = dindef

# Template against which to find pixel shift
	call clgstr ("template",tempfiles,SZ_PATHNAME)
	call clgstr ("tempdir",tempdir,SZ_PATHNAME)
	ldir = strlen (tempdir)
	if (ldir > 0 && tempdir[ldir] != '/' && tempdir[ldir] != '$') {
	    tempdir[ldir+1] = '/'
	    tempdir[ldir+2] = EOS
	    }

# Multispec template numbers (use only first if multiple files)
	tband = clgeti ("tempband")
	call clgstr ("tempnum",tempnums0,SZ_LINE)
	if (decode_ranges (tempnums0, tspec_range, MAX_RANGES, ntspec0) == ERR)
	    call error (1, "SHIFTPIX: Illegal multispec list")

# Optional integer pixel tshift between spectra and and template
	tshift = clgetd ("tshift")

# Allocate pixel shift and velocity vectors
	npts2 = npts * 2
	if (xcor == NULL) {
	    call countemp(tempfiles,tempdir,tempnums0,mspec,echelle,debug,ntmp)
	    ntpts = npts2 * ntmp
	    call malloc (xcor, ntpts, TY_REAL)
	    call malloc (xvel, ntpts, TY_REAL)
	    }

	if (debug) {
	    call printf ("SHIFTPIX: %d points in spectrum, up to %d points in transforms\n")
		call pargi (npts)
		call pargi (npts2)
	    call printf ("SHIFTPIX: %d templates being cross-correlated\n")
		call pargi (ntmp)
	    call flush (STDOUT)
	    }

# Optional object velocity to be used to shift templates
	czguess = 0.d0
	zinit = FALSE
	switch (vinit) {
	    case IVZERO:
		czguess = 0.d0
		zinit = FALSE
	    case IVGUESS:
		czguess = clgetd ("czguess")
		zinit = TRUE
	    case IZGUESS:
		czguess = c0 * clgetd ("czguess")
		zinit = TRUE
	    case IVXC:
		if (spxvel != dindef) {
		    czguess = spxvel
		    zinit = TRUE
		    }
	    case IVEM:
		if (spevel != dindef) {
		    czguess = spevel
		    zinit = TRUE
		    }
	    case IVCOMB:
		if (spvel != dindef) {
		    czguess = spvel
		    zinit = TRUE
		    }
	    default:
	    }

# Number of times to shift template
	nzpass = 0
	nzpass = clgeti ("nzpass")
	nzpass = nzpass + 1

# Allocate rebinned spectra and wavelength vectors
	if (shspec == NULL) {
	    call malloc (shspec, npts2, TY_REAL)
	    call malloc (shtemp, npts2, TY_REAL)
	    call malloc (wltemp, npts2, TY_REAL)
	    }

# If not computing shift for spectrum, skip to result plotting
	if (correlate == COR_NO) {
	    if (ntemp < 1)
		ntemp = 1
	    ntx = ntemp
	    if (itmax < 1) itmax = 1
	    minwav = W0(specsh)
	    if (minwav0 != dindef && minwav0 > minwav) minwav = minwav0
	    maxwav = W1(specsh)
	    if (maxwav0 != dindef && maxwav0 < maxwav) maxwav = maxwav0
	    do itemp = 1, ntemp {
		twl1[itemp] = minwav
		twl2[itemp] = maxwav
		}
	    go to endcorr_
	    }

# For each template, compute shift
refit_
	rmax = 0.d0
	itemp = 0
	spvqual = 0

	ntx = 0
	ntemp = 0
	templist = imtopen (tempfiles)

# Initialize new template file
newtemp_
	tspec = -1

# Get next template spectrum file name from the list
	if (imtgetim (templist, tempfile, SZ_PATHNAME) == EOF)
	    go to endspec_

# Check for specified section of template spectrum
	tempwl1[itemp+1] = dindef
	tempwl2[itemp+1] = dindef
	tempsec[1] = EOS
	ip = stridx (colon, tempfile)
	if (ip > 0) {
	    lfile = strlen (tempfile)
	    tempfile[ip] = EOS
	    ip = ip + 1
	    iw = stridxs (cdash, tempfile[ip])
	    if (iw > 0) {
		tempfile[ip+iw-1] = EOS
		call sscan (tempfile[ip])
		    call gargd (tempwl1[itemp+1])
		call sscan (tempfile[ip+iw])
		    call gargd (tempwl2[itemp+1])
		    call gargstr (tempsec, SZ_LINE)
		}
	    }

# Check for specified apertures in multispec template file
	ilb = stridxs (lbracket,tempfile)
	if (echelle)
	    ntspec = 1
	else if (ilb > 0) {
	    lfile = strlen (tempfile)
	    tempfile[ilb] = EOS
	    ilb = ilb + 1
	    ip = ilb
	    jp = 1
	    while (stridx (rbracket, tempfile[ip]) == 0 && ip <= lfile) {
		tempnums[jp] = tempfile[ip]
		tempfile[ip] = EOS
		ip = ip + 1
		jp = jp + 1
		}
	    tempnums[jp] = EOS
	    if (decode_ranges(tempnums,tspec_range,MAX_RANGES,ntspec)==ERR)
		call error (1, "SHIFTPIX: Illegal template multispec list")
	    }
	else {
	    call strcpy (tempnums0,tempnums,SZ_LINE)
	    if (decode_ranges (tempnums,tspec_range,MAX_RANGES,ntspec) == ERR)
		call error (1, "SHIFTPIX: Illegal template multispec list")
	    }
#	if (debug) {
#	    call printf ("SHIFTPIX: next template file is %s [%s] = %d aps\n")
#		call pargstr (tempfile)
#		call pargstr (tempnums)
#		call pargi (ntspec)
#	    }

# Check for readability of template spectrum
	if (stridx (cslash,tempfile) == 1 || stridx (cdollar,tempfile) > 0 ||
	    stridx (cdot,tempfile) == 1)
	    call strcpy (tempfile,temppath,SZ_PATHNAME)
	else {
	    call strcpy (tempdir,temppath,SZ_PATHNAME)
	    call strcat (tempfile,temppath,SZ_PATHNAME)
	    }
	if (imaccess (temppath, READ_ONLY) == NO) {
	    call eprintf ("SHIFTPIX: cannot read template file %s \n")
		call pargstr (temppath)
	    go to newtemp_
	    }

# Get next multispec aperture number from list
newtap_
	if (echelle)
	    tspec = mspec
	else if (get_next_number (tspec_range, tspec) == EOF)
	    go to newtemp_

# Start new template here if it has been cached
nextemp_
	itemp = itemp + 1

# Load template spectrum
	tempvel[itemp] = 0.d0
	tempshift[itemp] = 0.d0
	temphcv[itemp] = 0.d0
	call gettemp (tempfile,tspec,tband,tempdir,tempim,tempsh,itemp)
	tpixshift = wcs_getshift()
	if (tempim == ERR)
	    go to newtap_

# Set up identifying name for template
	if (tspec > 0) {
	    if (tempsec[1] != EOS) {
		call sprintf (tempid[1,itemp], SZ_PATHNAME,"%s:%s[%d]")
		    call pargstr (tempsec)
		    call pargstr (tempfile)
		    call pargi (tspec)
		}
	    else {
		call sprintf (tempid[1,itemp], SZ_PATHNAME,"%s[%d]")
		    call pargstr (tempfile)
		    call pargi (tspec)
		}
	    }
	else if (tempsec[1] != EOS) {
	    call sprintf (tempid[1,itemp], SZ_PATHNAME, "%s:%s")
		call pargstr (tempsec)
		call pargstr (tempfile)
	    }
	else {
	    call sprintf (tempid[1,itemp], SZ_PATHNAME, "%s")
		call pargstr (tempfile)
	    }
	ntx = ntx + 1
	ntemp = itemp

# Initialize shift results for this template
	tz = 0.d0
	cz[itemp] = 0.d0
	czerr[itemp] = 0.d0
	czr[itemp] = 0.d0
	zvel[itemp] = 0.d0
	nptemp = NP2(tempsh) - NP1(tempsh) + 1

# Set pixel limits of template WCS
	if (plttemp || pixfill) {
	    tempwcs = SY(tempsh) + NP1(tempsh) - 1
#	    twlwcs = SX(tempsh) + NP1(tempsh) - 1
	    call malloc (twlwcs, nptemp, TY_REAL)
	    do i = NP1(tempsh), NP2(tempsh), 1 {
		Memr[twlwcs+i-1] = real (wcs_p2w (double(i)))
		}
	    }
	if (plttemp) {
	    call plotspec (nptemp,Memr[tempwcs],tempid[1,itemp],
			   Memr[twlwcs],xlab,nsmooth)
	    }

# Spectrum emission line chopping flag
	call setschop (specim, tempim, itemp)

# Turn off template continuum subtraction, if not needed
	tcont = TRUE
	call imgbpar (tempim,"SUBCONT",tcont)
	if (debug) {
	    if (tcont)
		call printf ("SHIFTPIX: Template continuum will be subtracted\n")
	    else
		call printf ("SHIFTPIX: Template continuum already subtracted\n")
		call flush (STDOUT)
	    }
	tscont[itemp] = tcont

# Eliminate bad pixels from the template spectrum
	if (pixfill) {
	    if (echelle)
		iord = tspec
	    else
		iord = 0
	    call filllist (nptemp,Memr[tempwcs],Memr[twlwcs],iord,pixfill,debug)
	    if (plttemp) {
		call plotspec (nptemp,Memr[tempwcs],tempid[1,itemp],
			       Memr[twlwcs],xlab,nsmooth)
		}
	    }
	if (plttemp || pixfill) {
	    call mfree (twlwcs, TY_REAL)
	    }

# Remove continuum from spectrum and template
	divcon = FALSE
	call imgbpar (tempim,"DIVCONT",divcon)
	if (divcon)
	    conproc = ZEROCONT
	else
	    conproc = SUBCONT
	if (debug) {
	    if (divcon)
		call printf ("SHIFTPIX: Continuum divided, not subtracted\n")
	    else
		call printf ("SHIFTPIX: Continuum subtracted, not divided\n")
		call flush (STDOUT)
	    }
	tconproc[itemp] = conproc

# Ignore wavelength limit parameters and use file limits
	overlap = FALSE
	call imgbpar (tempim,"OVERLAP",overlap)
	if (minwav0 == dindef && maxwav0 == dindef)
	    overlap = FALSE
	if (debug) {
	    if (overlap)
		call printf ("SHIFTPIX: Use entire overlapping wavelength space\n")
	    else
		call printf ("SHIFTPIX: Use wavelength limit parameters, if set\n")
		call flush (STDOUT)
	    }
	toverlap[itemp] = overlap

# Zero-padding of transforms before shift computation
	zpad = FALSE
	switch (tzpad) {
	    case ZPAD:
		zpad = TRUE
	    case NOZPAD:
		zpad = FALSE
	    case TEMPFILE:
		call imgbpar (tempim,"ZEROPAD",zpad)
	    }

# Set redshift velocity limits
	minvel = clgetd ("minvel")
	maxvel = clgetd ("maxvel")
	if (maxvel != dindef && minvel != dindef && maxvel < minvel) {
	    tmpvel = minvel
	    minvel = maxvel
	    maxvel = tmpvel
	    }

# Add known velocity components as 1+z's
	if (correlate == COR_PIX) {
	    tz = 1.d0
	    ispix = TRUE
	    }
	else {
	    ispix = FALSE
	    tz = (1.d0 + (spechcv / c0)) *
	 	(1.d0 + (tempvel[itemp] / c0)) *
	 	(1.d0 + (tempshift[itemp] / c0)) *
	 	(1.d0 + (tshift / c0)) /
	 	(1.d0 + (temphcv[itemp] / c0))
	    if (debug) {
		call printf ("SHIFTPIX: spechcv= %.4f, tempvel= %.4f, temphcv= %.4f\n")
		    call pargd (spechcv)
		    call pargd (tempvel[itemp])
		    call pargd (temphcv[itemp])
		call printf ("SHIFTPIX: tshift= %.4f, tempshift= %.4f\n")
		    call pargd (tshift)
		    call pargd (tempshift[itemp])
		call printf ("SHIFTPIX: tz = %f\n")
		    call pargd (tz)
		call flush (STDOUT)
		}
	    }

# Set z offset for initial pass through loop
	if (zinit) {
	    zshift = (1.d0 + (czguess / c0)) / tz
	    if (debug) {
		call printf ("SHIFTPIX: zshift is %f\n")
		    call pargd (zshift)
		call flush (STDOUT)
		}
	    }
	else
	    zshift = 1.d0

# If shifting template on second pass, start here
	do izpass = 1, nzpass {

	    tmstrt = testrt[itemp]
	    tmfnsh = tefnsh[itemp]

	# Modify wavelength limits of template for velocity guess
	    if (zinit || izpass > 1) {
		tmshft = dlog10 (zshift)
		tmstrt = tmstrt + tmshft
		tmfnsh = tmfnsh + tmshft
		}
	    else
		tmshft = 0.d0

	    if (correlate == COR_PIX) {
		spwl1 = 1.d0
		spwl2 = double (npspec)
		}
	    else if (spstrt < spfnsh) {
		spwl1 = 10.d0 ** spstrt
		spwl2 = 10.d0 ** spfnsh
		}
	    else {
		spwl2 = 10.d0 ** spstrt
		spwl1 = 10.d0 ** spfnsh
		}
	    if (correlate == COR_PIX) {
		tmwl1 = 1.d0
		tmwl2 = double (nptemp)
		}
	    else if (tmstrt < tmfnsh) {
		tmwl1 = 10.d0 ** tmstrt
		tmwl2 = 10.d0 ** tmfnsh
		}
	    else {
		tmwl2 = 10.d0 ** tmstrt
		tmwl1 = 10.d0 ** tmfnsh
		}

	# Check to make sure that template and spectrum overlap
	    if (spwl1 > tmwl2) {
		call printf ("SHIFTPIX: Spectrum starts after template %s ends: %9.2f > %9.2f\n")
		    call pargstr (tempfile)
		    call pargd (spwl1)
		    call pargd (tmwl2)
		go to endtemp_
		}
	    if (spwl2 < tmwl1) {
		call printf ("SHIFTPIX: Spectrum ends before template %s starts %9.2f < %9.2f\n")
		    call pargstr (tempfile)
		    call pargd (spwl2)
		    call pargd (tmwl1)
		go to endtemp_
		}

	# Figure lambda and log-lambda limits and increment
	    minwav1 = max (spwl1, tmwl1)
	    maxwav1 = min (spwl2, tmwl2)
	    if (tempwl1[itemp] != dindef)
		minwav = tempwl1[itemp]
	    else if (overlap)
		minwav = minwav1
	    else if ((minwav0==dindef)||(minwav0<minwav1))
		minwav = minwav1
	    else
		minwav = minwav0

	    if (tempwl2[itemp] != dindef)
		maxwav = tempwl2[itemp]
	    else if (overlap)
		maxwav = maxwav1
	    else if ((maxwav0==dindef)||(maxwav0>maxwav1))
		maxwav = maxwav1
	    else
		maxwav = maxwav0

	    if (correlate == COR_PIX) {
		call wcs_set (specsh)
		call wcs_pixshift (pixshift)
		minwav = double (np1)
		maxwav = double (np2)
	 	}

	    twl1[itemp] = minwav
	    twl2[itemp] = maxwav

	# Fill in globally available values of wavelength
	    logmin = dlog10 (minwav)
	    logmax = dlog10 (maxwav)
	    logw0 = logmin
	    dlogw = (logmax - logmin) / double (npts - 1)
	    wave0 = minwav
	    delwav = (maxwav - minwav) / double (npts - 1)

	    if (debug) {
		call printf("SHIFTPIX: overlap %9.3f to %9.3f by %9.3f, npts = %d \n")
		    call pargd (minwav)
		    call pargd (maxwav)
		if (correlate == COR_VEL)
		    call pargd (10.d0 ** dlogw)
		else
		    call pargd (delwav)
		    call pargi (npts)
		call flush (STDOUT)
		}

	#  Rebin spectrum to overlap this template
	    if (debug) {
		call printf("SHIFTPIX: spectrum %9.3f to %9.3f = %d %d -> %d\n")
		    call pargd (spwl1)
		    call pargd (spwl2)
		    call pargi (SN(specsh))
		    call pargi (DC(specsh))
		    call pargi (npts)
		call flush (STDOUT)
		}
	    spshft = 0.d0
	    if (correlate == COR_WAV || correlate == COR_PIX) {
		if (debug) {
		    call printf ("SHIFTPIX: Rebin to %d pixels from %10.5f by %10.5f\n")
			call pargi (npts)
			call pargd (wave0)
			call pargd (delwav)
			}
		call rebin (Memr[spectrum],specsh, spshft,
			 Memr[shspec], npts, wave0, delwav, pixshift, ispix)
		}
	    else {
		if (debug) {
		    call printf ("SHIFTPIX: Rebin to %d pixels from %10.5f by %10.5f\n")
			call pargi (npts)
			call pargd (logw0)
			call pargr (dlogw)
			}
		call rebinl (Memr[spectrum],specsh, spshft,
			 Memr[shspec], npts, logw0, dlogw, pixshift)
		}
	    if (debug) {
		call printf ("SHIFTPIX: Spectrum before %9.3f %9.3f %9.3f\n")
		    call pargr (Memr[spectrum])
		    call pargr (Memr[spectrum+1])
		    call pargr (Memr[spectrum+2])
		call printf ("SHIFTPIX: Spectrum after  %9.3f %9.3f %9.3f\n")
		    call pargr (Memr[shspec])
		    call pargr (Memr[shspec+1])
		    call pargr (Memr[shspec+2])
		call flush (STDOUT)
		}

	#  Rebin template to overlap this spectrum
	    if (debug) {
		call printf("SHIFTPIX: template %9.3f to %9.3f = %d %d -> %d\n")
		    call pargd (tmwl1)
		    call pargd (tmwl2)
		    call pargi (SN(tempsh))
		    call pargi (DC(tempsh))
		    call pargi (npts)
		call flush (STDOUT)
		}
	    if (correlate == COR_WAV || correlate == COR_PIX)
		call rebin (Memr[SY(tempsh)], tempsh, tmshft,
			     Memr[shtemp],npts,wave0,delwav,tpixshift,ispix)
	    else
		call rebinl (Memr[SY(tempsh)], tempsh, tmshft,
			     Memr[shtemp], npts, logw0, dlogw, tpixshift)
	    if (debug) {
		call printf ("SHIFTPIX: Template before %9.3f %9.3f %9.3f\n")
		    call pargr (Memr[SY(tempsh)])
		    call pargr (Memr[SY(tempsh)+1])
		    call pargr (Memr[SY(tempsh)+2])
		call printf ("SHIFTPIX: Template after  %9.3f %9.3f %9.3f\n")
		    call pargr (Memr[shtemp])
		    call pargr (Memr[shtemp+1])
		    call pargr (Memr[shtemp+2])
		call flush (STDOUT)
		}

	# Set up wavelength vector for intermediate graphs
	    if (correlate == COR_PIX || correlate == COR_WAV) {
		do j = 1, npts {
		    Memr[wltemp+j-1] = wave0 + ((j-1) * delwav)
		    }
		}
	    else {
		do j = 1, npts {
		    Memr[wltemp+j-1] = 10 ** (logw0 + (j-1) * dlogw)
		    }
		}

	# Set known portion of velocity
	    if (zinit || izpass > 1)
		tvel = ((tz * zshift) - 1.d0) * c0
	    else
		tvel = (tz - 1.d0) * c0


	    if (debug) {
		call printf ("SHIFTPIX: tvel = %f\n")
		    call pargd (tvel)
		call printf ("SHIFTPIX: filter is %d %d %d %d\n")
		    call pargi (lo)
		    call pargi (toplo)
		    call pargi (topnrn)
		    call pargi (nrun)
		call flush (STDOUT)
		}

	# Cross-correlate this template/spectrum combination
	    ixcor = xcor + ((itemp - 1) * npts2)
	    ixvel = xvel + ((itemp - 1) * npts2)
	    call xcorfit (Memr[shspec],Memr[shtemp],Memr[wltemp],
			  ncor, Memr[ixcor],Memr[ixvel],itemp,
			  tcz,tczerr,tczr)

	# Save velocity and R-value
	    cz[itemp] = tcz
	    czerr[itemp] = tczerr
	    czr[itemp] = tczr
	    if (correlate == COR_PIX || correlate == COR_WAV) {
		zvel[itemp] = tcz
		if (zinit || izpass > 1) {
		    zvel[itemp] = zvel[itemp] + zshift
		    tvshift[itemp] = zshift
		    }
		else
		    tvshift[itemp] = 0.d0
		}
	    else {
		zvel[itemp] = (1.d0 + (tcz / c0)) * tz
		if (zinit || izpass > 1) {
		    zvel[itemp] = zvel[itemp] * zshift
		    tvshift[itemp] = (zshift - 1.d0) * c0
		    }
		else
		    tvshift[itemp] = 0.d0
		zvel[itemp] = (zvel[itemp] - 1.d0) * c0
		}

	    if (nzpass > 1 && izpass < nzpass) { 
		if (correlate == COR_PIX || correlate == COR_WAV)
		    zshift = zshift + tcz
		else
		    zshift = zshift * (1.d0 + (tcz / c0))
		if (debug) { 
		    call printf ("SHIFTPIX: pass %d: v= %9.3f, 1+z= %7.5f")
			call pargi (izpass)
			call pargd (tcz)
			call pargd (zshift)
		    call printf (" r= %6.3f h= %6.4f arms= %8.6f\n")
			call pargd (tczr)
			call pargd (thght[itemp])
			call pargd (tarms[itemp])
		    call flush (STDOUT)
		    }
		}
	    }
	# End of cz iteration

	if (debug) {
      	    call printf("SHIFTPIX: End of cz iteration for template %d, %d points\n")
		call pargi (itemp)
		call pargi (ncor)
	    call flush (STDOUT)
	    }

# Save velocity and cross-correlation vectors for summary plot
	ncvel[itemp] = ncor
	if (ncor < npts2) {
	    do i = ncor+1, npts2 {
		Memr[ixcor+i] = 0.0
		Memr[ixvel+i] = 0.0
		}
	    }
	if (debug) {
      	    call printf("SHIFTPIX: Cross-correlation saved\n")
	    call flush (STDOUT)
	    }

# if R is new max, set up velocity vector for correlation plot
	if (tczr >= rmax) {
	    rmax = tczr
	    npmax = tnpfit[itemp]
	    xcrmax = zvel[itemp]
	    itmax = itemp
	    if (debug) {
      		call printf("SHIFTPIX: Best correlation template: %s %5.2f\n")
		    call pargstr(tempname[1,itemp])
		    call pargd (rmax)
		call flush (STDOUT)
		}
	    }
	newresults = TRUE

# Print the filtered cross-correlation function if requested
	filexcor = FALSE
	filexcor = clgetb ("xcor_file")
	if (filexcor)
            call xcfile (specfile,specim,mspec,tempfile,tspec,itemp,ncor,
			 Memr[ixvel],Memr[ixcor],txfile[1,itemp])
 
# Close the template image
endtemp_
	call close_image (tempim,tempsh)
	if (debug) {
      	    call printf("SHIFTPIX: Template %d closed\n")
		call pargi (itemp)
	    call flush (STDOUT)
	    }

	ntspec = ntspec - 1
	if (ntspec > 0)
	    go to newtap_
	go to newtemp_

# Close the template list for this object spectrum
endspec_
	call imtclose (templist)
endcorr_

# Save results for this object spectrum if at least one template correlated
	IM_UPDATE(specim) = NO
	if (ntx > 0) {
	    spxvel = zvel[itmax]
	    spxerr = czerr[itmax]
	    spxr = czr[itmax]

	# Combine correlation and emission velocities for this object spectrum
	    if (newresults) {
	        call vcombine (spxvel,spxerr,spxr,spevel,speerr,spnlf,
			       spvel,sperr,debug)

	# Save correlation information to archive file for this object
		arcflag = FALSE
		arcflag = clgetb ("archive")
		if (arcflag)
		    call xcarch (specpath)
		}

	    rmode = 2
	    rmode = clgeti ("report_mode")
	    call xcrslts (specfile,mspec,specim,ncor,rmode,filexcor,txfile)

	# Plot results
#	    if (debug) {
#		ixvel = xvel + ((itmax - 1) * npts2)
#		call printf ("SHIFTPIX:  template %d: vel= %.3f - %.3f\n")
#		    call pargi (itmax)
#		    call pargr (Memr[ixvel])
#		    call pargr (Memr[ixvel+ncvel[itmax]-1])
#		ixcor = xcor + ((itmax - 1) * npts2)
#		call printf ("SHIFTPIX:  template %d: xcor= %.3f - %.3f\n")
#		    call pargi (itmax)
#		    call pargr (Memr[ixcor])
#		    call pargr (Memr[ixcor+ncvel[itmax]-1])
#		}
	    itm0 = itmax
	    call xcplot (Memr[wlwcs],Memr[specwcs],npspec,Memr[xvel],
			 Memr[xcor],ncvel,npts2,specfile,mspec,specim,newfit)
	    if (itm0 != itmax)
		newresults = TRUE
	    if (newfit) {
		newfit = FALSE
		if (waverest != 0.d0) {
		    czguess = spvel
		    zinit = TRUE
		    }
		go to refit_
		}

	# Save best Cz, error, and R to image header, if requested
	# after checking for writability
	    if ((correlate != COR_NO && savevel) || qplot) {
		if (imaccess (specpath, READ_WRITE) == NO) {
		call eprintf ("XCSAO: cannot write to %s; not saving results\n")
		    call pargstr (specpath)
		    }
		else if (savevel) {
		    if (debug) call printf ("SHIFTPIX: saving velocity\n")
		    call vwhead (mspec,specim)
		    if (debug) call printf ("SHIFTPIX: saving correlation\n")
		    call xcwhead (mspec,specim)
		    if (debug) call printf ("SHIFTPIX: results saved\n")
		    }
	        else if (qplot && newresults) {
		    if (debug) call printf ("SHIFTPIX: saving quality and results\n")
		    call qwhead (mspec,specim,"xcsao")
		    call xcwhead (mspec,specim)
		    }
	        else if (qplot) {
		    if (debug) call printf ("SHIFTPIX: saving quality\n")
		    call qwhead (mspec,specim,"xcsao")
		    }
		}
	    }

	else {
	    call printf ("SHIFTPIX: No templates; no cross-correlation\n")
	    }


# Close the object spectrum image
endxc_
	if (debug) call printf ("SHIFTPIX: closing image\n")
	if (!wenable)
	    IM_UPDATE(specim) = NO
	call close_image (specim, specsh)
	call mfree (wlwcs, TY_REAL)
	return

end

# Jun  9 2008	New subroutine based on xcfit.x
