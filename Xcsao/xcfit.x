# File rvsao/Xcsao/xcfit.x
# May 1, 2014
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss, Johns Hopkins University and others

# Copyright(c) 1994-2014 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# XCFIT finds the best template and velocity for one object spectrum
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include <smw.h>
include "rvsao.h"
include	"xcv.h"
include	"contin.h"

procedure xcfit (specfile,specdir,mspec,mband, oshift)

char	specfile[ARB]	# Object spectrum file name
char	specdir[ARB]	# Directory for object spectra
int	mspec		# Object aperture to read from multispec file
int	mband		# Object band to read from multispec file
int	oshift		# Shift between echelle orders cross-correlated

pointer spectrum	# Object spectrum
pointer specim		# Object image header structure

bool	arcflag		# archive record writing flag
bool	newfit		# true to redo cross-correlation
bool	ispix		# TRUE if cross-correlating in pixel space, else FALSE

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
char	zero_pad[SZ_LINE]	# Zero-pad before correlating or from file
				# (no | yes | tempfile)

pointer	specwcs			# Pointer to first data point with WCS
pointer	wlwcs			# Pointer to first wavelength with WCS
int	npspec			# Number of spectrum points with WCS
pointer	tempwcs			# Pointer to first template data point with WCS
pointer	twlwcs			# Pointer to first template wavelength with WCS
int	nptemp			# Number of template points with WCS

pointer	ixcor		# cross-correlation returned from xcorfit
pointer	ixvel		# cross-correlation velocities from xcorfit
int	ncvel[MAXTEMPS]	# Number of points in each correlation

double	zshift		# guessed Z for template offset
double	czguess		# initial velocity or Z for template offset
double	tmstrt,tmfnsh	# template limits after zguess applied
int	izpass		# Count for template shifts
int	nzpass		# Number of template shifts
bool	zinit		# TRUE if initial value for z is set
int	tspec		# Template spectrum to read from multispec file
int	tzpad		# Flag for zero-padding for one template
int	ntpts		# ntemp * npts2
bool	divcon		# TRUE if template header says to divide by continuum

pointer	templist	# List of template files

double	tcz		# Velocity returned from xcorfit
double	tczerr		# Velocity error returned from xcorfit
double	tczr		# R returned from xcorfit
double	rmax		# Maximum R value for spectrum and best template
double 	minwav0,maxwav0	# Starting,ending wavelengths from parameter file
double 	minwav1,maxwav1	# Overlap wavelength limits from spectrum and template
double 	minwav,maxwav	# Wavelength limits for cross-correlation
double 	logmin,logmax	# Log-10 of wavelength limits for cross-correlation
bool	filexcor	# Flag to save correlation to a file
double	spwl1,spwl2	# Spectrum first and last wavelengths in angstroms
double	tmwl1,tmwl2	# Template first and last wavelengths in angstroms
int	ntx		# Number of templates actually cross-correlated
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
	wlwcs = NULL
	twlwcs = NULL

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
	    call eprintf ("XCFIT: cannot read spectrum path %s \n")
		call pargstr (specpath)
	    return
	    }

	call getspec (specpath, mspec, mband, spectrum, specim, wenable)
	if (specim == ERR) {
	#   call eprintf ("XCFIT: Error reading spectrum %s\n")
	# 	call pargstr (specpath)
	    wenable = FALSE
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
	    call printf ("XCFIT: %d-pixel spectrum from %.3fA(%d) - %.3fA(%d)\n")
		call pargi (specpix)
		call pargd (wcs_p2w (double (np1)))
		call pargi (np1)
		call pargr (Memr[SX(specsh)+NP2(specsh)-1])
		call pargd (wcs_p2w (double (np2)))
	    call flush (STDOUT)
	    }

# Turn off spectrum continuum subtraction, if not needed
	tcont = TRUE
	call imgbpar (specim,"SUBCONT",tcont)
	divcon = FALSE
	call imgbpar (specim,"DIVCONT",divcon)
	if (tcont) {
	    if (divcon) {
		contfit = TRUE
		conproc = ZEROCONT
		}
	    else {
		contfit = TRUE
		conproc = SUBCONT
		}
	    }
	else {
	    contfit = FALSE
	    conproc = NOCONT
	    }
#	if (debug) {
#	    if (contfit)
#		call printf ("XCFIT: Spectrum continuum will be subtracted\n")
#	    else
#		call printf ("XCFIT: Spectrum continuum already subtracted\n")
#	    }

# Eliminate zeroes at left end of spectrum
	while (Memr[spectrum+np1-1] == 0.0) {
	    np1 = np1 + 1
	    }
#	if (debug) {
#	    call printf ("XCFIT: %d pixels were zero at blue end of spectrum\n")
#		call pargd (np1 - NP1(specsh))
#	    call flush (STDOUT)
#	    }

# Eliminate zeroes at right end of spectrum
	while (Memr[spectrum+np2-1] == 0.0)
	    np2 = np2 - 1
#	if (debug) {
#	    call printf ("XCFIT: %d pixels were zero at red end of spectrum\n")
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
	    call printf ("XCFIT: %8g < counts < %8g from %.3fA(%d) - %.3fA(%d)\n")
		call pargd (spmin)
		call pargd (spmax)
		call pargr (Memr[wlwcs+np1-1])
		call pargi (np1)
		call pargr (Memr[wlwcs+np2-1])
		call pargi (np2)
	    call flush (STDOUT)
	    }
	if (debug) {
	    call printf ("XCFIT: maximum count in spectrum is %f\n")
		call pargd (spmax)
	    call printf ("         minimum count in spectrum is %f\n")
		call pargd (spmin)
	    call flush (STDOUT)
	    }

# Renormalize spectrum if requested
	renorm = clgetb ("renormalize")

# If RENORM is not set, set it if maximum counts are less than 1
	if (!renorm) {
	    if (spmax < 1.0 && spmin > -1.0)
		renorm = TRUE
	    else
		renorm = FALSE
	    }

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
		    call printf ("XCFIT: mean count for %d pixels is %f\n")
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
		    call printf ("XCFIT: minimum added count is %f\n")
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
		call eprintf ("XCFIT: Spectrum is all zeroes\n")
		go to endxc_
		}

#	If plot enabled, show the renormalized object spectrum.
	    if (pltspec)
		call plotspec (npspec, Memr[specwcs], specid,
			       Memr[wlwcs], xlab, nsmooth)
	    }
 
# Optional wavelength limits for cross-correlation
	minwav0 = clgetd ("st_lambda")
	if (minwav0 == 0.d0)
	    minwav0 = dindef
	maxwav0 = clgetd ("end_lambda")
	if (maxwav0 == 0.d0)
	    maxwav0 = dindef

# Templates against which to correlate spectra
	call clgstr ("templates",tempfiles,SZ_PATHNAME)
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
	    call error (1, "XCFIT: Illegal multispec list")

# Zero padding on, off, or read from template file
	call clgstr ("zeropad",zero_pad,SZ_LINE)
	tzpad = strdic (zero_pad,zero_pad,SZ_LINE, XC_ZTYPES)

# Fraction of spectrum to apodize at each end
	han = clgetr ("bell_window")

# Filter type and filter constants [always cosine-bell]
	filter = COSBELL
	lo0     = clgeti ("low_bin")
	toplo0  = clgeti ("top_low")
	nrun0   = clgeti ("nrun")
	tpnrn0 = clgeti ("top_nrun")
	if (debug) {
	    call printf ("XCFIT: filter is %d %d %d %d\n")
		call pargi (lo0)
		call pargi (toplo0)
		call pargi (tpnrn0)
		call pargi (nrun0)
	    }

# Optional tshift to be used for all spectra and templates
# tshift can also be read from the TSHIFT parameter in a template header
	tshift = clgetd ("tshift")

# If npts not a power of two, determine next largest power of two:
	npts = clgeti ("ncols")
	if (npts == INDEF)
	    npts = alog (real (specpix)) / LN2 +  1.
	else
	    npts = alog (real (npts)) / LN2 +  1.
	npts = 2 ** npts

# Check for unreasonable filtering
	if (lo0 < 1) lo0 = 1
	if (lo0 > npts) lo0 = 1
	if (toplo0  < 1) toplo0 = 1
	if (toplo0 > npts) toplo0 = 1
	if (nrun0 < 1) nrun0 = npts
	if (nrun0 > npts) nrun0 = npts
	if (tpnrn0 < 1) tpnrn0 = npts
	if (tpnrn0 > npts) tpnrn0 = npts

# Allocate correlation and velocity vectors
	npts2 = npts * 2
	if (xcor == NULL) {
	    call countemp(tempfiles,tempdir,tempnums0,mspec,echelle,debug,ntmp)
	    ntpts = npts2 * ntmp
	    call malloc (xcor, ntpts, TY_REAL)
	    if (xvel != NULL)
		call mfree (xvel, TY_REAL)
	    call malloc (xvel, ntpts, TY_REAL)
	    }

	if (debug) {
	    call printf ("XCFIT: %d points in spectrum, up to %d points in transforms\n")
		call pargi (npts)
		call pargi (npts2)
	    call printf ("XCFIT: %d templates being cross-correlated\n")
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
	    if (shtemp != NULL)
		call mfree (shtemp, TY_REAL)
	    call malloc (shtemp, npts2, TY_REAL)
	    if (wltemp != NULL)
		call mfree (wltemp, TY_REAL)
	    call malloc (wltemp, npts2, TY_REAL)
	    }

# If not cross-correlating spectrum skip to result plotting
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

# For each template, compute cross-correlation
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
	if (echelle) {
	    ntspec = 1
	    }
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
		call error (1, "XCFIT: Illegal template multispec list")
	    }
	else {
	    call strcpy (tempnums0,tempnums,SZ_LINE)
	    if (decode_ranges (tempnums,tspec_range,MAX_RANGES,ntspec) == ERR)
		call error (1, "XCFIT: Illegal template multispec list")
	    }
#	if (debug) {
#	    call printf ("XCFIT: next template file is %s [%s] = %d aps\n")
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
	    call eprintf ("XCFIT: cannot read template file %s \n")
		call pargstr (temppath)
	    go to newtemp_
	    }

# Get next multispec aperture number from list
newtap_
	if (echelle) {
	    tspec = mspec + oshift
	    if (debug) {
		call printf ("XCFIT: Spectrum order %d, template order %d (%d)\n")
		    call pargi (mspec)
		    call pargi (tspec)
		    call pargi (oshift)
		}
	    }
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

# Initialize cross-correlation results for this template
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

# Template emission line chopping flag
	call settchop (specim, tempim)

# Turn off template continuum subtraction, if not needed
	tcont = TRUE
	call imgbpar (tempim,"SUBCONT",tcont)
	if (debug) {
	    if (tcont)
		call printf ("XCFIT: Template continuum will be subtracted\n")
	    else
		call printf ("XCFIT: Template continuum already subtracted\n")
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
	if (twlwcs != NULL) {
	    call mfree (twlwcs, TY_REAL)
	    twlwcs = NULL
	    }

# Remove continuum from spectrum and template
	if (!divcon) {
	    call imgbpar (tempim,"DIVCONT",divcon)
	    if (divcon)
		conproc = ZEROCONT
	    else
		conproc = SUBCONT
	    }
	if (debug) {
	    if (divcon) {
		if (tcont)
		    call printf ("XCFIT: Object + template continua divided, not subtracted\n")
		else
		    call printf ("XCFIT: Object continuum divided, not subtracted\n")
		}
	    else if (tcont)
		call printf ("XCFIT: Object +template continua subtracted, not divided\n")
	    else
		call printf ("XCFIT: Object continuum subtracted, not divided\n")
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
		call printf ("XCFIT: Use entire overlapping wavelength space\n")
	    else
		call printf ("XCFIT: Use wavelength limit parameters, if set\n")
		call flush (STDOUT)
	    }
	toverlap[itemp] = overlap

# Zero-padding of transforms before cross-correlation
	zpad = FALSE
	switch (tzpad) {
	    case ZPAD:
		zpad = TRUE
	    case NOZPAD:
		zpad = FALSE
	    case TEMPFILE:
		call imgbpar (tempim,"ZEROPAD",zpad)
	    }

# Set number of points in correlation depending on whether zero-padding
	    if (zpad) {
		ncor = 2 * npts
		if (debug)
		    call printf ("XCFIT: Correlation is zero-padded\n")
		call flush (STDOUT)
		if (lo0 > 1)
		    lo = 2 * lo0
		else
		    lo = lo0
		if (toplo0 > 1)
		    toplo = 2 * toplo0
		else
		    toplo = toplo0
		topnrn = 2 * tpnrn0
		if (topnrn > ncor)
		    topnrn = ncor
		nrun = 2 * nrun0
		if (nrun > ncor)
		    nrun = ncor
		}
	    else {
		ncor = npts
		lo = lo0
		toplo = toplo0
		topnrn = tpnrn0
		nrun = nrun0
		if (debug)
		    call printf ("XCFIT: Correlation is not zero-padded\n")
		call flush (STDOUT)
		}

# Set fraction of correlation peak to fit
	pkfrac = clgetd ("pkfrac")
	call imgdpar (tempim,"PEAKFRAC",pkfrac)

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
		call printf ("XCFIT: spechcv= %.4f, tempvel= %.4f, temphcv= %.4f\n")
		    call pargd (spechcv)
		    call pargd (tempvel[itemp])
		    call pargd (temphcv[itemp])
		call printf ("XCFIT: tshift= %.4f, tempshift= %.4f\n")
		    call pargd (tshift)
		    call pargd (tempshift[itemp])
		call printf ("XCFIT: tz = %f\n")
		    call pargd (tz)
		call flush (STDOUT)
		}
	    }

# Set z offset for initial pass through loop
	if (zinit) {
	    zshift = (1.d0 + (czguess / c0)) / tz
	    if (debug) {
		call printf ("XCFIT: zshift is %f\n")
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
		call printf ("XCFIT: Spectrum starts after template %s ends: %9.2f > %9.2f\n")
		    call pargstr (tempfile)
		    call pargd (spwl1)
		    call pargd (tmwl2)
		go to endtemp_
		}
	    if (spwl2 < tmwl1) {
		call printf ("XCFIT: Spectrum ends before template %s starts %9.2f < %9.2f\n")
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
		call printf("XCFIT: overlap %9.3f to %9.3f by %9.3f (%9.6f), npts = %d \n")
		    call pargd (minwav)
		    call pargd (maxwav)
		    call pargd (delwav)
		    call pargd (dlogw)
		    call pargi (npts)
		call flush (STDOUT)
		}

	#  Rebin spectrum to overlap this template
	    if (debug) {
		call printf("XCFIT: spectrum %9.3f to %9.3f = %d %d -> %d\n")
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
		    call printf ("XCFIT: Rebin to %d pixels from %10.5f by %10.5f\n")
			call pargi (npts)
			call pargd (wave0)
			call pargd (delwav)
			}
		call rebin (Memr[spectrum],specsh, spshft,
			 Memr[shspec], npts, wave0, delwav, pixshift, ispix)
		}
	    else {
		if (debug) {
		    call printf ("XCFIT: Rebin to %d pixels from %10.6f to %10.6f by %10.6g (%10.6g)\n")
			call pargi (npts)
			call pargd (logw0)
			call pargd (logmax)
			call pargd (dlogw)
			call pargd (pixshift)
		    }
		call rebinl (Memr[spectrum],specsh, spshft,
			 Memr[shspec], npts, logw0, dlogw, pixshift)
		}
	    if (debug) {
		call printf ("XCFIT: Spectrum before %9.3f %9.3f %9.3f\n")
		    call pargr (Memr[spectrum])
		    call pargr (Memr[spectrum+1])
		    call pargr (Memr[spectrum+2])
		call printf ("XCFIT: Spectrum after  %9.3f %9.3f %9.3f\n")
		    call pargr (Memr[shspec])
		    call pargr (Memr[shspec+1])
		    call pargr (Memr[shspec+2])
		call flush (STDOUT)
		}

	#  Rebin template to overlap this spectrum
	    if (debug) {
		call printf("XCFIT: template %9.3f to %9.3f = %d %d -> %d\n")
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
		call printf ("XCFIT: Template before %9.3f %9.3f %9.3f\n")
		    call pargr (Memr[SY(tempsh)])
		    call pargr (Memr[SY(tempsh)+1])
		    call pargr (Memr[SY(tempsh)+2])
		call printf ("XCFIT: Template after  %9.3f %9.3f %9.3f\n")
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
		call printf ("XCFIT: tvel = %f\n")
		    call pargd (tvel)
		call printf ("XCFIT: filter is %d %d %d %d\n")
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
		    call printf ("XCFIT: pass %d: v= %9.3f, 1+z= %7.5f")
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
      	    call printf("XCFIT: End of cz iteration for template %d, %d points\n")
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
      	    call printf("XCFIT: Cross-correlation saved\n")
	    call flush (STDOUT)
	    }

# if R is new max, set up velocity vector for correlation plot
	if (tczr >= rmax) {
	    rmax = tczr
	    npmax = tnpfit[itemp]
	    xcrmax = zvel[itemp]
	    itmax = itemp
	    if (debug) {
      		call printf("XCFIT: Best correlation template: %s %5.2f\n")
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
      	    call printf("XCFIT: Template %d closed\n")
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
#		call printf ("XCFIT:  template %d: vel= %.3f - %.3f\n")
#		    call pargi (itmax)
#		    call pargr (Memr[ixvel])
#		    call pargr (Memr[ixvel+ncvel[itmax]-1])
#		ixcor = xcor + ((itmax - 1) * npts2)
#		call printf ("XCFIT:  template %d: xcor= %.3f - %.3f\n")
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
		else {
		    if (savevel) {
			if (debug) call printf ("XCFIT: saving velocity\n")
			call vwhead (mspec,specim)
			if (debug) call printf ("XCFIT: saving correlation\n")
			call xcwhead (mspec,specim)
			if (debug) call printf ("XCFIT: results saved\n")
			}
		    else if (qplot && newresults) {
			if (debug) call printf ("XCFIT: saving quality and results\n")
			call qwhead (mspec,specim,"xcsao")
			call xcwhead (mspec,specim)
			}
		    else if (qplot) {
			if (debug) call printf ("XCFIT: saving quality\n")
			call qwhead (mspec,specim,"xcsao")
			}
		    if (savespec) {
			call putspec (specim, mspec, spectrum, TY_REAL, debug)
			savespec = FALSE
			}
		    }
		}
	    }

	else {
	    call printf ("XCFIT: No templates; no cross-correlation\n")
	    }


# Close the object spectrum image
endxc_
	if (specim == ERR)
	    specim = NULL
	if (specim != NULL) {
	    if (!wenable)
		IM_UPDATE(specim) = NO
	    call close_image (specim, specsh)
	    }
	if (wlwcs != NULL) {
	    call mfree (wlwcs, TY_REAL)
	    wlwcs = NULL
	    }
	return
end


procedure countemp (tempfiles, tempdir, tempnums0, mspec, echelle, debug, ntemp)

char	tempfiles[ARB]		# List of template spectra
char	tempdir[ARB]		# Directory for template spectra
char	tempnums0[ARB]		# List of multispec spectra to read
int	mspec			# Object aperture to read from multispec file
bool	echelle			# If true, template apertures track object
bool	debug
int	ntemp			# Number of templates (returned)

char	tempfile[SZ_PATHNAME]	# Template spectrum file name
char	temppath[SZ_PATHNAME]	# Template spectrum path name
char	tempnums[SZ_LINE]	# List of multispec spectra to read
pointer templist		# List of template files
int     ntspec,ntspec0		# Number of template multispec spectra
int     tspec_range[3,MAX_RANGES]
int	tspec			# Template spectrum to read from multispec file
int	ip, jp, lfile
char	lbracket[3]	# left bracket, parenthesis, brace
char	rbracket[3]	# right bracket, parenthesis, brace
char	cdot, colon

int	stridx(), stridxs()
int	decode_ranges(),get_next_number()
int	imtgetim(), imaccess(), strlen()
pointer	imtopen()

define  newtemp_ 40
define  newtap_  50
define  endspec_ 60

begin
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	cdot = '.'
	colon = ':'

	templist = imtopen (tempfiles)
	if (decode_ranges (tempnums0, tspec_range, MAX_RANGES, ntspec0) == ERR)
	    call error (1, "COUNTEMP: Illegal multispec list")
	ntemp = 0

# Initialize new template file
newtemp_
	tspec = -1

# Get next template spectrum file name from the list
	if (imtgetim (templist, tempfile, SZ_PATHNAME) == EOF)
	    go to endspec_

# Check for specified section of template spectrum
	ip = stridx (colon, tempfile)
	if (ip > 0) {
	    lfile = strlen (tempfile)
	    tempfile[ip] = EOS
	    do jp = ip, lfile {
		tempfile[ip] = EOS
		}
	    }

# Check for specified apertures in multispec template file
	ip = stridxs (lbracket,tempfile)
	if (echelle)
	    ntspec = 1
	else if (ip > 0) {
	    lfile = strlen (tempfile)
	    tempfile[ip] = EOS
	    ip = ip + 1
	    jp = 1
	    while (stridx (tempfile[ip],rbracket) == 0 && ip <= lfile) {
		tempnums[jp] = tempfile[ip]
		tempfile[ip] = EOS
		ip = ip + 1
		jp = jp + 1
		}
	    tempnums[jp] = EOS
	    if (decode_ranges (tempnums,tspec_range,MAX_RANGES,ntspec) == ERR)
		call error (1, "COUNTEMP: Illegal template multispec list")
	    }
	else {
	    call strcpy (tempnums0,tempnums,SZ_LINE)
	    if (decode_ranges (tempnums,tspec_range,MAX_RANGES,ntspec) == ERR)
		call error (1, "COUNTEMP: Illegal template multispec list")
	    }
	if (debug) {
	    call printf ("COUNTEMP: template file %d is %s [%s] = %d aps\n")
		call pargi (ntemp+1)
		call pargstr (tempfile)
		call pargstr (tempnums)
		call pargi (ntspec)
	    call flush (STDOUT)
	    }

# Check for readability of template spectrum
	if (stridx ("/",tempfile) == 1 || stridx ("$",tempfile) > 0 ||
	    stridx (cdot,tempfile) == 1)
	    call strcpy (tempfile,temppath,SZ_PATHNAME)
	else {
	    call strcpy (tempdir,temppath,SZ_PATHNAME)
	    call strcat (tempfile,temppath,SZ_PATHNAME)
	    }
	if (imaccess (temppath, READ_ONLY) == NO) {
	    call eprintf ("COUNTEMP: cannot read template file %s \n")
		call pargstr (temppath)
	    go to newtemp_
	    }
	if (debug) {
	    call printf ("COUNTEMP: template path %d is %s\n")
		call pargi (ntemp+1)
		call pargstr (temppath)
	    call flush (STDOUT)
	    }

# Get next multispec aperture number from list
newtap_
	if (echelle)
	    tspec = mspec
	else if (get_next_number (tspec_range, tspec) == EOF)
	    go to newtemp_

# Add this template to counter
	ntemp = ntemp + 1
	if (debug) {
	    call printf ("COUNTEMP: %d templates\n")
		call pargi (ntemp)
	    call flush (STDOUT)
	    }

	ntspec = ntspec - 1
	if (ntspec > 0)
	    go to newtap_
	go to newtemp_

endspec_
	call imtclose (templist)
	if (debug) {
	    call printf ("COUNTEMP: counted %d templates\n")
		call pargi (ntemp)
	    call flush (STDOUT)
	    }
	return
end


# Set spectrum emission line chopping flag

procedure setschop (specim, tempim, itemp)

pointer specim		# Object image header structure
pointer tempim		# Template image header structure
int	itemp		# Index for this tmeplate in template list

real rejem,rejabs
bool	chop0		# TRUE if emission lines already removed from spectrum
char	chopstr[SZ_LINE]	# Chop out emission lines
				# (yes | no | tempfile | specfile)
real	clgetr()
int	strdic()

include "rvsao.com"
include "contin.com"
include "results.com"

begin

	schop = FALSE
	emrej[1] = clgetr ("s_em_reject")
	abrej[1] = clgetr ("s_abs_reject")
	call clgstr ("s_emchop",chopstr,SZ_LINE)
        sfchop = strdic (chopstr, chopstr, SZ_LINE, XC_CTYPES)

	switch (sfchop) {

	    case CHOP:
		schop = TRUE

	    case NOCHOP:
		schop = FALSE

	    case TEMPFILE:
		schop = FALSE
		call imgbpar (tempim,"CHOPEM",schop)

	# Read spectrum chopping limits from template header, if present
		if (schop) {
		    rejem = emrej[1]
		    call imgrpar (tempim,"EMREJ",rejem)
		    emrej[1] = rejem
		    rejabs = abrej[1]
		    call imgrpar (tempim,"ABSREJ",rejabs)
		    abrej[1] = rejabs
		    if (debug) {
			call printf ("SETSCHOP: Keep spectrum data from -%f to %f sigma\n")
			    call pargr (abrej[1])
			    call pargr (emrej[1])
			}
		    }

	    case SPECFILE:
		schop = TRUE
		call imgbpar (specim,"CHOPEM",schop)

	    default:
		schop = FALSE
	    }
	chop0 = FALSE
	call imgbpar (specim,"EMCHOP",chop0)
	if (chop0)
	    schop = FALSE
	if (abrej[1] < emrej[1]) {
	    tachop[itemp] = schop
	    tschop[itemp] = FALSE
	    }
	else {
	    tschop[itemp] = schop
	    tachop[itemp] = FALSE
	    }

	if (debug) {
	    if (schop)
		call printf ("SETSCHOP: Spectrum emission lines removed\n")
	    else
		call printf ("SETSCHOP: Spectrum emission lines left in\n")
	    }
	return
end


# Set template emission line chopping flag

procedure settchop (specim, tempim)

pointer specim		# Object image header structure
pointer tempim		# Template image header structure

real rejem,rejabs
bool	chop0		# TRUE if emission lines already removed from spectrum
char	chopstr[SZ_LINE]	# Chop out emission lines
				# (yes | no | tempfile | specfile)
real	clgetr()
int	strdic()

include "rvsao.com"
include "contin.com"
include "results.com"

begin
	tchop = FALSE
	emrej[2] = clgetr ("t_em_reject")
	abrej[2] = clgetr ("t_abs_reject")
	call clgstr ("t_emchop",chopstr,SZ_LINE)
        tfchop = strdic (chopstr, chopstr, SZ_LINE, XC_CTYPES)
	switch (tfchop) {
	    case CHOP:
		tchop = TRUE
	    case NOCHOP:
		tchop = FALSE
	    case TEMPFILE:
		tchop = TRUE
		call imgbpar (tempim,"CHOPEM",tchop)

	# Read template chopping limits from template header, if present
		if (tchop) {
		    rejem = emrej[2]
		    call imgrpar (tempim,"EMREJ",rejem)
		    emrej[2] = rejem
		    rejabs = abrej[2]
		    call imgrpar (tempim,"ABSREJ",rejabs)
		    abrej[2] = rejabs
		    if (debug) {
			call printf ("SETTCHOP: Keep template data from -%f to %f sigma\n")
			    call pargr (abrej[2])
			    call pargr (emrej[2])
			}
		    }
	    case SPECFILE:
		tchop = TRUE
		call imgbpar (specim,"CHOPEM",tchop)
	    }
	chop0 = FALSE
	call imgbpar (tempim,"EMCHOP",chop0)
	if (chop0)
	    tchop = FALSE
	if (debug) {
	    if (tchop)
		call printf ("SETTCHOP: Template emission lines removed\n")
	    else
		call printf ("SETTCHOP: Template emission lines left in\n")
	    }
	return
end

# Jul 13 1994	New subroutine
# Jul 29 1994	Move header writing to xchead
# Aug  3 1994	Fix erroneous subroutine calls
# Aug  3 1994	Set emission line chopping switch for each spectrum template
# Aug  8 1994	Add optional transform zero-padding
# Aug  8 1994	Add spectrum directory
# Aug 10 1994	Drop wavelength limits from calls to XCPLOT and XCRSLTS
# Aug 12 1994	Change name of header-writing subrooutine
# Aug 17 1994	Move ZPAD to labelled common
# Aug 18 1994	Fix VCOMBINE call
# Dec  8 1994	Move mwcs initialization to GETSPEC
# Dec 19 1994	Read filter parameters here instead of in main program

# Jan 31 1995	Change lengths of file and directory names to sz_pathname
# Feb 15 1995	Print multi-pass results only if DEBUG is on
# Feb 27 1995	Deal correctly with null spectrum and template directories
# Mar 10 1995	Add option to read zero-padding flag from template
# Mar 13 1995	Use CLOSE_IMAGE to make sure allocated space is freed
# Mar 13 1995	Allow template-directed line chopping
# Mar 13 1995	Rename zero_pad parameter for ease of relearning
# Mar 15 1995	Allow template-directed template continuum removal
# Mar 29 1995	Read rejection limits from xcsao.par, not contpars
# Apr  4 1995	Move stack pointer setting to start of subroutine
# Apr  6 1995	Change options for initial velocity setting
# Apr 11 1995	Chop lines by default if template or spectrum header not set
# May  8 1995	Don't add / to template directory if it ends in $
# May 15 1995	Change all sz_fname to sz_line, which is 100 chars longer
# May 16 1995	Add template-driven option to divide, not subtract, continuum
# May 18 1995	Zero best correlation before setting it
# Jun 21 1995	Fix bug initializing velocity
# Jun 26 1995	Add overlap header parameter to override wavelength limits
# Jul 13 1995	Add debugging argument to VCOMBINE
# Jul 19 1995	Call VWHEAD as well as XCWHEAD
# Aug 11 1995	Set peak fraction from template header or parameter file
# Sep 20 1995	Add CORRELATE parameter to bypass correlation for qplot
# Oct  3 1995	If not CORRELATE, don't save results in header unless QPLOT
# Oct  4 1995	Improve error message if no templates
# Oct 13 1995	Always call SFREE before returning
# Oct 13 1995	Print limits in angstroms; fix bug setting limits from .par

# Feb 16 1996	Check filter limits to keep them within spectrum
# Feb 22 1996	Set quality flag to 0 if correlating
# Aug  7 1996	Use smw.h

# Jan 31 1997	Add code to eliminate bad night sky lines
# Feb  5 1997	If initial velocity is INDEF, start at zero
# Mar 14 1997	Drop SPECSH from GETSPEC call; pass it in common
# Mar 14 1997	Label spectrum dispersion axis from spectrum header
# Apr  7 1997	Save ALL cross-correlations
# Apr 14 1997	Add option to plot template spectra
# Apr 25 1997	Add initial velocity option izguess for Z instead of velocity
# Apr 30 1997	Pass NCOR to XCFILE so proper number of points is plotted
# Aug 27 1997	Add band arguments for spectrum and template files
# Dec  4 1997	Add spectrum and template aperture number to XCFILE call
# Dec 22 1997	Make x-axis label work even if units aren't set in header

# Feb  9 1998	Make maximum number of correlation points 16384, not 8192
# Feb 11 1998	Allocate correlation and velocity per template arrays
# Feb 11 1998	Add code to count number of templates being correlated
# May 15 1998	Allow minvel and maxvel to be INDEFD
# Jun 12 1998	Use pixel limits of spectrum dispersion from header
# Dec 29 1998	Fix bug when using full template pathnames in list

# Mar 11 1999	Add write argument to getspec()
# Mar 17 1999	Multiply renormalized fluxes by 1000
# Apr  6 1999	Save templates between object spectra
# Jun 30 1999	Fix bug which zeroed template velocity and HCV
# Jul 16 1999	Fix bug so QPLOT will update headers
# Jul 29 1999	Fix problems with paths in template names
# Aug 18 1999	Use int conproc and tconproc instead of bool divcon and tdivcon
# Sep 15 1999	Read old results from header only if not correlating

# Jan 25 2000	Add option to not remove object spectrum continuum (SUBCONT=F)
# Jun 30 2000	Set renormalization flag if maximum counts/pixel is less than 1
# Aug  2 2000	Turn off template caching if echelle flag is on
# Aug  2 2000	Print template spectrum name, not file, on plot
# Sep 19 2000	Make correlate multi-valued, not just boolean
# Sep 27 2000	Read template wavelength/pixel limits from list file

# Jan 24 2001	Switch limits if spectrum is reversed
# Jul 31 2001	Add delpix so pixel correlations work

# Feb  5 2002	If ncol is INDEF, set to next higher power of 2 past npix
# Sep 17 2002	Allocate large vectors only once, on first call of xcfit()
# Sep 18 2002	Comment out all template caching code
# Dec  2 2002	Rewrite renormalization to match EMSAO
# Dec  2 2002	Plot original and renormalized spectrum if plot option on
# Dec  2 2002	Fix bad pixels before renormalization
# Dec  2 2002	Compute renormalization factor over selected spectrum portion

# Jun  3 2003	If renormalizing and minimum is <0, add pedestal pf spmin+1

# May 25 2005	Deal with pixel shift when rebinning and plotting
# Nov  3 2005	Force renormalization only if -1 < spectrum values < 1
# Dec 20 2005	Free wlwcs and twlwcs, which have been malloc'ed but not ever freed

# Jan 12 2006	If not write-enabled, do not attempt to update spectrum header
# Aug 17 2006	Drop pixels with a value of zero from the ends of the spectrum
# Sep  5 2006	Drop trailing spaces from header keyword queries
# Sep 27 2006	Fix bug in normalization

# Jan 26 2007	Do not change data structure when dropping zeroes at ends
# Jan 31 2007	Fix wavelength vector computation
# Jan 31 2007	Move line rejection code into subroutines
# Feb  1 2007	Reset initial czguess after return from xcplot if requested
# Jun 25 2007   Allocate all transform vectors instead of fixing dimension
# Aug  2 2007	Drop no-longer-used tmp* template caching buffers

# Mar  5 2008	Add pixfill to fill bad lines in templates if pixels used
# Mar  6 2008	Plot template if plttemp and pixfill
# Mar  7 2008	Drop use of minpix, maxpix, delpix
# Mar 10 2008	Set spwl1, spwl2, tmwl1, tmwl2 to pixel limits if COR_PIX
# Mar 10 2008	Add ispix flag for rebinning in pixels, not wavelength
# May  9 2008	Set contproc true if removing continuum from object spectrum
# May 30 2008	Fix filtering so that zero-padded spectra are not filtered

# Mar  4 2009	Improve log rebinning diagnostics
# Mar 13 2009	Add oshift argument to cross-correlate shifted echelle orders
# Jul 14 2009	Set spectrum continuum division from spectrum keyword DIVCONT

# Jul 09 2010	End properly if spectrum error
# Jul 22 2010	Set allocated array pointers to NULL when not allocated

# May  1 2014	Write out spectrum if it has been changed
