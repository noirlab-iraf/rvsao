# File rvsao/Xcor/t_pxcsao.x
# May 1, 2014
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss, Johns Hopkins University and others

# Copyright(c) 2007-2014 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# PXCSAO is an IRAF task for obtaining redshifts and velocity dispersions
# using cross correlation methods.  For arguments, see parameter file xcsao.par.
# Information is shared in common blocks defined in "rvsao.com".
# It differs from XCSAO only in that it saves the fit velocity as a parameter
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include	<smw.h>
include "rvsao.h"
include	"xcv.h"
include	"emv.h"

procedure t_pxcsao ()

int	i
char	specfile[SZ_PATHNAME]	# Object spectrum file name
char	specpath[SZ_PATHNAME]	# Object spectrum path name
char	tempfiles[SZ_PATHNAME]	# List of template spectra
char	specdir[SZ_PATHNAME]	# Directory for object spectra
#char	filt_type[SZ_LINE]	# Filter for transform
				# (ramp | welch | hanning | cos-bell)
char	svel_corr[SZ_LINE]	# Type of velocity correction for spectrum
				# (none | file | heliocentric | barycentric)
char	tvel_corr[SZ_LINE]	# Type of velocity correction for template
				# (none | file | heliocentric | barycentric)
bool	savevel0		# Save velocity, error, and R in data file header
int	rmode		# Report format (1=normal,2=one-line)
int	logfiles	# List of log files
char	logfile[SZ_PATHNAME] # Log file name
char	wtitle[20]	# Title for wavelength plots of spectrum

int	mspec		# Object aperture to read from multispec file
int	mband		# Object band to read from multispec file

pointer	speclist	# List of spectrum files
char	str[SZ_LINE]
int	fd
char	vel_plot[SZ_LINE]	# type of velocity for redshifting plot
				# correlation|emission|combination|search

int	nmspec0		# Number of object multispec spectra
int	mspec_range[3,MAX_RANGES]
int	ip,jp,lfile	# Limits for multispec aperture decoding
char	lbracket[3]	# "[({"
char	rbracket[3]	# "])}"
double	sumvel, sumerr, sumr, avgvel, avgerr, avgr, dnap
double	medvel, qvel1, qvel2
int	nap
pointer	evel
char	vel_init[SZ_LINE]	# type of velocity for initial value
bool	echelle         # If true, template multispec numbers track object
char	tempnums[SZ_LINE]       # List of multispec spectra to read
int	tspec_range[3,MAX_RANGES]
int	tspec
int	ntspec
int	oshift  
 
bool	clgetb()
int	clpopnu(), clgeti(), clgfil(), open()
int	strdic(), stridx(), stridxs()
real	clgetr()
int	decode_ranges(),get_next_number()
int	imtgetim(), imaccess(), strlen(), ldir, clscan()
pointer	imtopenp()

define	newspec_ 10
define	newap_	 20
define	endxc_	 90

include	"rvsao.com"
include "results.com"
include	"emv.com"
include	"xcor.com"
include	"xcorf.com"
include	"xplt.com"
 
begin
	c0 = 299792.5
	qplot = FALSE
	nfound = 0
	maxpix = 0
	maxpts4 = 0
	call sprintf (lbracket,3,"[({")
	call sprintf (rbracket,3,"])}")
	call sprintf (wtitle,20,"Wavelength")
	ntmp = 0
	waverest = 0.d0
	specref = 0

# Initialize various utility vector pointers
	xcor = NULL
	xvel = NULL
	shspec = NULL
	shtemp = NULL
	wltemp = NULL
	xind = NULL
	xifft = NULL
	ft1 = NULL
	ft2 = NULL
	ftcfn = NULL
	tft = NULL
	pft = NULL
	spexp = NULL
	xcont = NULL
	scont = NULL
	smspec = NULL
	cspec = NULL
	smcspec = NULL
	fraclev = NULL
	xlev = NULL
	evel = NULL

# Get task parameters.

# Spectra to cross-correlate
	speclist = imtopenp ("spectra")

# Multispec spectrum numbers (use only first if multiple files)
	call clgstr ("specnum",specnums,SZ_LINE)
	if (decode_ranges (specnums, mspec_range, MAX_RANGES, nmspec0) == ERR){
	    call sprintf (str, SZ_LINE, "T_XCSAO: Illegal multispec list <%s>")
		call pargstr (specnums)
	    call error (1, str)
	    }
	call clgstr ("specdir",specdir,SZ_PATHNAME)
	ldir = strlen (specdir)
	if (specdir[1] != EOS && specdir[ldir] != '/') {
	    specdir[ldir+1] = '/'
	    specdir[ldir+2] = EOS
	    }
	mband = clgeti ("specband")
	mext = clgeti ("specext")

# Templates against which to correlate spectra
	call clgstr ("templates",tempfiles,SZ_PATHNAME)

# Optional correlation plot, where peak may be selected by cursor
	pltcor  = clgetb ("xcor_plot")

# Optional intermediate data plot switches
	pltspec = clgetb ("obj_plot")
	plttemp = clgetb ("temp_plot")
	pltcon  = clgetb ("contsub_plot")
	pltapo  = clgetb ("apodize_plot")
	pltfft  = clgetb ("fft_plot")
	pltuc   = clgetb ("uxcor_plot")
	plttft  = clgetb ("tfft_plot")

# Print processing information
	debug  = clgetb ("debug")

# Continuum fit parameter pset
	call cont_get_pars()

# Number of times to smooth (1-2-1) final data plot
	nsmooth = clgeti ("nsmooth")

# Velocity center and width of summary page cross-correlation plot
	xcr0 = clgetr ("cvel")
	xcrdif = clgetr ("dvel")

# Type of fit for correlation peak and fraction of peak to fit
	pkmode0 = clgeti ("pkmode")

# Type of heliocentric velocity correction to be used
	call clgstr ("svel_corr",svel_corr,SZ_LINE)
	svcor = strdic (svel_corr,svel_corr,SZ_LINE, HC_VTYPES)
	call clgstr ("tvel_corr",tvel_corr,SZ_LINE)
	tvcor = strdic (tvel_corr,tvel_corr,SZ_LINE, HC_VTYPES)

# Type of velocity for initial redshift
	call clgstr ("vel_init",vel_init,SZ_LINE)
	vinit = strdic (vel_init,vel_init,SZ_LINE,XC_VTYPES)

# Image header result flag
	savevel0 = FALSE
	savevel0 = clgetb ("save_vel")

# Report mode for log file
	rmode = 1
	rmode = clgeti ("report_mode")

# Maximum length for filename in output
	maxlength = 0
	maxlength = clgeti ("maxlength")

# Initialize emission and absorption lines for labelling
	call eminit (FALSE)
 
# Open log files and write a header.
	logfiles = clpopnu ("logfiles")
	call fseti (STDOUT, F_FLUSHNL, YES)
	i = 0
	call strcpy ("rvsao.xcsao",taskname,SZ_LINE)
	while (clgfil (logfiles, logfile, SZ_PATHNAME) != EOF) {
	    fd = open (logfile, APPEND, TEXT_FILE)
	    if (fd == ERR) break
	    if (rmode == 1) {
		call loghead (taskname,str)
		call fprintf (fd, "%s\n")
		    call pargstr (str)
		}
	    i = i + 1
	    logfd[i] = fd
	    }
	nlogfd = i
	call clpcls (logfiles)

# Type of velocity for plotting emission and absorption lines
	if (clscan("vel_plot") != EOF) {
	    call clgstr ("vel_plot",vel_plot,SZ_LINE)
	    vplot = strdic (vel_plot,vel_plot,SZ_LINE,PL_VTYPES)
	    }
	else
	    vplot = VCORREL

# Set echelle order shift, if any
	echelle = clgetb ("echelle")
	oshift = 0
	if (echelle) {
	    call clgstr ("tempnum",tempnums,SZ_LINE)
	    if (strlen (tempnums) > 0) {
		if (decode_ranges (tempnums, tspec_range, MAX_RANGES, ntspec) == ERR)
		call error (1, "T_XCSAO: Illegal template multispec list")
		if (get_next_number (tspec_range, tspec) != EOF) {
		    if (get_next_number (mspec_range, mspec) != EOF) {
			oshift = tspec - mspec
			if (debug) {
			    call printf ("T_XCSAO: Shifting by %d orders\n")
				call pargi (oshift)
			    }
			}
		    }
		}
	    call malloc (evel, nmspec0+5, TY_DOUBLE)
	    }

# Print tab table headers
	if (rmode < 0)
	    call xcrshead (rmode)

# Get next object spectrum file name from the list
newspec_
	if (imtgetim (speclist, specfile, SZ_PATHNAME) == EOF)
	   go to endxc_

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
	    call printf ("XCSAO: next file is %s [%s] = %d aps\n")
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
	sumvel = 0.d0
	sumerr = 0.d0
	sumr = 0.d0
	dnap = 0.d0
	nap = 0
newap_
	savevel = savevel0
	savespec = FALSE
	edited = FALSE
	if (nmspec <= 0)
	    go to newspec_
	if (get_next_number (mspec_range, mspec) == EOF)
	    go to newspec_

	call xcfit (specfile, specdir, mspec, mband, oshift)

	if (echelle) {
	    nap = nap + 1
	    Memd[evel+nap-1] = zvel[itmax]
	    sumvel = sumvel + zvel[itmax]
	    sumerr = sumerr + czerr[itmax]
	    sumr = sumr + czr[itmax]
	    dnap = dnap + 1.d0
	    }

# Move on to next aperture or next image
	nmspec = nmspec - 1
	if (nmspec > 0)
	    go to newap_
	go to newspec_
 
# Close the log files
endxc_	do i = 1, nlogfd {
	    call close (logfd[i])
	    }

	if (debug) {
	    call printf ("Best Template(%d) %s: %.3f (%.3f) R = %.4f\n")
		call pargi (itmax)
		call pargstr (tempid[1,itmax])
		call pargd (zvel[itmax])
		call pargd (czerr[itmax])
		call pargd (czr[itmax])
	    }

# Save some results in parameter file
	call clpstr ("besttemp", tempid[1,itmax])
	if (echelle && dnap > 0.d0) {
	    avgvel = sumvel / dnap
	    call clputd ("velocity", avgvel)
	    avgerr = sumerr / dnap
	    call clputd ("velerr", avgerr)
	    avgr = sumr / dnap
	    call clputd ("r",avgr)
	    if (debug) {
		call printf ("PXCSAO: Median of %d shifts\n")
		    call pargi (nap)
		call flush (STDOUT)
		}
	    call median (Memd[evel], nap, medvel, qvel1, qvel2)
	    if (debug) {
		call printf ("PXCSAO: %d shifts mean=%.3f, median=%.3f\n")
		    call pargi (nap)
		    call pargd (avgvel)
		    call pargd (medvel)
		call flush (STDOUT)
		}
	    call clputd ("velmed",medvel)
	    call clputd ("velq1",qvel1)
	    call clputd ("velq2",qvel2)
	    }
	else {
	    call clputd ("velocity", zvel[itmax])
	    call clputd ("velerr", czerr[itmax])
	    call clputd ("r",czr[itmax])
	    }
	if (debug) {
	    call printf ("PXCSAO: About to close speclist\n")
	    call flush (STDOUT)
	    }

#  Close spectrum list
	call imtclose (speclist)

# Free processing vectors allocated in xcfit()
	if (debug) {
	    call printf ("PXCSAO: About to free xcfit vectors\n")
	    call flush (STDOUT)
	    }
	if (xcor != NULL) {
	    call mfree (xcor, TY_REAL)
	    xcor = NULL
	    }
	if (xcor != NULL) {
	    call mfree (xvel, TY_REAL)
	    xvel = NULL
	    }
	if (shspec != NULL) {
	    call mfree (shspec, TY_REAL)
	    shspec = NULL
	    }
	if (shtemp != NULL) {
	    call mfree (shtemp, TY_REAL)
	    shtemp = NULL
	    }
	if (wltemp != NULL) {
	    call mfree (wltemp, TY_REAL)
	    wltemp = NULL
	    }

# Free processing vectors allocated in xcorfit()
	if (debug) {
	    call printf ("PXCSAO: About to free xcorfit vectors\n")
	    call flush (STDOUT)
	    }
	if (xind != NULL) {
	    call mfree (xind, TY_REAL)
	    xind = NULL
	    }
	if (xifft != NULL) {
	    call mfree (xifft, TY_REAL)
	    xifft = NULL
	    }
	if (pft != NULL) {
	    call mfree (pft, TY_REAL)
	    pft = NULL
	    }
	if (tft != NULL) {
	    call mfree (tft, TY_COMPLEX)
	    tft = NULL
	    }
	if (ftcfn != NULL) {
	    call mfree (ftcfn, TY_COMPLEX)
	    ftcfn = NULL
	    }
	if (ft1 != NULL) {
	    call mfree (ft1, TY_COMPLEX)
	    ft1 = NULL
	    }
	if (ft2 != NULL) {
	    call mfree (ft2, TY_COMPLEX)
	    ft2 = NULL
	    }
	if (spexp != NULL) {
	    call mfree (spexp, TY_REAL)
	    spexp = NULL
	    }
	if (xcont != NULL) {
	    call mfree (xcont, TY_REAL)
	    xcont = NULL
	    }

# Free plotting vectors allocated in xcplot()
	if (debug) {
	    call printf ("PXCSAO: About to free xcplot vectors\n")
	    call flush (STDOUT)
	    }
	if (scont != NULL) {
	    call mfree (scont, TY_REAL)
	    scont = NULL
	    }
	if (smspec != NULL) {
	    call mfree (smspec, TY_REAL)
	    smspec = NULL
	    }
	if (cspec != NULL) {
	    call mfree (cspec, TY_REAL)
	    cspec = NULL
	    }
	if (smcspec != NULL) {
	    call mfree (smcspec, TY_REAL)
	    smcspec = NULL
	    }

# Free plotting vectors allocated in xcorplot()
	if (debug) {
	    call printf ("PXCSAO: About to free xcorplot vectors\n")
	    call flush (STDOUT)
	    }
	if (xlev != NULL) {
	    call mfree (xlev, TY_REAL)
	    xlev = NULL
	    }
	if (fraclev != NULL) {
	    call mfree (fraclev, TY_REAL)
	    fraclev = NULL
	    }

# Free vector used to compute echelle median shift in t_pxcsao
	if (debug) {
	    call printf ("PXCSAO: About to free echelle result vector\n")
	    call flush (STDOUT)
	    }
	if (evel != NULL) {
	    call mfree (evel, TY_DOUBLE)
	    evel = NULL
	    }
	if (debug) {
	    call printf ("PXCSAO: Exiting\n")
	    call flush (STDOUT)
	    }

end

# Jul 31 2007	New task based on xcsao which writes results to parameters
# Aug  2 2007	Drop nctemp as templates are no longer cached
# Aug 13 2007	Free xcont and spexp at end of program
# Aug 13 2007	Set maxpix = 0 so plot buffers are initialized when needed
# Aug 13 2007	Include xplt.com to allocate buffers only when needed
# Oct 29 2007	Save mean of aperture velocities in parameter file if echelle

# Apr 23 2008	Add median and quartile velocities to ouput parameters
# Aug  5 2008	Move median subroutine to Util/ directory as separate file

# Mar 13 2009	Compute order shift and send to xcfit 

# Jul 08 2010	Read parameter specext into mext in rvsao.com
# Jul 22 2010	Initialize all pointers and free any which are not set to NULL

# May  1 2014	Initialize edited and savespec to FALSE for each spectrum
