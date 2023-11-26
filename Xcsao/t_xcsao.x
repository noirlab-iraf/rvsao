# File rvsao/Xcor/t_xcsao.x
# May 1, 2014
# By Jessica Mink, Harvard-Smithsonian Center for Astrophysics
# After Gerard Kriss, Johns Hopkins University and others

# Copyright(c) 1989-2014 Smithsonian Astrophysical Observatory
# You may do anything you like with this file except remove this copyright.
# The Smithsonian Astrophysical Observatory makes no representations about
# the suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
 
# XCSAO is an IRAF task for obtaining redshifts and velocity dispersions
# using cross correlation methods.  For arguments, see parameter file xcsao.par.
# Information is shared in common blocks defined in "rvsao.com".
 
include	<imhdr.h>
include	<imio.h>
include	<fset.h>
include	<smw.h>
include "rvsao.h"
include	"xcv.h"
include	"emv.h"

procedure t_xcsao ()

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
bool	echelle		# If true, template multispec numbers track object
char	tempnums[SZ_LINE]	# List of multispec spectra to read
int	tspec_range[3,MAX_RANGES]
int	tspec
int	ntspec
int	oshift
 
bool	clgetb()
int	clpopnu(), clgeti(), clgfil(), open()
int	strdic(), stridx(), stridxs()
real	clgetr()
int	decode_ranges(),get_next_number()
char	vel_init[SZ_LINE]	# type of velocity for initial value
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
	maxlength = 0

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
	spexp = NULL
	xcont = NULL
	scont = NULL
	smspec = NULL
	cspec = NULL
	smcspec = NULL
	fraclev = NULL
	xlev = NULL

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
newap_
	savevel = savevel0
	savespec = FALSE
	edited = FALSE
	if (nmspec <= 0)
	    go to newspec_
	if (get_next_number (mspec_range, mspec) == EOF)
	    go to newspec_

	call xcfit (specfile, specdir, mspec, mband, oshift)

# Move on to next aperture or next image
	nmspec = nmspec - 1
	if (nmspec > 0)
	    go to newap_
	go to newspec_
 
# Close the log files
endxc_	do i = 1, nlogfd {
	    call close (logfd[i])
	    }                                              

#  Close spectrum list
	call imtclose (speclist)

# Free processing vectors allocated in xcfit()
	if (debug) {
	    call printf ("XCSAO: About to free xcfit vectors\n")
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
	    call printf ("XCSAO: About to free xcorfit vectors\n")
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
	    call printf ("XCSAO: About to free xcplot vectors\n")
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
	    call printf ("XCSAO: About to free xcorplot vectors\n")
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

end
 
# Mid 1970's	Paul Schechter -- Data windowing routines.  Some functions
#				  from the IBM SSP library.
# June	1987	Gerard Kriss -- Wrote major portions in SPP.
# June	1988	Stephen Levine -- major revisions
# Sept	1989	Doug Mink--add parameters for 2-d echelle
# April	1990	add more fft filters and bcv
#		add option to resize velocity graph
# June	1990	add emission line chopping
# July	1990	rebin spectrum for each template
#		add John Tonry's fft
#		add Guillermo Torres' peak fit
# Sept	1990	add optional archive record output
# Oct	1990	add pixel limits for redshift
#		add normalization of spectra
#		add smoothing to spectrum display
# Dec	1990	add archive flag
# Jan	1991	write velocity to image header
#		add velocities or multiply 1+v/c's
# Mar	1991	plot final velocity correlation
# June	1991	make all velocities double
#		add IRAF continuum fitting
# Sept	1991	add second pass to align template
#		add plot options at end
# Nov 14 1991	Change label strings to char variables instead of literals
#		move vcombine and write combined velocities to header
# Nov 15 1991	Clean up debug formats
# Nov 18 1991	Change lsrplot calls to spplot and xcplot
# Nov 20 1991	Free str and logfile if necessary
# Dec  5 1991	Set velxc, vxerr, and vr for absorbtion line plot
# Dec 12 1991	Free object and template spectra when done
# Dec 16 1991	Set wavelength vector for region being cross-correlated

# Feb 18 1992	If multispec, write results to aperture-dependent strings
# Mar 27 1992	Get shift limits in velocity; get velocity vector from xcorfit
# Apr 20 1992	Pass mspec as argument
# Apr 22 1992	Put all cursor interaction and plotting into xcplot
# May 22 1992	Initialize mspec and tspec to -1
# Aug 12 1992	Add HISTORY line about XCSAO; drop polynomial continuum option
# Oct  9 1992	Read "echelle" parameter
# Nov 30 1992	Move spectrum smoothing to xcplot

# Feb  1 1993	Allow INDEF for wavelength limit parameters
# Feb  2 1993	Exit if a zero spectrum is found when renormalizing
# May  4 1993	Add option to save correlation to a file
# May 10 1993	Fix correlation writing
# May 18 1993	Move correlation file output to subroutine xcprint
# Jun  2 1993	Move wavelength <-> pixel conversion to subroutines
# Jun 14 1993	Use run-time allocation for all spectrum and wavelength vectors
# Jun 16 1993	Set version to 1.1; use shdr for wavelength conversions
# Jul  1 1993	Pass spectrum header structure to rebinl
# Jul  7 1993	Add spectrum header to getspec and gettemp
# Jul  8 1993	Close sh properly
# Aug  4 1993	Print error message and exit if there is no overlap
# Aug 20 1993	Deal with double wavelength vector
# Dec  2 1993	Only write history line once for multispec files
# Dec  2 1993	Pass mspec to xc_rslts
# Dec  3 1993	Update to version 1.2

# Feb  3 1994	Proceed gracefully if files are not useable
# Feb 11 1994	Add correlation file name list.
# Mar 23 1994	Add multispec aperture specification for each image
# Mar 23 1994	Update to version 1.3
# Apr  6 1994	Change smooth parameter to nsmooth for consistency
# Apr 12 1994	Return MWCS header pointer from getspec and gettemp
# Apr 13 1994	Drop unused variable imax
# Apr 15 1994	Initialize world coordinate system for spectrum
# Apr 19 1994	Pass object spectrum label as variable, not literal
# Apr 20 1994	Drop 2nd argument from IMTOPEN call
# Apr 26 1994	Use spectrum header wavelength vector directly
# Apr 26 1994	Use arrays rather than pointers when possible
# May  9 1994	Add number of times to smooth spectrum argument to PLOTSPEC call
# May 23 1994	Fix error message when not able to write to file
# May 23 1994	Keep smoothed spectrum local to XCPLOT
# Jun 15 1994	Move filtered-template flag from parameter list to image file
# Jun 15 1994	Ignore template directory for template names with / in them
# Jun 23 1994	Keep MWCS pointer in SHDR structure
# Jun 24 1994	Set spectrum velocity in common
# Jun 24 1994	Update to version 1.4
# Jul 29 1994	Read version from fquot.h
# Aug  3 1994	Change emission line chopping flags
# Aug  3 1994	Change common and header from fquot to rvsao
# Aug  4 1994	Keep nmspec in labelled common
# Aug  8 1994	Add specdir as parameter
# Aug 10 1994	Write log heading in one line instead of two
# Aug 24 1994	Fix bug so SPECDIR can be null
# Dec 19 1994	Move renorm and filter parameters to XCFIT

# Jan 31 1995	Change lengths of file and directory names to sz_pathname
# Mar  9 1995	Move zero-padding to XCFIT
# Apr  5 1995	Change solar system velocity correction flags
# May 15 1995	Change all sz_fname to sz_line, which is 100 chars longer
# May 17 1995	Fix bug in bracket -> mspec code
# Aug 11 1995	Move PKFRAC to XCFIT
# Oct  4 1995	Default SPECNUMS to 0 if null string is input

# Aug  7 1996	Use smw.h

# Jan 24 1997	Drop phase plot; add transformed transform plot
# Apr  9 1997	Deal with null aperture list in DECODE_RANGES, not here
# Apr 14 1997	Add option to plot template spectra
# Aug 27 1997	Add parameter for spectrum multspec band

# Apr  6 1999	Initialize number of cached templates to zero
# Jul 16 1999	Initialize qplot to FALSE
# Jul 28 1999	Print tab table headers, if appropriate
# Aug 19 1999	Set vplot from vel_plot parameter

# Sep 17 2002	Initialize and free processing vectors (allocated in xcfit)

# Dec  1 2006	Fix bug so apertures can be in brackets in file name

# Jun 20 2007	Initialize and free transform vectors in xcorf.com
# Jun 25 2007	Add xcont and spexp pointer initialization
# Aug  2 2007	Drop nctemp as templates are no longer cached
# Aug 13 2007	Free xcont and spexp at end of program
# Aug 13 2007	Set maxpix = 0 so plot buffers are initialized when needed
# Aug 13 2007	Include xplt.com to allocate buffers only when needed

# Mar 13 2009	Compute order shift and send to xcfit 

# Jul 22 2010   Initialize all pointers and free any which are not set to NULL

# May  1 2014	Initialize edited and savespec to FALSE for each spectrum
